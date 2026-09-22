# TASK: 結合解決の統合エントリポイント(`AtomGroup::resolve_bonds`)

> 本タスクは `RUST_PORT_SPEC.md` §3.12(2026-09-22追記)の計画を実装するものです。背景・設計方針・スコープの詳細は同節を参照してください。

## ブランチ運用(MUST)

- `develop` から `feature/unified-bond-resolution` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。
- **重要な依存関係(MUST)**: 本タスクは `fix/bond-setup-duplicate-bonds`(`Bond::setup()`が既存結合を重複登録するバグの修正)が`develop`にマージされていることを前提とする。**そのタスクがまだ`develop`にマージされていない場合は、本タスクに着手せず、マージされるまで待つこと。** マージ済みであることを`git log develop`等で確認してから`develop`から本タスクのブランチを切ること。

## 着手前に必ず読むこと

- `RUST_PORT_SPEC.md` §3.8(結合情報の優先順位方針)・§3.9(CCDテンプレートDB)・§3.12(本タスクの計画)。
- `rust/crates/proteindf-bridge/src/bond.rs`の`Bond::setup()`(修正済みの重複防止ロジックを含む)。
- `rust/crates/proteindf-bridge/src/atom_group.rs`の`apply_ccd_bond_templates`(既存結合を上書きしないロジック)。

## 対象

### 1. 統合エントリポイントの追加

- `AtomGroup::resolve_bonds(&mut self, db: &CcdTemplateDb) -> Result<()>`のようなメソッドを`atom_group.rs`に追加する(名前は実装者の判断でよいが、`RUST_PORT_SPEC.md` §3.12の意図が伝わる名前にすること)。
- 内部処理:
  1. `self.apply_ccd_bond_templates(db)`を呼ぶ。
  2. `Bond::setup()`を呼ぶ(`fix/bond-setup-duplicate-bonds`の修正により、ステップ1で登録済みの結合は重複せず、残りの原子ペアのみヒューリスティックで補完される前提)。
- ファイル由来の明示的結合(CONECT等、呼び出し前から`self`に既に存在する結合)は、`apply_ccd_bond_templates`・`Bond::setup()`いずれも上書きしないため、自然に最優先で尊重される。追加の特別処理は不要なはずだが、実装時に動作を確認すること。

### 2. `db`引数による拡張性

- 引数として`&CcdTemplateDb`を受け取ること(内部で`CcdTemplateDb::global()`を決め打ちしない)。呼び出し側が組み込みデフォルトDBをそのまま渡すことも、§3.11の`insert`/`merge`で拡張した独自DBを渡すこともできるようにする。

### 3. `RUST_PORT_SPEC.md` §3.8の更新

- 「呼び出し側(「結 (YUI)」等)の推奨利用パターン」のコード例を、新しい`resolve_bonds()`を使う形に更新する(既存の`Bond::setup()`単体の説明は残しつつ、統合エントリポイントを推奨する形にすること)。

## 完了の定義(Definition of Done)

1. 実PDBフィクスチャ(`1hls.pdb`)で、`resolve_bonds()`を1回呼ぶだけで、CCDテンプレートによる結合次数(GLU側鎖CD=OE1等の二重結合)とヒューリスティックによる補完(ペプチド結合等、CCDテンプレートでは対象にならない残基間結合)の両方が、**重複なく**正しく得られることを検証する回帰テストを追加すること。
2. ファイル由来の結合(CONECT等)がある`AtomGroup`に対して`resolve_bonds()`を呼んでも、その結合が上書き・重複されないことを検証すること。
3. `Bond::setup()`単体・`apply_ccd_bond_templates()`単体を直接呼ぶ既存のテスト(§3.9フェーズA・§3.10)が、実装追加後も引き続き全てパスすることを確認すること(既存APIの意味を変えていないことの回帰確認)。
4. `cargo clippy` / `cargo fmt` を通すこと。
5. `RUST_PORT_SPEC.md` §3.12・§3.8に実施内容・完了ステータスを追記すること。

## スコープ外

- `Bond::setup()`自体の変更(純粋なヒューリスティックのままとする)。
- `apply_ccd_bond_templates()`自体の変更。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- `Bond::setup()`を「CCDテンプレートも内部で使う」ように変更しない(§3.12の背景で説明した通り、既存のベースラインテストが壊れるため)。
- `fix/bond-setup-duplicate-bonds`がマージされる前に本タスクに着手しない。
