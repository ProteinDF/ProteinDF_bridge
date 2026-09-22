# TASK: `AtomGroup::setup()`を賢いデフォルトの結合解決エントリポイントにする

> ユーザーから「呼び出し側からすると`Bond::setup()`を呼べばデフォルトで良きに計らってほしい。ヒューリスティックの方を別名にして`setup()`一発で済むようにしてほしい」との要望があった。さらに「AtomGroupの内容更新時に暗黙的に走ってほしい」との要望もあったが、性能・正確性・実装コストの観点から「ファイルローダーの読み込み完了」という区切りに限定した暗黙実行に落とし込んだ。詳細は`RUST_PORT_SPEC.md` §3.14を参照。

## ブランチ運用(MUST)

- `develop` から `feature/setup-default-entrypoint` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。
- 依存関係: §3.13(`AtomGroup::resolve_bonds`)がマージ済みであることを前提にする（マージ済み: `develop`コミット`37b0c6d`）。

## 背景

現状、結合解決には2つのAPIが存在し、どちらを呼ぶべきか呼び出し側から見て分かりにくい:
- `Bond::setup(&mut self, mol: &mut AtomGroup) -> Result<()>`（`bond.rs`、純粋な共有結合半径ヒューリスティックのみ）
- `AtomGroup::resolve_bonds(&mut self, db: &CcdTemplateDb) -> Result<()>`（`atom_group.rs`、CCDテンプレート→ヒューリスティックの賢い版）

「`setup`という名前」＝「呼べば良きに計らってくれる賢いデフォルト」という直感を実現するため、以下の方針でリネーム・整理する。**`Bond::setup()`自体をCCDテンプレート対応にする（`bond.rs`が`ccd_templates.rs`に依存する形にする）ことは行わない**（`bond.rs`は基礎プリミティブとして`ccd_templates.rs`に依存しないレイヤリングを維持する、という§3.13の既存方針を踏襲）。

## 対象

### 1. `bond.rs`のリネーム

- `Bond::setup()` を `Bond::setup_heuristic()` にリネームする。ロジックは一切変更しない（純粋な共有結合半径ヒューリスティックのみ、という意味を保つ）。
- `rust/crates/proteindf-bridge-py/src/bond.rs`（Pythonバインディング）内の呼び出し箇所（`self.inner.setup(&mut mol.inner)`）も新名称に追従させる。Python側の公開メソッド名（PyO3の`#[pyo3(name = ...)]`等）をどうするかは、既存のPython API命名規則（`proteindf_bridge`の既存メソッド名踏襲方針、`RUST_PORT_SPEC.md` §4参照）に従って判断してよいが、**Python側の公開名を変更する場合は変更内容をユーザー経由でClaudeに報告すること**（互換性への影響があるため）。
- `bond.rs`内のテスト・ベンチマーク的記述（`bond.setup(&mut ag)`呼び出し、`"Bond::setup() time for..."`のようなログ文言）も新名称に更新する。

### 2. `AtomGroup`側APIの整備(`atom_group.rs`)

- `AtomGroup::resolve_bonds(&mut self, db: &CcdTemplateDb) -> Result<()>` を `AtomGroup::setup_with_db(&mut self, db: &CcdTemplateDb) -> Result<()>` にリネームする。ロジックは変更しない（§3.11の拡張DBを使う場合の明示的エントリポイントとして残す）。
- 新規に `AtomGroup::setup(&mut self) -> Result<()>` を追加する。内部で `self.setup_with_db(CcdTemplateDb::global())` を呼ぶだけの薄いラッパーとする。

### 3. 各フォーマットローダーでの暗黙実行

`get_atomgroup()`（またはパース処理の完了直前）で、`ag.get_bond_list().is_empty()` の場合のみ `ag.setup()?` を呼んでから返すようにする。対象:

- `format/pdb.rs`の`get_atomgroup()`
- `format/mmcif.rs`の`get_atomgroup()`
- `format/amber_prmtop.rs`の`get_atomgroup()`
- `format/gro.rs`の`get_atomgroup()`
- `format/mol2.rs`は`get_atomgroup()`が`&AtomGroup`（参照）を返す設計のため、`parse_str()`内で`self.set_by_atomgroup(&ag)`を呼ぶ直前に同様のチェックを行う。

いずれも「ファイル由来の結合が既にある場合は何もしない」ため、MOL2/PRMTOP/PDB(CONECT)等、明示的結合情報を持つフォーマットでは実質的にno-opとなり、§3.8の優先順位方針（ファイル由来結合は上書きしない）を壊さないこと。**この前提が崩れていないか、明示的結合情報を持つフォーマットの既存テストが引き続きパスすることで必ず確認すること。**

### 4. ドキュメント更新

- `RUST_PORT_SPEC.md` §3.8の呼び出し側推奨パターンは既にこのタスクの完成形を前提とした記述に更新済み（`ag.setup()`/`ag.setup_with_db()`/`Bond::setup_heuristic()`を参照する形）。実装がこの記述と整合していることを確認すること。
- 各ローダーファイル内の`Bond::setup()`を参照するdocコメント（`mol2.rs`・`pdb.rs`・`amber_prmtop.rs`）も新名称（`Bond::setup_heuristic()`または`AtomGroup::setup()`、文脈に応じて適切な方）に更新する。
- `RUST_PORT_SPEC.md` §3.14に実施内容・完了ステータスを追記すること。

## 完了の定義(Definition of Done)

- [x] 1. `Bond::setup()`への参照が名称`Bond::setup_heuristic()`に統一され、`cargo build --workspace`（Pythonバインディング含む）が通ること。
- [x] 2. `AtomGroup::setup()`（引数なし）を呼ぶだけで、CCDテンプレート＋ヒューリスティックによる結合解決が行われることを検証する回帰テストを追加すること（既存の`test_resolve_bonds_*`系テストを新API名に追従させる形でよい）。
- [x] 3. 各ローダー（PDB・mmCIF・PRMTOP・GRO・MOL2）について、明示的結合情報を持たない入力に対して`get_atomgroup()`を呼んだだけで結合が自動解決されていることを検証する回帰テストを、フォーマットごとに最低1件追加すること。
- [x] 4. 明示的結合情報を持つ入力（PDBのCONECT、MOL2のBONDセクション、PRMTOPのBONDS等）に対して`get_atomgroup()`を呼んでも、ファイル由来の結合が上書き・重複されないことを確認する既存テストが引き続きパスすること。
- [x] 5. `cargo clippy --workspace --all-targets -- -D warnings` / `cargo fmt --check` を通すこと。
- [x] 6. `RUST_PORT_SPEC.md` §3.14に実施内容・完了ステータスを追記すること。

## スコープ外

- CCDテンプレート自体・共有結合半径テーブル自体の変更。
- 「AtomGroupの内容更新のたびに暗黙的に再実行する」方式の実装（§3.14の背景で説明した理由により、本タスクでは採用しない。ローダー完了時点に限定した暗黙実行のみを行うこと）。
- 既存Pythonコード(`proteindf_bridge/`)は変更しない。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- `Bond::setup_heuristic()`・`apply_ccd_bond_templates()`・`AtomGroup::setup_with_db()`自体のロジック（結合検出アルゴリズム、既存結合の上書き防止ロジック）を変更しない。本タスクは名称変更と呼び出しタイミングの整理のみを対象とする。
- ローダーの暗黙`setup()`呼び出しを、「結合情報の有無を確認せず常に呼ぶ」形にしない（`ag.get_bond_list().is_empty()`のチェックを必ず経由すること。これを怠るとファイル由来結合を上書きする重大なバグになる）。
