# TASK: CCD結合テンプレートDBの実行時拡張(ユーザー提供の外部CCDデータ)

> 本タスクは `RUST_PORT_SPEC.md` §3.11(2026-09-22追記)の計画を実装するものです。背景・設計方針・スコープの詳細は同節を参照してください。§3.9フェーズA(組み込み29残基DB)・§3.10(共有結合半径ベース判定)とは独立です。

## ブランチ運用(MUST)

- `develop` から `feature/ccd-runtime-extension` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。
- 依存関係: なし。`format/mmcif.rs`のCCDパース(2026-09-22、`AROM`/`QUAD`結合次数対応済み)・`ccd_templates.rs`(フェーズA、`CcdBondTemplate`/`CcdTemplateDb`定義済み)を土台にする。

## 着手前に必ず読むこと

- `RUST_PORT_SPEC.md` §3.9(組み込みDBの設計)・§3.11(本タスクの計画)。
- `rust/crates/proteindf-bridge/src/ccd_templates.rs`(`CcdBondTemplate`/`CcdTemplateDb`の既存定義)。
- `rust/crates/proteindf-bridge/src/format/mmcif.rs`の`SimpleMmcif::get_atomgroup`(CCD形式のパース、特に`_chem_comp_bond.value_order`変換ロジックと`extract_atoms_and_name`)。この結合次数変換ロジックを重複実装せず再利用すること。

## 対象

### 1. `SimpleMmcif`データブロックから`CcdBondTemplate`を組み立てる変換関数

- `CcdBondTemplate::from_mmcif_block(block: &MmcifDataBlock, comp_id: &str) -> Result<CcdBondTemplate>`のようなAPIを新設する(配置場所は`ccd_templates.rs`が自然だが、`format/mmcif.rs`とのモジュール依存関係を見て判断してよい)。
- 内部実装は`format/mmcif.rs`の`_chem_comp_atom`/`_chem_comp_bond`パースロジック(結合次数変換テーブル、`AROM => 1`/`QUAD => 4`含む)を再利用すること。同じ変換ロジックを2箇所に書かないこと。
- CCDのCCD形式ブロック(`_atom_site`を持たない、`has_atom_site() == false`のブロック)のみを対象とする。`_atom_site`ブロックを渡された場合はエラーを返すこと。

### 2. `CcdTemplateDb`への動的登録・マージ手段

- `CcdTemplateDb::insert(&mut self, template: CcdBondTemplate)`(1件登録)を追加する。
- `CcdTemplateDb::merge(&mut self, other: &CcdTemplateDb)`(複数DBの合成)を追加する。同一`comp_id`が両方に存在する場合の優先順位(後勝ち・先勝ち等)を決めてdocコメントに明記すること。
- 組み込みのデフォルトDB(`CcdTemplateDb::global()`)自体は不変のままとする。ユーザーは`CcdTemplateDb::default()`(既存の`Default`実装、`global()`のクローン)を起点に、独自データを`insert`/`merge`で追加登録する使い方を想定する。

### 3. 利用パターンのドキュメント化

- 「ユーザーが`components.cif`(または個別コンポーネントのCCDファイル)をダウンロード → `SimpleMmcif`でロード → 対象コンポーネントを`CcdBondTemplate::from_mmcif_block`で変換 → 自分の`CcdTemplateDb`インスタンスに`insert` → `AtomGroup::apply_ccd_bond_templates`にそのDBを渡す」という一連の流れを、`ccd_templates.rs`のモジュールレベルdocコメントに実例コードとして示すこと。

## 完了の定義(Definition of Done)

1. [x] 既存の`tests/data/ALA.cif`(単一コンポーネントのCCDフィクスチャ)を`SimpleMmcif`でロードし、`CcdBondTemplate::from_mmcif_block`で変換した結果が、組み込みDBの`ALA`エントリ(12結合・C=O二重結合、§3.9フェーズAで実データ検証済み)と一致することを検証する回帰テストを追加すること。(`test_from_mmcif_block_ala_cif_matches_embedded`で確認)
2. [x] 組み込みDBに存在しない架空の合成コンポーネント(CCD形式の合成データ、既存のmmCIFテストで使われている手法と同様のもの)を変換・`insert`し、`apply_ccd_bond_templates`がそれを正しく適用できることを検証する回帰テストを追加すること。(`test_from_mmcif_block_synthetic_custom_ligand`で確認)
3. [x] `_atom_site`ブロックを`from_mmcif_block`に渡した場合に適切にエラーになることを検証すること。(`test_from_mmcif_block_rejects_atom_site`で確認)
4. [x] `merge`の優先順位(同一`comp_id`が重複する場合の挙動)を検証するテストを追加すること。(`test_ccd_template_db_merge_precedence`で確認)
5. [x] `cargo clippy` / `cargo fmt` を通すこと。(ワークスペース全ターゲット警告0・fmtパス確認)
6. [x] `RUST_PORT_SPEC.md` §3.11に実施内容・完了ステータスを追記すること。(完了追記済み)

## スコープ外

- CCD全件の自動ダウンロード・キャッシュ機構(クレートにネットワーク取得を組み込まない。あくまでユーザーが自分でファイルを用意する前提)。
- 組み込みデフォルトDB(29残基)自体の拡張。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない(この機能はPython版に存在しないため)。
- `format/mmcif.rs`の結合次数変換ロジック(`AROM`/`QUAD`対応含む)を重複実装しない。既存のロジックを呼び出す・共有する形にすること。
- クレートにネットワーク取得機能を追加しない(ユーザーが用意したバイト列・ファイルを受け取るだけの設計を維持すること)。
