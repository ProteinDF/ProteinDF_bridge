# TASK: CCD結合テンプレートデータベース フェーズA(結合次数補完)

> 本タスクは `RUST_PORT_SPEC.md` §3.9(2026-09-22追記)の計画のフェーズAを実装するものです。背景・設計方針・スコープの詳細は同節を参照してください。フェーズB(水素付加)は本タスクの対象外です。

## ブランチ運用(MUST)

- `develop` から `feature/ccd-bond-templates` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。
- 依存関係: なし。`format/mmcif.rs`のCCDパース(2026-09-22、`AROM`/`QUAD`結合次数対応済み)を土台にする。

## 着手前に必ず読むこと

- `RUST_PORT_SPEC.md` §3.8(ファイル由来結合とVDWヒューリスティックの優先順位方針)・§3.9(本タスクの計画)。
- `rust/crates/proteindf-bridge/src/format/mmcif.rs`の`SimpleMmcif::get_atomgroup`(CCD形式のパース、`_chem_comp_bond.value_order`の変換ロジック含む)。
- `rust/crates/proteindf-bridge/src/brd.rs`(既存のMessagePack往復フォーマットインフラ、`rmp-serde`の使い方)。

## 対象

### 1. CCDサブセットデータの抽出

- **対象コンポーネント**: 標準アミノ酸20種(ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL)、標準核酸8種(DA, DC, DG, DT, A, C, G, U)、水(HOH)、および既存テストフィクスチャ(`1hls.pdb`/`2MGO.pdb`/`3i3zH.pdb`等)に登場する非標準残基があれば追加で含める。
- **データソース**: wwPDBは各コンポーネントを個別のCCD形式ファイルとして配布している(例: `https://files.rcsb.org/ligands/view/ALA.cif`)。ネットワークアクセスが利用可能であれば、これらを実際に取得して使うこと。**取得したCCDファイルの内容(原子名・結合トポロジー・結合次数)を手作業で書き起こしたり、記憶や推測で捏造したりしないこと。** 化学構造データは正確性が最優先であり、誤った結合次数を埋め込むと下流の全ての利用箇所に静かに伝播する。ネットワークアクセスができない、または実データの入手方法が分からない場合は、実装を止めてユーザー経由でClaudeに相談すること。
- **抽出スクリプト**: `scripts/`配下に1回限りの抽出スクリプト(Python、既存`scripts/`ディレクトリの慣習に合わせる。例: `scripts/build_ccd_bond_templates.py`)を作成し、取得したCCDファイル群から`_chem_comp.id`・`_chem_comp_atom.atom_id`・`_chem_comp_bond.{atom_id_1,atom_id_2,value_order}`を抽出してMessagePack形式にシリアライズする。このスクリプトは1回実行して出力データファイルを生成するためのものであり、実行時には呼ばれない(生成物だけをリポジトリにコミットする)。

### 2. テンプレートデータの同梱

- 生成したMessagePackファイルを`rust/crates/proteindf-bridge/src/data/ccd_bond_templates.msgpack`のような場所に配置し、`include_bytes!`でバイナリに埋め込むこと(§3.9で確定した方針。外部ファイルを実行時にパスで探す方式は採用しないこと — wasm32ターゲット・PyO3 wheel配布のポータビリティを損なうため)。
- 初回参照時に`std::sync::OnceLock`等で遅延デシリアライズし、静的なルックアップテーブル(`HashMap<String, CcdBondTemplate>`等)として保持すること。

### 3. テンプレート構造体とAPI設計

- 新規モジュール`rust/crates/proteindf-bridge/src/ccd_templates.rs`を作成する。
- `CcdBondTemplate { comp_id: String, atoms: Vec<String>, bonds: Vec<(String, String, usize)> }`(原子名ペア+結合次数)のような構造体を定義する。
- `CcdTemplateDb::lookup(comp_id: &str) -> Option<&CcdBondTemplate>`のようなルックアップAPIを提供する。

### 4. `AtomGroup`への適用

- `AtomGroup::apply_ccd_bond_templates(&mut self, db: &CcdTemplateDb)`のようなメソッド(配置場所は`atom_group.rs`または`ccd_templates.rs`、設計時に判断してよい)を追加する。各residueレベルのグループについて、その`name`でテンプレートDBを引き、原子名の対応が取れる結合ペアについて結合次数を設定する。
- **§3.8の優先順位に第3階層として追加すること**: (1) ファイル由来の明示的結合(既存の`ag.get_bond_list().is_empty()`判定を通過済みのもの) > (2) CCDテンプレートによる結合(本タスク) > (3) `Bond::setup()`のVDWヒューリスティック。
- **既存の結合を上書きしないこと**(テンプレートは「無い結合の補完」であり、ファイル由来の結合が既にある場合はそれを尊重する)。同一原子ペアに矛盾する結合次数がある場合の扱い(ファイル由来を優先する等)は、実装前に方針をユーザー経由でClaudeに確認すること。

## 完了の定義(Definition of Done)

1. [x] 標準アミノ酸(例: グルタミン酸のCOOH側鎖、アルギニンのグアニジノ基等、二重結合や複数の結合次数パターンを含む残基)を含む実PDBフィクスチャ(`1hls.pdb`等)で、`Bond::setup()`単独では次数1にしかならない結合が、`apply_ccd_bond_templates`適用後に正しい次数になることを検証する回帰テストを追加すること。(`test_apply_ccd_bond_templates_1hls_real_pdb`でパス確認)
2. [x] ファイル由来の明示的結合(CONECT等)がある場合、テンプレート適用がそれを上書きしないことを検証するテストを追加すること。(`test_apply_ccd_bond_templates_does_not_overwrite_existing_bonds`でパス確認)
3. [x] テンプレートDBに存在しない残基名(非標準・見つからない)の場合、エラーにならず何もせず処理が続行されることを検証すること。(`test_apply_ccd_bond_templates_unknown_component_safe_skip`でパス確認)
4. [x] `cargo clippy` / `cargo fmt` を通すこと。(全ターゲット `-- -D warnings` パス確認)
5. [x] `RUST_PORT_SPEC.md` §3.9のフェーズA該当箇所に実施内容・完了ステータスを追記すること。(完了追記済み)

## スコープ外

- フェーズB(水素付加)。
- CCD全件(`components.cif`全体)のデータベース化。
- プロトネーション状態・互変異性体の推定。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない(この機能はPython版に存在しないため)。
- **化学構造データ(原子名・結合トポロジー・結合次数)を手作業で捏造・推測しない。** 実際のCCDデータソースから取得すること。入手方法に迷ったら実装を止めてClaudeに相談すること。
- 実行時に外部データファイルをパスで探す方式を採用しない(`include_bytes!`によるバイナリ埋め込み方式で確定済み)。
- ファイル由来の既存の結合情報をテンプレートで上書きしない。
