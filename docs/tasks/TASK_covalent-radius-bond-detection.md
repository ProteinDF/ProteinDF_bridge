# TASK: `Bond::setup()`のVDWヒューリスティックを共有結合半径ベースに変更

> 本タスクは `RUST_PORT_SPEC.md` §3.10(2026-09-22追記)の計画を実装するものです。背景・設計方針・スコープの詳細は同節を参照してください。

## ブランチ運用(MUST)

- `develop` から `feature/covalent-radius-bond-detection` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。
- 依存関係: なし。§3.9のCCDテンプレートDB(フェーズA、2026-09-22完了・マージ済み)とは独立(結合次数の判定はCCDテンプレート、結合の有無の判定は本タスクの担当で、役割が異なる)。

## 着手前に必ず読むこと

- `RUST_PORT_SPEC.md` §3.10(本タスクの計画、背景・数値根拠を含む)。
- `rust/crates/proteindf-bridge/src/bond.rs`の`Bond::setup()`(現行のVDW半径ヒューリスティック、PR#35で導入した`CellList`によるO(N)近傍探索の実装含む)。
- `rust/crates/proteindf-bridge/src/periodic_table.rs`の既存`VDW`テーブル(共有結合半径テーブルはこれと並行する形で追加する。削除・置き換えはしないこと)。

## 対象

### 1. 共有結合半径テーブルの追加

- `periodic_table.rs`に`COVALENT_RADIUS`(または同等の名前)テーブルを新設する。
- **値は信頼できる文献から取得すること。** 推奨: Cordero et al., "Covalent radii revisited", *Dalton Trans.*, 2008, 2832–2838(pymatgen・ASE等の主要ツールが採用している現代的な標準参照値)。値の出典をコード中のdocコメントに明記すること。
- **値を記憶や推測で埋めないこと。** 既存の`VDW`テーブルと同様、実際の文献値を確認して使うこと。出典に自信が持てない元素がある場合は、実装を止めてユーザー経由でClaudeに相談すること。
- `PeriodicTable::covalent_radius(atom: impl IntoAtomId) -> Result<f64>`のような、既存の`vdw()`と対になるAPIを追加する。

### 2. 許容値(トレランス)の決定

- 文献的に一般的な加算マージン(例: 0.4Å)または乗算係数(例: 1.2倍)のいずれかを採用し、名前付き定数として`bond.rs`に定義する(現行の`MAX_DENSE_MATRIX_ATOMS`と同様の扱い)。
- どちらの方式・どの値を採るかは実装者の判断でよいが、根拠(参考にした文献・ツールの慣習、例えばOpenBabelは加算0.45Å、ASEの`natural_cutoffs`は乗算1.2倍を既定にしている等)をdocコメントに残すこと。

### 3. `Bond::setup()`の判定式変更

- `Bond::setup()`内の判定式`distance(p, q) <= vdw(p) + vdw(q) + 0.4`を、`distance(p, q) <= covalent_radius(p) + covalent_radius(q) + <新しい許容値>`に置き換える。
- **VDW半径テーブル自体・`PeriodicTable::vdw()`は削除しないこと**(将来的に非結合接触判定等、別用途で使われる可能性があるため残す)。
- PR#35で実装した`CellList`によるO(N)近傍探索の仕組み(動的セルサイズの計算含む)はそのまま踏襲し、セルサイズの計算に使う半径を`vdw`から`covalent_radius`に差し替えること。
- **Python版`bond.py`からの意図的な乖離であることをdocコメント・PR説明に明記すること**(Phase 3 PR#11で`superposer_quaternion.py`の実バグを模倣しなかった前例と同種の判断)。

## 完了の定義(Definition of Done)

1. [x] 実PDBフィクスチャ(`1hls.pdb`等)で、既知の共有結合(主鎖のペプチド結合、既知のジスルフィド結合等)が変更後も正しく検出され続けることを確認する回帰テストを追加すること。(`test_covalent_bond_detection::test_known_covalent_bonds_1hls`で確認)
2. [x] 変更前のカットオフでは「結合あり」と誤検出されるが実際には共有結合ではない原子ペアを合成データで構築し(例: 新カットオフの外・旧カットオフの内に収まる距離に2原子を配置)、変更後は「結合なし」と正しく判定されることを検証する回帰テストを追加すること。(`test_covalent_bond_detection::test_non_covalent_contact_excluded`で確認)
3. [x] 既存の`test_bond_setup`系テスト(想定結合数・結合距離を検証しているもの)を、新しい判定式に合わせて確認・必要なら更新すること。想定結合数が変わる場合は、その理由をテストのコメントに明記すること。(`bond.rs`内のテスト更新とコメント追記完了)
4. [x] PR#35のベンチマークテスト(`test_benchmark_1m_atoms`、`#[ignore]`)が新しい判定式でも問題なく動作することを確認すること(結果の結合数が変わるのは想定内、クラッシュ・著しい性能劣化がないことを確認)。(6.36秒で正常動作確認)
5. [x] `cargo clippy` / `cargo fmt` を通すこと。(ワークスペース全ターゲット警告0・フォーマット正常確認)
6. [x] `RUST_PORT_SPEC.md` §3.10に実施内容・完了ステータス・採用したトレランス値とその根拠を追記すること。(完了追記済み)

## スコープ外

- 結合次数の判定(引き続き§3.9のCCDテンプレート、テンプレートが無い場合は次数1のまま)。
- VDW半径テーブル自体の削除・置き換え。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- **共有結合半径の値を記憶・推測で捏造しない。** 実際の文献値を確認して使うこと。
- VDW半径テーブル・`PeriodicTable::vdw()`を削除・変更しない(共有結合半径テーブルを追加するのみ)。
