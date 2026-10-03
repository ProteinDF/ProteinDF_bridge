# TASK: Pythonバインディングの拡充 (PR#45〜48)

> 2026.10.0時点の`proteindf_bridge_rs`は、Phase 6以降の機能(`.brd`・`Modeling`・`Neutralize`、水素結合・DSSP・CH-π・`InteractionSet`、結合解決`AtomGroup::setup`・CCDテンプレートDB、水素付加)をPythonに公開していない。ユーザー確認(2026-10-03)のうえ、A(基盤・結合解決)→B(水素付加)→C(解析)→D(既存Python機能の移行)の順に公開する。設計方針は`RUST_PORT_SPEC.md` §4.3にある。**着手前に§4.3を必ず精読すること。** 本ファイルは実行チェックリストである。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaude(`/code-review`)が担当する。
- 着手順はPR#45 → PR#46 → PR#47 → PR#48。各PRは前のPRがレビューで承認され、developにマージされてから始める。
- `develop`から`feature/py-bindings-prNN`(名前はagyの判断でよい)を切って作業する。
- **`develop`へは自分でマージしない。** 完了したら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認を待つこと。
- **完了報告では、コマンドの出力(`test result`行、`Ran N tests`、`OK`など)を要約・再構成せず、実際の出力行をそのまま貼ること。** 過去に、報告されたテスト件数やテストファイル一覧が実際と食い違ったことがある(`docs/rust-port-handoff.md`の教訓)。

## 全PR共通の完了の定義

1. `cargo test --workspace`、`cargo clippy --workspace --all-targets -- -D warnings`、`cargo fmt --all -- --check`が通る。
2. 既存のPythonテスト(`python -m unittest discover -s tests`、純Python版のテストと`tests/test_rs_*.py`のすべて)と、そのPRで追加したテストが通る。実行方法は§4.3「共通の設計方針」5。
3. 新しい`#[pyclass]`には`from_py_object`/`skip_from_py_object`が明示されている。
4. コアクレート(`proteindf-bridge`)を変更していない。変更した場合は、その内容と理由を完了報告に書く。

## PR#45: A. 基盤・結合解決

### 対象

1. `AtomGroup.setup()`・`AtomGroup.setup_with_db(db)`(§3.13〜3.15の結合解決の標準の入口)。既存の`Bond.setup()`(ヒューリスティックのみ)はそのまま残す。
2. `CcdTemplateDb`クラス。
   - 組み込みの29テンプレートを持つDBを得る方法(例: `CcdTemplateDb.builtin()`)。
   - ユーザー提供のCCDファイル(mmCIF、複数のデータブロックを含みうる)を読み込んで追加する方法(§3.11)。
   - `lookup(comp_id)`で、テンプレートの原子(名前・元素・理想座標)と結合(原子の組・次数)を読み取り専用で参照できること。`len()`・`merge()`など、Rust側の公開APIに対応するもの。
3. 階層規約の検証: `AtomGroup.validate_schema()`(違反の一覧を返す。各違反の種類・パス・説明が取れること)、`is_model_level()`・`is_chain_level()`・`is_residue_level()`。
4. 二次構造フィールドの読み書き: 残基レベルの`AtomGroup`の`secondary_structure`(`"H"`/`"E"`/`"-"`または`None`)。
5. mmCIF: `SimpleMmcif.get_structure_atomgroup_with_report()`(および`_for_block_with_report`)。解決できなかった`_struct_conn`の一覧(結合ID・種類・どちらの相手のどこが見つからなかったか)が取れること。

### 完了の定義

1. 1HLS・2FB4・1WCTなど既存の実データで、`setup()`後の結合数や結合の組が、Rust側のテスト(`test_bond_resolution.rs`・`test_ccd_templates.rs`・`test_mmcif_writer.rs`のPR#43部分など)と同じ基準値になることを、Pythonから確認する。基準値をどのRustテストから引用したかをコメントに書く。
2. ユーザー提供CCDファイル(既存フィクスチャ`ALA.cif`など)を読み込んで`setup_with_db()`に使えることを確認する。
3. 部分木のコピーに`setup()`を呼んでも元の木が変わらないこと(§4.3「共通の設計方針」2)をテストで確認し、docstringに書く。
4. `validate_schema()`・`with_report`の結果を、違反あり/なしの両方のケースで確認する。

## PR#46: B. 水素付加

### 対象

1. `AtomGroup.add_missing_hydrogens(db=None)`(`db`を省略したら組み込みDB)。呼び出したオブジェクト自身に水素を追加する。
2. 結果レポート`OverallHydrogenationReport`(`total_added_hydrogens`・`total_removed_hydrogens`・`hydrogenated_residues`・`skipped_residues`・`step_errors`・`residue_reports`)と、残基ごとの`HydrogenationReport`を読み取り専用で公開する。

### 完了の定義

1. `1hls.pdb`を読み込み → `setup()` → `add_missing_hydrogens()` → `SimpleMmcif.save_structure()` → 再読み込み、という一連の流れがPythonから通り、追加された水素の数などがRust側のテスト(`test_orchestrator.rs`・`test_mmcif_writer.rs`の`test_hydrogenation_pipeline_roundtrip`)と同じ基準値になることを確認する。
2. 結晶水やテンプレートのない残基を含む構造で、`skipped_residues`・`step_errors`がPythonから確認できることをテストする。

## PR#47: C. 解析(Phase 8・9)

### 対象

1. 主鎖の水素結合(`calc_backbone_hbonds`)と、側鎖の水素結合(`calc_sidechain_hbonds`・`_with_options`)。
2. DSSP: `calc_secondary_structure(chain)`(残基ごとの結果を返す)と`apply_secondary_structure(chain)`(`AtomGroup`に書き戻す。PR#45の`secondary_structure`フィールドで読める)。
3. CH-π相互作用の検出。距離・角度の閾値は引数で変えられること(§3.3「決め打ちにしない」)。
4. `InteractionSet`の構築と、MessagePack・YAMLでの書き出し・読み込み。

### 完了の定義

1. `1hls.pdb`などで、各解析の結果(件数・代表的な組・エネルギーや角度など)がRust側のテスト(`test_hydrogen_bond.rs`・`test_sidechain_hydrogen_bond.rs`・`test_secondary_structure.rs`・`test_ch_pi.rs`・`test_interaction_set.rs`)と同じ基準値になることを、Pythonから確認する。
2. `InteractionSet`をPythonで書き出して読み直し、内容が一致することを確認する。

## PR#48: D. 既存Python機能の移行(Phase 6)

### 対象

1. `.brd`の読み書き: 純Python版`functions.py`の`load_atomgroup`・`save_atomgroup`に相当するもの(プレーンなMessagePack)と、YUIヘッダー形式(Magic + Version + zstd、§1・Phase 6)の読み書き。
2. `Modeling`: 純Python版`modeling.py`の公開メソッド(`get_ACE`・`get_NME`・`get_ACE_simple`・`get_NME_simple`・`add_methyl`・`get_NH3`・`select_residues`・`get_last_index`・`neutralize_*`など)。
3. `Neutralize`: 純Python版`neutralize.py`の`Neutralize`クラス。

### 完了の定義

1. 既存のphase1〜7のテストと同じ方式で、純Python版と結果(原子数・原子名・座標など)を1対1で比較するテストを追加する。Rust版が意図して純Python版と異なる点(`docs/rust-port-handoff.md`の教訓、例えば`neutralize.py`の`_exempt_list`)がある場合は、その理由をテストのコメントに書く。
2. 純Python版で書いた`.brd`をRust版で読めること、およびその逆を確認する。

## 全PR完了後

- `RUST_PORT_SPEC.md` §4.3に「実施内容・検証」と既知の限界を追記する。
- `docs/rust-port-handoff.md`に記録する。
