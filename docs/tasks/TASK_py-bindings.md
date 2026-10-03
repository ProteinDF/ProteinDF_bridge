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

## PR#45 レビュー結果(1回目、2026-10-03、要修正)

`feature/py-bindings-pr45`(`495eeee`・`a8c260b`・`042e1e8`)をレビューした。コアのクレートは変更されていない。`cargo test --workspace`(318 passed)、clippy、fmt、Pythonテスト一式(`unittest discover`、155件)が通ることをClaudeが確認した。基準値はRust側のテストから引用されており、引用元も明記されている。

### 修正依頼

1. **【不具合】`MmcifStructureReport.atomgroup`がアクセスのたびに構造全体の新しいコピーを返す。** Claudeが1WCTで確認したところ、`rep.atomgroup is rep.atomgroup`は`False`で、`rep.atomgroup.setup()`を呼んだ後に`rep.atomgroup`を見ても結合は10本(ファイル由来のみ)のままだった。`setup()`はその場限りのコピーにかかって黙って失われる。大きな構造ではアクセスのたびに全体を複製する性能上の問題もある。レポートが`AtomGroup`をPythonオブジェクト(`Py<PyAtomGroup>`)として1つだけ持ち、毎回同じオブジェクトを返すようにする。`rep.atomgroup is rep.atomgroup`であることと、`rep.atomgroup.setup()`の結果が`rep.atomgroup`に残ることをテストする。
2. **【不具合】`CcdTemplateDb.add_from_file()`が、壊れたデータブロックを黙って無視する。** 1つでも読み込めたブロックがあると、読み込めなかったブロックのエラーを捨てて成功扱いにする。Claudeが、正常な`ALA`と矛盾する重複原子を持つ`BAD`の2ブロックを含むファイルで確認したところ、例外は出ず、戻り値1で`BAD`はDBに入っていなかった。すべてのブロックを先に解析し、1つでも失敗したら、DBを一切変更せずに例外を出す(どのブロックがなぜ失敗したかをメッセージに含める)。このケースをテストする。
3. **【軽微】ファイル読み込みなどRust側の処理に由来するエラーを`PyValueError`で出している箇所がある(`add_from_file`)。** §4.3「共通の設計方針」4に従い、`BrError`系(`to_py_err`)に揃える。引数の値そのものが不正な場合(`secondary_structure`への不正なコードなど)は`ValueError`のままでよい。
4. **【軽微】テストファイル名`tests/test_rs_pr45.py`とクラス名`TestRsPr45...`を、内容に合った名前(例: `tests/test_rs_bond_resolution.py`)に変える。** PR#44と同じ理由(PR番号を名前に使わない)。PR#46以降も同様にする。

### 完了の定義(修正後)

1. 上記1〜4に対応し、同じブランチに追加コミットする。
2. 全PR共通の完了の定義1〜4を満たす。

## PR#45 レビュー結果(2回目、2026-10-03、収束・マージ済み)

修正コミット`9ad0369`を確認した。Claudeが実際に試し、次を確認した。

1. `rep.atomgroup is rep.atomgroup`が`True`になり、`rep.atomgroup.setup()`の結果が残る(1WCTで結合10本→224本)。
2. `ALA`と壊れた`BAD`を含むファイルで`add_from_file()`が`BrInputError`を出し、DBは変更されない。構造データのブロック(`_atom_site`を含むもの)が混ざったファイルも、黙って読み飛ばさずエラーにするようになった。
3. エラーは`BrError`系に揃った。
4. テストファイルは`tests/test_rs_bond_resolution.py`に改名された。

コアのクレートは変更なし。`cargo test --workspace`(318 passed)、clippy、fmt、Pythonテスト一式(157件)をClaudeが確認した。

**残っている軽微な点**: `repr(MmcifStructureReport)`が常に`atoms=0`と表示される(最上位に直接ある原子の数を作成時点で数えているため)。→ PR#46で、全原子数をその時点の値で表示するよう修正する。

ユーザー承認のうえ、2026-10-03にdevelopへマージした(`5211af2`)。**PR#45は完了。** 次はPR#46。

## PR#46 レビュー結果(1回目、2026-10-03、要修正・軽微)

`feature/py-bindings-pr46`(`5a1e80e`・`0d42d03`・`9a9551d`)をレビューした。agyは最終検証の途中で利用上限に達し、完了報告は出なかったため、Claudeが直接検証した。実バグはない。

- Pythonから`1hls.pdb` → `setup()` → `add_missing_hydrogens()` → mmCIF書き出し → 再読み込みが通り、水素の数(元379 → 15追加 → 394)が`test_orchestrator.rs`の基準値と一致する(Rust側に該当のアサートがあることを確認した)。
- `repr(MmcifStructureReport)`が全原子数(1WCTで`atoms=218`)を表示するようになった。
- コアのクレートは変更なし。`cargo test --workspace`(318 passed)、clippy、fmt、Pythonテスト一式(164件)をClaudeが確認した。

### 修正依頼(ユーザー判断により、agyが対応する)

1. **【軽微】`residue_reports`の並び順が実行ごとに変わる。** Rust側の`HashMap`の順序をそのまま使っているため(`docs/rust-port-handoff.md`教訓5「コレクションの順序決定性」)。キー(残基パス)で決定的にソートした順で辞書を作る。2回実行して順序が同じであることをテストする。
2. **【軽微】`skipped_residues`・`step_errors`・`residue_reports`が外から書き換えられる。** Claudeが`rep.skipped_residues.append(...)`を試したところ、レポートの中身が変わった。§4.3「共通の設計方針」3(読み取り専用)に従い、`skipped_residues`・`step_errors`は`(path, reason)`のタプルのタプル(または毎回新しいリスト)、`residue_reports`は毎回新しい辞書(中の`HydrogenationReport`は共有してよい)を返す。外から変更してもレポートが変わらないことをテストする。なお「大きなオブジェクトはコピーせず同じものを返す」というPR#45の方針は`AtomGroup`のような大きな構造が対象で、これらの小さな一覧には当てはまらない。
3. **【軽微】`tests/test_rs_hydrogenation.py`のコメントが、`test_mmcif_writer.rs`の306行目に394という値があると書いているが、実際にはない**(394は`test_orchestrator.rs`の値。`test_mmcif_writer.rs`は再読み込み前後の水素数が等しいことを確認している)。引用元を正しく書き直す。

### 完了の定義(修正後)

1. 上記1〜3に対応し、同じブランチに追加コミットする。
2. 全PR共通の完了の定義1〜4を満たす。

## PR#46 レビュー結果(2回目、2026-10-03、収束・マージ済み)

修正コミット`14eda0b`を確認した。`residue_reports`は残基パスの自然順(`sort_nicely`)で返り、2回実行して順序が同じことをClaudeが確認した。`skipped_residues`・`step_errors`・`residue_reports`はアクセスのたびに新しい一覧を返し、外から変更してもレポートは変わらない。テストのコメントの引用元も修正された。コアのクレートは変更なし。`cargo test --workspace`(318 passed)、clippy、fmt、Pythonテスト一式(165件)をClaudeが確認した。

**PR#47以降で踏襲する方針**: `AtomGroup`のような大きな構造は、アクセスのたびにコピーせず同じオブジェクトを返す。結果の一覧(リスト・辞書)は、外からの変更でレポートが変わらないよう、アクセスのたびに新しいコンテナを返す(中の個々の結果オブジェクトは共有してよい)。

ユーザー承認のうえ、2026-10-03にdevelopへマージした(`7557c84`)。**PR#46は完了。** 次はPR#47。

## PR#47 レビュー結果(1回目、2026-10-03、収束・マージ済み)

`feature/py-bindings-pr47`(`cb9f1eb`・`d42394c`)をレビューした。指摘なし。

- 基準値の引用元(`test_hydrogen_bond.rs`・`test_interaction_set.rs`・`test_ch_pi.rs`・`test_secondary_structure.rs`・`test_sidechain_hydrogen_bond.rs`)の該当行に実際にその値があることを、Claudeが抜き取りで確認した。
- 結果の一覧はアクセスのたびに新しいコンテナを返し、主鎖の水素結合の順序が決定的であることもテストされている。`apply_secondary_structure`は呼び出したオブジェクト自身を変更し、部分木のコピーでは元の木が変わらない。
- 1hlsで`InteractionSet.detect_all`が`total=38(disulfide 3, salt_bridge 0, hydrogen_bond 28, ch_pi 7)`になること、CH-πの閾値を省略した場合と既定値(4.5Å・40°)を明示した場合の結果が一致することをClaudeが確認した。
- コアのクレートは変更なし。`cargo test --workspace`(318 passed)、clippy、fmt、Pythonテスト一式(180件)をClaudeが確認した。

**参考(コア側の既知の限界、対応不要)**: CH-πの閾値に負の距離などを渡してもエラーにならず、結果が空になる(コア側に閾値の検証がない)。

ユーザー承認のうえ、2026-10-03にdevelopへマージした(`5e0fb4c`)。**PR#47は完了。** 次はPR#48。

## PR#48 レビュー結果(1回目、2026-10-03、収束・マージ済み)

`feature/py-bindings-pr48`(`63c101a`・`9733c9c`・`790def4`)をレビューした。指摘なし。`get_ACE`・`get_NME`は原子数と座標(小数4桁)を、`Neutralize`は1hls(782→792原子)の原子名と座標を純Python版と照合している。`.brd`は純Python版との相互読み書きを確認している。意図した違い(`Neutralize`の`_exempt_list`、2026.9.3で修正済みの`arbitary_rotate_matrix`)は根拠とともにテストに書かれている。`Neutralize.neutralized`は同じオブジェクトを返し、存在しないファイルの読み込みは`BrError`になる。`modeling.rs`の`#![allow(non_snake_case)]`は純Python版のメソッド名(`get_ACE`など)に合わせるためのもので妥当。コアのクレートは変更なし。`cargo test --workspace`(318 passed)、clippy、fmt、Pythonテスト一式(196件)をClaudeが確認した。ユーザー承認のうえ、2026-10-03にdevelopへマージした(`6df79ae`)。

**PR#45〜48はすべて完了。** `RUST_PORT_SPEC.md` §4.3に実施内容・既知の限界を、`docs/rust-port-handoff.md`に記録した。
