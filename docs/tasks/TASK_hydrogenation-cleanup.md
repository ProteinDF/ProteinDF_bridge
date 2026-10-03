# TASK: 水素付加オーケストレーターのエラー記録の一元化と後片付け (PR#41)

> §3.16(PR#37〜40)の水素付加はdevelopにマージ済み。PR#40のレビュー2〜5回目で、エラー記録が2系統に分かれていることが「片方を直すともう片方が壊れる」系の不具合を繰り返し生んだと分かった(`RUST_PORT_SPEC.md` §3.16 PR#40「設計上の技術的負債」、`docs/tasks/TASK_ccd-hydrogenation.md`のPR#40レビュー5回目)。mmCIF書き出し(§3.17、`TASK_mmcif-writer.md`)に着手する前に、PR#40の経緯を覚えているうちに整理する。ユーザー確認済み(2026-10-03)。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaude(`/code-review`)が担当する。
- `develop`から`feature/hydrogenation-cleanup-pr41`(名前はagyの判断でよい)を切って作業する。
- **`develop`へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`の確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認を待つこと。
- **共有ワークツリー注意(MUST)**: 作業開始前に`git status`と`git branch --show-current`を確認する(`docs/rust-port-handoff.md` §1.4)。

## 対象

### 1. 挙動を固定するテストを先に追加する(リファクタリング前に必須)

`orchestrator.rs`の`hydrogenate_single_residue`について、次の組み合わせごとに「`residue_reports`・`skipped_residues`・`step_errors`のどれに、何回記録されるか」を**現在の挙動のまま**固定するテストを`test_orchestrator.rs`に追加する。

- Step1(主鎖): 対象外(アミノ酸でない) / 成功 / 失敗
- Step2(CCD): 成功 / 失敗(テンプレートはあるが処理に失敗) / テンプレート未検出
- 変更の有無: あり / なし

到達できない組み合わせはテストのコメントにその旨を書けばよい。合成データで各組み合わせを作るのが難しい場合は、`test_partially_modified_unknown_residue_recorded_in_reports_only`などの既存テストの作り方を参考にする。

**注意**: 現在の実装では、Step1とStep2がどちらも`Err`で変更がない場合、同じ残基が`step_errors`と`skipped_residues`の両方に載る。これも「現在の挙動」としてそのまま固定する。この挙動がおかしいと判断した場合は、変更せずにユーザー経由でClaudeに報告すること(このPRは挙動を変えないリファクタリングである)。

### 2. エラー記録を関数末尾の1か所にまとめる

- Step1・Step2の`Err`分岐での`report.record_error(...)`の即時呼び出しをやめ、エラー情報を変数に保持するだけにする。
- `step2_err: Option<String>`は「テンプレート未検出(想定内)」と「処理失敗(想定外)」という意味の異なる2つのケースを1つの変数で表しており、PR#40で不具合が繰り返した根本原因になっている。小さな列挙型(例: `enum StepOutcome { NotApplicable, Ok, Failed(String), TemplateMissing(String) }`。名前・形はagyの判断でよい)で区別する。
- 関数末尾の判定ブロック(`has_modifications` × エラーの種類)だけで`residue_reports`・`skipped_residues`・`step_errors`への記録先を決める。
- これにより`step_errors`の線形スキャンによる重複防止(`!report.step_errors.iter().any(...)`、PR#40レビュー5回目の指摘1、O(n²))は不要になるので削除する。
- `OverallHydrogenationReport`の公開フィールドとdocコメントの意味は変えない。

### 3. ドキュメントの更新

- `docs/rust-port-handoff.md` §5の「CCD参照構造による水素付加」の項目が「計画済み・未着手」のまま残っている。完了済み(PR#37〜40、2026-09-27にdevelopにマージ)に更新し、§3のフェーズ実績にも水素付加(PR#37〜40)を1項目として追加する。
- `RUST_PORT_SPEC.md` §3.16 PR#40の「設計上の技術的負債」の末尾にある「PR#41で対応予定」を、完了の記録に書き換える。

## 完了の定義

1. 対象1で追加したテストが、**リファクタリングの前と後の両方で**同じ結果になる(先にテストだけをコミットし、次にリファクタリングをコミットする2段階にすると、レビュー時に挙動が変わっていないことを確認しやすい)。
2. 既存の`test_orchestrator.rs`の9テストを含め、`cargo test --workspace`がすべて成功する。
3. `cargo clippy --workspace --all-targets -- -D warnings`と`cargo fmt --all -- --check`で警告・エラーが出ない。
4. 完了したら、ブランチ名と完了内容をユーザー経由でClaudeに報告する。

## PR#41 レビュー結果(1回目、2026-10-03、要修正・軽微)

`feature/hydrogenation-cleanup-pr41`(`bc8e189`・`37fc18e`・`c638937`)をレビューした。**リファクタリング本体に実バグはない。**

- 追加された組み合わせテスト(`test_orchestrator.rs`、計24件)を、リファクタリング前のコミット`bc8e189`で実行しても全件成功することをClaudeが実際に確認した。コードを読んで、`step_errors`への記録順(Step1 → Step2)、記録先、`skipped_residues`の理由の選び方が元と同じであることも確認した。削除した重複防止ガードが効いていたのは到達不能な組み合わせだけである。
- `cargo test --workspace`・`cargo fmt --check`の成功もClaudeが確認した。

### 修正依頼

1. **`proteindf-bridge-py/src/lib.rs`の`#![allow(deprecated)]`を削除すること(スコープ外の変更)。** developでは元々、pyo3 0.29の非推奨警告(`FromPyObject`のopt-in化、27件)で`cargo clippy --workspace --all-targets -- -D warnings`が失敗している。これはClaudeがTASKの完了条件に`--workspace`を書いたことによる見落としである。クレート全体の非推奨警告をまとめて黙らせると、今後の非推奨にも気づけなくなる。pyo3の警告は、別タスク`docs/tasks/TASK_pyo3-from-py-object.md`で、`#[pyclass(from_py_object)]`などを付ける正式な方法で直す(ユーザー判断、2026-10-03)。
2. **`test_unreachable_step_outcome_combinations`の`#[test]`を外すこと。** 中身がコメントだけで何も検証していない。到達不能な4パターンの説明は、テストモジュール内のコメントとして残す。

### 完了の定義(修正後)

1. 上記2点を同じブランチに追加コミットする。
2. `cargo test --workspace`と`cargo fmt --all -- --check`が通ること。clippyは**PR#41に限り** `cargo clippy -p proteindf-bridge --all-targets -- -D warnings`(コアのクレートのみ)で確認する。`--workspace`での失敗はdevelopに元からあるもので、上記の別タスクで解消する。

## PR#41 レビュー結果(2回目、2026-10-03、収束・マージ済み)

修正コミット`8ccd53b`を確認した。`proteindf-bridge-py`への変更はなくなり、空のテストは説明コメントに置き換わった。`cargo test --workspace`(全件成功)、`cargo clippy -p proteindf-bridge --all-targets -- -D warnings`、`cargo fmt --check`をClaudeが確認した。ユーザー承認のうえ、2026-10-03にdevelopへマージした(`77af040`)。**PR#41は完了。** 次はpyo3の警告対応(`TASK_pyo3-from-py-object.md`)とmmCIF書き出し(`TASK_mmcif-writer.md`のPR#42)。
