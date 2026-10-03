# TASK: pyo3 0.29の`FromPyObject`非推奨警告の解消

> developでは`cargo clippy --workspace --all-targets -- -D warnings`が、`proteindf-bridge-py`でのpyo3 0.29の非推奨警告27件により失敗している(`Clone`を実装する`#[pyclass]`型の`FromPyObject`自動実装がopt-inに変わる件)。PR#41のレビューで判明し、ユーザー判断(2026-10-03)で、クレート全体の`#![allow(deprecated)]`で黙らせるのではなく、正式な方法で直すことにした。mmCIF書き出しのPythonバインディング(PR#44、`TASK_mmcif-writer.md`)より前に片付ける。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaude(`/code-review`)が担当する。
- **PR#41(`TASK_hydrogenation-cleanup.md`)がdevelopにマージされてから着手する。** PR#42(mmCIF書き出し)とは独立しているので、どちらを先にしてもよい。
- `develop`から`fix/pyo3-from-py-object`(名前はagyの判断でよい)を切って作業する。
- **`develop`へは自分でマージしない。** 完了したら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認を待つこと。

## 対象

1. 警告が出ている各`#[pyclass]`について、Pythonから引数として値で受け取る(`FromPyObject`が必要な)型には`#[pyclass(from_py_object)]`を、不要な型には`#[pyclass(skip_from_py_object)]`を付ける。**どちらにするかは型ごとに、Rust側の関数シグネチャでその型を値(`T`)や`PyRef<T>`以外の形で受け取っている箇所があるかを調べて決める。** 迷ったら、従来の挙動を保つ`from_py_object`にする。
2. `#![allow(deprecated)]`や`#[allow(deprecated)]`で警告を黙らせないこと。

## 完了の定義

1. `cargo clippy --workspace --all-targets -- -D warnings`が通る。
2. `cargo test --workspace`と`cargo fmt --all -- --check`が通る。
3. Pythonバインディングの既存テスト(Python版との比較テスト等、Phase 5で作られたもの)が通る。実行方法(`maturin develop`等)と結果を完了報告に書く。
4. `skip_from_py_object`にした型の一覧と、そう判断した理由を完了報告に書く。
