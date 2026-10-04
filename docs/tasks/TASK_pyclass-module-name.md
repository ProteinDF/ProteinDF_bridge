# TASK: Rust版クラスの`__module__`を`proteindf_bridge.rs`に直す

> PR#49(`TASK_packaging.md`)でRust版の拡張モジュールを`proteindf_bridge.rs`としてimportする形にしたが、`rust/crates/proteindf-bridge-py/src`の各`#[pyclass(module = "proteindf_bridge_rs")]`(47か所)が古いままである。PR#51のレビュー(2026-10-04)で判明し、ユーザー承認(2026-10-05)のうえ別タスクとした。

## 影響

importや動作には影響しないが、`rs.AtomGroup.__module__`が存在しないモジュール名`proteindf_bridge_rs`を返すため、`help()`の表示、Sphinxの自動APIドキュメント、`repr`やエラーメッセージ中の型名に古い名前が出る。将来`pickle`に対応する場合は正しいモジュール名が必須になる。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaude(`/code-review`)が担当する。
- `develop`から専用のブランチを切って作業する。**`develop`へは自分でマージしない。**

## 対象

1. `#[pyclass(module = ...)]`をすべて`proteindf_bridge.rs`に直す。`#[pyfunction]`や例外型(`create_exception!`など)にもモジュール名を指定している箇所があれば、同様に直す。
2. 古い名前`proteindf_bridge_rs`が、Rustのソース・Pythonのテスト・文書(`docs/`の`api`を含む)に残っていないか確認し、残っていれば直す。過去の経緯として記録している文書(`RUST_PORT_SPEC.md`の過去の節、`docs/tasks/`の過去のレビュー記録)は書き換えない。

## 完了の定義

1. Rust版のすべてのクラス・例外型の`__module__`が`proteindf_bridge.rs`になることを、Pythonのテストで確認する(モジュール内の型を列挙して一括で確認する)。
2. `python -m unittest discover -s tests`、`cargo test --workspace`、`cargo clippy --workspace --all-targets -- -D warnings`、`cargo fmt --all -- --check`が通る。
