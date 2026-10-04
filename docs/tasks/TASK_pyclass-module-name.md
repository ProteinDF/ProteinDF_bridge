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

## レビュー結果(1回目、2026-10-05、収束・マージ済み)

`fix/pyclass-module-name`(`66d9bf7`・`58f504e`)をレビューした。指摘なし。Claudeがこのブランチを新しい環境にインストールして確認したところ、Rust版のすべてのクラス(例外型を含む)の`__module__`が`proteindf_bridge.rs`になった(`<class 'proteindf_bridge.rs.AtomGroup'>`)。すべての型を一括で確認するテスト`tests/test_rs_module_name.py`が追加されている。古い名前はRustのソース・テスト・API文書に残っておらず、残っているのは`RUST_PORT_SPEC.md`・`docs/rust-port-handoff.md`の過去の経緯の記述のみ。そのうち、現在の読者を誤解させる§4.2の利用指針とhandoffの命名規約には、Claudeが統合後の名前への読み替えの注記を加えた。Pythonテスト一式(198件)、`cargo test --workspace`(318件)、clippy、fmtをClaudeが確認した。ユーザー承認のうえ、2026-10-05にdevelopへマージした(`d1b7388`)。**本タスクは完了。**
