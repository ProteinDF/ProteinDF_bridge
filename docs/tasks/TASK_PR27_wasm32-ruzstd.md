# TASK_PR27: wasm32ターゲットの`ruzstd`自動選択

> 本タスクは `docs/rust-port-handoff.md` の「Phase 10」節(RUST_PORT_SPEC.md §9対応)から抽出したものです。全体の背景・優先順位・他PRとの関係は同ドキュメントを参照してください。

## ブランチ運用(MUST)

- `develop` から `feature/phase10-pr27` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。承認が出るまで次のタスクの実装に着手しない。
- 依存関係: なし。Phase 10の中で最も変更が軽いタスクなので最初に着手してよい。

## 背景

`RUST_PORT_SPEC.md` §9(「結」リポジトリのフェーズ6e調査で見つかった統合ブロッカー)の「中」優先度項目。現状 `rust/crates/proteindf-bridge/Cargo.toml` は `default = ["zstd"]` のみで、wasm32ビルド時は消費側(YUI)が毎回 `--no-default-features --features ruzstd` を指定する必要がある。これをCargo.toml側の設定だけで解決してほしい、というのがYUI側の要望。

## 対象

`rust/crates/proteindf-bridge/Cargo.toml` に `[target.'cfg(target_arch = "wasm32")'.dependencies]` セクションを追加し、wasm32ターゲットでは `zstd` ではなく `ruzstd` が自動的に有効になるようにする。YUI自身の `core/Cargo.toml` が同じパターンを既に採用しているので、参照可能であれば実際のコードを確認して揃えること。

## 完了の定義(Definition of Done)

1. `cargo check --target wasm32-unknown-unknown` が追加のフラグなしで成功すること。
2. 通常のネイティブビルド(`cargo build --workspace`)は引き続き既定で `zstd`(Cバインディング版)を使うこと(既定の挙動を壊さない)。
3. `cargo test --workspace` が既存同様パスすること。
4. `cargo clippy` / `cargo fmt` を通すこと。

## スコープ外

- `brd.rs` のAPI変更は不要(どのfeatureが選ばれるかが変わるだけで、公開APIは変わらない)。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- Phase 10の他タスク(schema検証・ファイル由来結合・`Bond::setup()`のスケーラビリティ等)には手を出さない。
