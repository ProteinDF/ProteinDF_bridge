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

## レビュー指摘 (2026-09-19, code-review)

実装(コミット `ad55de7`, `aeab5bc`)をレビューした結果、正当性の問題は見つかりませんでした(ネイティブ/wasm32双方・`zstd`/`ruzstd`双方のfeature組み合わせでビルド・`cargo clippy -D warnings`・`cargo fmt --check`・既存テスト・クロスバックエンド往復テストが全てパス)。

以下、Minorな指摘が1件あります。修正必須ではありませんが、対応する場合はご確認ください。

- **`rust/crates/proteindf-bridge/Cargo.toml:21` と `:24`** — `ruzstd`のバージョン指定(`"0.9"`)が、`[target.'cfg(not(target_arch = "wasm32"))'.dependencies]`と`[target.'cfg(target_arch = "wasm32")'.dependencies]`の2箇所に重複して書かれています。
  - **懸念**: 将来どちらか一方だけバージョンを上げた場合、`cargo`は重複について警告を出さないため、ネイティブビルドとwasm32ビルドで異なる`ruzstd`バージョンが解決される可能性があります。両バックエンドが生成する zstd フレームの相互互換性は同一バージョン前提で確認されているため、バージョンがずれると微妙な非互換のリスクがあります。
  - **提案**: target非依存の`ruzstd`をoptional依存として一本化し(例: `[dependencies]`に`ruzstd = { version = "0.9", optional = true }`を置く)、wasm32側はtarget-gatedな`[target.'cfg(target_arch = "wasm32")'.dependencies.ruzstd] optional = false`相当、もしくはfeatureのtarget別デフォルト指定で有効化する形にまとめると重複を解消できます。対応方針はagyの判断に委ねます。
