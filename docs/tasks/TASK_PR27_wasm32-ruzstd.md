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

## レビュー結果(2026-09-19、要修正)

`feature/phase10-pr27` をClaudeがレビューした。`cargo check --target wasm32-unknown-unknown` は成功し(DoDの1は満たす)、`cargo fmt`/`cargo clippy`/`cargo test`も単体では通るが、**実際にビルドして検証したところ、ネイティブターゲット側の既存機能を壊す回帰が2件見つかった。マージ前に修正すること。**

### 実バグ(要修正)

1. **`rust/crates/proteindf-bridge/src/brd.rs`(29行目付近): `zstd_compress`/`zstd_decompress`の選択がCargoの`zstd`/`ruzstd`featureではなく`target_arch`のcfgだけで決まるようになっている。** そのため非wasm32ターゲットでは常に`zstd`(Cバインディング版)が使われ、featureによる明示選択が効かなくなった。
   - **再現手順**: `cargo build --no-default-features --features ruzstd`をネイティブ(非wasm32)で実行する。ビルドは成功するが、ビルドログ上は`zstd`/`zstd-sys`がコンパイルされ、`ruzstd`は一切使われない。エラーも警告も出ないため気づけない。
   - これは`--no-default-features --features ruzstd`(Cツールチェイン不要ビルドの既存エスケープハッチ)を壊す回帰である。修正前のdocコメントにはこの用途が明記されていた。
2. **`rust/crates/proteindf-bridge/Cargo.toml`(22行目付近): `zstd`/`ruzstd`featureが`dep:`指定のない空リスト(`zstd = []`, `ruzstd = []`)になっており、実体を失っている。** 一方`brd.rs`の相互排他`compile_error!`ガード(24〜27行目付近)はそのまま残っている。
   - **再現手順**: `cargo build --all-features`を実行すると、もはや実体のない矛盾のために`compile_error!("features \`zstd\` and \`ruzstd\` are mutually exclusive...")`で失敗する。

### 修正方針

「wasm32では既定で`ruzstd`が自動選択される」という当初の目的は達成しつつ、**ネイティブターゲットでは従来通りCargoの`zstd`/`ruzstd`featureで明示的に選択できる状態を復元すること。** 具体的には、`brd.rs`の`zstd_compress`/`zstd_decompress`の選択をfeatureベースのcfg(`#[cfg(feature = "zstd")]`/`#[cfg(feature = "ruzstd")]`)に戻し、`Cargo.toml`側は`[target.'cfg(target_arch = "wasm32")'.dependencies]`セクションで**wasm32ターゲットの既定featureセットにのみ`ruzstd`を含める**(ネイティブの既定は`zstd`のまま)という形で両立させること。`zstd`/`ruzstd`featureの`dep:`指定・相互排他ガードは実体を伴う形に戻すこと。

### 完了の定義(修正後、再レビュー依頼前に確認すること)

1. `cargo check --target wasm32-unknown-unknown`が追加のフラグなしで成功すること(既存DoD、再確認)。
2. ネイティブで`cargo build --no-default-features --features ruzstd`を実行すると、実際に`ruzstd`がリンクされ`zstd`/`zstd-sys`はコンパイルされないこと(ビルドログで確認)。
3. `cargo build --all-features`が成功すること、または`zstd`/`ruzstd`の相互排他が依然として意味を持つ形で`compile_error!`が発生すること(矛盾のための失敗ではなく)。
4. `cargo test --workspace`・`cargo clippy`・`cargo fmt`が通ること。
5. 修正後、同じ`feature/phase10-pr27`ブランチに追加コミットし、再度ユーザー経由でClaudeにレビュー依頼すること。
