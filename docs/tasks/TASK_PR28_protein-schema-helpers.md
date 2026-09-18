# TASK_PR28: protein schemaの形式化と検証ヘルパー

> 本タスクは `docs/rust-port-handoff.md` の「Phase 10」節(RUST_PORT_SPEC.md §9対応)から抽出したものです。全体の背景・優先順位・他PRとの関係は同ドキュメントを参照してください。

## ブランチ運用(MUST)

- `develop` から `feature/phase10-pr28` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。承認が出るまで次のタスクの実装に着手しない。
- 依存関係: なし。着手推奨順序ではPR#27の次。

## 着手前に必ず読むこと(重要、二重実装を避けるため)

**`rust/crates/proteindf-bridge/src/format/mod.rs` の `Format` 構造体(Phase 2 PR#4で実装済み)を必ず読むこと。** `Format::is_residue`/`is_chain`/`is_protein`/`is_models` が「直下に原子を持たない」「サブグループが次階層の条件を満たす」という構造的判定を既に提供しており、§9が要求する検証機能の大部分をカバーしている。ゼロから設計しないこと。本タスクはこれを土台にする。

## 背景

`RUST_PORT_SPEC.md` §9の「高」優先度項目。「結(YUI)」側の調査で、`/model_N/chain_id/res_key/atom_key` というパス深さによるmodel/chain/residueの区別は、現状 `biopdb.py` 等の実装コードにのみ暗黙的に存在し、Rust版 `atom_group.rs` にはこれを検証する関数(`is_model_level()`/`is_chain_level()`/`is_residue_level()` 相当)が無いと指摘されている。`AtomGroup` 自体はスキーマレスな汎用木のため、この規約を破るデータ(残基ラッパーなしでchain直下に置かれるHETATM/水分子等)が来ても検出できない。

## 対象

1. `/model_N/chain_id/res_key/atom_key` というパス階層規約を、`RUST_PORT_SPEC.md`(または `atom_group.rs` のモジュールdocコメント)に明文化する。
2. `AtomGroup` に `is_model_level()`/`is_chain_level()`/`is_residue_level()` を追加する。これは `path()`(例: `/model_1/A/6/`)の**パス深さ**に基づく位置的判定とする。`Format::is_chain` 等の**構造的**判定(原子を直接持たない・サブグループが条件を満たす)とは軸が異なることをdocコメントで明記すること——正常データでは両者は一致するが、規約違反データ(残基ラッパーなしでchain直下に置かれたHETATM/水分子等)ではパス深さは「chainレベル」のままなのに構造判定は崩れる、という乖離が生じる。この乖離こそが検出したい違反である。
3. 規約違反を列挙するヘルパー(例: `AtomGroup::validate_schema() -> Vec<SchemaViolation>`)を新設する。

## 完了の定義(Definition of Done)

1. 正常な階層構造(model→chain→residue→atom)で3つの `is_*_level()` 全てが期待通りの値を返すことをテストすること。
2. 規約違反データ(chain直下にHETATM/水分子を直接配置したもの)を合成し、`validate_schema()` がこれを検出することをテストすること。
3. `Format` の既存メソッドとの役割の違いをdocコメントで明記すること。
4. `cargo clippy` / `cargo fmt` を通すこと。

## スコープ外

- mmCIF/PDB/PRMTOPパーサ自体の変更(規約違反データを検出できるようにするだけで、弾く・直すのは対象外)。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- **`Format`(`format/mod.rs`)の既存実装を必ず読んでから設計すること。** 車輪の再発明をしない。
- Phase 10の他タスク(ファイル由来結合・`Bond::setup()`のスケーラビリティ・二次構造書き戻し等)には手を出さない。
