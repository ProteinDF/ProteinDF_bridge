# TASK_PR29: 二次構造情報の`AtomGroup`への書き戻し

> 本タスクは `docs/rust-port-handoff.md` の「Phase 10」節(RUST_PORT_SPEC.md §9対応)から抽出したものです。全体の背景・優先順位・他PRとの関係は同ドキュメントを参照してください。

## ブランチ運用(MUST)

- `develop` から `feature/phase10-pr29` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。承認が出るまで次のタスクの実装に着手しない。
- 依存関係: なし(PR#28とは独立)。着手推奨順序ではPR#28の次。

## 背景

2026-09-19に `RUST_PORT_SPEC.md` §9へ追記された新規「高」優先度項目。`calc_secondary_structure(chain: &AtomGroup) -> Vec<SecondaryStructure>`(Phase 8 PR#23で実装済み、`secondary_structure.rs`)は結果を別のVecとして返すのみで、`AtomGroup` ツリー自体には反映されない(`AtomGroup` に汎用メタデータフィールドが無いため)。一方 `Bond::setup()` は `mol.add_bond(...)` で結果を `AtomGroup` 自体に書き戻す設計になっており、一貫していない。YUI側は現状、residueのpath文字列をキーとする一時的なサイドマップで代替している(bridge側の対応までの暫定措置)。

## 対象

1. `bonds: Vec<BondRecord>` と同格の、residueレベルの `AtomGroup` が持つ専用フィールド `secondary_structure: Option<SsCode>` を追加する(汎用メタデータ袋ではなく、`bonds` と同じ「specific typed field」パターンをYUI側は希望している)。
2. `calc_secondary_structure` と対になる `apply_secondary_structure(chain: &mut AtomGroup)`(算出結果を対応するresidueグループの `secondary_structure` フィールドに書き戻す関数)を追加する。

## 完了の定義(Definition of Done)

1. `apply_secondary_structure` を呼んだ後、chain内の各residueグループの `secondary_structure` フィールドが、`calc_secondary_structure` が返すVecの対応するエントリと一致することをテストすること(Phase 8の基準値データ、`1hls.pdb` のchain A/Bで検証)。
2. `AtomGroup` の `Clone`・マージ演算(`merge`/`BitAnd`/`BitOr`/`BitXor`)が新フィールドを正しく扱う(消えない・上書きロジックが妥当)ことを確認すること。**Phase 1是正事項5で `bonds` フィールドのマージ漏れが実際にバグとして見つかった前例があるため、同じ轍を踏まないこと。**
3. `cargo clippy` / `cargo fmt` を通すこと。

## スコープ外

- 8状態DSSP分類への拡張(Phase 8のスコープ外のまま)。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない(この機能はPython版に存在しないため)。
- **`AtomGroup` のマージ・集合演算全て(`merge`/`BitAnd`/`BitOr`/`BitXor`/`Clone`)を洗い出してから着手すること。** 新フィールド追加のたびに一部の演算だけ対応漏れするパターンがPhase 1是正事項5で実際に起きている。
- Phase 10の他タスク(schema検証・ファイル由来結合・`Bond::setup()`のスケーラビリティ等)には手を出さない。
