# TASK_PR31: MOL2読み込み(`format/mol2.rs`)

> 本タスクは `docs/rust-port-handoff.md` の「Phase 10」節(RUST_PORT_SPEC.md §9対応)から抽出したものです。全体の背景・優先順位・他PRとの関係は同ドキュメントを参照してください。

## ブランチ運用(MUST)

- `develop` から `feature/phase10-pr31` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。承認が出るまで次のタスクの実装に着手しない。
- 依存関係: なし。TASK_PR32(PRMTOP)・TASK_PR33(PDB CONECT)と並行作業可能。

## 背景

`RUST_PORT_SPEC.md` §9の「高」優先度項目(ファイル由来の明示的な結合トポロジーの読み込み)の一部。**Python版 `mol2.py` の `SimpleMol2` は書き込み専用で、`load`/読み込みに相当するメソッドが一切存在しない。** 「Python版との1:1比較」という受け入れ基準が使えない完全新規機能のため、独立検証(ラウンドトリップ・合成フィクスチャでの手動検証)で正しさを担保すること。

## 対象

`SimpleMol2` に `load`/`from_str`(既存の `save`/`get_text` と対になる読み込み)を追加し、`@<TRIPOS>ATOM`・`@<TRIPOS>BOND` セクションをパースして結合情報付きの `AtomGroup` を構築する `get_atomgroup()` 相当のメソッドを追加する。

## 完了の定義(Definition of Done)

1. 自身の `save()` が生成したMOL2テキストを `load` でラウンドトリップし、原子数・座標・結合(原子ペア・結合次数)が保持されることを検証すること。
2. 小さな合成MOL2フィクスチャ(新規追加)でパース結果を手動検証すること。
3. `cargo clippy` / `cargo fmt` を通すこと。

## 全フォーマット共通の後続作業(このPR単体の完了条件ではない)

TASK_PR31(本タスク)・TASK_PR32(PRMTOP)・TASK_PR33(PDB CONECT)の**3つ全てが完了・マージされた後**、以下の方針を `RUST_PORT_SPEC.md`(§2表の備考、または新規節)に明文化する後続タスクがある: 「ファイルに明示的な結合情報があればそれを使い、`Bond::setup()`(VDW半径ヒューリスティック)は呼ばない。ファイルに結合情報がない場合のみ `Bond::setup()` にフォールバックする」という優先順位の確立。**本PRの完了条件には含めない。** 3つのうち最後に着手する担当者(またはユーザー)が別途対応すること。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない(この機能はPython版に存在しないため)。
- Phase 10の他タスク(schema検証・`Bond::setup()`のスケーラビリティ等)には手を出さない。
