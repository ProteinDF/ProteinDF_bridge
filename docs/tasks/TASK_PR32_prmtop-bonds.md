# TASK_PR32: PRMTOP `BONDS_*`セクションのパース(`format/amber_prmtop.rs`)

> 本タスクは `docs/rust-port-handoff.md` の「Phase 10」節(RUST_PORT_SPEC.md §9対応)から抽出したものです。全体の背景・優先順位・他PRとの関係は同ドキュメントを参照してください。

## ブランチ運用(MUST)

- `develop` から `feature/phase10-pr32` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。承認が出るまで次のタスクの実装に着手しない。
- 依存関係: なし。TASK_PR31(MOL2)・TASK_PR33(PDB CONECT)と並行作業可能。

## 背景

`RUST_PORT_SPEC.md` §9の「高」優先度項目(ファイル由来の明示的な結合トポロジーの読み込み)の一部。**Python版 `amber_prmtop.py` は `ATOM_NAME`/`CHARGE`/`ATOMIC_NUMBER` セクションのみ読み込み、BOND関連セクションは未対応。** 「Python版との1:1比較」という受け入れ基準が使えない完全新規機能のため、Amber公式フォーマット仕様に基づく独立検証で正しさを担保すること。

## 対象

`%FLAG BONDS_WITHOUT_HYDROGEN`/`%FLAG BONDS_INC_HYDROGEN` セクションをパースする。Amber PRMTOP形式ではこれらは `(atom1_idx*3, atom2_idx*3, bond_type_idx)` の3つ組のフラットな整数配列(`%FORMAT(10I8)`)であり、**原子インデックスは0-basedで3倍された値(座標配列オフセット)である点に注意すること。** パース結果を `get_atomgroup()` の結合情報として追加する。

## 完了の定義(Definition of Done)

1. Amber公式フォーマット仕様に基づき、小さな合成PRMTOPフィクスチャ(数原子・数結合)を用意し、期待される結合ペアが正しくパースされることを検証すること。
2. インデックス変換(`/3`、0-based→内部表現)の境界値(最初/最後の原子)を検証すること。
3. `cargo clippy` / `cargo fmt` を通すこと。

## 全フォーマット共通の後続作業(このPR単体の完了条件ではない)

TASK_PR31(MOL2)・TASK_PR32(本タスク)・TASK_PR33(PDB CONECT)の**3つ全てが完了・マージされた後**、以下の方針を `RUST_PORT_SPEC.md`(§2表の備考、または新規節)に明文化する後続タスクがある: 「ファイルに明示的な結合情報があればそれを使い、`Bond::setup()`(VDW半径ヒューリスティック)は呼ばない。ファイルに結合情報がない場合のみ `Bond::setup()` にフォールバックする」という優先順位の確立。**本PRの完了条件には含めない。** 3つのうち最後に着手する担当者(またはユーザー)が別途対応すること。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない(この機能はPython版に存在しないため)。
- Phase 10の他タスク(schema検証・`Bond::setup()`のスケーラビリティ等)には手を出さない。
