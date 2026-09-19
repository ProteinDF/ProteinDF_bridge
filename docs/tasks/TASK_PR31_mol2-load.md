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

## 実施内容と検証結果 (2026-09-19)

### 1. `SimpleMol2` への読み込み API 追加 (`format/mol2.rs`)
- `from_file`, `from_str`, `load`, `parse_str` メソッドを追加。
- `@<TRIPOS>MOLECULE`, `@<TRIPOS>ATOM`, `@<TRIPOS>BOND` の各セクションをパース。
- SYBYL atom type（例: `"C.3"`, `"N.pl3"`, `"O.3"`, `"C.ar"`）および原子名から元素記号を正確に推定（`deduce_symbol_from_mol2`）。
- 座標（x, y, z）、部分電荷（charge）をパースし、各原子を格納。
- 結合種別（1, 2, 3, ar, am, du, un 等）を解釈し、`AtomGroup::add_bond` を用いて結合トポロジーを構築。
- `get_atomgroup(&self) -> &AtomGroup` を追加。

### 2. ラウンドトリップ検証 (`format/mol2.rs` インラインテスト)
- `test_mol2_roundtrip`: `SimpleMol2::save` / `get_text` で出力した MOL2 を `from_str` で読み戻し、分子名・原子数・座標・結合（ペア・次数）および再出力テキストの完全一致を検証。
- `test_mol2_from_file_roundtrip`: 一時ファイルを経由したファイル保存・読み込みのラウンドトリップを検証。
- `test_mol2_load_synthetic_fixture`: 芳香族結合 (`ar`)、部分電荷、SYBYL 型を含む合成フィクスチャのパースを検証。

### 3. 合成フィクスチャと統合テストの追加
- `tests/data/sample.mol2`: エタノール（9原子、8結合、部分電荷付き）の合成フィクスチャを追加。
- `tests/test_mol2.rs`: `sample.mol2` の読み込み・原子座標・電荷・結合トポロジーの手動検証テストを追加。

