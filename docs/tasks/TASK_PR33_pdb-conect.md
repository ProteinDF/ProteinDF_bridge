# TASK_PR33: PDB CONECTレコードの読み込み(`format/pdb.rs`)

> 本タスクは `docs/rust-port-handoff.md` の「Phase 10」節(RUST_PORT_SPEC.md §9対応)から抽出したものです。全体の背景・優先順位・他PRとの関係は同ドキュメントを参照してください。

## ブランチ運用(MUST)

- `develop` から `feature/phase10-pr33` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。承認が出るまで次のタスクの実装に着手しない。
- 依存関係: なし。TASK_PR31(MOL2)・TASK_PR32(PRMTOP)と並行作業可能。

## 背景

`RUST_PORT_SPEC.md` §9の「高」優先度項目(ファイル由来の明示的な結合トポロジーの読み込み)の一部。**Python版 `biopdb.py` はSSBOND(ジスルフィド)のみを結合として取り込み、CONECTレコードは未対応。** 「Python版との1:1比較」という受け入れ基準が使えない完全新規機能のため、独立検証(合成/実PDBフィクスチャでの手動検証)で正しさを担保すること。

## 対象

`CONECT` レコード(serial番号1つ+最大4つの結合相手serial番号)をパースし、`Pdb::get_atomgroup()` にSSBOND同様の方法で結合情報として追加する。serial番号からAtomGroup内の実際の原子への対応付けが必要(既存のSSBOND実装がserial→pathマッピングを持っていれば再利用すること)。

## 完了の定義(Definition of Done)

1. CONECTレコードを含む実PDBフィクスチャ(既存フィクスチャになければ新規に小さな合成PDBフィクスチャを用意する)でパース結果を検証すること。
2. 1つのCONECT行に複数の結合相手(最大4つ)が書かれているケースを検証すること。
3. `cargo clippy` / `cargo fmt` を通すこと。

## 全フォーマット共通の後続作業(3つのPR全て完了後、いずれかの担当者が対応)

TASK_PR31(MOL2)・TASK_PR32(PRMTOP)・TASK_PR33(本タスク)の**3つ全てが完了・マージされた後**、以下の方針を `RUST_PORT_SPEC.md`(§2表の備考、または新規節)に明文化すること:

> ファイルに明示的な結合情報があればそれを使い、`Bond::setup()`(VDW半径ヒューリスティック)は呼ばない。ファイルに結合情報がない場合のみ `Bond::setup()` にフォールバックする。

各 `get_atomgroup()` はこの3PR完了後、結合情報が取得できればそれを設定済みの状態で `AtomGroup` を返す想定なので、呼び出し側(YUI)が `AtomGroup::get_number_of_bonds() > 0` 等で判定してから `Bond::setup()` を呼ぶかどうかを決める、という利用パターンをドキュメント化する。**本PR単体の完了条件には含めない**——3つのうち最後にマージされたブランチの担当者(またはユーザー)が、3つ全てのマージを確認してから対応すること。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない(この機能はPython版に存在しないため)。
- Phase 10の他タスク(schema検証・`Bond::setup()`のスケーラビリティ等)には手を出さない。
