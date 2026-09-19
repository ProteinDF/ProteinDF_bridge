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

## 実施内容 (feature/phase10-pr33)

1. **`Pdb` への CONECT レコードサポート追加 (`format/pdb.rs`)**:
   - `Pdb` 構造体に `conects: Vec<(usize, usize)>` フィールドを追加。
   - ゲッター `pub fn conects(&self) -> &[(usize, usize)]` を追加。
   - `parse_str` に `CONECT` レコードパーサーを実装。中心原子 serial（カラム 7-11）と最大4つの相手原子 serial（カラム 12-16, 17-21, 22-26, 27-31）をパース。
   - 自己結合（`serial == partner`）を除外し、`(min(s1, s2), max(s1, s2))` で一意化して PDB 特有の双方向冗長記述や重複を排除して登録。
   - `get_atomgroup()` において、原子登録時に `serial -> Atom` マップを構築し、SSBOND 同様に `model.add_bond(a1, a2, 1)` を呼び出して結合を登録。
   - SSBOND と CONECT の両方で同一結合（例: ジスルフィド結合）が記述されている場合に二重登録されないよう、原子ペアパスを用いた deduplication を導入。
2. **単体テスト・結合テストの実装 (`tests/test_pdb.rs`)**:
   - `test_2mgo_real_pdb_conect`: 実PDBフィクスチャ `2MGO.pdb`（行2872に `CONECT 6 89` が存在）のパースおよび SSBOND との重複排除（全20モデルで各1結合、計20結合）を検証。
   - `test_conect_multiple_partners_synthetic`: 1行に最大4つの結合相手を持つケース（中心炭素 1 に結合する水素 2, 3, 4, 5）および双方向冗長行を含む合成 PDB で、4本の結合が正しくパースされ、`ag.resolve_bond()` で両端原子が解決されることを検証。
   - `test_conect_invalid_error`: CONECT レコード内の serial が整数でない場合のエラーハンドリングを検証。
3. **品質検証**:
   - `cargo clippy --workspace --all-targets -- -D warnings`: 警告ゼロでパス
   - `cargo fmt --check`: 差分なしでパス
   - `cargo test --workspace`: 全テスト PASS (85 unit tests, 90 integration tests = 計175テストすべてパス)

