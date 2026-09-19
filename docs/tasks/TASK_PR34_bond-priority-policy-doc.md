# TASK_PR34: ファイル由来結合とVDWヒューリスティックの優先順位確立(ドキュメントのみ)

> 本タスクは `docs/rust-port-handoff.md` の「Phase 10」節(RUST_PORT_SPEC.md §9対応)から抽出したものです。全体の背景・優先順位・他PRとの関係は同ドキュメントを参照してください。TASK_PR31(MOL2)・TASK_PR32(PRMTOP)・TASK_PR33(PDB CONECT)の3つ全てがマージされたことを受けて起票された、ドキュメントのみのフォローアップタスクです。

## ブランチ運用(MUST)

- `develop` から `feature/phase10-pr34` を作成して作業する。
- **`develop` へは自分でマージしない。** 作業・`cargo clippy`/`cargo fmt`確認(コード変更は無いが、`RUST_PORT_SPEC.md`にコード例を含める場合は整合性を確認すること)が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。
- 依存関係: TASK_PR31・TASK_PR32・TASK_PR33が全てマージ済みであること(2026-09-19時点で3つとも`develop`にマージ済み)。TASK_PR35+(`Bond::setup()`のスケーラビリティ改善)とは独立で、順序を問わない。

## 背景

TASK_PR31〜33により、`SimpleMol2::get_atomgroup()`・`AmberPrmtop::get_atomgroup()`・`Pdb::get_atomgroup()`のいずれもが、ファイルに明示的な結合情報(MOL2の`@<TRIPOS>BOND`、PRMTOPの`BONDS_*`、PDBの`CONECT`/`SSBOND`)があればそれをパースして`AtomGroup`の結合トポロジーとして設定するようになった。一方、`Bond::setup()`(VDW半径ヒューリスティックによる距離ベースの結合推定)は既存のまま残っており、どちらを優先すべきかという方針がコード上どこにも明文化されていない。「結(YUI)」側がこの3フォーマットのローダーを使う際、`Bond::setup()`と重複して呼んでしまう・あるいはどちらを信頼すべきか迷う、という混乱を避けるため、優先順位をドキュメントとして確立する。

## 対象

1. `RUST_PORT_SPEC.md`(§2表の備考、または新規節)に以下の方針を明文化すること:
   > ファイルに明示的な結合情報があればそれを使い、`Bond::setup()`(VDW半径ヒューリスティック)は呼ばない。ファイルに結合情報がない場合のみ`Bond::setup()`にフォールバックする。
2. 各`get_atomgroup()`(`SimpleMol2`/`AmberPrmtop`/`Pdb`)は、結合情報が取得できればそれを設定済みの状態で`AtomGroup`を返す、という現状の実装が上記方針と一致していることを確認し、その旨を明記すること。
3. 呼び出し側(YUI等)向けの利用パターンとして、「`AtomGroup::get_number_of_bonds() > 0`等で判定してから`Bond::setup()`を呼ぶかどうかを決める」という指針をドキュメント化すること。

## 完了の定義(Definition of Done)

1. `RUST_PORT_SPEC.md`に上記の優先順位方針が明文化されていること。
2. `SimpleMol2`/`AmberPrmtop`/`Pdb`の`get_atomgroup()`のdocコメントが、ファイル由来の結合情報を優先する旨(またはそれに反する現状があればその旨)を反映していること。矛盾が見つかった場合は、コード修正ではなくClaudeへの報告を優先すること(本タスクはドキュメントのみが完了条件であり、コード修正が必要と判断した場合は別タスクとして切り出す)。
3. `docs/rust-port-handoff.md`のPhase 10節(またはその完了報告)に、本タスクの完了を記録すること。

## スコープ外

- `Bond::setup()`自体の実装変更(スケーラビリティ改善はTASK_PR35+の担当)。
- MOL2/PRMTOP/PDB以外のフォーマット(例: mmCIF)への同様の結合優先順位の適用(将来必要になれば別タスク)。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- コード変更を伴う作業(`get_atomgroup()`の挙動変更など)が必要だと判断した場合、独断で実装せず、まずユーザー経由でClaudeに報告し指示を仰ぐこと。
- Phase 10の他タスク(`Bond::setup()`のスケーラビリティ等、TASK_PR35+のスコープ)には手を出さない。

## 実施内容 (feature/phase10-pr34)

1. **`RUST_PORT_SPEC.md` の更新**:
   - 新規節「3.8 結合情報の優先順位方針（ファイル由来結合 vs VDWヒューリスティック）」を追加。
   - 「ファイルに明示的な結合情報があればそれを優先して使い、`Bond::setup()`は呼ばない。ファイルに結合情報がない場合のみ`Bond::setup()`にフォールバックする」という基本方針を明文化。
   - 各フォーマット（MOL2, PRMTOP, PDB）の結合パース実装状況を明記。
   - 呼び出し側（「結 (YUI)」等）における推奨利用パターン（`ag.get_bond_list().is_empty()` による分岐コード例）を記述。
   - §2 モジュール対応表および §9 の該当項目に完了ステータス・3.8節への参照を追記。
2. **各フォーマットの doc コメント整備**:
   - `SimpleMol2::get_atomgroup()` (`format/mol2.rs`)
   - `AmberPrmtop::get_atomgroup()` (`format/amber_prmtop.rs`)
   - `Pdb::get_atomgroup()` (`format/pdb.rs`)
   それぞれの doc コメントに、ファイル由来の明示的結合情報が設定された `AtomGroup` を返す旨、および `RUST_PORT_SPEC.md` §3.8 の優先順位ポリシーへの参照を追記。
3. **引き継ぎドキュメントの更新**:
   - `docs/rust-port-handoff.md` の Phase 10 PR#34 節に完了記録を追記。
4. **品質検証**:
   - `cargo clippy --workspace --all-targets -- -D warnings`: 警告ゼロでパス
   - `cargo fmt --check`: 差分なしでパス
   - `cargo test --workspace`: 全テスト PASS (85 unit tests, 90 integration tests = 計175テストすべてパス)

