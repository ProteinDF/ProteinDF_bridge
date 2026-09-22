# TASK: ローダーからの暗黙`setup()`呼び出しの撤廃

> §3.14で各フォーマットローダーに追加した「明示的結合情報がない場合、内部で自動的に`AtomGroup::setup()`を呼ぶ」仕様を、ユーザーの判断により撤廃する。理由は`RUST_PORT_SPEC.md` §3.15を参照（パースと結合解決という責務の混在、テストの摩擦、柔軟性の低下）。`AtomGroup::setup()`/`AtomGroup::setup_with_db()`という明示APIは変更せずそのまま維持する。

## ブランチ運用(MUST)

- `develop` から `chore/remove-implicit-loader-setup` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。
- 依存関係: §3.14(`feature/setup-default-entrypoint`)がマージ済みであることを前提にする（マージ済み: `develop`コミット`9c65eae`）。

## 対象

1. 以下のローダーから、§3.14で追加した「`get_bond_list().is_empty()`の場合に`setup()`を呼ぶ」処理を削除し、パースした生の結果（ファイル由来結合があればそれを含む、無ければ結合ゼロ）をそのまま返すようにする:
   - `format/pdb.rs`の`get_atomgroup()`
   - `format/mmcif.rs`の`get_atomgroup()` / `get_structure_atomgroup_for_block()`
   - `format/amber_prmtop.rs`の`get_atomgroup()`
   - `format/gro.rs`の`get_atomgroup()`
   - `format/mol2.rs`の`parse_str()`
2. `AtomGroup::setup()` / `AtomGroup::setup_with_db()` / `Bond::setup_heuristic()` 自体（§3.14・§3.12のロジック）は変更しない。
3. §3.14で追加した`AtomGroup::clear_bonds()`は、ローダーの暗黙実行を前提としたテストのために追加されたものである。暗黙実行の撤廃に伴い、これに依存していたテスト(`test_implicit_setup_loaders_without_bonds`・`test_implicit_setup_preserves_explicit_file_bonds`等)を、明示的に`ag.setup()`を呼ぶ形に書き換える。書き換えた結果`clear_bonds()`が他のテストからも使われなくなった場合は削除して構わない(既存の`test_apply_ccd_bond_templates_1hls_real_pdb`等で引き続き使われている場合は残してよい)。
4. 各ローダーファイル内のdocコメント（`AtomGroup::setup`が自動実行される旨の記述、§3.14で更新したもの）を、明示呼び出しが必要である旨に戻す。
5. `RUST_PORT_SPEC.md` §3.8の呼び出し側推奨パターンは、本タスクの内容に合わせて明示呼び出し形式に先行して更新済み（`ag.setup()`を呼び出し側が明示的に呼ぶ形）。実装がこの記述と整合していることを確認すること。

## 完了の定義(Definition of Done)

1. 上記5ローダーが、明示的結合情報を持たない入力に対して`get_atomgroup()`を呼んだ場合に、結合ゼロの`AtomGroup`を返すことを検証する回帰テストを追加・更新すること(§3.14で追加した暗黙実行の回帰テストを、明示呼び出し前提のテストへ書き換える形でよい)。
2. `ag.setup()`を明示的に呼んだ場合には、従来通りCCDテンプレート＋ヒューリスティックによる結合解決が行われることを確認する既存テスト(`test_atomgroup_setup_1hls_real_pdb`等)が引き続きパスすること。
3. 明示的結合情報を持つ入力（PDBのCONECT、MOL2のBONDセクション、PRMTOPのBONDS等）に対する既存テストが引き続きパスすること(この変更による影響はないはずだが回帰確認)。
4. `cargo clippy --workspace --all-targets -- -D warnings` / `cargo fmt --check` を通すこと。
5. `RUST_PORT_SPEC.md` §3.15に実施内容・完了ステータスを追記すること。

## スコープ外

- CCDテンプレート自体・共有結合半径テーブル自体の変更。
- `AtomGroup::setup()` / `AtomGroup::setup_with_db()` / `Bond::setup_heuristic()`自体のロジック変更。本タスクはローダー側の呼び出し削除とテスト・ドキュメントの追従のみを対象とする。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- ファイル由来の明示的結合情報（CONECT、MOL2 BOND、PRMTOP BONDS等）のパース処理自体は一切変更しない（本タスクは「結合情報が無い場合の暗黙フォールバック」の削除のみが対象であり、明示的結合情報のパースはこれまで通り維持すること）。
