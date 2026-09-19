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

## レビュー結果(2026-09-19、要修正)

`feature/phase10-pr28` をClaudeがレビューした(独立検証エージェントによる再現確認込み)。正常系のテスト・`Format`との役割分担のdocコメント・`cargo clippy`/`cargo fmt`/`cargo test`はいずれも問題ない。**しかし`validate_schema()`が導入目的(スキーマ違反の検出)を果たせない実バグが見つかった。マージ前に修正すること。**

### 実バグ(要修正、最重要)

1. **`path_depth()`(`atom_group.rs`、163行目付近)が、木構造を辿った実際のネスト深さではなく、`self.path`文字列中の`/`の個数を数えているだけである。** `set_group`/`set_path`が`format!("{}{}/", self.path, key)`という単純な文字列連結でpathを構築しており、`key`自体に`/`が含まれないことを検証していない。
   - **再現手順(独立検証エージェントが実際にテストコードを追加して確認済み)**: キー`"A/B"`を持つグループをmodel直下にchainとして付け、直下に原子を置く(=chain直下に原子があるという明確な規約違反)。`path_depth()`が本来の2ではなく3と計算され、`is_residue_level()`が`true`を返してしまう。結果、`validate_schema()`はこの構造全体に対して**違反0件**を返す——本来検出すべき違反を完全に見逃す。
   - これは「パス文字列の見た目」ではなく「木の実際の深さ」を使うべき典型的なバグである。`validate_schema()`はまさに「信頼できない入力から規約違反を検出する」ための機能であり、その検証ロジック自体が細工されたキー名(または`/`をたまたま含む正当なキー名)で無力化されてしまうのは本末転倒。
2. **同じ`path_depth()`が空文字列キーの連鎖でも深さを誤カウントする。** `path_depth()`は空セグメントを`.filter(|s| !s.is_empty())`で除外する一方、`set_group`等は空文字列キーそのものを拒否しない。
   - **再現手順**: root→`""`→`""`→`""`→葉(原子あり)という実際のネスト深さ4の構造を作ると、計算上の深さは0になり、違反自体は報告されるものの`path`/`depth`が`DirectAtomsAtNonResidueLevel { path: "////", depth: 0 }`のような無意味な値になる(本来は`ExcessiveDepth`(depth 4)も併せて検出されるべきところ、誤カウントのためルートノードに見えて検出されない)。キーの組み合わせ次第では完全な見逃しも起こりうる。

### 修正方針

`path_depth()`を、pathの文字列表現から逆算するのではなく、**木を実際に辿って(親から子へ再帰呼び出し時にカウンタを渡す、またはグループ挿入時に深さを記録する等)正しいネスト深さを計算する**方式に変更すること。あわせて、`set_group`(および`set_path`/`update_paths`が使う経路)で`key`が空文字列または`/`を含む場合にエラーを返す(または`BrError`系で拒否する)ガードを追加することを推奨する(規約違反データの「検出」だけでなく「そもそも壊れたpathを作らせない」という二段構えにできる)。

### 副次的な指摘(必須ではないが対応推奨)

3. **`atom_group.rs`のSPDXヘッダーが本PRで失われている。** 今回追加した`//!`モジュールdocコメントに置き換わる形で、他の全ファイル(`atom.rs`等)が持つ`SPDX-FileCopyrightText`/`SPDX-License-Identifier`ヘッダーが消えている。復元すること。
4. **インラインテスト`test_schema_violations_subgroups_in_residue_and_excessive_depth`(334行目付近)が`matches!`でvariantの種類だけを検証しており、`path`/`depth`/`group_keys`の値まで検証していない。** 上記の実バグ1・2はこのテストでは検出できない。`path`/`depth`等の値まで厳密に検証するテストに強化するか、既存の`tests/test_schema.rs`側の厳密なテスト(`test_schema_violation_subgroup_in_residue_and_depth`)と同等の厳密さに揃えること。

### 完了の定義(修正後、再レビュー依頼前に確認すること)

1. 上記「再現手順」の2ケース(キーに`/`を含む場合、空文字列キーを連鎖させる場合)を回帰テストとして追加し、正しい深さ・正しい違反が検出されることを確認すること。
2. 通常のPDB/mmCIF由来の正常系データでは既存の挙動が変わらないことを確認すること(既存テストが全てパスすること)。
3. 可能であればSPDXヘッダーを復元すること。
4. `cargo clippy` / `cargo fmt` / `cargo test`を通すこと。
5. 修正後、同じ`feature/phase10-pr28`ブランチに追加コミットし、再度ユーザー経由でClaudeにレビュー依頼すること。
