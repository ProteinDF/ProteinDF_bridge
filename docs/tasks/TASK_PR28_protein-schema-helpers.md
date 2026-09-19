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

## 再レビュー結果(2026-09-19、修正コミット`c4b8b16`確認・マージブロッカーなし)

上記の実バグ1・2、SPDXヘッダー、テストの厳密化は全て修正・確認済み(回帰テスト`test_schema_regression_key_with_slash`/`test_schema_regression_empty_string_keys`を含め、ビルド・`cargo clippy -D warnings`・`cargo fmt --check`・`cargo test --workspace`は全てパス)。**この修正はマージ可能。**

以下、Minorな残存懸念が1件あります。マージをブロックするものではありませんが、対応する場合はご確認ください。

5. **`atom_group.rs`の`set_path()`(261〜263行目付近)に、`self.depth == 0 && self.path != "/"`のときだけ旧来の文字列分割方式(`path.split('/').filter(...).count()`)で`depth`を再計算するfallbackが残っている。** これは`set_group`経由(=`load_atomgroup`/`load_brd_yui`など外部データ読み込み経路)では`depth`が事前に非ゼロにセットされるため発火せず、**実際の読み込み経路には影響しないことを確認済み**。しかし`set_path`自体は`pub fn`であり、`modeling.rs`のACE/NMEキャップ生成(`answer.set_path("/ACE")`等)や`tests/test_ssbond.rs`で、未アタッチのfreshなgroupに直接パスを与える用途に実際に使われている。
   - **再現手順(検証済み、一時テストで確認後revert)**: `AtomGroup::new()`(本来depth=0)に対して直接`g.set_path("/A/B".to_string())`を呼ぶと、`g.path_depth()`が本来の0ではなく2と誤計算される。
   - 現状の呼び出し箇所(`modeling.rs`の"ACE"/"NME"、`test_ssbond.rs`の"model_1")はいずれも単一セグメント・スラッシュなしのリテラルなので実害はないが、「文字列ではなく真の木構造深さを使う」という本PRの設計意図・docコメントの説明("ensuring robustness against keys containing slashes or empty segments")と矛盾する経路が残っている。将来この関数を攻撃者制御パスに対して直接呼ぶコードが追加されると、修正したはずのバグが別の入口から再発しうる。
   - **提案(任意)**: このfallbackを削除し、「`depth`は常に`set_group`/`update_paths`経由でのみ設定され、ルートの初期値0のみが正」という不変条件に統一する。`modeling.rs`側で直接`set_path`を呼んでいる箇所は、深さ管理が必要なら`set_group`で組み立てるよう見直すか、不要なら現状維持でよい。対応方針はagyの判断に委ねる。

## 再々レビュー結果(2026-09-19、修正コミット`cf1b68c`確認・新たな要修正あり)

上記5番の指摘どおり`set_path()`のfallbackは削除されました(ビルド・`cargo clippy -D warnings`・`cargo fmt --check`・`cargo test --workspace`は全てパス)。しかしこの削除の仕方により、**`modeling.rs`の既存呼び出しパターンで新たな不整合が生じる**ことを確認しました。マージ前の対応を推奨します(前回までの実バグ1・2ほど深刻ではありませんが、5番より一段重い"要検討"扱いとしています)。

### 新たに判明した問題(要検討)

6. **fallback削除により、「別の場所からclone/抽出したサブツリーに`set_path()`で新しいパスを直接与える」パターンで`depth`が古い(クローン元の)値のまま取り残され、新しい`path()`と矛盾する状態になる。** `set_path()`はもはや`depth`を一切更新しないため、`update_paths()`は誤った`self.depth`を基準に子孫の`depth`を再配布し続ける。
   - **再現手順(検証済み、一時テストで確認後revert)**: 木構造中で深さ2に位置するグループ(`/A/1/`)を`clone()`し、`set_path("/ACE".to_string())`で単一セグメントのパス(本来は深さ1)に付け替えると、`path()`は`/ACE/`になるが`path_depth()`は2のまま(本来期待される1にならない)。
   - **実際に該当するコード**: `rust/crates/proteindf-bridge/src/modeling.rs`の`get_ACE`(101〜102行目付近)・`get_NME`(170〜171行目付近)が、まさにこのパターン(`best.get_group(...)`で取得した既存木の一部を`.clone()`し、`set_path("/ACE")`/`set_path("/NME")`で付け替え)を使っている。現状これらの戻り値に対して`path_depth()`/`is_residue_level()`/`validate_schema()`を呼ぶ本番コード・テストはまだ存在しないため即座の実害は顕在化していないが、これらのキャップ残基を将来ペプチドモデルに組み込んで検証する処理を書いた瞬間に、`depth`と`path`が矛盾したデータとして混入する。
   - なお、`load_atomgroup`/`load_brd_yui`などの外部データ読み込み経路(`set_group`のみ使用)には影響なし。
3. **修正方針(提案)**: `set_path()`を「常に`depth`を`path`から独立に扱う(呼び出し元が`set_group`で正しく再アタッチするまでは`depth`は未定義扱い)」という設計にするなら、`set_path()`内で`depth`を（文字列パースではなく)**確実に無効化する**——例えば「この経路で`path`が変更されたら`depth`もその場で`path`のセグメント数から算出し直す」のではなく、**`set_path`を直接公開APIとして深さ整合性が必要な用途に使わせない**方針が良い。具体的には:
   - (a) `modeling.rs`のget_ACE/get_NMEを、`.clone()` + `set_path()`ではなく、新規`AtomGroup`を`with_name`で作って親に`set_group`で正しく再アタッチする(あるいは深さをリセットする専用メソッド、例えば`detach_and_set_path()`が`depth`もpathのトップレベル呼び出し時点で明示的に0または1にリセットする)形に直す。
   - (b) 最小限の対応として、`set_path()`のドキュメントコメントに「このメソッドは`depth`を更新しない。クローンしたサブツリーを付け替える場合は呼び出し元が`depth`の整合性に責任を持つこと(理想的には`set_group`経由で親に再アタッチすること)」という注意書きを追加する。
   - いずれの対応も必須ではありませんが、`modeling.rs`のget_ACE/get_NMEの戻り値に対して将来`validate_schema()`等を呼ぶ計画があるなら、(a)を強く推奨します。

### 完了の定義(対応する場合)

1. 上記6番の再現手順を回帰テストとして追加する(clone後の`set_path`でdepthが新パスと矛盾しないことを確認)。
2. `modeling.rs`のget_ACE/get_NMEを見直すか、`set_path()`のdocコメントに注意書きを追加する。
3. `cargo clippy` / `cargo fmt` / `cargo test`を通すこと。
4. 対応してもしなくても、対応方針をユーザー経由でClaudeに報告すること。**現時点でこの1件のみを理由にマージをブロックする必要はありません**(実際に踏むコードパスが現状存在しないため)が、agyの判断で対応するか、既知の制約としてドキュメント化するかを決めてください。

## 指摘6対応結果 (2026-09-19)

提案(a)および(b)の両方を採用して対応を完了しました。

1. **`AtomGroup::set_path_with_depth` の新設**:
   - `pub fn set_path_with_depth(&mut self, mut new_path: String, depth: usize)` を追加。サブツリーの切り離しやスタンドアロンのルート作成時に、任意の基準パスと深さを整合させて再設定可能とし、子孫の `depth`/`path` も正しく再帰更新する。
2. **`set_path` の仕様明確化と doc コメント拡充**:
   - `set_path` は `self.depth` を維持する仕様であることをドキュメントに明記し、深さを変更してデタッチ・再配置する場合は `set_path_with_depth` または `set_group` を使用すべき注意書きを追加。
3. **`modeling.rs` の ACE / NME 生成処理の更新**:
   - `get_ACE` および `get_NME` で `answer.set_path_with_depth("/ACE".to_string(), 1)` / `answer.set_path_with_depth("/NME".to_string(), 1)` を使用するよう修正。
4. **回帰テスト追加**:
   - `tests/test_schema.rs` に `test_schema_regression_cloned_subtree_set_path_depth` を追加し、クローンしたサブツリーに `set_path_with_depth` を呼んだ際の `path_depth()`、子グループ/原子の path・depth 整合性を検証。

