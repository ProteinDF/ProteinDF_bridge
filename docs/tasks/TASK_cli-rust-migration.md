# TASK: コマンドのRust版への移行 (PR#52〜55)

> コマンド32本のうち、構造を扱う26本を純Python版からRust版(`proteindf_bridge.rs`)に切り替える。設計の根拠・スコープの判断は`RUST_PORT_SPEC.md` §4.5にある(ユーザー確認済み、2026-10-11)。**着手前に§4.5を必ず精読すること。** 本ファイルは実行チェックリストである。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaude(`/code-review`)が担当する。
- 着手順はPR#52 → PR#53 → PR#54 → PR#55。各PRは前のPRがレビューで承認され、developにマージされてから始める。
- `develop`から`feature/cli-migration-prNN`を切って作業する。**`develop`へは自分でマージしない。** タグやReleaseは作らない。
- **完了報告は、起動時に指示されたパス(リポジトリの外)に、決まった見出しの短いファイルとして書く。** コマンドの出力は要約・再構成せず、実際の出力行をそのまま貼る。

## 全PR共通のルール

1. コマンド名・引数・オプションは変えない。
2. 純Python版のパッケージ(`proteindf_bridge`の各モジュール)は削除・変更しない。コマンド(`proteindf_bridge/cli/`)だけを書き換える。
3. 移行に必要なRust版のAPIが足りない場合は、Pythonバインディング(`rust/crates/proteindf-bridge-py`)に追加する。コアのクレート(`rust/crates/proteindf-bridge`)の変更が必要な場合は、変更せずに報告する。
4. Pythonバインディングに追加する場合は、これまでの方針(`RUST_PORT_SPEC.md` §4.3「共通の設計方針」、`docs/tasks/TASK_py-bindings.md`の確立した方針)に従う。
5. テストは、新しい一時環境にパッケージを`pip install`して、`python -P -m unittest discover -s tests`で実行する。Rustは`cargo nextest run --workspace --status-level fail --final-status-level fail`と`cargo test --doc`、`cargo clippy --workspace --all-targets -- -D warnings`、`cargo fmt --all -- --check`。

## PR#52: 移行前の出力の固定と、APIの対応表

### 対象

1. 26本のコマンドそれぞれについて、リポジトリのフィクスチャ(`proteindf_bridge/data/`・`rust/crates/proteindf-bridge/tests/data/`)を入力にして**純Python版のまま**実行し、出力(書き出したファイルと標準出力)を基準ファイルとして保存し、現在の出力と一致することを確かめるテストを作る(例: `tests/test_cli_golden.py`と`tests/cli_golden/`)。各コマンドで、主な引数・オプションの組み合わせを少なくとも1つずつ含める。
   - 出力に実行ごとに変わる値(日時など)や浮動小数点の表記揺れがある場合は、比較の前に正規化する方法を決めて、その理由をテストに書く。
   - `.brd`(MessagePack)の比較は、バイト列ではなく、読み込んだ内容(原子名・元素・座標・結合など)で行う。
2. 対応表`docs/cli-migration-notes.md`を作る。各コマンドについて、使っている純Python版のAPIと、それに対応するRust版のAPI(`proteindf_bridge.rs`)の有無、移行するPR(PR#53〜55)を書く。足りないAPIがあれば一覧にする。

### 完了の定義

1. 26本すべてに基準ファイルとテストがあり、純Python版のままで通る。対象外の6本はテストに含めない(理由を書く)。
2. 対応表に26本すべてが載っており、足りないRust版のAPIの一覧がある。
3. 全PR共通のルール5のテストが通る。

## PR#53: フォーマット変換(10本)

`pdb2brd.py`・`brd2pdb.py`・`brd2gro.py`・`gro2brd.py`・`brd2xyz.py`・`xyz2brd.py`・`brd2txt.py`・`mmcif2txt.py`・`mmcif2mol2.py`・`read_amber_prmtop.py`

### 方針(Claude決定、2026-10-11)

1. 対応表(`docs/cli-migration-notes.md` §2)に従い、構造の読み書きは`proteindf_bridge.rs`の`load_atomgroup`・`save_atomgroup`・`Pdb`・`SimpleGro`・`Xyz`・`SimpleMmcif`・`SimpleMol2`・`AmberPrmtop`を使う。
2. **`AtomGroup`ではない汎用データ**(`mmcif2txt.py`・`mmcif2mol2.py`・`read_amber_prmtop.py`が書き出す辞書)のMessagePack読み書きは、Pythonの`msgpack`パッケージ(既に`pyproject.toml`の依存にある)をコマンドの中で直接使う。`proteindf_bridge.rs`にMessagePackの汎用関数は追加しない。構造データではないため。
3. **`AtomGroup.save_csv`**(`brd2txt.py`)は、Pythonバインディング(`rust/crates/proteindf-bridge-py`)に追加する。列・書式は純Python版(`proteindf_bridge/atomgroup.py`の`save_csv`)と同じにし、Pythonのテストを追加する。コアのクレートは変更しない(必要なら報告)。書き出しは途中で失敗しても壊れたファイルを残さない(一時ファイルに書いて名前を変える。PR#42と同じ)。
4. `read_amber_prmtop.py`の出力先が`output.brd`に固定されている点は、コマンドの引数を変えないルールに従い、そのままにする。

### 完了の定義(PR#53)

1. 「PR#53〜55共通の完了の定義」1〜3を満たす。
2. 10本それぞれについて、PR#52の基準ファイルとの比較結果(一致/違いあり)を報告する。違いがあったものは、違いの具体的な中身(例: どのキーがどう変わったか、原子数・結合数の変化)と、「意図した違い」か「不具合」かの分類と根拠を書く。
3. `.brd`の基準ファイルを更新した場合は、更新前後の原子数・結合数を示す。

## PR#54: 構造の操作・表示(14本)

`brd-box.py`・`brd-density.py`・`brd-formula.py`・`brd-show-res.py`・`brd-show-bonds.py`・`brd-setup-bond.py`・`brd-select.py`・`brd-select-path.py`・`brd-divide.py`・`brd-divide-mainchain.py`・`brd-renumber-resid.py`・`brd-restructure.py`・`reorder.py`・`crystallize.py`

### 方針(Claude決定、2026-10-11)

1. **`AVOGADRO_CONST`**: `proteindf_bridge.rs`のモジュール定数として、純Python版(`proteindf_bridge/common.py`)と同じ値`6.02214076e+23`(CODATA 2018の定義値)を公開する。
2. **`AtomGroup.restructure(reference, range=1.0e-5)`**: Pythonバインディング(`rust/crates/proteindf-bridge-py`)に追加する。アルゴリズムは純Python版(`proteindf_bridge/atomgroup.py`の`restructure`と、そこから呼ばれる`_copy_attributes`・`_get_rest_of_frame_molecule`・`_assign_rest_molecule`)と同じにし、バインディングに既にある`Select_AtomGroup`・`Select_Atom`・`select`・`^`(`__xor__`)に相当するコアの機能を組み合わせて実装する。コアのクレートは変更しない(必要なら報告)。純Python版と同じ入力で結果(`.brd`の内容)が一致することをテストする。
3. **`brd-setup-bond.py`**: 純Python版は`AttributeError`で動かない(PR#52の記録)。Rust版では`AtomGroup.setup()`で結合を解決し`save_atomgroup`で保存する。基準ファイルは新しい出力で置き換えてよい(「不具合の修正」に分類)が、スナップショットとの一致だけでなく、`docs/cli-migration-notes.md` §3の5に書いた正しさの基準(1hlsでGLU 4の主鎖`C=O`がちょうど1本で次数2、鎖BのARG 22に次数2の結合が2本以上、ペプチド結合とCYS6–CYS11のジスルフィド結合があること、結合の重複がないこと)を、コマンドの出力`.brd`に対してアサートするテストを追加する。
4. 上の2・3以外のコマンドは、対応表のとおり既存のRust版APIで置き換える。

### 完了の定義(PR#54)

1. 「PR#53〜55共通の完了の定義」1〜3を満たす。
2. 14本それぞれについて、PR#52の基準ファイルとの比較結果(一致/違いあり)を報告する。違いがあったものは、違いの具体的な中身と、「意図した違い」か「不具合」かの分類と根拠(実際のコードの行や仕様の節)を書く。
3. `.brd`の基準ファイルを更新した場合は、更新前後の原子数・結合数を示す。

## PR#55: モデリング(2本)

`neutralize.py`・`superposer.py`

### 方針(Claude決定、2026-10-12)

1. 対応表のとおり、`Neutralize`・`Superposer`・`Superposer_quaternion`・`load_atomgroup`・`save_atomgroup`の既存のRust版APIで置き換える。
2. **`neutralize.py`**: 追加されるイオンの数・種類・位置が基準ファイルと違う場合は、純Python版とRust版のアルゴリズムの違い(実際のコードの行)を調べて分類する。「意図した違い」とする場合は、出力`.brd`について、全体の電荷が0になること(または純Python版と同じ目標電荷)、追加したイオンの数と元素、イオンと既存の原子・イオン同士が近づきすぎていないこと(最短距離)をアサートするテストを追加する。乱数を使う処理があれば、実行ごとに結果が変わらないかを確かめて報告する。
3. **`superposer.py`**: 数値(回転行列・RMSDなど)が浮動小数点の誤差の範囲で違う場合は、許容誤差を決めて比べ、その値と理由をテストに書く。許容誤差を広げて違いを隠さないこと。

### 完了の定義(PR#55)

1. 「PR#53〜55共通の完了の定義」1〜3を満たす。
2. 2本それぞれについて、PR#52の基準ファイルとの比較結果(一致/違いあり)を報告する。違いがあったものは、違いの具体的な中身と分類・根拠(実際のコードの行や仕様の節)を書く。
3. 26本すべての移行が終わったことを確かめる: `proteindf_bridge/cli/`の対象26本で、純Python版の`proteindf_bridge`(`proteindf_bridge.rs`以外)をimportしているものがないことを、`grep`の実際の出力で示す。

## PR#53〜55共通の完了の定義

1. 対象のコマンドが`proteindf_bridge.rs`だけを使う(純Python版の`AtomGroup`などをimportしない)。
2. PR#52の基準ファイルとの比較で違いが出たものを、すべて「意図した違い」か「不具合」に分類する。意図した違いは、根拠(`RUST_PORT_SPEC.md`の節や`docs/rust-port-handoff.md`の教訓)とともに`docs/cli-migration-notes.md`に記録し、基準ファイルを更新する。不具合は直す。**基準ファイルを、理由なく新しい出力で上書きしないこと。**
3. 全PR共通のルール5のテストが通る。

## PR#52 レビュー結果(1回目、2026-10-11、要修正)

`feature/cli-migration-pr52`(`14ccf6f`・`ed646ac`)をレビューした。26本すべてにテストがあり、コマンドの`main()`を直接呼んで、書き出したファイル(`.brd`は読み込んだ内容で)と標準出力を基準ファイルと比べる作りは妥当。正規化の方法(`.brd`はキーをソートし浮動小数点を6桁に丸めたJSON、テキストは改行と行末空白を揃える、`crystallize.py`の一時パスを置き換える)も妥当。純Python版のパッケージとRustのコードは変更されていない。Pythonテスト一式(224件)と`cargo nextest`(317件、2件スキップ)をClaudeが確認した。

### 修正依頼

1. **【不十分】基準ファイルの一部が空または意味のない出力になっており、移行後の確認に使えない。** 多くのテストの入力が、階層を持たない小さなキャッピングペプチド`ACE_ALA_NME.brd`(22原子、結合なし)であるため。Claudeが基準ファイルを確認したところ、次のとおりだった。
   - `brd2gro`: 原子0件のGROしか出ていない。
   - `brd-show-res`: `unknown bridge file format.`しか出ていない。
   - `brd-show-bonds`: 出力が空(入力に結合がない)。
   - `brd-renumber-resid`: 原子0件の`.brd`。
   - `neutralize`: 電荷を持たないペプチドなのでイオンが追加されず、入力と同じ22原子。
   
   各コマンドの主な処理が実際に働く入力に変える。具体的には、モデル・鎖・残基の階層を持つタンパク質の`.brd`(例: `1hls.pdb`から純Python版で作ったもの)を入力にし、`neutralize`は電荷を持つ残基のあるタンパク質で実行してイオンが追加されることを基準に含める。`brd-show-bonds`は結合を持つ入力にする。テストの中で純Python版を使って入力を作るか、作った入力をフィクスチャとして追加するかは任せる。**どのコマンドも、出力が空・0件・エラーメッセージだけの基準ファイルにしないこと。** 純Python版の既存の不具合で意味のある出力が得られない場合は、その旨と原因を`docs/cli-migration-notes.md`に書く(次の2を参照)。
2. **【記録】`brd-setup-bond.py`は純Python版で`AttributeError`になる(`bridge.Bond`が公開されていない)ため、移行前の基準として使える出力がない。** `docs/cli-migration-notes.md`に、PR#54で移行する際の正しさの基準(Rust側のテストの基準値、例えば`test_bond_resolution.rs`・`test_ccd_templates.rs`のどの値を使うか)を書いておく。
3. **【報告】`cargo nextest`の結果は、`Summary`の行を貼る**(1回目の報告では`Finished`の行しかなかった。Claudeが実行し、`317 tests run: 317 passed, 2 skipped`を確認した)。

### 完了の定義(修正後)

1. 上記1〜3に対応し、同じブランチに追加コミットする。
2. 修正後の基準ファイルについて、ファイルごとの大きさと、`.brd`の基準ファイルごとの原子数を一覧で示す。
3. 全PR共通のルール5のテストが通る。

## PR#52 レビュー結果(2回目、2026-10-11、要修正・文書のみ)

修正コミット`f75798b`を確認した。入力を階層を持つタンパク質(`1hls.pdb`・`2MGO.pdb`から作った`.brd`)に変え、`brd2gro`(782原子)・`brd-show-res`(51残基)・`brd-show-bonds`(2MGOのSSBOND)・`brd-renumber-resid`(782原子)・`neutralize`(782→792原子)など、どの基準ファイルも主な処理が働く中身になった(修正依頼1は対応済み)。

### 修正依頼

1. **【誤り】`docs/cli-migration-notes.md`の`brd-setup-bond.py`の正しさの基準(修正依頼2への対応)が、存在しないテストを引用している。** `tests/test_rs_bond_resolution.rs`の`test_resolve_bonds_1hls_real_pdb`、`tests/test_rs_phase2.rs`の`test_bond_setup_covalent_tolerance_and_exclusion`は、どちらもリポジトリに存在しない(Claudeが`grep`で確認)。また「GLU 4 の CD=OE1」は誤りで、実際のテストが確認しているのはGLU 4の主鎖`C=O`である。次の実在するテストと、その中で実際にアサートしている内容に書き直す。
   - `rust/crates/proteindf-bridge/tests/test_ccd_templates.rs`の`test_atomgroup_setup_1hls_real_pdb`(710行目〜): GLU 4の`C=O`結合がちょうど1本で次数2、鎖BのARG 22に次数2の結合が2本以上(`C=O`と`CZ=NH2`)、など。
   - 同ファイルの`test_resolve_bonds_preserves_existing_file_bonds`(805行目〜): ファイル由来の結合が保持され、重複しない。
   - `rust/crates/proteindf-bridge/tests/test_covalent_bond_detection.rs`の`test_non_covalent_contact_excluded`・`test_known_covalent_bonds_1hls`: 共有結合半径による判定で、非共有の接触を結合とみなさない。
   - **引用するテスト名・ファイル名・アサートの内容は、必ず実際のファイルを開いて確認してから書くこと。記憶や推測で書かないこと。**

### 完了の定義(修正後)

1. 上記1に対応し、同じブランチに追加コミットする。引用した各テストについて、ファイル名・行番号・該当するアサート行をそのまま報告に貼る。

## PR#52 完了(2026-10-11)

`60bbcef`で修正依頼に対応済み。引用された4つのテスト(`test_atomgroup_setup_1hls_real_pdb` 710行目、`test_resolve_bonds_preserves_existing_file_bonds` 805行目、`test_non_covalent_contact_excluded` 45行目、`test_known_covalent_bonds_1hls` 91行目)が実在し、アサートの内容が一致することをClaudeが確認した。Pythonテスト一式224件OK(新しい環境)。純Python版のパッケージとRustのコードに変更なし。`--no-ff`でdevelopにマージした(`1be62f7`)。

## PR#53 レビュー結果(1回目、2026-10-11、要修正)

`feature/cli-migration-pr53`(`c0a7927`・`8a7a53c`)をレビューした。10本すべてが`proteindf_bridge.rs`と`msgpack`だけをimportしている。基準ファイルとの違い3件のうち、`mmcif2txt`(純Python版トークナイザが引用符付きの値を分割していた不具合の修正)と`brd2gro`(ボックス行をGROMACSと同じ`%10.5f`にした)は意図した違いとして妥当。Pythonテスト一式227件、`cargo nextest`(`317 tests run: 317 passed, 2 skipped`)、clippy・fmtをClaudeが新しい環境で確認した。

### 修正依頼

1. **【不具合・テストが空振り】`tests/test_rs_phase6.py`の`test_simple_mmcif_data_property`(740行目)は、`ALA.cif`のパスが`rust/proteindf_bridge/data/ALA.cif`(存在しない)に解決されるため、`if os.path.exists(...)`の中のアサートが一度も実行されずに成功している。** パスを`proteindf_bridge/data/ALA.cif`に正しく直し、`if`でスキップせず、ファイルがなければ失敗するようにする。
2. **【テスト不足】`SimpleMmcif.from_data`のテストがない。また、`mmcif2mol2.py`の`-w`(スクラッチ書き出し)・`-r`(スクラッチ読み込み)の経路が基準テストで通っていない。** `data` → `msgpack` → `from_data`の往復で内容が一致するテストと、`mmcif2mol2.py`を`-w`で実行してから`-r`で実行し、mol2の出力が通常の実行と同じになるテストを追加する。不正な入力(辞書でない値、ブロックの形が違うもの)を`from_data`に渡したときに`BrError`の系統のエラーになることも確かめる。
3. **【書式の違いを記録していない】`AtomGroup.save_csv`の改行がLF(`\n`)だが、純Python版(`csv.writer`の既定)はCRLF(`\r\n`)。** 基準テストは改行を正規化し、`test_atomgroup_save_csv_parity`も行ごとに`strip()`するため、この違いが検出されていない。TASKの方針3(列・書式は純Python版と同じ)とRFC 4180に合わせてCRLFにし、パリティテストは改行を含めたバイト列で比較する。
4. **【書式の違い】`save_csv`の浮動小数点の書式がPythonの`repr`と一致しない場合がある。** Rustの`{}`は指数表記を使わないため、例えば電荷`1e-05`はPythonでは`1e-05`、Rust版では`0.00001`になる(`1e16`以上も同様)。Pythonの`repr(float)`と同じ規則(指数が-4未満または16以上で指数表記)にし、`1e-05`・`-2.5e-07`・`1e16`・`123.0`・`0.1`などの値でPython版と一致することをテストする。
5. **【APIの方針】`Pdb.debug`は何もしないプロパティ(getterは常に`False`、setterは値を捨てる)であり、公開APIに意味のない属性を足している。** `Pdb.debug`をバインディングから削除し、`pdb2brd.py`側で`pdb_obj.debug = debug`の行を削除する(`-d/--debug`オプションはルール1に従い受け付けたまま残す)。この扱いを`docs/cli-migration-notes.md`に書く。
6. **【根拠の誤り】`brd2pdb.py`のB-factorの違い(1.00→0.00)の根拠に挙げている`RUST_PORT_SPEC.md` §3.17は、mmCIF書き出しの仕様であり、PDB書き出しの根拠ではない。** 分類(意図した違い)は妥当なので、根拠を実際のコードに直す: 純Python版は読み込み時の既定値が`0.0`(`proteindf_bridge/biopdb.py`の216行目)なのに、`set_by_atomgroup`の経路でだけ`1.0`(同474行目の`setdefault("temp_factor", 1.0)`)になる不整合があり、Rust版は読み込み・書き出しとも`0.0`で一貫している(`rust/crates/proteindf-bridge/src/format/pdb.rs`)。行番号は実際のファイルで確認してから書く。

### 完了の定義(修正後)

1. 上記1〜6に対応し、同じブランチに追加コミットする。
2. 追加・変更したテストについて、ファイル名・テスト名・行番号を報告に書く。修正依頼1については、修正後のテストのアサートが実行されることの根拠(例: ファイルを一時的に別名にすると失敗すること)を示す。
3. 全PR共通のルール5のテストが通る。

## PR#53 完了(2026-10-11)

修正コミット(`8554dc8`・`d3c9590`・`3a046fc`・`3878ad3`・`649ce56`)で修正依頼1〜6に対応済み。`test_simple_mmcif_data_property`のアサートが実行されること、`from_data`の往復と不正入力、`mmcif2mol2`の`-w`/`-r`経路、`save_csv`のCRLFとPythonの`repr`と同じ浮動小数点の書式(バイト列で比較)をテストで確認した。引用された行番号(`biopdb.py` 216・474行目、`pdb.rs` 266・559行目)をClaudeが実物で確認した。Pythonテスト一式231件OK(新しい環境)、`cargo nextest` `317 tests run: 317 passed, 2 skipped`、clippy・fmt OK。コアのクレートと純Python版のパッケージ(`cli/`以外)に変更なし。`--no-ff`でdevelopにマージした(`9b7857f`)。

## PR#54 レビュー結果(1回目、2026-10-12、要修正)

`feature/cli-migration-pr54`(`d7797c7`・`3fb7806`・`8b42ce9`・`e32b1df`)をレビューした。14本のコマンドは`proteindf_bridge.rs`だけを使っている(他のimportは標準ライブラリのみ)。`AVOGADRO_CONST`は純Python版と同じ値。`AtomGroup.restructure`は純Python版(`atomgroup.py`の`restructure`・`_copy_attributes`・`_get_rest_of_frame_molecule`・`_assign_rest_molecule`)と同じ手順(`Select_Atom`の既定の距離0.1も同じ)で、基準ファイルと一致した。`brd-setup-bond.py`の出力は782原子・793結合で、正しさの基準(GLU 4の`C=O`、ARG 22、ペプチド結合、CYS6–CYS11、重複なし)をアサートするテストがある。Pythonテスト一式235件、`cargo nextest`(`317 tests run: 317 passed, 2 skipped`)、clippy・fmtをClaudeが新しい環境で確認した。コアのクレートと純Python版のパッケージ(`cli/`以外)に変更なし。

### 修正依頼

1. **【APIの方針】`AtomGroup.erase_atom`をバインディングに追加しているが、純Python版の`erase_atom`は非推奨(`proteindf_bridge/atomgroup.py`の`erase_atom`は`DeprecationWarning`を出し、`remove_atom`を使うよう案内している)。** 新しいAPIに非推奨の名前を持ち込まない。`erase_atom`をバインディングから削除し、`reorder.py`は`remove_atom`を使う。`test_atomgroup_remove_atom_and_erase_atom`も合わせて直す(`remove_atom`だけのテストにし、テスト名も変える)。
2. **【記録の誤り】`brd-formula.py`の違いの記録(`docs/cli-migration-notes.md`の表と報告)に誤りがある。**
   - 報告では純Python版の出力を`N379C254O76H67S6`と書いているが、実際の基準ファイルは`N67C254O76H379S6`。
   - Rust版の実装の場所は`rust/crates/proteindf-bridge/src/core/atom_group.rs`ではなく`rust/crates/proteindf-bridge/src/atom_group.rs`(547行目の`get_formula`)。
   - 「化学式表現として標準的な仕様」は誤り。化学式の標準的な書き方はHill方式(C、H、その後アルファベット順。1hlsなら`C254H379N67O76S6`)であり、Rust版の原子番号順はHill方式ではない。記録は「純Python版は出現順で入力によって変わる。Rust版は原子番号順で入力によらず決まる。どちらもHill方式ではない」と事実どおりに書き、分類は「意図した違い(決定的になった)」のままでよい。Hill方式にするかどうかは、コアの変更になるため本PRでは扱わない。
3. **【記録】`brd-setup-bond.py`の出力`.brd`は、実行ごとにバイト列が変わる。** Claudeが同じ入力で3回実行したところ、3つとも`md5sum`が異なり、違いは鎖レベル(`/groups/model_1/groups/A/bonds`(21本)・`.../B/bonds`(29本))の結合リストの並び順だけだった。基準テストで結合リストを並べ替えて比べる対応は、テストとしては妥当。ただし原因はコアのクレート側にあり、本PRでは直さない(ルール3)。`docs/cli-migration-notes.md`に「既知の問題: `setup()`後の鎖レベルの結合の並び順が実行ごとに変わる。コアで対応する(別タスク)」と書く。報告にある原因の推測(`spatial.rs`の`HashMap`)は、確かめていないなら「未確認」と書く。

### 完了の定義(修正後)

1. 上記1〜3に対応し、同じブランチに追加コミットする。
2. 全PR共通のルール5のテストが通る。

## PR#54 完了(2026-10-12)

修正コミット(`e72486d`・`886fa39`)で修正依頼1〜3に対応済み。`erase_atom`はバインディング・コマンド・テストから削除され、`reorder.py`は`remove_atom`を使う。`brd-formula.py`の記録(純Python版`N67C254O76H379S6`、Rust版`src/atom_group.rs`の`get_formula`、どちらもHill方式ではない)と、`brd-setup-bond.py`の結合の並び順の既知の問題(原因未確認、コアで別タスク)をClaudeが確認した。Pythonテスト一式235件OK(新しい環境)、`cargo nextest` `317 tests run: 317 passed, 2 skipped`、clippy・fmt OK。`--no-ff`でdevelopにマージした(`1e63ffc`)。マージ後、`tests/test_rs_phase6.py`に残っていた`erase_atom`に触れるコメントをClaudeが直した。

### 別タスクの候補(PR#55の後)

1. `AtomGroup.setup()`後の鎖レベルの結合の並び順を決定的にする(コアの修正)。
2. 化学式(`get_formula`)をHill方式にするかどうか(コアの変更、ユーザー判断待ち)。

## PR#55 レビュー結果(1回目、2026-10-12、要修正・文書のみ)

`feature/cli-migration-pr55`(`dbd716c`・`dc6991f`・`77aac4c`)をレビューした。`neutralize.py`は基準ファイルと一致し、782→792原子、Na+ 6個・Cl- 4個、全体の電荷0、イオンと原子・イオン同士の最短距離、2回の実行で出力が同一であることをアサートしている。`superposer.py`のRMSDの差(約6e-16)は浮動小数点の丸めの範囲で、許容誤差`1e-14`は妥当。26本すべてが`proteindf_bridge.rs`だけを使い、純Python版を使うのは対象外の6本だけであることを`git grep`の出力で確認した。

### 修正依頼

1. **【根拠の誤り】`superposer.py`のRMSDの違いの根拠(`docs/cli-migration-notes.md`と報告)が実物と違う。**
   - 「numpy SVD」「nalgebra Kabsch SVD」は誤り。コアのクレートは`nalgebra`に依存していない(`rust/crates/proteindf-bridge/Cargo.toml`)。純Python版も Rust版も、対称行列の固有値分解を使う同じ手順である(純Python版`proteindf_bridge/superposer.py`の`trr.eig()`、Rust版`rust/crates/proteindf-bridge/src/superposer.rs`の`trr_sym.eig()`)。違いは固有値分解の実装(純Python版はnumpy、Rust版は自前の`Matrix::eig`)による丸めの差である。
   - 引用している行番号(`superposer.py` 82〜88行目は`center1`のプロパティ、`superposer.rs` 24〜30行目は構造体のフィールド)は、該当する処理を指していない。実際の行を開いて直す。
2. **【分類の書き直し】`superposer.py`の標準出力(`print(ag)`)の書式の違いを「Rust版のフォーマット規格に従う変更」としているが、Rust版の`Display for AtomGroup`(`rust/crates/proteindf-bridge/src/atom_group.rs`)は入れ子の深さに応じた字下げをしないため、階層の見通しが悪くなっている。** Claudeが1hls(モデル→鎖→残基の3階層)で確かめたところ、`group key=model_1`・`group key=A`・`group key=1`の見出しがすべて同じ字下げで並んだ。純Python版の`_get_str`は階層ごとに字下げする。基準ファイルの更新はこのままでよい(Rust版の出力を標準にする方針)が、記録は「意図した違い」ではなく「既知の問題: Rust版の`Display`は入れ子の字下げをしない。コアで対応する(別タスク)」と事実どおりに書く。

### 完了の定義(修正後)

1. 上記1・2に対応し、同じブランチに追加コミットする。引用する行番号は実際のファイルで確かめ、該当行をそのまま報告に貼る。

## PR#55 完了(2026-10-12)

修正コミット`7b4ddac`で修正依頼1・2に対応済み(文書とテストのコメントのみ)。RMSDの違いの根拠(純Python版`superposer.py` 272行目`trr.eig()`→`matrix.py` 501行目`numpy.linalg.eigh`、Rust版`superposer.rs` 305行目`trr_sym.eig()`→`matrix.rs` 593行目のJacobi法)と、`Display`の字下げを「既知の問題」とした記録をClaudeが実物で確認した。`--no-ff`でdevelopにマージした(`fd78b84`)。

**これで構造を扱う26本すべての移行が完了した。**

### 次の作業(ユーザー決定、2026-10-12)

1. コアの修正(`docs/tasks/TASK_core-determinism-display.md`): `setup()`後の結合の並び順を決定的にする/`Display for AtomGroup`を階層ごとに字下げする。
2. リリース(2026.10.3)。
3. その後、化学式(`get_formula`)をHill方式にする。
