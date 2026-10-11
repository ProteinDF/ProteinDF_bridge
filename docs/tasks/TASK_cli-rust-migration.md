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

## PR#54: 構造の操作・表示(14本)

`brd-box.py`・`brd-density.py`・`brd-formula.py`・`brd-show-res.py`・`brd-show-bonds.py`・`brd-setup-bond.py`・`brd-select.py`・`brd-select-path.py`・`brd-divide.py`・`brd-divide-mainchain.py`・`brd-renumber-resid.py`・`brd-restructure.py`・`reorder.py`・`crystallize.py`

## PR#55: モデリング(2本)

`neutralize.py`・`superposer.py`

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
