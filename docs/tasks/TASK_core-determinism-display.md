# TASK: コアの修正 — 結合の並び順の決定性と、AtomGroupの表示の字下げ (PR#56)

> CLI移行(`docs/tasks/TASK_cli-rust-migration.md`、PR#52〜55)のレビューで見つかった、コアのクレート(`rust/crates/proteindf-bridge`)の2つの問題を直す。ユーザー確認済み(2026-10-12)。このあとリリース(2026.10.3)を行い、化学式のHill方式はその後の別タスクとする。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaudeが担当する。
- `develop`から`feature/core-determinism-display`を切って作業する。**`develop`へは自分でマージしない。** タグやReleaseは作らない。
- 完了報告は、起動時に指示されたパス(リポジトリの外)に、決まった見出しの短いファイルとして書く。コマンドの出力は要約・再構成せず、実際の出力行をそのまま貼る。
- テスト名・ファイル名・行番号を報告や文書に書くときは、必ず実際のファイルを開いて確認してから書く。確かめていない推測を事実として書かない。

## 対象

### 1. `AtomGroup.setup()`後の結合の並び順を決定的にする

- **現象**: `brd-setup-bond.py`を同じ入力(1hls)で実行するたびに、出力`.brd`のバイト列が変わる。違いは鎖レベル(`/groups/model_1/groups/A/bonds`・`.../B/bonds`)の結合リストの並び順だけ(Claude確認、2026-10-12)。
- **原因(Claude確認)**: `src/bond.rs`の`Bond::setup_heuristic`が、`src/spatial.rs`の`CellList::for_each_neighbor_pair`のコールバックの順に`bonds_to_add`へ積み、その順に`mol.add_bond`している。`CellList`はセルを`HashMap<(i32, i32, i32), Vec<usize>>`(`spatial.rs`の`grid`)で持つため、走査順がハッシュのシードに依存する。
- **直し方**: `setup_heuristic`で`mol.add_bond`する前に、`bonds_to_add`を原子の添字`(p, q)`の順に並べ替える。`CellList`自体の走査順を変える必要はない(ほかの利用者はテストだけ)。
- **監査**: コアのクレートで、`HashMap`・`HashSet`の走査順が出力(`.brd`・PDB・mmCIFなどのファイル、Pythonに返すリストや辞書、標準出力)の順序に影響する箇所が他にないかを調べ、あれば同じく決定的にする。調べた箇所と結論(影響あり/なし、その理由)を報告に一覧で書く。
- **テスト**:
  - Rust: 1hlsで`setup()`を同じプロセス内で2回以上(別々の`AtomGroup`に対して)実行し、`get_bond_list()`の結果が順序も含めて完全に一致すること。`HashMap`のシードはインスタンスごとに変わるので、プロセス内の繰り返しで違いを検出できることを、修正前のコードで失敗することによって確かめて報告する。
  - Python: `brd-setup-bond.py`を2回実行した出力`.brd`のバイト列が一致すること(`tests/test_cli_golden.py`)。基準テストの結合リストの並べ替え(`_canonicalize_brd`)は、必要がなくなれば外して、並び順も含めて基準ファイルと比べる。外せない場合は理由を書く。

### 2. `Display for AtomGroup`を階層ごとに字下げする

- **現象**: `src/atom_group.rs`の`impl fmt::Display for AtomGroup`は入れ子の深さを持たないため、`print(ag)`でモデル・鎖・残基の見出しがすべて同じ字下げで並ぶ。また原子の行の先頭に原子のキーが重複して付く(行末にパスがあるため)。
- **直し方**: 純Python版の`proteindf_bridge/atomgroup.py`の`_get_str`と同じ構造にする。
  - 深さ`d`の字下げは`"  " * d`。
  - グループの見出しは`{字下げ}# group key={key} name={name}`(最上位のkeyは空)。純Python版の` parent=`は出さなくてよい(Rust版の`AtomGroup`は親を持たないため)。
  - 子グループを再帰的に(深さ+1で)出し、続いて原子を`{字下げ}{atom}`(`Display for Atom`をそのまま使う。原子のキーを前に付けない)、続いて結合を`{字下げ}bond {path1} <-{order}-> {path2}`で出す。
  - 出力順は、これまでどおり`groups()`・`atoms()`の順(既に決定的)に従う。
- Python側の`str(ag)`(`rust/crates/proteindf-bridge-py/src/atom_group.rs`の`__str__`)はこの`Display`を使っているので、そのまま反映される。
- **テスト**: 3階層(モデル→鎖→残基)の小さな構造で、出力の全文を期待値と比べるRustのテストを追加する。結合を持つ構造で`bond`行が出ることも確かめる。

### 3. 基準ファイルの更新

1と2により、`tests/cli_golden/`のうち`brd_setup_bond`と`superposer.stdout.txt`(ほか、`print(ag)`を使うコマンドがあればそれも)が変わる。変わったファイルごとに、違いの中身を確かめ、`docs/cli-migration-notes.md`の「既知の問題」の記録(PR#54の`brd-setup-bond.py`、PR#55の`superposer.py`)を「解決済み(本タスク)」に更新する。**基準ファイルは、上記1・2による違いだけであることを確かめてから更新する。**

## 守ること

- `RUST_PORT_SPEC.md`・`docs/tasks/`の過去の記録は書き換えない(本ファイルへの追記もしない)。
- 純Python版のパッケージは変更しない。
- 意味のまとまりごとにコミットする。

## 完了の定義

1. 対象1〜3に対応し、テストを追加している。
2. 次のテストが通る: 新しい一時環境にパッケージを`pip install`して`python -P -m unittest discover -s tests`、`cargo nextest run --workspace --status-level fail --final-status-level fail`、`cargo test --doc`、`cargo clippy --workspace --all-targets -- -D warnings`、`cargo fmt --all -- --check`。
3. `brd-setup-bond.py`を3回実行し、出力`.brd`の`md5sum`がすべて同じであることを、実際の出力で示す。
