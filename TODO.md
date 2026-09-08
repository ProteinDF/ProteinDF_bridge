# TODO

`SPEC.md` の解析(§9 既知の不具合一覧)から起票した対応事項。

**2026-09-08 追記**: `SPEC.md` を実ソースと再照合する監査を行い、下記の「優先度:
高・中・低」の項目(mail.py の注記を除く)はすべて実装側で修正済みであることを
ソースコード読解で確認した。あわせて `SPEC.md` 側にもこれらの修正が反映されていない
記述(既に直った不具合が「既知の不具合」として残っていた)が見つかったため、
`SPEC.md` を実装に合わせて修正した。その過程で新たに1件、未修正の実バグ
(`functions.py` `save_yaml()`)を発見したので下に追加した。

その後、壊れていた `.venv`(`pyyaml`/`msgpack`/`numpy`/`pip` すら未導入)を
`uv pip install --python .venv/bin/python -e .` で復旧し、
`python -m unittest discover -s tests -v` を実際に実行して確認した:
86件中84件がpass。上記「優先度:高・中・低」の修正済み項目に対応するテストは
すべてpassし、実行レベルでも修正が有効であることを確認できた。残る2件は
今回のSPEC.md監査とは無関係な、以前から存在していたdoctestの不具合
(下記「新規発見」参照)。

## 優先度: 高(実行時に確実に例外で落ちる)

- [x] `functions.py:62-68` `save_yaml()`: `get_yaml()` が返す `str` を
      `open(yaml_path, "wb")`(バイナリモード)へ書き込んでおり、呼び出すと必ず
      `TypeError: a bytes-like object is required, not 'str'` になっていた。
      `open(path, "w")` に修正し、回帰テスト(`tests/test_functions.py`)を追加した。
- [x] `superposer_quaternion.py`: `Superposer_quaternion` に `superimpose()` を実装した。
      `scripts/superposer.py FILE1 FILE2 -q` が正常に動作するようになった。
- [x] `superposer_quaternion.py:143-163`: `calc()` メソッドを修正し、
      一括計算を実行して `rmsd` を返す正常なメソッドとして実装した。

## 優先度: 中(特定条件下で例外/誤動作)

- [x] `dbmanager.py` `get_results()`(304-323行): ループのインデントを直し、
      複数行結果が全部返るように修正した。`data` が0件の場合の挙動も修正。
- [x] `dbmanager.py` `create_table()`(90-93行): `field_names` が辞書の場合に
      `.items()` を使用し、フォーマットを修正した。
- [x] `xyz.py` `Xyz.__init__`(47-48行): 文字列パス指定時の引数参照(`rhs`)を修正した。
- [x] `matrix.py:534` `SymmetricMatrix.get_raw_data()`: `range(dim * (dim + 1) // 2)` に修正した。
- [x] `bond.py`: ファイル冒頭に `from .matrix import SymmetricMatrix` と
      `from .atomgroup import AtomGroup` を追加した。

## 優先度: 低(実害は小さいが直しておきたい)

- [x] デバッグ用 `print()` の削除・`logging` への置き換え:
  - `atomgroup.py:551` `assign_charges()`
  - `amber_prmtop.py` の `_check_data()`, `_read_atom_name()`, `_read_charges()`,
    `_read_atomic_number()`
  - `superposer.py` `_get_rotation_matrix()` ほか
  - `superposer_quaternion.py`
- [x] `neutralize.py` `_exempt_list()`: `model` を引数で受け取り `IonPair(model)` を呼ぶように修正。
- [x] `setup.cfg` と `proteindf_bridge/_version.py` のバージョン表記を `2026.8.0` に一致させた。
- [x] `SPEC.md`: 上記の修正が反映されず「既知の不具合」として残っていた10件の記述
      (matrix/bond/xyz/atomgroup/amber_prmtop/neutralize/superposer/superposer_quaternion×2/
      dbmanager×2/バージョン乖離)を実装と付き合わせて修正・削除した。あわせて
      CLIスクリプト節の記載漏れ(`relax_protein.py`)や誤り(`brd-select.py`と
      `brd-select-path.py`を同一視、`load_atomgroup`/`save_atomgroup`を「一貫して使用」
      としていた記述)も修正した。
- [ ] `mail.py`: `smtp_password` が設定ファイルに平文保存される点を、
      運用ドキュメント([[pdf-dev-proteindf-bridge]])に注意書きとして残すか、
      keyring 等への移行を検討する。

## 環境関連

- [x] `.venv` に `pyyaml`/`msgpack`/`numpy` 等の依存パッケージ(`pip`自体も)が
      インストールされておらず、`proteindf_bridge` の import 自体ができない状態
      だった。`uv pip install --python .venv/bin/python -e .`(`setup.cfg` の
      `install_requires` を使用)で復旧し、`python -m unittest discover -s tests`
      が実行できる状態にした。

## 新規発見(2026-09-08、テスト実行時に判明。SPEC.md §9とは無関係)

- [x] `position.py` のdoctest(35-91行付近)が2件失敗していたのを修正した:
      - `p.norm()`(戻り値 `self`)を docstring 側で `>>> _ = p.norm()` と受けて
        `repr` が表示されないようにした。
      - `dot()` が `numpy.dot()` の戻り値(`numpy.float64`)をそのまま返しており
        numpy 2.x の repr 変更で doctest が壊れていたため、`vector.py` の同種
        メソッドに合わせて `float(...)` で包んで plain `float` を返すようにした。
- [ ] `ssbond.py` のdoctest(13-20行付近)が失敗する: `Pdb('./data/1hls.pdb')` が
      テスト実行時のカレントディレクトリに依存しており、`FileNotFoundError` になる
      (doctestが相対パスに依存していて自己完結していない)。

## ドキュメント関連(`docs/TODO.md` から再掲・関連)

- [ ] GitHub Pages有効化、APIリファレンス日本語訳など(詳細は `docs/TODO.md` 参照)
