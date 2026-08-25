# TODO

`SPEC.md` の解析(§9 既知の不具合一覧)から起票した対応事項。優先度は影響範囲の見立てで、
実行テストによる再現確認はまだ行っていない。

## 優先度: 高(実行時に確実に例外で落ちる)

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
- [ ] `mail.py`: `smtp_password` が設定ファイルに平文保存される点を、
      運用ドキュメント([[pdf-dev-proteindf-bridge]])に注意書きとして残すか、
      keyring 等への移行を検討する。

## ドキュメント関連(`docs/TODO.md` から再掲・関連)

- [ ] GitHub Pages有効化、APIリファレンス日本語訳など(詳細は `docs/TODO.md` 参照)
