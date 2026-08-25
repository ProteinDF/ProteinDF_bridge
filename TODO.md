# TODO

`SPEC.md` の解析(§9 既知の不具合一覧)から起票した対応事項。優先度は影響範囲の見立てで、
実行テストによる再現確認はまだ行っていない。

## 優先度: 高(実行時に確実に例外で落ちる)

- [ ] `superposer_quaternion.py`: `Superposer_quaternion` に `superimpose()` を実装する。
      `scripts/superposer.py FILE1 FILE2 -q` が RMSD 表示直後に `AttributeError` で
      必ず落ちる(`Superposer.superimpose()` 相当の処理が丸ごと欠けている)。
- [ ] `superposer_quaternion.py:143-163`: `calc()` メソッドを正しく実装する(後で実装予定)。
      現状は未定義変数(`atom_group1`, `atom_group2`, `position1`, `position2`)と
      存在しないメソッド名(`self.match_positions` 等、正しくは `self._match_positions`)
      を参照しているため、呼び出し時の動作を整理して実装し直す。

## 優先度: 中(特定条件下で例外/誤動作)

- [ ] `dbmanager.py` `get_results()`(304-323行): ループのインデントを直し、
      複数行結果が全部返るようにする(現状は最後の1行のみ)。`data` が0件の場合の
      `NameError` も合わせて修正。
- [ ] `dbmanager.py` `create_table()`(90-93行): `field_names` が辞書(型指定あり)の
      場合の処理を修正する。`.items()` の呼び忘れと、フォーマット文字列
      (`'{name} {type}'.format(k, v)`)の名前/位置引数不一致を直す。
- [ ] `xyz.py` `Xyz.__init__`(47-48行): `Xyz("path/to/file.xyz")` のように文字列を
      渡すコンストラクタ経路が `NameError` になる(`self.load(file_path)` の
      `file_path` が未定義、正しくは `rhs`)。
- [ ] `matrix.py:534` `SymmetricMatrix.get_raw_data()`: `range(dim * (dim + 1) / 2)`
      を `range(dim * (dim + 1) // 2)` に修正する(Python3では真division結果の
      `float` を `range()` に渡せず `TypeError`)。
- [ ] `bond.py`: ファイル冒頭に `from .matrix import SymmetricMatrix` と
      `from .atomgroup import AtomGroup` を追加する(現状は他モジュール経由の
      import に依存した偶然の動作)。

## 優先度: 低(実害は小さいが直しておきたい)

- [ ] デバッグ用 `print()` の削除・`logging` への置き換え:
  - `atomgroup.py:551` `assign_charges()`
  - `amber_prmtop.py` の `_check_data()`, `_read_atom_name()`, `_read_charges()`,
    `_read_atomic_number()`
  - `superposer.py` `_get_rotation_matrix()` ほか(呼ぶたびに標準出力を汚染する)
  - `superposer_quaternion.py` ほぼ全メソッド
- [ ] `neutralize.py` `_exempt_list()`: `self._model` 未設定のまま参照している。
      呼び出し自体は現在コメントアウトされているデッドコードだが、将来
      再度有効化する場合は先に直す必要がある(「対イオン追加の除外リスト」機能が
      実質無効化されたままになっている点も要検討)。
- [x] `setup.cfg` と `proteindf_bridge/_version.py` のバージョン表記を `2026.8.0` に一致させた。
- [ ] `mail.py`: `smtp_password` が設定ファイルに平文保存される点を、
      運用ドキュメント([[pdf-dev-proteindf-bridge]])に注意書きとして残すか、
      keyring 等への移行を検討する。

## ドキュメント関連(`docs/TODO.md` から再掲・関連)

- [ ] GitHub Pages有効化、APIリファレンス日本語訳など(詳細は `docs/TODO.md` 参照)
