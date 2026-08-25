# ProteinDF_bridge 仕様書

> 本ドキュメントは `proteindf_bridge/` 以下のソースコード(`proteindf_bridge/*.py`
> 全33ファイル、`scripts/*.py` 全CLI)を解析して作成した内部仕様書です。

作成状況(全項目完了):

- [x] 1. 概要・全体アーキテクチャ
- [x] 2. コア基盤モジュール (error, common, path, str_processing, functions, utils, `__init__`)
- [x] 3. データモデル (periodictable, vector, matrix, position, atom, bond, atomgroup)
- [x] 4. フォーマット変換 (format, xyz, gro, mol2, mmcif, amber_prmtop, biopdb)
- [x] 5. 構造操作 (select, modeling, superposer, superposer_quaternion, neutralize, ionpair, ssbond, aminoacid)
- [x] 6. その他インフラ (dbmanager, mail)
- [x] 7. CLIスクリプト仕様 (`scripts/*.py`)

---

## 1. 概要

`ProteinDF_bridge` は、量子化学計算エンジン [ProteinDF](https://proteindf.github.io/) と、
PDB / GROMACS / MOL2 / mmCIF / Amber prmtop / XYZ といった一般的な分子構造フォーマットとを
橋渡しする Python ライブラリである。中心となるのは独自の分子表現オブジェクトモデル
(`AtomGroup` / `Atom` / `Bond`) で、これをコンパクトなバイナリ形式 (msgpack、通称 **brd 形式**)
または YAML との間でシリアライズできる。ライブラリ本体に加えて、`scripts/` 配下に
30 個超の単機能 CLI ツール (`brd-*`, `*2brd`, `brd2*` など) が提供される。

パッケージ名: `proteindf_bridge` / バージョン: `2024.3.0`
(`proteindf_bridge/_version.py` の `__version__` は `2022.2.5` のままで、`setup.cfg` の
`version = 2024.3.0` と乖離している — 要確認)

## 2. 全体アーキテクチャ

```
proteindf_bridge/
├── コア基盤       error, common, path, str_processing, functions, utils
├── 数理データ型    vector, matrix, periodictable, position
├── 分子データモデル atom, bond, atomgroup, aminoacid
├── フォーマットI/O format, xyz, gro, mol2, mmcif, amber_prmtop, biopdb
├── 構造操作       select, modeling, superposer, superposer_quaternion,
│                  neutralize, ionpair, ssbond
└── インフラ       dbmanager, mail
```

典型的なデータフローは次の通り:

```
外部フォーマット (PDB/gro/mol2/mmCIF/prmtop/XYZ)
        │  *2brd 系スクリプト / biopdb, gro, mol2, mmcif, amber_prmtop モジュール
        ▼
   AtomGroup オブジェクトツリー (階層構造: モデル→チェイン→残基→原子)
        │  brd-* 系スクリプト (select, divide, renumber-resid, setup-bond, ...)
        ▼
   加工済み AtomGroup
        │  brd2* 系スクリプト / functions.save_atomgroup (msgpack) or YAML
        ▼
外部フォーマット (PDB/gro/XYZ) または .brd ファイル (ProteinDF本体の入力)
```

`AtomGroup` はキー(文字列)で子要素(`AtomGroup` または `Atom`)を保持する再帰的な
辞書ライクコンテナであり、パス文字列 (`/chain_id/res_id/atom_name` 形式、`path.py` 参照) で
任意の深さの要素にアクセスできる。

## 3. コア基盤モジュール

### 3.1 `error.py` — 例外クラス

| クラス | 継承 | 用途 |
| --- | --- | --- |
| `BrError` | `Exception` | 本パッケージの基底例外。`errmsg` 属性を持つ |
| `BrInputError` | `BrError` | 入力値のエラー。`__init__(expr, msg)` |
| `BrValueError` | `BrError` | 値のエラー。`__init__(expr, msg)` |

パッケージ全体でこの3クラスのみが例外として使われる(型ヒエラルキーは浅い)。

### 3.2 `common.py` — 定数

- `AVOGADRO_CONST = 6.02214076e+23`

### 3.3 `path.py` — パス文字列ユーティリティ

`AtomGroup` のパス表記 (`/chain/res/atom` 形式) を扱う静的メソッド群。

| メソッド | 説明 |
| --- | --- |
| `Path.split_path(path)` | `/` 区切りでパスをリストに分割(先頭の空要素は除去) |
| `Path.get_chain_id(path)` | パスの第1要素(チェインID)を取得 |
| `Path.get_res_id(path)` | パスの第2要素(残基ID)を取得 |

例: `Path.get_chain_id("/B/17/561_CA")` → `'B'`

### 3.4 `str_processing.py` — 文字列/エンコーディングユーティリティ (`StrUtils`)

Python 2/3 両対応(`unicode`/`str`/`bytes` の互換定義)を主目的とした文字列処理群。

| メソッド | 説明 |
| --- | --- |
| `sort_nicely(list)` | 自然順ソート(数字部分を数値として比較) |
| `add_spaces(s, n)` / `del_spaces(s, n)` / `num_spaces(s)` / `unindent_block(s)` | 行頭インデントの操作 |
| `get_common_str(s1, s2)` | 2文字列の先頭からの共通部分を返す |
| `to_unicode(s)` / `to_unicode_dict(d)` / `to_unicode_list(l)` | bytes を再帰的に UTF-8 文字列へ変換(msgpack読み込み後の後処理に使用) |
| `to_bytes(s)` | 文字列を UTF-8 bytes へ変換 |
| `str_to_bool(s)` | 文字列/数値を真偽値に変換( `"Y"`, `"T"` で始まる、または非0の数値で `True`) |
| `check_pickled(data)` | データが pickle 可能かを再帰的に検証(デバッグ用) |

### 3.5 `functions.py` — YAML / MsgPack I/O ヘルパー

パッケージ直下に公開される、ファイル入出力のためのモジュールレベル関数群。

| 関数 | 説明 |
| --- | --- |
| `locate()` | 呼び出し元の (ファイル名, 関数名, 行番号) を返すデバッグ用関数 |
| `load_yaml(path)` / `parse_yaml(text)` | YAML ファイル/文字列を読み込み、`yaml.safe_load_all` の結果をリストで返す(複数ドキュメント対応) |
| `get_yaml(data)` / `save_yaml(data, path)` | Python オブジェクトを YAML 文字列化/保存 |
| `load_msgpack(path)` / `save_msgpack(data, path)` | msgpack ファイルの読み書き。読み込み後は `StrUtils` で再帰的に unicode 化 |
| `load_atomgroup(path)` | msgpack (`.brd`) ファイルを読み込み `AtomGroup` を構築 |
| `save_atomgroup(atomgroup, path)` | `AtomGroup.get_raw_data()` を msgpack として保存 |

`.brd` ファイルは実体として **msgpack でシリアライズされた `AtomGroup` の内部データ構造**
であり、`load_atomgroup`/`save_atomgroup` がその読み書きの正式な入り口となる。

### 3.6 `utils.py` — AtomGroup 操作ユーティリティ (`Utils`)

| メソッド | 説明 |
| --- | --- |
| `Utils.remove_WAT(atomgroup)` | 残基名が `HOH`/`WAT` のグループを再帰的に除去した新しい `AtomGroup` を返す(非破壊) |
| `Utils.get_sequential_residue_id(protein, chain_id, res_id)` | チェイン・残基IDから、タンパク質全体を通した通し番号(1始まり)を算出。該当なしなら `0` |

### 3.7 `__init__.py` — パッケージ公開API

`proteindf_bridge` パッケージの `__init__.py` は主要クラス/関数を `from .module import *`
形式でパッケージ直下に再エクスポートしている。主な公開名:

- 例外: `BrError`, `BrInputError`, `BrValueError`
- 数理型: `Vector`, `Matrix`, `SymmetricMatrix`, `identity_matrix`, `Position`, `PeriodicTable`
- 分子データモデル: `Atom`, `AtomGroup`, `AminoAcid`
- ユーティリティ: `Path`, `Utils`, `Select` 系一式(`select.py` の内容を丸ごと import)
- フォーマットI/O: `Pdb` (`biopdb.py`), `SimpleMmcif`, `SimpleMol2`, `SimpleGro`, `Xyz`, `AmberPrmtop`
- 構造操作: `Modeling`, `IonPair`, `SSBond`, `Neutralize`, `Superposer`, `Superposer_quaternion`
- インフラ: `DbManager`, `Mail`

`__all__` はコメントアウトされており未使用 — `import *` 時の公開範囲はモジュール側の
実装に委ねられている(明示的なエクスポート制御はない)。

## 4. データモデル

### 4.1 `periodictable.py` — `PeriodicTable`

原子番号1〜118の元素記号・原子量・van der Waals半径を保持する静的テーブルクラス
(インスタンス化せず `@staticmethod` として利用する)。

| メソッド | 説明 |
| --- | --- |
| `get_num_of_atoms()` | テーブルに含まれる元素数(ダミー原子 `X` を含む) |
| `get_symbol(atomic_number)` | 原子番号→元素記号 |
| `get_atomic_number(symbol)` | 元素記号(大文字小文字は正規化)→原子番号 |
| `atomic_weight(atom)` | 原子番号 or 記号→原子量 |
| `vdw(atom)` | 原子番号 or 記号→van der Waals半径(Å)。未定義元素は `0.0` |
| `__contains__(symbol)` | インスタンスに対する `in` 演算子(テーブル存在チェック) |

原子番号 `0` は `'X'`(不明/ダミー原子)。`Atom` のデフォルトはこの `X`。

### 4.2 `vector.py` — `Vector`

`numpy.ndarray` をラップした1次元ベクトル型。`Matrix`/`SymmetricMatrix` の行・列ベクトル、
`Position` の内部計算などで使われる。

| 種別 | メンバ |
| --- | --- |
| 生成 | `Vector(size)` / `Vector(list または ndarray)` / `Vector(Vector)` |
| プロパティ | `max`, `min`, `data` (numpy配列参照) |
| 操作 | `size()`, `resize(n)`, `get(i)`/`set(i,v)`, `to_list()`, `abs()`, `argsort()`, `flip()` |
| バッファI/O | `get_buffer()`/`set_buffer()` (生バイト列、`get_ndarray()`はディープコピー) |
| 演算子 | `+`, `-`, `*`(スカラー or 内積), 単項 `-`, `len()`, `[]` |

### 4.3 `matrix.py` — `Matrix` / `SymmetricMatrix` / `identity_matrix`

`numpy.ndarray` ベースの一般行列 (`Matrix`, `type="GE"`) と対称行列
(`SymmetricMatrix`, `type="SY"`, `Matrix` を継承)。対称行列は下三角のみを実データとして
持ち、`get`/`set` 時に `row < col` なら自動的に転置して下三角側にアクセスする。

| メソッド (`Matrix`) | 説明 |
| --- | --- |
| `rows`/`cols`/`type`/`data` | 形状・種別・numpy配列への参照(プロパティ) |
| `get(r,c)`/`set(r,c,v)`/`add(r,c,v)` | 要素アクセス |
| `resize(rows,cols)` | サイズ変更(既存データは左上から保持) |
| `transpose()` | 転置(破壊的) |
| `select(r0,c0,r1,c1)` | 部分行列の抽出 |
| `get_row_vector(r)`/`get_col_vector(c)` | 行/列を `Vector` として取得 |
| `inverse()`/`pseudo_inverse()` | 逆行列/擬似逆行列 (numpy.linalg) |
| `get_symmetric_matrix()` | 対称性を検証しつつ `SymmetricMatrix` に変換 |
| `get_raw_data()` | `{"row","col","type":"GE","data":[...]}` (行優先の全要素配列) にシリアライズ |
| 演算子 | `+`, `-`, `*` (スカラー/行列/`Vector`との積), `==` |

| メソッド (`SymmetricMatrix` 追加分) | 説明 |
| --- | --- |
| `dim` | 次元数 (`rows == cols` が前提) |
| `get_general_matrix()` | 上三角も埋めた `Matrix` に変換 |
| `eig()` | `numpy.linalg.eigh` による固有値 `Vector` / 固有ベクトル `Matrix` (列ベクトルとして) |
| `get_raw_data()` | `{"row","col","type":"SP","data":[...]}` (上三角パック形式) |

`identity_matrix(dim)` はモジュールレベル関数で、`dim×dim` の単位行列 (`SymmetricMatrix`) を返す。

> **既知の不具合**: `SymmetricMatrix.get_raw_data()`(matrix.py:534)は
> `range(dim * (dim + 1) / 2)` を呼んでいるが、Python 3 では `/` が真division となり
> `float` を `range()` に渡すことになるため **`TypeError` で例外になる**(`//` にすべき)。
> `.brd`/msgpack への保存経路で対称行列を直接シリアライズする場合に影響する可能性がある。

### 4.4 `position.py` — `Position`

3次元座標 (x, y, z) を表すクラス。`list`/`tuple`/`numpy.ndarray`/`"x, y, z"` 形式の文字列/
他の `Position` から構築できる。

| 分類 | メンバ |
| --- | --- |
| プロパティ | `x`, `y`, `z`, `xyz` (リスト参照) |
| 距離 | `distance_from(other)`, `square_distance_from(other)` |
| 演算 | `+`, `-`, `*` (内積 or スカラー倍), `/`, 単項 `-`, `dot()`, `cross()`, `norm()`(正規化・破壊的) |
| 変換 | `move_to(pos)`, `rotate(3x3 Matrix)` |
| 比較 | `__eq__` は `distance_from < epsilon (1.0E-5)` で判定(厳密な浮動小数点一致ではない) |
| シリアライズ | `get_raw_data()` → `[x, y, z]` のリスト、`__getstate__`/`__setstate__` で pickle 対応 |

### 4.5 `atom.py` — `Atom`

分子構造の最小単位。原子番号・座標・力・電荷・名前・ラベル・パス・親への参照を持つ。

| 分類 | メンバ |
| --- | --- |
| 生成 | `Atom()` / `Atom(Atom)` (コピー) / `Atom(dict)` (`set_by_raw_data`) / `Atom(str)` (元素記号として) / `Atom(int)` (原子番号として)。キーワード引数 `symbol`, `position`/`xyz`, `force`, `name`, `label`, `charge`, `path`, `parent` にも対応 |
| プロパティ | `atomic_number`, `symbol` (相互変換), `xyz`/`position`, `force`, `name`, `label`, `charge`, `path`, `is_real` (原子番号>0), `vdw`, `weight()` |
| 移動 | `move_to()`, `shift_by()`, `rotate()`, `*=` |
| 比較 | `__eq__`: 原子番号と座標が一致すれば等しい(電荷は比較対象外、コメントアウト済み) |
| シリアライズ | `get_raw_data()` → `{"Z","name","Q","xyz","force"}` / `set_by_raw_data(dict)` で復元。未知キーは `logger.debug` で無視 |

`parent` は `AtomGroup` を指すが、循環import回避のためコンストラクタ内では型チェック
(`assert isinstance(..., AtomGroup)`)がコメントアウトされている。

### 4.6 `bond.py` — `Bond`

`AtomGroup` 内の全原子ペアについて、原子間距離と van der Waals 半径の和 (+0.4Å の許容値)
から結合の有無を推定し、`AtomGroup.add_bond()` で結合情報を追加するヘルパークラス。
`brd-setup-bond.py` (結合の自動推定)の実体。

| メソッド | 説明 |
| --- | --- |
| `setup(mol)` | `mol` (`AtomGroup`) 内の全原子を対象に距離行列・結合行列を作成し、`mol.add_bond()` を呼び出す |
| `_list_atoms(mol)` | 再帰的に全原子を `self._atoms` に収集 |
| `_make_distance_matrix()` | 原子間距離の `SymmetricMatrix` を作成 |
| `_make_bond_matrix()` | `distance <= vdw(p) + vdw(q) + 0.4` なら結合ありとする `SymmetricMatrix` (0/1) を作成 |

> **既知の不具合**: `bond.py` は `SymmetricMatrix` と `AtomGroup` を使用しているが、
> ファイル冒頭で **どちらも import していない**。`proteindf_bridge/__init__.py` で
> `from .matrix import ...` や `from .atomgroup import ...` が先に実行され、それらの名前が
> たまたま同一プロセスのグローバル名前空間に存在する場合のみ動作してしまう可能性があり、
> `bond.py` を単体で `import` した場合は `NameError` になる。

### 4.7 `atomgroup.py` — `AtomGroup`

パッケージの中核クラス。`OrderedDict` により子要素(サブ `AtomGroup` と `Atom`)を
キー付きで保持する再帰的コンテナで、モデル→チェイン→残基→原子のような階層構造を
そのまま表現できる。

**内部構造**

- `self._groups`: `{key: AtomGroup}` の `OrderedDict`
- `self._atoms`: `{key: Atom}` の `OrderedDict`
- `self._bonds`: `(atom1_relpath, atom2_relpath, order)` タプルのリスト(このグループを
  基準とした相対パスで保持される)
- `self._path`: このグループの絶対パス(`/` 終端)。子を追加するたびに `_update_path()` で
  子孫のパスを再計算する
- `self._parent`: 親 `AtomGroup` への参照
- `self._sort_atoms` / `self._sort_groups`: `"nice"` を指定すると `atoms()`/`groups()` の
  イテレーション順を `StrUtils.sort_nicely()`(自然順)でソートする

**パスアクセス**

`set_atom("/group1/subgroup1/C1", atom)` のように `/` 区切りのパスを直接キーとして渡すと、
存在しない中間グループを自動生成しながら再帰的に配置する(`set_atom` 内の分割ロジック)。
`__getitem__`/`get_group`/`get_atom` はキー一致 → 名前(`name`属性)一致の順で検索する。

**主なAPI一覧**

| カテゴリ | メソッド |
| --- | --- |
| 集計 | `get_number_of_groups()`, `get_number_of_atoms()`, `get_number_of_all_atoms()` (再帰), `sum_of_atomic_number()`, `get_atom_kinds()`, `get_atom_kinds_count()`, `get_formula()` (Hill式風の組成式), `formula()` (別実装、`get_atom_list()`ベース), `get_number_of_bonds()` |
| グループ操作 | `groups()`, `get_group()`, `set_group()`, `has_group[key\|name]()`, `remove_group()` (`erase_group()`は非推奨エイリアス), `get_group_list()` |
| 原子操作 | `atoms()`, `get_atom()`, `set_atom()`, `has_atom[key\|name]()`, `remove_atom()` (`erase_atom()`は非推奨), `pickup_atoms()`, `get_atom_list()`, `get_path_list()` |
| 電荷 | `charge` (全原子の`charge`合計), `real_nuclei_charge` (実原子のみの核電荷和), `nuclei_charge` (ダミー原子は`charge`で代替した核電荷和), `assign_charges(Vector)` (深さ優先順に電荷を割り当て) |
| 幾何 | `center()` (原子数重み付き重心), `box()` (バウンディングボックス min/max の `Position` タプル), `shift_by()`, `rotate()`, `*=` (スカラー倍、座標のみ) |
| 選択/再構成 | `select(Select)` (`select.py`のSelectorで部分木を抽出), `restructure(reference)` (フラットな座標リストを参照構造にはめ込み直す。PDB由来の座標をbrd構造に流し込む用途), `get_family(path)` (パスから共通祖先を辿って要素を検索) |
| 結合 | `add_bond(atom1, atom2, order)`, `get_bond_list()` (全結合を絶対パスのタプルリストで取得) |
| 集合演算 | `&`/`__and__` (共通部分), `\|`/`__or__` (`merge()`によるマージ), `^`/`__xor__` (対称差)。それぞれ `&=`, `\|=`, `^=` のインプレース版あり |
| フォーマット出力 | `get_xyz()` (XYZ形式文字列), `__str__()` (デバッグ用ツリー表示), `save_csv()`/`_get_csv_list()` |
| シリアライズ | `get_raw_data()` → `{"groups":{...}, "atoms":{...}, "name":..., "bonds":[...], "sort_atoms":..., "sort_groups":...}` / `set_by_dict_data(dict)` で復元。`__getstate__`/`__setstate__` で pickle・msgpack 両対応 |

`_bonds` は「このグループを基準とした相対パス」で正規化される
(`_add_bond_normalize` → `get_family(common_path)` で結合の両端に共通する最も近い祖先
グループを探し、そこに結合情報を格納する)。そのため `get_bond_list()` は再帰的に
子孫の相対パスを自分の絶対パスと連結して展開する。

> **既知の不具合**: `assign_charges()`(atomgroup.py:551)内に `print(index, len(charges))`
> というデバッグ出力がそのまま残っている。

## 5. フォーマット変換

### 5.1 `format.py` — `Format`

`AtomGroup` の階層構造が PDB 的な「モデル群→モデル→チェイン→残基→原子」という規約に
従っているかを検証する静的メソッド群。検証のみで例外は投げず、`logger` へ警告/デバッグ
ログを出しつつ `bool` を返す。

| メソッド | 期待する構造 |
| --- | --- |
| `is_residue(res)` | 直下にサブグループを持たず、原子のみを持つ |
| `is_chain(chain)` | 直下に残基グループのみを持ち、原子を直接持たない |
| `is_protein(model)` | 直下にチェイングループのみを持ち、原子を直接持たない |
| `is_models(models)` | 直下にモデルグループのみを持ち、原子を直接持たない |

いずれも `AtomGroup.__init__` で `from .format import *` されるため、パッケージ直下に
`Format` として公開される。

### 5.2 `xyz.py` — `Xyz`

XYZ形式(1行目:原子数、2行目:コメント、以降:`記号 x y z`)の読み書き。

| メソッド | 説明 |
| --- | --- |
| `load(file_path)` | XYZファイルを読み込み内部リストに格納 |
| `save(file_path)` / `get_text()` | XYZ形式文字列として書き出し |
| `get_atom_group()` | 読み込んだ内容から `AtomGroup` を構築(フラットな1階層) |
| `set_by_atomgroup(atomgroup)` | `AtomGroup`(モデル階層を再帰的に辿る)からXYZ内部データへ変換 |

> **既知の不具合**: コンストラクタ `Xyz(file_path_str)` (xyz.py:47-48) は
> `self.load(file_path)` を呼んでいるが、実際の引数は `rhs` に束縛されており
> `file_path` という変数は定義されていない。**`Xyz("path/to/file.xyz")` の形で
> 文字列を渡すコンストラクタ経由の読み込みは `NameError` になる**
> (`Xyz().load(path)` のように明示的に呼べば問題なく動作する)。

### 5.3 `gro.py` — `SimpleGro`

GROMACS `.gro` 形式(固定カラム幅のテキスト)の読み書き。座標はGROMACSのnm単位から
Å単位(`AtomGroup`/`Atom`の内部単位)へ ×10 して変換する(逆方向は ×0.1)。

| メソッド | 説明 |
| --- | --- |
| `load(path)` | 固定カラム位置([0:5]=残基番号, [5:10]=残基名, [10:15]=原子名, [15:20]=原子番号, 以降 x/y/z, vx/vy/vz)でパース |
| `get_atomgroup()` | `model(1) → chain("_") → residue(res_id) → atom` の3階層構造を構築。原子記号は原子名の先頭1〜2文字から `PeriodicTable` に問い合わせて推定 |
| `set_by_atomgroup(atomgroup)` | `model→chain→residue→atom` の4階層 `AtomGroup` から `.gro` 内部データへ変換(残基番号は出現順に1から振り直す) |
| `__str__()` | `.gro` 形式のテキストを生成。box vectors は `atomgroup.box()` の対角成分のみ(斜方晶近似) |

### 5.4 `mol2.py` — `SimpleMol2`

Tripos MOL2 形式(`@<TRIPOS>MOLECULE` / `ATOM` / `BOND` セクション)への書き出し専用
(読み込みは未実装)。

| メソッド | 説明 |
| --- | --- |
| `set_by_atomgroup(atomgroup)` | 対象の `AtomGroup` を保持し、原子順を記録したインデックステーブルを作成 |
| `save(path)` / `__str__()` | MOLECULE/ATOM/BONDの3セクションを組み立てて出力。BONDセクションは `AtomGroup.get_bond_list()` のパスから `Select_Path` で該当原子を検索し、インデックステーブル上の番号に変換して出力 |

### 5.5 `mmcif.py` — `SimpleMmcif`

mmCIF (`data_...` ブロック、`key value` 形式、`loop_` テーブル、`;...;` の複数行値) を
正規表現ベースで解析する簡易パーサ。主に化合物定義 (`_chem_comp_atom.*`,
`_chem_comp_bond.*`) の読み込みに使われる。

| メソッド | 説明 |
| --- | --- |
| `load(path)` | mmCIFファイル全体を `{data_block_name: (kv_dict, [loop_table, ...])}` にパース |
| `_get_line()` | 1論理行を取得。`;`で始まる複数行値ブロックを1つの引用符付き値に結合する下請け |
| `get_molecule_names()` | パース済みデータブロック名の一覧 |
| `get_atomgroup(name)` | 指定した化合物(データブロック)を `AtomGroup` に変換。`_chem_comp_atom.pdbx_model_Cartn_*_ideal` を優先し、無ければ `model_Cartn_*` を使用して座標を決定。`D`(重水素)は `H` 扱い |
| (bond) | `_chem_comp_bond.value_order` (`SING`/`DOUB`/`TRIP`) を結合次数 1/2/3 に変換して `add_bond()` |
| `load_msgpack()`/`save_msgpack()` | パース結果自体をmsgpackで保存/復元(brd形式とは別のキャッシュ用途) |

### 5.6 `amber_prmtop.py` — `AmberPrmtop`

Amber の `prmtop`(トポロジ)+`inpcrd`(座標)ペアを読み込み、`AtomGroup` に変換する。

| メソッド | 説明 |
| --- | --- |
| `_load_prmtop(path)` | `%FLAG ATOM_NAME` / `%FLAG CHARGE` / `%FLAG ATOMIC_NUMBER` セクションのみを固定長パースで抽出(他のFLAGは無視) |
| `_load_inpcrd(path)` | 1行目タイトル・2行目原子数の後、12桁固定幅×3列で座標を読み込む |
| `get_atomgroup()` | フラットな1階層の `AtomGroup` を構築 |

電荷は Amber 内部単位から電荷素量(e)単位へ `charge_amber / 18.2223` で変換される
(Amberのprmtop格納値は電荷を18.2223倍したスケール)。

> **実装メモ**: `_check_data()`, `_read_atom_name()`, `_read_charges()`,
> `_read_atomic_number()` にデバッグ用 `print()` がそのまま残っており、
> `AmberPrmtop(...)` を呼ぶだけで標準出力に件数やパース中の行が出力される。

### 5.7 `biopdb.py` — `Pdb`

PDBフォーマットの読み書きを担当する、本パッケージで最も作り込まれたフォーマットI/O
クラス。単純なパース/シリアライズに加えて、**Amber(AmberTools/`tleap`/`reduce`)と
一般(formal)PDBの命名規則の差異を吸収するリネームテーブル**を持つのが特徴。

| メソッド | 説明 |
| --- | --- |
| `load(path)` | `ATOM`/`HETATM`/`MODEL`/`TER`/`SSBOND` レコードを固定カラム位置でパース。チェインIDが空白の場合はチェイン出現順にA, B, C, ... を自動割当。元素記号が欠損している場合は原子名からのヒューリスティックで推定 |
| `get_atomgroup(select_model=None, select_altloc="A")` | `model_<N> → chain → residue(res_seq) → atom` の階層で `AtomGroup` を構築。`SSBOND` レコードはチェイン・残基から `SG` 原子を検索し `add_bond()` でジスルフィド結合として登録。altloc(異なる配座)は指定した1つ(既定 `"A"`)のみ採用 |
| `set_by_atomgroup(atomgroup, is_charge2tempfactor=False)` | `AtomGroup`→内部データへ逆変換。C末端(`OXT`を含む残基)や各チェインの末尾に `TER` レコードを自動挿入。`is_charge2tempfactor=True` で電荷をB-factor欄に書き出す(電荷の可視化用) |
| `__str__()` | PDB形式のテキストを生成(`ATOM`/`TER`行を固定カラムでフォーマット) |
| `renumber()` | 全モデルの `serial` 番号を1から振り直す |
| `get_modpdb_atomgroup(ag_protein)` | `mode`(`None`または`"AMBER"`)に応じて残基名・原子名を変換した複製を返す。`_modpdb_res`/`_modpdb_resatom`/`_modpdb_atom` の3段階で「HIS→HID/HIE/HIP」「NA→Na+/NA」「NME内の原子名差異」等を変換する |
| `_rename_to_amber_dialect(res)` | `HIS` 残基を、Hδ/Hε の有無から `HID`/`HIE`/`HIP` に自動判別してリネーム(プロトン化状態の推定) |

`AmberToolsVer` クラス変数(既定22)により、AmberToolsのバージョンによって異なる
`NME`(N-メチルアミド末端)原子名テーブルを切り替える。

`main()` はCLIエントリポイント(`optparse` 使用、Python標準の非推奨モジュール)で、
PDBファイルを読み込んでパース結果をそのまま標準出力に表示するのみ
(`scripts/`配下には対応する専用スクリプトは無く、`biopdb.py` 単体実行用)。

## 6. 構造操作

### 6.1 `select.py` — `Select` セレクタ群

`AtomGroup.select()` に渡す「条件オブジェクト」の集合。共通インターフェースは
`Select.is_match(obj)` (`obj` は `Atom` または `AtomGroup`)を実装すること。

| クラス | 条件 |
| --- | --- |
| `Select_Symbol(symbol)` | 元素記号が一致する原子 |
| `Select_Name(name)` | `name`属性(前後空白を除去)が完全一致 |
| `Select_Path(query, use_wildcard=True)` | パス文字列で選択。**非推奨** — `use_wildcard=True`(既定)だと `Select_Path_wildcard` 相当、`False` だと完全一致のみ。`use_wildcard=True` 使用時に `logger.warning` で非推奨警告 |
| `Select_Path_simple(path)` | パスの完全一致のみ |
| `Select_Path_wildcard(pattern)` | `*`→`.*`, `?`→`?` に変換した正規表現でパスにマッチ(`Select_Path`の後継、ワイルドカード専用) |
| `Select_PathRegex(regex)` | パスに対する任意の正規表現 (`re.search`) |
| `Select_Range(pos, d)` | 指定座標から半径 `d` 以内の原子(距離の二乗で比較、`Position`) |
| `Select_Atom(atom, distance=0.1)` | 原子番号が一致し、かつ座標が `distance` 以内の原子(同一性判定に使用) |
| `Select_AtomGroup(ref_atomgroup, range=1.0e-5)` | 参照 `AtomGroup` に含まれる原子(原子番号+距離)のいずれかに一致する原子。`AtomGroup.restructure()` の内部で使用 |

### 6.2 `aminoacid.py` — `AminoAcid`

標準/非標準の20+αアミノ酸3文字コード(`HIE`/`HIP`/`CYX`/`ASX`/`GLX`/`XAA`等の
プロトン化状態・曖昧表記を含む)のリストを保持し、`AtomGroup`(残基)がアミノ酸か
どうかを名前だけで判定する。実質的にコンストラクタ引数は使われないユーティリティクラス。

| メソッド | 説明 |
| --- | --- |
| `is_aminoacid(atomgroup)` | `atomgroup.name` が既知のアミノ酸コード一覧に含まれるか |

### 6.3 `ssbond.py` — `SSBond`

タンパク質モデル内のジスルフィド結合(S-S結合)を、`CYS`/`CYX` 残基の `SG` 原子間距離
から検出する。

| メソッド | 説明 |
| --- | --- |
| `get_bonds()` | (キャッシュしつつ)`(path1, path2)` のタプルのリストを返す |
| `_check()` | 全チェインから `CYS`/`CYX` 残基の `SG` 原子を収集し、`_check_SGs()` へ |
| `_check_SGs(SGs)` | 全ペアの距離を計算し、`_ss_bond_max_length` (2.1Å × 1.1 = 2.31Å) 未満なら結合と判定 |

`biopdb.Pdb.get_atomgroup()` は `PDB`ファイル中の `SSBOND` レコードを直接使うため、
このクラスは(座標のみから)**ジスルフィド結合を推定し直す**独立した経路として使われる
(例: `SSBOND`レコードを持たない構造ファイルに対する結合推定)。

### 6.4 `ionpair.py` — `IonPair`

タンパク質中の酸性/塩基性残基や末端基を検出し、4.0Å以内にある組を「イオン対」として
列挙する。

| 分類種別 | 対象 |
| --- | --- |
| アニオン | `GLU`(側鎖COO⁻の中心), `ASP`(同), C末端(`OXT`を持つ残基のCOO⁻中心) |
| カチオン | `LYS`(NZ位置), `ARG`(グアニジノ基の中心 or NH1側 or NH2側の3通り), N末端(`H3`を持つ残基のNH3⁺方向) |

| メソッド | 説明 |
| --- | --- |
| `get_ion_pairs()` | 全アニオン×カチオンの組み合わせで距離4.0Å未満のものを `(anion_path, cation_path, anion_type, cation_type)` のリストで返す |
| `_get_center_GLU/ASP/LYS/ARG/Nterm/Cterm()` | 各官能基の「代表座標」を計算(COO⁻は2つの酸素の中点、グアニジノ基は3原子の重心など) |

`Neutralize` から「既に電荷的に中和されている(=無視すべき)残基」の判定に利用される。

### 6.5 `neutralize.py` — `Neutralize`

タンパク質全体を電荷的に中性化するため、酸性/塩基性側鎖や末端基の近傍に対イオン
(Na⁺/Cl⁻)を配置する。`brd-*`系ではなく `scripts/neutralize.py` から呼ばれる。

| 対象残基 | 追加する対イオン | 配置ロジック |
| --- | --- | --- |
| N末端 (`H3`所持、`PRO`は別ロジック) | Cl⁻ | NH3型/NH2型の重心から外側 3.187Å |
| C末端 (`OXT`所持) | Na⁺ | COO型の中点から外側 2.521Å |
| `GLU`/`ASP` | Na⁺ | 側鎖カルボキシル基のCOO型判定 |
| `LYS` | Cl⁻ | NH3型(側鎖アミン)判定 |
| `ARG` | Cl⁻ | グアニジノ基(中央/NH1側/NH2側)判定 |
| `FAD` | Na⁺ ×2 | リン酸基(POO型)2箇所 |

`Modeling` クラス(6.6節)の `neutralize_*` 系メソッド(座標計算)を呼び出し、
`_add_ions()` で命名衝突を避けながら `AtomGroup` に対イオンを追加する。

> **既知の不具合(未使用コード)**: `_exempt_list()`(neutralize.py:20)は
> `self._model` を参照しているが、このクラスのどこにも `self._model` を設定する
> コードがない(`__init__` は `self._neutral_obj` のみ設定)。ただし
> `_neutralize()` 内でこのメソッドへの呼び出しはコメントアウトされている
> (`exempt_list = []  # self._exempt_list()`)ため、現状は実行されず影響はない
> (=イオン対による「既に中和済みなので対イオンを追加しない」除外ロジックは
> 実質的に無効化されている)。

### 6.6 `modeling.py` — `Modeling`

アミノ酸残基の欠損部分(末端保護基 ACE/NME、中性化イオンの座標など)を、
参照構造の重ね合わせ(`Superposer`)や幾何計算で補完するクラス群。

- **ACE/NME末端キャップの付加**: パッケージ同梱の参照構造
  `data/ACE_ALA_NME_{trans1,trans2,cis1,cis2}.brd`(4種のアミド結合配座)を
  `Modeling.__init__()` でロードしておき、`get_ACE(res, next_aa)` /
  `get_NME(res, next_aa)` で対象残基のCA/N/C/O等を参照構造にフィッティング
  (`Superposer`でRMSD最小の配座を選択)して末端キャップ原子を生成する。
  `get_ACE_simple`/`get_NME_simple` は隣接残基のCαをそのままメチル基として使う簡易版。
- **メチル基付加**: `add_methyl(C1, C2)` はエタン分子のテンプレート座標を
  `arbitary_rotate_matrix()`(任意軸回転行列の生成、ロドリゲスの回転公式に相当)で
  回転・並進させてC1側に水素3つを付加する。
- **NH3幾何生成**: `get_NH3(angle, length)` は正四面体角に基づくアンモニア型
  配置(N-H×3)を生成する(プロトン化アミン等のテンプレート用)。
- **中性化用イオン座標**: `neutralize_Nterm/Cterm/GLU/ASP/LYS/ARG/FAD` と、
  内部で使う `_get_neutralize_pos_NH3_type/NH2_type/COO_type/POO_type` が、
  各官能基の幾何中心から一定距離(官能基の種類ごとに固定値: NH3型3.187Å,
  COO型2.521Å, POO型2.748Å)だけ外側に対イオンの座標を計算する。`Neutralize`
  クラス(6.5節)から呼ばれる。
- `select_residues(chain, from_resid, to_resid)`: 連続した残基番号の範囲を
  `AtomGroup` として抽出。

### 6.7 `superposer.py` — `Superposer` (Kabschアルゴリズム)

2つの `AtomGroup` に共通するキーの原子ペアを抽出し、剛体変換(回転+並進)で
重ね合わせたときのRMSDと回転行列を求める。実装はKabschアルゴリズム
(共分散行列の特異値分解/固有値分解による最適回転の導出)。

| メソッド/プロパティ | 説明 |
| --- | --- |
| `Superposer(ag1, ag2)` | `ag1`(動かす側)と `ag2`(基準)を受け取り、共通原子キーの座標ペアを抽出 |
| `rmsd` | 最適重ね合わせ後のRMSD(遅延評価プロパティ、内部で`rotation_mat`等を計算) |
| `rotation_mat` | 最適回転行列(3x3 `Matrix`) |
| `superimpose(atomgroup)` | 任意の`AtomGroup`(`ag1`と同じ座標系にあるもの)を、求めた回転・並進で`ag2`の座標系に重ね合わせて返す |

回転行列の導出は、共分散行列 `R = Σ p2ᵢ p1ᵢᵀ` の `RᵀR` を対称行列として固有値分解し、
`_make_right_handed()` で右手系に補正した基底 `a`, `b` から `r_ij = Σ b_ki a_kj` として
回転行列を組み立てる、という手順(`_get_rotation_matrix()`)。

`Superposer.__init__` にはコメントアウトされた `_calc()` 呼び出しが残っており、
実際には各プロパティへの初回アクセス時に遅延計算される設計になっている。

> **実装メモ**: `_get_rotation_matrix()` 内に多数の `print()` デバッグ出力
> (`eigval`, `eigvec`, `eigvec2`, `make right handled`, `b` 等)が残っており、
> `Superposer` を使うたびに標準出力へ大量のログが出力される。

### 6.8 `superposer_quaternion.py` — `Superposer_quaternion`(四元数法、実験的)

`Superposer`(Kabsch法)とは別に、剛体変換の最適化を**四元数**で行う代替実装。
プロパティ経由の遅延評価チェーン(`center1/2` → `r_A/r_B`(重心补正後座標) →
`va/vb`(補助ベクトル) → `matB`(4x4対称行列) → `eigval`(最小/最大固有ベクトル
=最適四元数) → `matR`(四元数から回転行列) → `rmsd`)という設計は `Superposer` と
同様だが、`__init__` に対応する `superimpose()` 相当のメソッドは実装されていない
(回転行列 `matR` を取得した後の座標変換は呼び出し側の責務、または未実装のまま)。

> **既知の不具合**: `calc()` メソッド(superposer_quaternion.py:143-163)は
> `atom_group1`, `atom_group2`, `position1`, `position2` など**未定義の変数**を
> 参照しており、また `self.match_positions`/`self.calc_center`/`self.make_B`/
> `self.make_R` のように(実際の実装は `_match_positions` 等アンダースコア始まり
> の別名)**存在しないメソッド名**を呼び出している。呼び出せば確実に
> `NameError`/`AttributeError` になる、事実上の**壊れたデッドコード**。
> このクラスを使う場合は `calc()` を呼ばず、`rmsd`/`matR` 等のプロパティに
> 直接アクセスする必要がある。
>
> また `_shift_positions()`, `_make_va()`, `_make_vb()`, `_make_B()`,
> `_get_r_A()`/`_get_r_B()` などほぼ全メソッドに `print()` デバッグ出力が残る。
>
> `scripts/superposer.py` は `-q`/`--quaternion` オプションで `Superposer_quaternion`
> を選択できるが、`Superposer_quaternion` には `superimpose()` メソッドが
> **存在しない**。`scripts/superposer.py` は `-q` 指定の有無によらず必ず
> `sp.superimpose(atomgroup1)` を呼び出す(§8.6参照)ため、**`superposer.py -q`
> を実行すると `rmsd` の計算・表示までは成功するが、その後 `AttributeError` で
> 必ず異常終了する**。

## 7. その他インフラ

### 7.1 `dbmanager.py` — `DbManager` / `DbTable`

`sqlite3` を薄くラップしたシンプルなORM風ヘルパー。分子構造そのものではなく、
計算ジョブ管理・メタデータ管理などの用途を想定した汎用DBアクセス層。

| メソッド | 説明 |
| --- | --- |
| `DbManager(db=':memory:', sql_debugout=False)` | DBファイルパス(既定はインメモリ)を開く |
| `create_table(name, field_names, primary_key=None)` | テーブル作成。`field_names` はリスト(型指定なし)または `{name: type}` の辞書 |
| `get_table_names()` / `has_table(name)` | テーブル一覧・存在確認 |
| `get_field_names(table)` / `get_primary_keys(table)` | カラム名一覧・主キー一覧の取得(`PRAGMA table_info`使用) |
| `insert(table, contents)` / `update(table, contents, where)` / `delete(table, where)` | `contents`/`where` は `{field: value}` 辞書、プレースホルダ(`?`)を使ったパラメータ化クエリでSQLインジェクションを回避 |
| `select(table, fields=None, where=None)` | `where` は辞書(AND結合のみ)または生SQL文字列。結果を `[{field: value}, ...]` で返す |
| `execute(sql, parameters=None)` / `get_results(sql)` | 任意SQLの実行 |
| `set_user_version(v)` / `get_user_version()` | `PRAGMA user_version` の読み書き(スキーマバージョン管理用) |
| `pp_table(table)` / `pp(data)` | 結果をテーブル状に整形して文字列化するデバッグ用プリティプリント |
| `__getitem__(key)` | `db[table_name]` で `DbTable` オブジェクトを取得(存在しなければ `None`) |

> **既知の不具合**: `create_table()`(dbmanager.py:90-93)で `field_names` に
> 辞書(型指定あり)を渡した場合、`for k, v in field_names:` は辞書を素のまま
> イテレートしてしまうため(`.items()`の呼び忘れ)キー文字列を分解しようとして
> 失敗し、さらに続く `'{name} {type}'.format(k, v)` も名前付きプレースホルダに
> 対して位置引数を渡しているため `KeyError: 'name'` になる。**型指定付きの
> `create_table()` は事実上動作しない**(型指定なしのリスト渡しのみ動作する)。
>
> **既知の不具合**: `get_results()`(dbmanager.py:304-323)は、複数行分の
> `row_items` を組み立てるループの外側で1回だけ `answer.append(row_items)`
> しているため(インデントの誤り)、**クエリ結果が2行以上あっても最後の1行しか
> 返らない**。さらに `data` が0行の場合は `row_items` が未定義のまま参照され
> `NameError` になる。`get_user_version()` は常に1行しか返らない
> `PRAGMA user_version` にしか使っていないため、この不具合は表面化していない。
> (同種の集計処理である `select()` メソッドは `answer = [{}] * len(data)` で
> 事前確保する実装になっており、こちらにはこの不具合はない)

### 7.2 `mail.py` — `Mail`

`smtplib` を使ったメール送信ヘルパー。ジョブ完了通知など運用系スクリプトからの
利用を想定。

| メソッド | 説明 |
| --- | --- |
| `load_config(path)` / `save_config(path)` | `configparser` 形式の設定ファイル(`[mail]` セクション: `smtp_server`, `smtp_port`, `use_SSL`, `smtp_account`, `smtp_password`, `from_address`)の読み書き |
| `send()` | `MIMEText` でメールを組み立て、`use_SSL` に応じて `SMTP_SSL` または `SMTP`+`STARTTLS` で送信。`smtp_account`/`smtp_password` で認証 |
| プロパティ | `smtp_server`, `smtp_port`, `use_SSL`, `smtp_account`, `smtp_password`, `from_address`, `to_address`, `charset`(既定 `ISO-2022-JP`、日本語メール想定), `subject`, `text` |

`smtp_password` は設定ファイルに**平文**で保存/読み込みされる。設定ファイルの
パーミッション管理は呼び出し側の責任となる。

## 8. CLIスクリプト仕様 (`scripts/*.py`)

全スクリプトは `import proteindf_bridge as bridge` した上で `argparse` を用いて
CLI引数を処理する薄いラッパーであり、内部で本体モジュール(§2〜7)のクラス/関数を
1〜数回呼び出すだけの構成になっている。共通して `-v`/`--verbose` を持つものが多い。
`.brd` ファイルの読み書きには一貫して `bridge.load_atomgroup()` /
`bridge.save_atomgroup()`(=`functions.py`、msgpack)が使われる。

### 8.1 フォーマット変換系

| スクリプト | 位置引数 | 主なオプション | 処理内容 |
| --- | --- | --- | --- |
| `pdb2brd.py` | `PDB_FILE`, `BRD_FILE` | `-m/--model`(モデル番号選択), `-l/--alt_loc`(既定`A`), `-d/--debug` | `biopdb.Pdb.load()`→`get_atomgroup()`→ brd保存 |
| `brd2pdb.py` | `FILE`(.brd) | `-o/--output`, `-a/--amber`(Amber方言リネーム), `-c/--charge2tempfactor` | brd読込→`Pdb.set_by_atomgroup()`→PDBテキスト出力 |
| `gro2brd.py` | `GRO_FILE`, `BRD_FILE` | `-v` | `SimpleGro.load()`→`get_atomgroup()`→brd保存 |
| `brd2gro.py` | `FILE`(.brd) | `-v` | brd読込→`SimpleGro.set_by_atomgroup()`→標準出力 |
| `xyz2brd.py` | `XYZ_PATH`, `BRD_PATH` | `-v` | `Xyz.load()`→`get_atom_group()`→brd保存 |
| `brd2xyz.py` | `FILE`(.brd) | `-o/--output`, `-v` | brd読込→`Xyz.set_by_atomgroup()`→XYZ出力 |
| `mmcif2txt.py` | `mmCIF_FILE` | `-w/--write`(msgpackキャッシュ書き出し), `-v` | `SimpleMmcif.load()`→パース結果を人間可読形式で表示 |
| `mmcif2mol2.py` | `mmCIF_FILE` | `-o/--output`, `-w/--write`, `-r/--read`(msgpackキャッシュから読込) | mmCIF(または再利用キャッシュ)→`AtomGroup`→`SimpleMol2`→MOL2出力 |
| `read_amber_prmtop.py` | `prmtop_file`, `inpcrd_file` | — | `AmberPrmtop`で読み込み内容を表示(brd変換は行わない模様) |
| `mpac2yml.py` | `mpac_path`, `yaml_path` | — | msgpack→YAML変換(`functions.load_msgpack`/`get_yaml`) |
| `yml2mpac.py` | `YAML_FILE`, `MPAC_FILE` | — | YAML→msgpack変換 |
| `mpac2txt.py` | `FILE` | — | msgpackの内容を人間可読表示 |
| `brd2txt.py` | `FILE`(.brd) | `-c/--csv`, `-v` | brd読込→`AtomGroup.__str__()`または`save_csv()`で出力 |
| `db2txt.py` | `FILE`(sqlite3) | `-v` | `DbManager`でDBを開き `pp()` でテーブル内容を表示 |

### 8.2 構造編集系

| スクリプト | 位置引数 | 主なオプション | 処理内容 |
| --- | --- | --- | --- |
| `brd-select.py` / `brd-select-path.py` | `FILE`(.brd) | `-q/--query`(既定`"*"`), `-o/--output` | `Select_Path_wildcard`等で `AtomGroup.select()` を実行 |
| `brd-divide.py` | `brd_path` | `-o/--output` | 構造をサブグループ単位に分割出力 |
| `brd-divide-mainchain.py` | `brd_path` | `-o/--output` | 主鎖/側鎖単位での分割出力 |
| `brd-restructure.py` | `target_brd_path`, `ref_brd_path` | `-o/--output_path`, `-r/--range` | `AtomGroup.restructure()` の呼び出し(§4.7参照) |
| `brd-renumber-resid.py` | `FILE`, `increment` | `-o/--output`, `-q/--query` | 指定クエリに一致する残基番号を `increment` だけシフト |
| `brd-setup-bond.py` | `FILE` | `-o/--output` | `Bond().setup(atomgroup)` で距離ベースの結合を自動推定(§4.6の不具合の影響を受けうる) |
| `brd-show-bonds.py` | `FILE` | `-o/--output` | `AtomGroup.get_bond_list()` を整形表示 |
| `brd-show-res.py` | `FILE` | — | 残基一覧を表示 |
| `remove_wat.py` | `FILE` | `-o/--output` | `Utils.remove_WAT()` の呼び出し |
| `reorder.py` | `INPUT_FILE`, `OUTPUT_FILE` | `-d/--debug` | 原子順の並べ替え |
| `neutralize.py` | `INPUT_FILE`, `OUTPUT_FILE` | `-d/--debug` | `Neutralize` の呼び出し(§6.5) |
| `crystallize.py` | `INPUT_BRD_PATH`, `OUTPUT_BRD_PATH` | `--num_x`, `--num_y`, `--num_z` | 単位格子を指定数だけ複製して結晶構造を生成 |
| `brd-setup-bond.py` 以外の未収載スクリプト(`superposer.py`) | 下記8.6参照 | | |

### 8.3 解析系

| スクリプト | 位置引数 | 主なオプション | 処理内容 |
| --- | --- | --- | --- |
| `brd-box.py` | `FILE`, `num_of_molecules`, `density` | `--use-nm` | `AtomGroup.box()`/密度からセルサイズを逆算 |
| `brd-density.py` | `FILE`, `num_of_molecules` | `--box`, `--use-nm` | セルサイズと分子数から密度を算出 |
| `brd-formula.py` | `FILE` | — | `AtomGroup.get_formula()` を表示 |

### 8.4 開発用

| スクリプト | 説明 |
| --- | --- |
| `doctest_runner.py` | `proteindf_bridge` 配下の各モジュールに埋め込まれた doctest を一括実行 |
| `module_inspect.py` | モジュール内容の調査用デバッグツール |

### 8.5 その他

`db2txt.py`(§8.1) は解析/開発というより `DbManager`(§7.1)の中身をダンプする
デバッグツールに近い位置づけ。

### 8.6 `superposer.py`

`FILE1`, `FILE2`(いずれも `.brd`)を読み込み、`Superposer`(既定)または
`-q/--quaternion` 指定時は `Superposer_quaternion` で重ね合わせ、RMSDと
重ね合わせ後の構造(`FILE1`側)を標準出力に表示する。

> **既知の不具合(実行時エラー)**: `main()`(superposer.py:63-76)は
> `use_quaternion` の値によらず必ず `sp.superimpose(atomgroup1)` を呼ぶが、
> `Superposer_quaternion` (§6.8) には `superimpose()` が実装されていない。
> **`superposer.py FILE1 FILE2 -q` を実行すると、RMSD値の表示までは成功するが
> 直後に `AttributeError` で異常終了する。** `-q` オプションは事実上使用不可。

## 9. 既知の不具合一覧(まとめ)

解析中に見つかった、修正を要すると思われる箇所の一覧(詳細は各節を参照)。

| # | 箇所 | 症状 | 深刻度目安 |
| --- | --- | --- | --- |
| 1 | `matrix.py:534` `SymmetricMatrix.get_raw_data()` | `range(dim*(dim+1)/2)` がPython3では`float`を渡すことになり`TypeError` | 中(対称行列を直接brd/msgpack保存する経路でのみ発現) |
| 2 | `bond.py`(全体) | `SymmetricMatrix`/`AtomGroup` の import 漏れ。単体importで`NameError`の恐れ | 低〜中(パッケージ経由の間接import次第で発現しない場合あり) |
| 3 | `xyz.py:47-48` `Xyz.__init__` | `Xyz("path")`のような文字列引数コンストラクタが`NameError`(`file_path`が未定義) | 中(`Xyz().load(path)`なら回避可) |
| 4 | `atomgroup.py:551` `assign_charges()` | デバッグ用`print()`の消し忘れ | 低(実害なし、ログ出力ノイズ) |
| 5 | `amber_prmtop.py` 各`_read_*`/`_check_data` | デバッグ用`print()`の消し忘れが多数 | 低(実害なし、標準出力ノイズ) |
| 6 | `neutralize.py:20` `_exempt_list()` | 未設定の`self._model`を参照(ただし呼び出し自体がコメントアウトされ未使用) | 低(現状デッドコード、ただし将来有効化すると壊れる) |
| 7 | `superposer.py` `_get_rotation_matrix()`ほか | デバッグ用`print()`が多数残り、呼ぶたびに標準出力を汚染 | 低〜中(ライブラリとして呼ばれた際の副作用) |
| 8 | `superposer_quaternion.py:143-163` `calc()` | 未定義変数・存在しないメソッド名を参照。呼び出せば必ず例外 | 高(ただし`calc()`自体は他から呼ばれていないデッドコード) |
| 9 | `superposer_quaternion.py`(全体) | `superimpose()`が未実装 | 高(`scripts/superposer.py -q`が実行時に必ず`AttributeError`で落ちる。§8.6) |
| 10 | `dbmanager.py:90-93` `create_table()` | 型指定辞書渡し時、`.items()`呼び忘れ+フォーマット文字列の名前/位置引数不一致で`KeyError` | 中(型指定なしのリスト渡しでは問題なし) |
| 11 | `dbmanager.py:304-323` `get_results()` | インデント誤りで複数行結果のうち最後の1行しか返らない、0行時は`NameError` | 中(`get_user_version()`は常に1行なので表面化しない) |
| 12 | `mail.py` `smtp_password` | 設定ファイルに平文で保存される | 低(運用上の注意点、バグではない) |
| 13 | `setup.cfg` vs `_version.py` | `version = 2024.3.0` (setup.cfg) と `__version__ = "2022.2.5"` (`_version.py`) が乖離 | 低(表示上の不整合) |

いずれも本仕様書作成のための静的解析(コードリーディング)で発見したもので、
実行テストによる再現確認は行っていない項目を含む。修正の要否・優先度は
別途プロジェクト側で判断されたい。
