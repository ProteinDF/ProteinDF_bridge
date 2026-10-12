# ProteinDF 計算結果データ形式仕様書 (`docs/qc-result-format.md`)

本ドキュメントは、ProteinDFの計算結果データ（HDF5形式 `reference.h5` / `pdfresults.h5`、および外部バイナリ行列 `.mat` / ベクトル `.vtr`）の形式調査結果をまとめ、Rust移植版 `proteindf-bridge` における計算結果I/O（PR#60）および軌道別Mulliken寄与計算（PR#61）の設計基礎とするものである。

調査対象・参照ソース:
- テストデータ: `/home/toshi/work/dev/pdf-dev/ProteinDF_test/serial_dev/*/reference.h5`
- 既存ツール: `/home/toshi/work/dev/pdf-dev/ProteinDF_pytools`
  - `proteindf_tools/pdfparam_hdf5.py`
  - `proteindf_tools/matrix.py`
  - `proteindf_tools/vector.py`
  - `proteindf_tools/orbinfo.py`
  - `proteindf_tools/basisset.py`
  - `scripts/pdf-archive-h5.py`
- ProteinDF本体: `/home/toshi/work/dev/pdf-dev/ProteinDF`
  - `src/libpdf/DfObject.cpp`
  - `devtool/pyscf_scf_check.py`

---

## 1. HDF5ファイルの中身の要約

### 1.1 ルート属性 (Root Attributes)

`reference.h5` のルートグループ (`/`) には、計算条件および計算結果のサマリ情報が属性として格納されている。
(`pdfparam_hdf5.py` 190-218行目 `_load_model` で読み込み)

| 属性名 | 型 (実データ) | 例 (`H2O_HF`) | 説明 |
| --- | --- | --- | --- |
| `method` | 文字列 (`str`) | `"rks"` | 計算手法 (`"rks"`: スピン制限閉殻, `"uks"`: スピン非制限開殻) |
| `xc_functional` | 文字列 (`str`) | `"HF"` | 交換相関汎関数 / 手法名 (例: `"HF"`, `"B3LYP"`, `"BLYP"`, `"SVWN5"`) |
| `num_of_atoms` | 64bit整数 (`int64`) | `3` | 分子内の原子数 ($N_{\text{atom}}$) |
| `num_of_AOs` | 64bit整数 (`int64`) | `24` | 原子基底関数 (AO) の総数 ($N_{\text{AO}}$) |
| `num_of_MOs` | 64bit整数 (`int64`) | `24` | 分子軌道 (MO) の総数 ($N_{\text{MO}}$) |
| `iterations` | 64bit整数 (`int64`) | `49` | 実行されたSCFイテレーション回数 |
| `max_iterations` | 64bit整数 (`int64`) | `100` | 最大イテレーション設定値 |
| `scf_converged` | ブール値 (`bool`) | `True` | SCFが収束したか否か |
| `step_control` | 文字列 (`str`) | `"create integral guess scf force"` | 実行ステップ制御文字列 |
| `guess` | 文字列 (`str`) | `"core"` | 初期推測解の種類 |
| `cut_value` | 文字列 (`str`) | `"1.0e-10"` | 積分カットオフ閾値 (文字列表現) |
| `CDAM_tau` | 文字列 (`str`) | `"1.0E-10"` | Cholesky分解近似パラメータ |
| `CD_epsilon` | 文字列 (`str`) | `"1.0E-6"` | Cholesky分解閾値 |
| `scf_acceleration` | 文字列 (`str`) | `"damping"` | SCF加速手法 (例: `"damping"`, `"diis"`) |
| `scf_acceleration_damping_damping_factor` | 64bit浮動小数点 (`float64`) | `0.67` | ダンピング係数 |
| `comment` | 文字列 (`str`) | `""` | ユーザーコメント |

> [!NOTE]
> `cut_value`, `CDAM_tau`, `CD_epsilon` は浮動小数点数値ではなく文字列 (`str`) として保存されている（実データで確認済み）。

---

### 1.2 グループ・データセットの木構造一覧

`serial_dev` の `reference.h5` (回帰テスト用基準値)、および `pdf-archive-h5.py` が生成する `pdfresults.h5` (標準計算結果アーカイブ) を走査した結果、木構造の基本骨格は共通しているが、含まれるデータセットの範囲が異なる。

#### (A) 標準アーカイブ `pdfresults.h5` (`pdf-archive-h5.py` 出力、実例: Ne B3LYP)
`pdf-archive-h5.py` により生成される標準の計算結果ファイルでは、重なり行列 `s`、Fock行列 `F`、密度行列 `P`、ハミルトニアン `h`, `h2` を含め、全反復の行列・ベクトルが完全に保存される。
```text
/
├── C/                                  [Group] 分子軌道(MO)係数行列 (全反復 1..iterations)
│   ├── rks_1                           [Dataset] shape=(N_AO, N_MO), dtype=float64
│   └── ...                             (rks_1..rks_{iterations})
├── F/                                  [Group] Fock行列 (全反復 1..iterations)
│   ├── rks_1                           [Dataset] shape=(N_AO, N_AO), dtype=float64
│   └── ...                             (rks_1..rks_{iterations})
├── P/                                  [Group] 密度行列 (全反復 1..iterations)
│   ├── rks_1                           [Dataset] shape=(N_AO, N_AO), dtype=float64
│   └── ...                             (rks_1..rks_{iterations})
├── TEs                                 [Dataset] shape=(iterations,), dtype=float64 (各SCFの全エネルギー)
├── basisset/                           [Group] 基底関数系定義
│   ├── cgto                            [Dataset] shape=(N_cgto, 4), dtype=|S32 (固定長バイト文字列)
│   ├── name                            [Dataset] shape=(N_spec, 2), dtype=object または |S32
│   └── pgto                            [Dataset] shape=(N_pgto, 5), dtype=|S32 (原始ガウス関数)
├── control/                            [Group]
│   └── file_base_name                  [Group] 外部ファイル名テンプレート
├── energy_level/                       [Group] 軌道エネルギー (全反復 1..iterations)
│   ├── rks_1                           [Dataset] shape=(N_MO,), dtype=float64
│   └── ...                             (rks_1..rks_{iterations})
├── h                                   [Dataset] shape=(N_AO, N_AO), dtype=float64 (1電子ハミルトニアン)
├── h2                                  [Dataset] shape=(N_AO, N_AO), dtype=float64 (2電子ハミルトニアン成分)
├── molecule/                           [Group] 分子構造
│   └── xyzc                            [Dataset] shape=(N_atom, 6), dtype=|S32 (元素, ラベル, x, y, z, 電荷)
├── occ/                                [Group] 軌道占有数
│   └── rks                             [Dataset] shape=(N_MO,), dtype=float64
└── s                                   [Dataset] shape=(N_AO, N_AO), dtype=float64 (重なり行列 S)
```
> [!NOTE]
> `pdfresults.h5` には `s`, `h`, `h2`, `P`, `F` が含まれる一方、`pop/mulliken` は保存されない（`pdf-archive-h5.py` にMulliken電荷を保存する呼び出しがない）。また、`C`, `F`, `P`, `energy_level` は全SCF反復（`1..iterations`）分が格納される。

#### (B) テストデータ `reference.h5` (RKS計算例: `H2O_HF`, `Ne_HF`, `N2_RHF` 等)
ProteinDF_test の回帰テスト基準値用 `reference.h5` では、ファイルサイズ節約のため最終反復データや特定項目のみが保持されている。
```text
/
├── C/                                  [Group] 分子軌道(MO)係数行列 (最終反復のみ)
│   └── rks_{iterations}               [Dataset] shape=(N_AO, N_MO), dtype=float64
├── TEs                                 [Dataset] shape=(iterations,), dtype=float64 (各SCFの全エネルギー)
├── basisset/                           [Group] 基底関数系定義
│   ├── cgto                            [Dataset] shape=(N_cgto, 4), dtype=|S32
│   ├── name                            [Dataset] shape=(N_spec, 2), dtype=object または |S32
│   └── pgto                            [Dataset] shape=(N_pgto, 5), dtype=|S32
├── control/file_base_name              [Group] 外部ファイル名テンプレート
├── energy_level/                       [Group] 軌道エネルギー (各反復 1..iterations)
│   ├── rks_1                           [Dataset] shape=(N_MO,), dtype=float64
│   └── ...
├── molecule/                           [Group] 分子構造
│   └── xyzc                            [Dataset] shape=(N_atom, 6), dtype=|S32
├── occ/                                [Group] 軌道占有数
│   └── rks                             [Dataset] shape=(N_MO,), dtype=float64
└── pop/                                [Group] ポピュレーション解析結果 (最終反復のみ)
    └── mulliken/
        └── rks_{iterations}           [Dataset] shape=(N_atom,), dtype=float64 (最終反復のMulliken原子電荷)
```

#### (C) テストデータ `reference.h5` (UKS計算例: `O2_UHF`, `O2_UB3LYP` 等)
```text
/
├── P/                                  [Group] 密度行列 (スピンα, β別の最終反復)
│   ├── uks_alpha_{iterations}         [Dataset] shape=(N_AO, N_AO), dtype=float64
│   └── uks_beta_{iterations}          [Dataset] shape=(N_AO, N_AO), dtype=float64
├── TEs                                 [Dataset] shape=(iterations,), dtype=float64
├── basisset/                           [Group] (RKSと同形式)
├── control/file_base_name              [Group] (RKSと同形式)
├── energy_level/                       [Group] 軌道エネルギー (スピンα, β別)
│   ├── uks_alpha_{iterations}         [Dataset] shape=(N_MO,), dtype=float64
│   └── uks_beta_{iterations}          [Dataset] shape=(N_MO,), dtype=float64
├── molecule/xyzc                       [Dataset] shape=(N_atom, 6), dtype=|S32
└── occ/                                [Group] 軌道占有数 (スピンα, β別)
    ├── uks_alpha                       [Dataset] shape=(N_MO,), dtype=float64
    └── uks_beta                        [Dataset] shape=(N_MO,), dtype=float64
```

> [!IMPORTANT]
> `N2_UHF` などの一部のUKSテストデータでは、`TEs`, `basisset`, `control`, `molecule` の4グループのみが存在し、`C`, `energy_level`, `occ`, `pop` 等の収束結果データセットが保存されていない（最小限の回帰テスト用基準値として作成されたため）。

---

### 1.3 各グループ・データセットの詳細仕様

実データを直接開いて確認した各データセットのレイアウトと属性仕様は以下の通りである。

#### (1) `molecule/xyzc` (分子構造)
- **パス**: `/molecule/xyzc`
- **データ型**: `|S32` (32バイト固定長バイト文字列の2次元配列)
- **形状**: `(N_atom, 6)`
- **属性**: `header = ['symbol', 'label', 'x', 'y', 'z', 'charge']`
- **カラムの意味**:
  1. `symbol` (列0): 元素記号 (例: `b"O"`, `b"H"`, `b"N"`)
  2. `label` (列1): 原子ラベル (通常は空文字列 `b""`)
  3. `x` (列2): X座標。**単位は Bohr (a.u.)**。文字列として格納 (例: `b"1.8897261254578281"`)
  4. `y` (列3): Y座標 (Bohr)
  5. `z` (列4): Z座標 (Bohr)
  6. `charge` (列5): 原子番号 / 核電荷 $Z$。文字列として格納 (例: 酸素は `b"8.0"`, 水素は `b"1.0"`)
- **座標系の注意点**: `fl_Userinput` で `geometry/cartesian/unit = angstrom` と指定されていても、ProteinDF内部およびHDF5内では **Bohr (atomic units)** に変換されて記録されている (1.0 Å = 1.8897261254578281 Bohr)。

#### (2) `basisset` (基底関数系)
`basisset` グループには3つのデータセットが存在する (`pdfparam_hdf5.py` 289-320行目)。

1. **`basisset/name`**:
   - **形状**: `(N_spec, 2)`
   - **データ型**: `object` または `|S32`
   - **属性**: `header = ['name', 'atom_label']`
   - **データ例**: `[[b'O-DZVP2.H', b'H'], [b'O-DZVP2.O', b'O']]` (原子種ごとの基底関数名対応)

2. **`basisset/cgto`** (縮約ガウス関数: CGTO):
   - **形状**: `(N_cgto, 4)`
   - **データ型**: `|S32` (文字列配列)
   - **属性**: `header = ['name', 'cgto_id', 'shell_type', 'scale_factor']`
   - **カラムの意味**:
     - `name`: 基底セット名 (例: `b"O-DZVP2.O"`)
     - `cgto_id`: 縮約GTO番号 (0-index, 例: `b"0"`, `b"1"`)
     - `shell_type`: 軌道殻種別 (`b"s"`, `b"p"`, `b"d"`, `b"f"`)
     - `scale_factor`: スケール因子 (例: `b"1.0"`)

3. **`basisset/pgto`** (原始ガウス関数: PGTO):
   - **形状**: `(N_pgto, 5)`
   - **データ型**: `|S32` (文字列配列)
   - **属性**: `header = ['name', 'cgto_id', 'pgto_id', 'coef', 'exp']`
   - **カラムの意味**:
     - `name`: 基底セット名 (例: `b"O-DZVP2.O"`)
     - `cgto_id`: 属するCGTOのID
     - `pgto_id`: CGTO内の原始GTO番号
     - `coef`: 縮約係数 (例: `b"0.0737289"`)
     - `exp`: ガウス関数指数 $\alpha$ (例: `b"7.483218"`)

#### (3) `C` (分子軌道係数行列)
- **パス**: `C/{run_type}_{iteration}` (例: `C/rks_49`)
- **データ型**: `float64`
- **形状**: `(N_AO, N_MO)` (2次元配列)
- **属性**:
  - `type`: `"GE"` (General Matrix: 一般行列)
  - `row`: 行数 (= $N_{\text{AO}}$)
  - `col`: 列数 (= $N_{\text{MO}}$)
  - `iteration`: イテレーション番号
  - `name`: `"C"`
- **行列の配置・並び順**:
  - C-contiguous (行優先) で格納。
  - **行 $p$ ($0 \le p < N_{\text{AO}}$)**: 原子基底関数 (AO) インデックス
  - **列 $i$ ($0 \le i < N_{\text{MO}}$)**: 分子軌道 (MO) インデックス
  - すなわち `C[p, i]` が AO $p$ の MO $i$ に対する展開係数 $C_{pi}$ である (実データ `H2O_HF` の $1s$ 軌道成分およびノルムで確認済み)。

#### (4) `energy_level` (軌道エネルギー)
- **パス**: `energy_level/{run_type}_{iteration}` (例: `energy_level/rks_49`)
- **データ型**: `float64`
- **形状**: `(N_MO,)` (1次元ベクトル)
- **属性**: `dim = N_MO`, `name = "energy_level"`
- **単位**: Hartree (a.u.)
- 例 (`H2O_HF` 最終イテレーション49):
  - 軌道 0: -20.5603 Hartree (O 1s)
  - 軌道 4 (HOMO): -0.5113 Hartree
  - 軌道 5 (LUMO): +0.1983 Hartree

#### (5) `occ` (軌道占有数)
- **パス**: `occ/{run_type}` (例: `occ/rks`)
- **データ型**: `float64`
- **形状**: `(N_MO,)` (1次元ベクトル)
- **属性**: `dim = N_MO`, `name = "occ"`
- 例 (`H2O_HF`): `[2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0, ...]` (10電子系、上位5軌道が占有数2.0)

#### (6) `pop/mulliken` (Mulliken原子電荷)
- **パス**: `pop/mulliken/{run_type}_{iteration}` (例: `pop/mulliken/rks_49`)
- **データ型**: `float64`
- **形状**: `(N_atom,)` (1次元ベクトル)
- **属性**: `dim = N_atom`
- **内容**: 各原子の正味のMulliken電荷 (net charge = $Z_A - \sum_{p \in A} (PS)_{pp}$)。
- 例 (`H2O_HF` 最終イテレーション49):
  - 原子0 (O): `-0.59339845`
  - 原子1 (H): `+0.29669921`
  - 原子2 (H): `+0.29669924`
  - (合計: `0.00000000`)

#### (7) `TEs` (全エネルギー推移)
- **パス**: `/TEs`
- **データ型**: `float64`
- **形状**: `(iterations,)` (1次元ベクトル)
- **属性**: `dim = iterations`
- **内容**: イテレーション 1 から iterations までの各イテレーションのSCF全エネルギー (Hartree)。

#### (8) `control/file_base_name` (外部ファイル名対応表)
- **パス**: `/control/file_base_name` (Group)
- **内容**: ProteinDFが作業ディレクトリ (`fl_Work/`) で入出力するバイナリファイル名のテンプレートが属性として記録されている。
  - `Spq_matrix = "Spq.mat"` (重なり行列 S)
  - `C_matrix = "C.%s.mat"` (MO係数行列 C)
  - `Ppq_matrix = "Ppq.%s.mat"` (密度行列 P)
  - `eigenvalues = "eigenvalues.%s.vtr"` (軌道エネルギー)
  - `occupation_vtr = "occupation.%s.vtr"` (占有数)
  - `pop/mulliken = "pop.mulliken.%s.vtr"` (Mulliken電荷)

#### (9) `s` (重なり行列: S行列)
- **パス**: `/s`
- **データ型**: `float64`
- **形状**: `(N_AO, N_AO)` (対称行列)
- **属性**: `type: "SY"`, `row: N_AO`, `col: N_AO`
- **内容**: 原子基底関数間の重なり積分行列 $S_{\mu\nu} = \langle \chi_\mu | \chi_\nu \rangle$。
- **実データでの確認**: `/tmp/test_compare/with_xc/pdfresults.h5` において shape=(14, 14), float64, 属性 `row=14, col=14, type="SY"` を確認。
- **保存処理**: `scripts/pdf-archive-h5.py` 128-133行目 `archive_s` (`entry.save_s_matrix(S)`)。

#### (10) `P` (密度行列)・`F` (Fock行列)
- **パス**: `P/{run_type}_{iteration}`, `F/{run_type}_{iteration}`
- **データ型**: `float64`
- **形状**: `(N_AO, N_AO)`
- **属性**: `type: "SY"`, `row: N_AO`, `col: N_AO`, `iteration: itr`
- **実データでの確認**: `/tmp/test_compare/with_xc/pdfresults.h5` に全反復 (1..27) 分が存在。
- **保存処理**: `scripts/pdf-archive-h5.py` 162-170行目 `archive_p`, 172-180行目 `archive_f`。

#### (11) `h`, `h2` (1電子・2電子ハミルトニアン行列)
- **パス**: `/h`, `/h2`
- **データ型**: `float64`
- **形状**: `(N_AO, N_AO)`
- **属性**: `type: "SY"`, `row: N_AO`, `col: N_AO`
- **実データでの確認**: `/tmp/test_compare/with_xc/pdfresults.h5` において shape=(14, 14), float64 を確認。
- **保存処理**: `scripts/pdf-archive-h5.py` 136-141行目 `archive_h`, 144-149行目 `archive_h2`。

---

## 2. 必要なデータの所在 (HDF5内か別ファイルか)

`reference.h5`(ProteinDF_testの比較用の参照ファイル)と、`pdf-archive-h5.py`が作る`pdfresults.h5`(標準計算結果アーカイブ)は、同じ形式の木構造であるが、含まれるデータセットが異なる。後者には `s`・`h`・`h2`・`P`・`F` が常に保存され、前者には保存されていない。

### `reference.h5` と `pdfresults.h5` のデータセット比較表

| データ項目 | HDF5パス / データセット | テスト参照ファイル `reference.h5` | 標準アーカイブ `pdfresults.h5` | `pdf-archive-h5.py` 該当行と処理 | 外部バイナリファイル名 (`control/file_base_name`) | 備考 |
| --- | --- | --- | --- | --- | --- | --- |
| **分子構造 (原子・座標)** | `/molecule/xyzc` | **あり** | **あり** | 65行目: `entry.save_model()` | なし | 全HDF5に存在 (Bohr単位) |
| **基底関数系** | `/basisset/{cgto, name, pgto}` | **あり** | **あり** | 65行目: `entry.save_model()` | なし | 全HDF5に存在 |
| **外部ファイル名対応** | `/control/file_base_name` | **あり** | **あり** | 65行目: `entry.save_model()` | なし | 全HDF5に存在 |
| **全エネルギー推移** | `/TEs` | **あり** | **あり** | 68, 106-113行目: `archive_TEs` (`entry.save_TEs`) | なし | 全SCF反復の全エネルギー |
| **重なり行列 (S)** | `/s` | **なし** (存在せず) | **あり** (14×14) | 71, 128-133行目: `archive_s` (`entry.save_s_matrix(S)`) | `Spq_matrix = "Spq.mat"` | `pdf-archive-h5.py` は常にSを保存 |
| **1電子ハミルトニアン (h)** | `/h` | **なし** | **あり** (14×14) | 74, 136-141行目: `archive_h` (`entry.save_h_matrix(h)`) | 未確認 | 対称行列 (SY) |
| **2電子ハミルトニアン (h2)** | `/h2` | **なし** | **あり** (14×14) | 75, 144-149行目: `archive_h2` (`entry.save_h2_matrix(h2)`) | 未確認 | 対称行列 (SY) |
| **軌道占有数** | `/occ/{run_type}` | **あり** | **あり** | 81, 96-103行目: `archive_occ` (`entry.save_occ`) | `occupation_vtr = "occupation.%s.vtr"` | RKSおよび一部UKS |
| **軌道エネルギー** | `/energy_level/{run_type}_{itr}` | **あり** (全反復) | **あり** (全反復) | 84, 116-126行目: `archive_energy_level` (`entry.save_energy_level`) | `eigenvalues = "eigenvalues.%s.vtr"` | 1..iterations反復分 |
| **MO係数行列 (C)** | `/C/{run_type}_{itr}` | **あり** (最終反復のみ) | **あり** (全反復 1..iterations) | 87, 152-160行目: `archive_c` (`entry.save_c_matrix`) | `C_matrix = "C.%s.mat"` | 一般行列 (GE) |
| **密度行列 (P)** | `/P/{run_type}_{itr}` | UKSのみ最終反復あり (RKSなし) | **あり** (全反復 1..iterations) | 90, 162-170行目: `archive_p` (`entry.save_density_matrix`) | `Ppq_matrix = "Ppq.%s.mat"` | 対称行列 (SY) |
| **Fock行列 (F)** | `/F/{run_type}_{itr}` | **なし** | **あり** (全反復 1..iterations) | 93, 172-180行目: `archive_f` (`entry.save_f_matrix`) | 未確認 | 対称行列 (SY) |
| **Mulliken原子電荷** | `/pop/mulliken/{run_type}_{itr}` | **あり** (RKS最終反復のみ) | **なし** (保存処理なし) | なし (`archive` 内に呼び出しなし) | `pop/mulliken = "pop.mulliken.%s.vtr"` | `reference.h5` はテスト用基準値として保持 |

### 外部バイナリファイルの形式 (`matrix.py`, `vector.py`)

ProteinDF本体が出力し、`ProteinDF_pytools` が読み書きするバイナリファイル（`.mat`, `.vtr`）の形式仕様は以下の通りである。

#### (1) バイナリ行列 (`.mat`: `matrix.py` 30-290行目)
- **ヘッダー (12バイト前後)**:
  - 構造: `(matrix_type, row, col)`
  - エンディアン: ネイティブ `=`、リトル `<`、ビッグ `>` (先頭バイトから自動判定)
  - C++出力の標準フォーマット: `=bii` (1バイトの種別 + 4バイト整数2つ = 9バイト) または `=iii` (4バイト整数3つ = 12バイト)
  - `matrix_type`:
    - `0` (RSFD: Row-Major 一般行列): ボディサイズ $= \text{row} \times \text{col} \times 8$ バイト
    - `1` (CSFD: Column-Major 一般行列): ボディサイズ $= \text{row} \times \text{col} \times 8$ バイト
    - `2` (SY: 対称行列): ボディサイズ $= \text{row} \times (\text{row} + 1) / 2 \times 8$ バイト
- **ボディ**:
  - 倍精度浮動小数点数 (`f64`, 8バイト) の連続。
  - SY (対称行列、重なり行列 $S$ など) の場合: **下三角成分**が `r in 0..dim`, `c in 0..=r` の順序で格納されている (`matrix.py` 262-265行目)。
  - ファイル全体のサイズ検証: `filesize == size_of_header + size_of_data` が厳密に成立することを検証してヘッダー形式を特定する (`matrix.py` 68-84行目)。

#### (2) バイナリベクトル (`.vtr`: `vector.py` 30-121行目)
- **ヘッダー (4または8バイト)**:
  - 構造: `(size,)` (`=i` または `=l`)
  - サイズ検証: `filesize == size_of_header + size * 8`
- **ボディ**:
  - `size` 個の倍精度浮動小数点数 (`f64`, 8バイト) の連続。

---

## 3. 軌道ごとのMulliken寄与に必要なデータの充足性判断

### 3.1 計算に必要なデータ要件

RUST_PORT_SPEC.md §3.5:
> YUIのMOブラウザ機能向けに、MO係数ベクトルと重なり行列(S)から軌道ごとの原子別Mulliken寄与を計算する関数を `qc_result/population.rs` に新設し、`OrbitalContribution { orbital_index, atom_contributions: Vec<(AtomId, f32)> }` として出力する。

Mullikenポピュレーション解析の定式化:
- MO $i$ の係数ベクトル: $\mathbf{C}_i = (C_{1i}, C_{2i}, \dots, C_{N_{\text{AO}} i})^T$
- 重なり行列: $\mathbf{S}$ ($N_{\text{AO}} \times N_{\text{AO}}$ の実対称行列)
- 軌道 $i$ に対する基底関数 $p$ の寄与 (gross population):
  $$q_p^{(i)} = C_{pi} (\mathbf{S} \mathbf{C}_i)_p = C_{pi} \sum_{q=1}^{N_{\text{AO}}} S_{pq} C_{qi}$$
- 原子 $A$ に属する基底関数の集合を $\mathcal{B}_A$ とすると、軌道 $i$ の原子 $A$ へのMulliken寄与:
  $$Q_A^{(i)} = \sum_{p \in \mathcal{B}_A} q_p^{(i)} = \sum_{p \in \mathcal{B}_A} C_{pi} (\mathbf{S} \mathbf{C}_i)_p$$
- MOの規格化条件より、全原子についての寄与の和は 1 となる: $\sum_A Q_A^{(i)} = \mathbf{C}_i^T \mathbf{S} \mathbf{C}_i = 1$。

したがって、計算には以下の3要素が必要である:
1. **MO係数行列 $C$** ($N_{\text{AO}} \times N_{\text{MO}}$)
2. **重なり行列 $S$** ($N_{\text{AO}} \times N_{\text{AO}}$ 対称行列)
3. **基底関数 (AO) と原子の対応付け** ($p \mapsto A$)

### 3.2 充足性判断と根拠

**判断**:
- **テストデータ (`reference.h5`) 単体**: 重なり行列 $S$ が含まれていないためそろわない。
- **標準計算結果アーカイブ (`pdfresults.h5`)**: `pdf-archive-h5.py` が生成するHDF5であれば、重なり行列 $S$ (`/s`)、MO係数行列 $C$ (`/C/...`)、基底関数系 (`/basisset/...`)、分子構造 (`/molecule/xyzc`) がすべて含まれているため、**完全にそろう**。

**根拠**:
1. **MO係数行列 $C$**:
   - `reference.h5` では `H2O_HF`, `Ne_HF`, `N2_RHF` などのRKSテストデータにおいて、`C/rks_{iterations}` として最終反復分がそろう。
   - `pdfresults.h5` では全反復 (`rks_1`〜`rks_{iterations}`) が `C` 配下に保存され、完全にそろう (`pdf-archive-h5.py` 152-160行目 `archive_c`)。
2. **重なり行列 $S$**:
   - `reference.h5` (ProteinDF_test の回帰テスト用基準値) にはデータセット `s` が保存されていない（サイズ節約等のため除外）。
   - 一方、ProteinDFの計算結果をHDF5にまとめる標準スクリプト `pdf-archive-h5.py` は、`archive_s`(128-133行目、`logger.debug("save S")` → `filepath = pdfparam.get_s_mat_path()` → `S = pdf.SymmetricMatrix(); S.load(filepath)` → `entry.save_s_matrix(S)`) により、**常に重なり行列をデータセット `s` として書き出す**。
   - 実データ (`/tmp/test_compare/with_xc/pdfresults.h5`、Ne B3LYP計算) においても、データセット `s` (shape=(14, 14), float64, type="SY") が存在することを確認済み。
3. **基底関数と原子の対応付け ($p \mapsto A$)**:
   - `molecule/xyzc` と `basisset` (`name`, `cgto`, `pgto`) から **完全に導出可能**。
   - ProteinDFの規約 (`orbinfo.py` 71-97行目 `OrbInfo._setup_by_param`):
     - 各原子 $A$ (0-index) について、`basisset/name` で指定された基底セット内の各 CGTO を走査。
     - 各 CGTO の `shell_type` (`s`, `p`, `d`) に応じて、$2 \times l + 1$ 個の基底関数 (s=1, p=3, d=5) を順番に展開して AO インデックスを割り振る。
     - これにより、各 AO $p \in \{0, \dots, N_{\text{AO}}-1\}$ がどの原子 $A$ に属するかが一意に確定する。

---

## 4. PR#60・PR#61の範囲の案

調査結果を踏まえ、PR#60 および PR#61 のスコープ案を以下のように提案する。

### PR#60: ProteinDF HDF5計算結果の読み込み (スコープ案)

1. **機能フラグ**:
   - コアクレート (`rust/crates/proteindf-bridge`) に `hdf5` オプション機能を新設 (`Cargo.toml`)。依存は `hdf5-metno`。
   - デフォルトビルドおよび wasm32 ターゲットでは無効のまま維持。
2. **データ構造体の設計 (`qc_result/`)**:
   - `QcResult`: 計算結果の統合コンテナ (ルート属性・分子・基底関数・軌道情報・電荷)。
   - `QcMolecule`: 原子のリスト、座標 (Å換算およびBohr保持)、核電荷。`molecule/xyzc` から生成。`AtomGroup` への変換メソッドを提供。
   - `QcBasisSet`: 各原子の基底関数系定義 (`cgto`, `pgto`)。
   - `OrbitalInfo`: AO総数、MO総数、AO-原子対応表、各AOの殻種別 (s, p, d) の管理。
   - `OrbitalEnergy`: 軌道エネルギーベクトル (`energy_level`)、占有数ベクトル (`occ`)、HOMO/LUMOインデックス。
   - `MoCoefficients`: MO係数行列 $C$ (行優先または列優先アクセス可能な2次元行列ラッパー)。
   - `MullikenPopulation`: 原子別Mulliken電荷ベクトル (`pop/mulliken`)。
3. **HDF5ローダー実装 (`qc_result/pdfparam_hdf5.rs`)**:
   - `H2O_HF_reference.h5`, `Ne_HF_reference.h5`, `N2_UHF_reference.h5` を読み込めるローダー。
   - 文字列データセット (`|S32`, `object`) の柔軟なデコード (固定長バイト配列からの文字列パース、数値文字列の `parse::<f64>()` 変換)。
   - 存在しないデータセット (`N2_UHF` での `C` や `occ` の欠損) を安全に `Option<T>` として扱う防御的設計。
4. **テスト**:
   - `tests/test_qc_hdf5.rs`: 今回コピーした3つのテストフィクスチャを用いた単体テスト。

### PR#61: 軌道別Mulliken寄与の計算 (スコープ案)

`pdf-archive-h5.py` が生成する標準HDF5では重なり行列 $S$ (`/s`) が常に提供されるため、S行列の入力はHDF5からの読み込みで足りる。外部バイナリファイル `Spq.mat` の読み込みは範囲から外す（必要になったら別途検討）。また、テストは合成のS行列ではなく、実際の計算結果で行う（テストデータとして `s` を含む多原子分子のHDF5をユーザーが用意する、ユーザー確認済み、2026-10-12）。

1. **計算コア (`qc_result/population.rs`)**:
   - `OrbitalContribution { orbital_index: usize, atom_contributions: Vec<(usize, f64)> }` の定義。
   - 計算関数: `calc_orbital_mulliken_contributions(mo_coeffs: &MoCoefficients, s_matrix: &SymmetricMatrix, orb_info: &OrbitalInfo) -> Result<Vec<OrbitalContribution>>`
   - 各軌道 $i$ について $q_p^{(i)} = C_{pi} (\mathbf{S}\mathbf{C}_i)_p$ を計算し、原子ごとに総和をとる。
   - 規格化検証: 各軌道の原子寄与総和が $1.0 \pm 10^{-5}$ になることをアサート。
2. **S行列の入力**:
   - HDF5データセット `/s` (対称行列) から重なり行列 $S$ をロードする（PR#60のHDF5ローダー機能を利用）。
   - 外部バイナリファイル `Spq.mat` のパーサは本PRの対象外とする。
3. **テスト戦略**:
   - ユーザーから提供される実計算結果（データセット `s` を含む多原子分子のHDF5ファイル）をテストフィクスチャとして使用する。
   - 実際の計算結果に対して、各軌道の原子寄与の和が 1.0 になること（MOの規格化条件）を検証する。

---

## 5. テストデータの選定と個人情報確認結果

### 5.1 選定ファイル

`ProteinDF_test/serial_dev/` から以下の3ファイルを選定し、`rust/crates/proteindf-bridge/tests/data/` にコピーした。

| コピー元 | コピー先ファイル名 | サイズ | 計算手法 | 特徴 |
| --- | --- | --- | --- | --- |
| `serial_dev/H2O_HF/reference.h5` | `H2O_HF_reference.h5` | 62,360 B (61 KB) | RKS (HF) | 3原子、24 AO、C行列・pop・エネルギー完備の多原子標準例 |
| `serial_dev/Ne_HF/reference.h5` | `Ne_HF_reference.h5` | 43,416 B (43 KB) | RKS (HF) | 1原子、14 AO、最小構成のRKS単原子例 |
| `serial_dev/N2_UHF/reference.h5` | `N2_UHF_reference.h5` | 20,784 B (21 KB) | UKS (HF) | 2原子、28 AO、最小構成のUKS例 (欠損グループへの耐性テスト用) |

合計サイズ: 126,560 B (約124 KB)。リポジトリを肥大化させないコンパクトな構成である。

### 5.2 個人データの確認方法と結果

コピー前に、リポジトリ外の一時環境 (`/tmp/survey_hdf5_venv`) において以下の2段階の検査を実施した。

1. **HDF5構造化データの走査 (`/tmp/check_privacy.py`)**:
   - 全属性名・属性値、および文字列データセット (`molecule/xyzc`, `basisset/*`) の全要素に対し、正規表現パターン (`/home/`, `/Users/`, ユーザー名 `toshi`, `hirano`, メールアドレスパターン, プライベートIP `192.168.*`, `10.*`, `127.0.0.1`, `localhost`) で走査。
   - **結果**: 3ファイルすべてにおいて検出件数 **0件**。
2. **バイナリ全体の文字列走査 (`strings` + `grep`)**:
   - ファイル全体のバイナリに対し `strings` を実行し、環境依存パス・ユーザー名・ホスト名等を検索。
   - **結果**: ローカルパス、ユーザー名、ホスト名、メールアドレス等の個人データは一切含まれていないことを確認 (**CLEAN**)。
