# proteindf-bridge — Rust移植仕様書 (RUST_PORT_SPEC.md)

本ドキュメントは、本リポジトリ（`ProteinDF_bridge`, Python, GPLv3）をRustに1:1移植した新プロジェクト **`proteindf-bridge`** の仕様を定義する。既存の [`SPEC.md`](./SPEC.md) は現行Python実装の仕様書であり、本ドキュメントはそれを"正"としてRust版が満たすべき要件・拡張範囲・新規実装項目を規定する立場にある。

> **命名について(2026-09-15)**: 当初は `pdf-bridge` という名称を検討していたが、「PDF」(Portable Document Format)との混同を避けるため `proteindf-bridge` に変更した。Pythonバインディングのパッケージ名は既存の純Python版 `proteindf_bridge` と共存できるよう `proteindf_bridge_rs` とする(§4参照)。

## 0. 背景・目的

- タンパク質正準分子軌道計算プログラム群（ProteinDF/QCLObot）の可視化・モデリング研究プラットフォーム「結 (YUI)」が、構造データI/O・モデリング前処理基盤として本ライブラリのRust版を利用する（YUI側の設計は `yui` リポジトリの `architecture.md` 10章・10.1・10.2・11章・12.2節を参照）。
- YUI専用にせず、**C/C++・Pythonからも利用可能な多言語ライブラリ**として設計し、既存 `ProteinDF_bridge`/`ProteinDF_pytools` ユーザーの移行コストを抑える。
- 既存Python実装は変更・破棄せず、当面は並行して保守する。Rust版が既存資産と同等以上の正当性を持つことを、既存テストスイートの移植によって担保する。

## 1. 全体方針

1. **1:1移植を基本方針とする。** 既存モジュールのクラス・メソッド・引数の意味をできる限りそのままRustの型・関数に対応させる。内部実装（例: `Vector`/`Matrix`を自前実装するか`nalgebra`等に委ねるか）は移植後に最適化してよいが、まず動作を一致させることを優先する。
2. **既存テストを移植の正当性チェックリストとする。** `tests/test_*.py` は1モジュール1ファイルの粒度で存在するため、対応するRustユニットテストを同粒度で作成し、同じ入出力を検証する。
3. **ネイティブ往復フォーマット（`.brd`）の互換性を維持する。** MessagePackベースのシリアライズは、YUI側の `core` が持つ `MessagePack + zstd` ヘッダー設計（`[Magic: "YUI\0"(4B)] + [Version(1B)] + [Compression Flag(1B)] + [Payload]`）と相互運用できるエンコーディング規約を採用する。

## 2. モジュール対応表（1:1移植範囲）

既存Pythonモジュール（`proteindf_bridge/*.py`）に対応するRustモジュールを以下のように設計する。

| Python (既存) | Rust (新規) | 備考 |
| --- | --- | --- |
| `error.py` (`BrError`等) | `error.rs` | `thiserror`等でエラー型を定義 |
| `periodictable.py` | `periodic_table.rs` | |
| `vector.py` | `vector.rs` | |
| `matrix.py` (`Matrix`, `SymmetricMatrix`) | `matrix.rs` | |
| `position.py` | `position.rs` | |
| `atom.py` | `atom.rs` | |
| `bond.py` | `bond.rs` | |
| `atomgroup.py` | `atom_group.rs` | 鎖→残基→原子の再帰木構造 |
| `select.py` (`Select_*`) | `selector.rs` | 6.2節参照。既存にない `And`/`Or`/`Not` ブール合成コンビネータを新設 |
| `aminoacid.py` | `amino_acid.rs` | |
| `ssbond.py` | `ssbond.rs` | ジスルフィド結合検出（SG-SG距離 < 2.31Å） |
| `ionpair.py` | `ion_pair.rs` | 塩橋検出（4.0Å未満） |
| `neutralize.py` | `neutralize.rs` | |
| `modeling.py` | `modeling.rs` | ACE/NME等の末端キャッピング |
| `superposer.py` | `superposer.rs` | Kabschアルゴリズムによる構造重ね合わせ |
| `superposer_quaternion.py` | `superposer_quaternion.rs` | 四元数法（実験的） |
| `xyz.py` | `format/xyz.rs` | |
| `gro.py` (`SimpleGro`) | `format/gro.rs` | |
| `mol2.py` (`SimpleMol2`) | `format/mol2.rs` | |
| `mmcif.py` (`SimpleMmcif`) | `format/mmcif.rs` | **堅牢化が必要**（3章参照） |
| `amber_prmtop.py` | `format/amber_prmtop.rs` | |
| `biopdb.py` (`Pdb`) | `format/pdb.rs` | |
| `functions.py`（YAML/MsgPack I/Oヘルパー） | `brd.rs` | ネイティブ `.brd` 往復フォーマット |

`dbmanager.py`/`mail.py`（DB・メール送信のインフラ機能）は可視化・構造I/Oと無関係なため、Rust移植のスコープ外とする。

## 3. 新規実装項目（既存Pythonライブラリに存在しない機能）

### 3.1 PDBx/mmCIFの堅牢化

レガシーPDB形式はカラム桁数制約により原子数99,999件・チェーンID表現等にハード上限があり、YUI側の「最大約100万原子」表示目標と両立しない。既存 `SimpleMmcif` は名前の通り簡易実装であるため、大規模構造を確実に往復できるレベルまで堅牢化した実装を `format/mmcif.rs` に新規作成する。**PDBx/mmCIFをRust版の主力フォーマットとする。**

### 3.2 二次構造推定（DSSP相当）

現行 `ProteinDF_bridge` には二次構造推定機能がない。Kabsch-Sanderの静電モデル（`E = q1*q2*(1/r(ON) + 1/r(CH) - 1/r(OH) - 1/r(CN)) * 332`、`E < -0.5 kcal/mol` で水素結合と判定）による主鎖 N-H...O=C 水素結合パターンからヘリックス/シートを推定する、オリジナルDSSPと同様のアルゴリズムを `secondary_structure.rs` に新規実装する。YUI側のリボン（カートゥーン）表示で使用する。この主鎖水素結合検出ロジックは3.3節の水素結合検出と共有する。

### 3.3 非共有結合相互作用の検出（`InteractionSet`）

分子構造から非共有結合的な相互作用を検出し、構造データとは独立した **`InteractionSet`**（`{ kind, atoms, distance, angle, donor/acceptor role }` のリスト）として出力する。MessagePack/YAMLで往復可能にする。

| 種別 | 実装方針 | モジュール |
| --- | --- | --- |
| ジスルフィド結合 | 既存 `ssbond.py` を1:1移植 | `ssbond.rs`（2章） |
| 塩橋 (Salt Bridge) | 既存 `ionpair.py` を1:1移植 | `ion_pair.rs`（2章） |
| 水素結合 (Hydrogen Bond) | **新規。** 主鎖分は3.2節のDSSPロジックを流用。側鎖分（Ser/Thr/Tyr水酸基、Asn/Gln/Hisアミド/イミダゾール等）はドナー/アクセプター原子タイプ表を新設し、距離(D...A < 3.5Å目安)・角度(D-H...A > 120°目安)で判定。水素なし構造向けの重原子簡易判定モードも用意。 | `hydrogen_bond.rs` |
| CH-π相互作用 | **新規。** 芳香環を持つ残基(PHE/TYR/TRP/HIS)の環構成原子テーブルを新設し、環の重心・法線ベクトルを算出。C-H側原子（または水素）と環重心の距離・環法線とのなす角を幾何基準とする。既定値（例: 距離4.5Å以内・角度40°以内）は**チューニング可能なパラメータ**として公開し、決め打ちにしない。 | `ch_pi.rs` |

近傍探索は全原子対の総当たり(O(n²))を避け、空間分割（BVH/Octree等）を用いる。

### 3.4 ProteinDF計算結果I/Oの統合（`ProteinDF_pytools`相当）

構造データ（本リポジトリ由来）とは別に、ProteinDFの計算結果（`pdfparam.h5`: 基底関数系・MO係数行列・軌道エネルギー・Mulliken電荷等のポピュレーション解析）を扱う姉妹プロジェクト **`ProteinDF_pytools`**（`orbinfo.py`, `matrix.py`/`vector.py`, `basisset.py`/`basis2.py`, `pdfparam_hdf5.py`, `poputils.py` 等）も `proteindf-bridge` に統合する。

- HDF5読み込みは `hdf5-metno` クレート（オリジナルの `hdf5` crateは保守停止のため）を使用する。
- モジュール: `qc_result/pdfparam_hdf5.rs`, `qc_result/basis_set.rs`, `qc_result/orbital_info.rs`, `qc_result/population.rs`。

### 3.5 軌道別Mulliken寄与計算（`OrbitalContribution`）

YUIのMOブラウザ機能（`architecture.md` 12.2節）向けに、MO係数ベクトルと重なり行列(S)から**軌道ごとの原子別Mulliken寄与**を計算する関数を `qc_result/population.rs` に新設し、`OrbitalContribution { orbital_index, atom_contributions: Vec<(AtomId, f32)> }` として出力する。

### 3.6 Gaussian CUBEファイルの読み込みパーサ

ProteinDF本体のC++ツール（`pdf-mkfld-dens`/`pdf-mkfld-mo`/`pdf-mkfld-esp`）が標準Gaussian CUBE形式で既にファイル出力するため、**生成処理（GTOのグリッド展開）は不要**。標準Gaussian CUBEファイルの読み込みパーサのみを `format/cube.rs` に新規実装すればよい。

### 3.7 QCLObotプレイブック連携のためのパス構文共有

`*.QCLO.yaml` の `brd_select: /model_1/A/6/` のようなフラグメント選択パス文字列は、2章の `selector.rs`（`Select_Path_wildcard`/`Select_PathRegex`相当）でそのまま解決できる構文とする。QCLObot側の変更は不要。

## 4. 多言語バインディング方針

- **Rust:** コアライブラリ本体。ネイティブクレートとしてYUIの `core`/`renderer-native` から直接利用する。
- **C/C++:** `cdylib` + `cbindgen` によるヘッダー生成でC ABIを公開する。
- **Python:** `PyO3` + `maturin` によるバインディングを提供し、既存 `ProteinDF_bridge`/`ProteinDF_pytools` ユーザーが最小コストで移行できるようにする（可能な限り既存Python APIの関数・クラス名を踏襲する）。パッケージ名は `proteindf_bridge_rs` とし、既存の純Python版 `proteindf_bridge` と共存インストールできるようにする。

## 5. ライセンス

GPLv3を継続する（本リポジトリと同一ライセンス）。

## 6. テスト移植方針

`tests/test_*.py`（`test_atom.py`, `test_atomgroup.py`, `test_select.py`, `test_mmcif.py`, `test_superposer.py` 等、全モジュール分）を、対応するRustの `#[cfg(test)]` モジュールとして1:1移植する。新規実装項目（3章）については、文献値や既知構造（例: ジスルフィド結合を持つ既知PDB構造）を用いたテストケースを新規に作成する。

## 7. 想定利用者

- 「結 (YUI)」可視化・モデリング研究プラットフォーム（`core`/`renderer-native` から依存）
- 既存 `ProteinDF`/`QCLObot` エコシステムのPython/C++ツール群（Pythonバインディング経由での段階的移行を想定）

## 8. 未決事項（今後この文書を更新して確定する）

- ~~Rustクレート/リポジトリの配置場所（新規別リポジトリとするか、本リポジトリ内に `rust/` として同居させるか）~~
  → **決定（2026-09-14）**: 本リポジトリ内に `rust/` として同居させる。Python版との1:1移植であることを踏まえ、テスト移植時の参照・差分レビューを同一リポジトリ内で完結させる。
- CH-π・水素結合検出の閾値パラメータのデフォルト値の最終決定（文献レビューが必要）
- `hdf5-metno` クレートの実運用検証（ProteinDFが生成するHDF5ファイルとの互換性確認）
- 内部数値計算（`Vector`/`Matrix`）を自前実装のまま保つか、`nalgebra`等の既存クレートに置き換えるかの最終判断（1:1移植完了後に検討）
