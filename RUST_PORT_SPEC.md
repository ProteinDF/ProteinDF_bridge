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
| `mol2.py` (`SimpleMol2`) | `format/mol2.rs` | 読み込み対応・明示的結合情報をパース（3.8節参照） |
| `mmcif.py` (`SimpleMmcif`) | `format/mmcif.rs` | **堅牢化が必要**（3章参照） |
| `amber_prmtop.py` | `format/amber_prmtop.rs` | 明示的結合情報をパース（3.8節参照） |
| `biopdb.py` (`Pdb`) | `format/pdb.rs` | SSBONDに加えCONECTレコードをパース（3.8節参照） |
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

### 3.8 結合情報の優先順位方針（ファイル由来結合 vs VDWヒューリスティック）

構造データが持つ明示的結合情報と、距離ベースの幾何学的結合推定（`Bond::setup()`）の優先順位について、以下の方針を確立する。

> **方針**: ファイルに明示的な結合情報があればそれを優先して使用し、`Bond::setup()`（VDW半径ヒューリスティック）は呼ばない。ファイルに結合情報が存在しない場合のみ、`Bond::setup()` にフォールバックする。

#### 各フォーマットの対応状況
TASK_PR31〜PR33の実装により、明示的結合情報を持つ主要フォーマットのローダーは、パース時に結合トポロジーを構築して`AtomGroup`に格納した状態で返す:
- **Tripos Mol2 (`format/mol2.rs`)**: `@<TRIPOS>BOND` セクションの結合ペア・結合次数をパースして`AtomGroup`に登録。
- **Amber PRMTOP (`format/amber_prmtop.rs`)**: `%FLAG BONDS_WITHOUT_HYDROGEN` および `%FLAG BONDS_INC_HYDROGEN` の座標オフセット配列から原子インデックスを算出して`AtomGroup`に登録。
- **PDB (`format/pdb.rs`)**: `SSBOND`（ジスルフィド結合）および `CONECT` レコードをパースして`AtomGroup`に登録。PDB仕様に基づく双方向冗長記述やSSBONDとの同一結合重複は自動的に排除（deduplication）される。

#### 呼び出し側（「結 (YUI)」等）の推奨利用パターン
各ローダーの `get_atomgroup()` は、明示的結合情報が存在する場合は設定済みの `AtomGroup` を返す。呼び出し側は以下のように結合情報の有無を判定し、存在しない場合のみ `Bond::setup()` を呼ぶ設計とする:

```rust
let mut ag = loader.get_atomgroup()?;

// 構造全体で結合情報が1件以上存在するか判定
// ※ PDB等の階層構造（root -> model -> chain...）を含め、
//    ag.get_bond_list().is_empty() で構造全体の結合有無を確実に判定できる。
if ag.get_bond_list().is_empty() {
    // 明示的結合情報がないフォーマット（例: XYZ、GRO、結合未定義のPDB等）のみ
    // VDW半径ヒューリスティックによる距離ベース結合推定にフォールバック
    Bond::setup(&mut ag)?;
}
```

この方針により、MOL2/PRMTOP/PDB由来の正確な結合トポロジーがヒューリスティック判定で上書き・二重定義されることを防止し、かつ結合情報を持たないフォーマットに対しても自動補完を提供する。

## 4. 多言語バインディング方針

- **Rust:** コアライブラリ本体。ネイティブクレートとしてYUIの `core`/`renderer-native` から直接利用する。
- **C/C++:** `cdylib` + `cbindgen` によるヘッダー生成でC ABIを公開する。
- **Python:** `PyO3` + `maturin` によるバインディングを提供し、既存 `ProteinDF_bridge`/`ProteinDF_pytools` ユーザーが最小コストで移行できるようにする（可能な限り既存Python APIの関数・クラス名を踏襲する）。パッケージ名は `proteindf_bridge_rs` とし、既存の純Python版 `proteindf_bridge` と共存インストールできるようにする。

### 4.1 クレート配布方式の指針 (PR#36)

現時点では crates.io への一般公開やプライベートレジストリでのバージョン管理は対応不要である。現状、「結 (YUI)」側はローカル開発環境における相対パス依存（`path = "../ProteinDF_bridge/rust/crates/proteindf-bridge"`）を前提にしている。

#### 現実的な次の選択肢: GitHub git依存
相対パス依存と crates.io 公開の中間の現実的な選択肢として、Cargoのgit依存機能が利用できる:
```toml
# リリース・特定タグを指定する場合
proteindf-bridge = { git = "https://github.com/<org>/ProteinDF_bridge", tag = "v0.1.0" }

# 特定のブランチやコミットを指定する場合
# proteindf-bridge = { git = "https://github.com/<org>/ProteinDF_bridge", branch = "main" }
```

**git依存のトレードオフ:**
- **メリット**: リポジトリ外部のプロジェクト（別PC環境やCI/CD）からも同一ディレクトリ配置を強制されることなくクレートを取得・ビルドできる。
- **デメリット・制約**:
  - crates.io のような SemVer に基づく柔軟なバージョン範囲解決やクレート単位の中央キャッシュ共有が効かず、指定したコミットやタグ単位での固定となる。
  - プライベートリポジトリの場合、ローカルビルド環境やCI runnerにおいてGitHubアクセストークン（PAT）やSSH鍵等のgit認証設定が必要となる。

**将来の移行判断基準:**
外部の共同研究者・サードパーティ利用者が増加した場合、またはマルチリポジトリ構成のCI/CDパイプラインにおいて相対パス/git依存の管理コストが増大した段階で、crates.io公開（パブリッククレート化）またはプライベートCargoレジストリ（Cloudsmith、JFrog等）の導入を再検討する。

### 4.2 Pythonバインディングの役割分担指針 (PR#36)

将来的に `proteindf-bridge-py` と YUI 独自の `core-py`（`yui` パッケージ）が共存し得るため、その役割分担と使い分け指針を以下のように定める。

- **`proteindf-bridge-py` (`proteindf_bridge_rs` パッケージ)**:
  - **対象・目的**: 既存の純Python版 `ProteinDF_bridge` / `ProteinDF_pytools` を利用しているユーザー向けの移行パス、およびスクリプトベースのバッチ解析・量子化学計算前処理パイプライン。
  - **責務**: 生体分子ファイルの高速パース/書き出し（PDB, mmCIF, Amber PRMTOP, MOL2等）、階層データモデル（`AtomGroup`）、結合判定、幾何重ね合わせ等、従来の `ProteinDF_bridge` の提供機能を高速化して提供する。既存Python APIの関数・クラス名を踏襲し、ユーザーが最小限の移行コストで高速化の恩恵を受けられるようにする。
- **`core-py` (`yui` パッケージ)**:
  - **対象・目的**: 「結 (YUI)」研究プラットフォーム向けの機能拡張、GUI/可視化連動、統合モデリングワークフロー。
  - **責務**: YUIプラットフォーム独自の機能（レンダリングシーン制御、ビューア状態同期、対話的操作イベントハンドリング、GUIプラグイン機能等）を提供する。
- **利用者の使い分け判断基準**:
  - 既存のPythonスクリプトや計算バッチ処理の高速化・移行が目的の場合は **`proteindf_bridge_rs`** を使用する。
  - YUIの可視化機能やUI・レンダラーと連動するアプリケーションやプラグインを開発する場合は **`yui`** を使用する。

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

## 9. 「結 (YUI)」からの要求リスト（2026-09-19、`yui`リポジトリ フェーズ6e調査より）

「結 (YUI)」が本クレートの`AtomGroup`/`Bond`を実際に統合する（`ROADMAP.md`フェーズ6e）にあたり、
YUI側の調査で見つかった、bridge側で対応してほしい項目。優先度順。

- **[高] protein schemaの形式化と検証ヘルパー (PR#28対応)**:
  タンパク質構造の標準パス階層規約を以下のように形式化する:
  ```text
  /model_N/chain_id/res_key/atom_key
  ```
  - **Level 0 (Root / Models)**: パス深さ0（例: `"/"`）。モデルグループ群を保持。直接の原子は不許可。
  - **Level 1 (Model)**: パス深さ1（例: `"/model_1/"`）。`is_model_level()`で判定。チェイングループ群を保持。直接の原子は不許可。
  - **Level 2 (Chain)**: パス深さ2（例: `"/model_1/A/"`）。`is_chain_level()`で判定。残基グループ群を保持。直接の原子は不許可。
  - **Level 3 (Residue)**: パス深さ3（例: `"/model_1/A/6/"`）。`is_residue_level()`で判定。直接の原子を保持。サブグループは不許可。
  - **Level 4 (Atom)**: パス深さ4（例: `"/model_1/A/6/CA"`）。葉ノード（原子）。

  **位置的判定と構造的判定の区別**:
  `AtomGroup`の`is_*_level()`は`path()`の深さに基づく「位置的」判定であり、`Format::is_chain`等の「構造的」判定（直下に原子がない・サブグループが要件を満たす）とは判定軸が異なる。
  規約違反データ（残基ラッパーなしでchain直下に配置されたHETATMや水分子等）では、パス深さはchainレベル（深さ2）のまま構造的判定が失敗するという乖離が生じる。
  この乖離および規約違反（非残基レベルの直接原子、残基内のサブグループ、深さ超過）を走査・検出するヘルパー`AtomGroup::validate_schema() -> Vec<SchemaViolation>`を提供する。
- **[高] ファイル由来の明示的な結合トポロジーの読み込み (PR#31〜34対応済み、3.8節参照)**:
  MOL2読み込み（PR#31）、PRMTOP `BONDS_*`パース（PR#32）、PDB `CONECT`レコードパース（PR#33）を
  実装し、各ローダーが明示的結合情報付きの`AtomGroup`を返すよう拡張した。また、ファイル由来結合を
  優先し、結合情報がない場合のみ`Bond::setup()`（VDW半径ヒューリスティック）にフォールバックする
  利用優先順位を3.8節に確立した。
- **[高] `Bond::setup()`のスケーラビリティ (PR#35対応済み)**: 従来は全原子ペアのO(n²)距離行列
  （`SymmetricMatrix`）だったが、`spatial.rs`に一様セルリスト`CellList`を新設し、`Bond::setup()`を
  動的セルサイズによるO(N)近傍探索に置き換えた。100万原子規模の合成構造で約3.4秒での完走を実測済み。
  外部から効率的な近傍ペアクエリ（半径内の原子ペア列挙）を投げられる低レベルAPI
  （`CellList::for_each_neighbor_pair`/`query_pairs_within`）も公開した。
  2,000原子超では`distmat`/`bondmat`（従来のO(n²)密行列フィールド）は`None`になる
  （`MAX_DENSE_MATRIX_ATOMS`定数、メモリ保護のため）。
- **[高] 二次構造情報の`AtomGroup`への書き戻し**: 現状`calc_secondary_structure(chain: &AtomGroup) -> Vec<SecondaryStructure>`は結果を別のVecとして返すのみで、`AtomGroup`ツリー自体には反映されない
  （`AtomGroup`に汎用メタデータフィールドが無いため）。一方`Bond::setup()`は`mol.add_bond(...)`で
  結果を`AtomGroup`自体に書き戻す設計になっており、一貫していない。`bonds: Vec<BondRecord>`と
  同格の、生物学的に意味の明確な専用フィールド（例: 各residueレベルの`AtomGroup`が持つ
  `secondary_structure: Option<SsCode>`）を追加し、`calc_secondary_structure`と対になる
  `apply_secondary_structure(chain: &mut AtomGroup)`のような書き戻し関数を提供してほしい
  （汎用メタデータ袋ではなく、`bonds`と同じ「specific typed field」パターンを希望）。
  YUI側はこれが無い間、residueのpath文字列をキーとする一時的なサイドマップで代替する
  （フェーズ6e-ii、`atom_group.rs`/`selector.rs`と同様「bridge実装までの一時代替」と明記）。
- **[中] パスベース`BondRecord`の効率的な解決**: `BondRecord`の`atom1_path`/`atom2_path`が
  文字列パスのため、大規模構造でこれを原子への参照へ解決するコストを確認したい。
  パス文字列→原子への効率的なルックアップAPI（O(1)またはO(log n)）が既にあるか、
  なければ追加してほしい。
- **[中] wasm32ターゲット向けのデフォルト設定**: `cargo check --target wasm32-unknown-unknown`は
  `--no-default-features --features ruzstd`を指定すれば成功することを確認したが、これを消費側が
  毎回指定するのではなく、`Cargo.toml`側で`[target.'cfg(target_arch = "wasm32")'.dependencies]`を
  使い、wasm32ターゲットでは自動的に`ruzstd`が使われるよう構成してほしい（YUI自身の
  `core/Cargo.toml`が`zstd`に対して既に行っているのと同じパターン）。
- **[低・将来] クレート配布方式 (PR#36対応完了、4.1節参照)**: 現時点では対応不要。YUI側は相対パス依存を前提とし、次の現実的な選択肢としてGitHub git依存の指針・トレードオフを4.1節に明記した。
- **[低・将来] Pythonバインディングの名前空間整理 (PR#36対応完了、4.2節参照)**: `proteindf-bridge-py`（`proteindf_bridge_rs`）とYUI独自の`core-py`（`yui`）の責務と使い分け判断基準を4.2節に明記した。
- 内部数値計算（`Vector`/`Matrix`）を自前実装のまま保つか、`nalgebra`等の既存クレートに置き換えるかの最終判断（1:1移植完了後に検討）
