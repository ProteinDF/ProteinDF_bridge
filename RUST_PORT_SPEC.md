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
§3.14で確立した方針により、各ローダーの `get_atomgroup()` は、明示的結合情報が存在する場合はそれを尊重し、存在しない場合は内部で自動的に結合解決（CCDテンプレート→ヒューリスティック）まで完了させた状態で `AtomGroup` を返す。呼び出し側は追加の呼び出しなしにそのまま使ってよい:

```rust
let ag = loader.get_atomgroup()?;
// ag は既に結合解決済み（ファイル由来 > CCDテンプレート > ヒューリスティックの優先順位で）。
```

原子を手動で追加・変更した後など、明示的に結合を再解決したい場合は `ag.setup()?`（組み込みCCD DBを使用）または `ag.setup_with_db(&db)?`（§3.11の拡張DBを使用）を呼び出す。純粋な幾何ヒューリスティックのみを行いたい場合（CCDテンプレートを適用しない場合）は、低レベルAPIの `Bond::setup_heuristic(&mut ag)?` を直接呼び出すことも可能である（詳細は§3.14）。

この方針により、MOL2/PRMTOP/PDB由来の正確な結合トポロジーがヒューリスティック判定で上書き・二重定義されることを防止し、かつ結合情報を持たないフォーマットに対しても自動補完を提供する。

### 3.9 CCD結合テンプレートデータベースによる結合情報補完 (フェーズA完了: 2026-09-22)

#### 背景

§3.8の方針で、ファイルに明示的結合情報(CONECT/BONDS_*/@<TRIPOS>BOND等)があればそれを優先する仕組みは確立した。しかし、**標準アミノ酸・核酸・水などの「よくある」残基を含む通常のPDB/mmCIFファイルは、多くの場合これらの明示的結合情報を持たない**(レガシーPDBのCONECTは通常HETATM分にしか付与されず、標準残基間のペプチド結合や側鎖内結合はファイルに記載されない)。この場合、現状は`Bond::setup()`(VDW半径ヒューリスティック、距離のみで判定し**結合次数は常に1として登録される**)にフォールバックするしかなく、二重結合(C=O等)や芳香環の結合次数情報が失われる。

wwPDB Chemical Component Dictionary (CCD, https://www.wwpdb.org/data/ccd) は、標準・非標準を問わずほぼ全ての残基/リガンドについて、正準原子名・結合トポロジー・結合次数(SING/DOUB/TRIP/QUAD/AROM)を定義済みである(`format/mmcif.rs`のCCDパース対応は2026-09-22時点で実装済み、§3.8参照)。これを**残基名をキーとするテンプレートデータベース**として整備し、構造データの結合情報を補完・是正するのに使う。RDKit・OpenBabel・PyMOL等の主要ツールが採用している標準的な「テンプレートベース結合推定」手法である。

#### スコープ: フェーズA(結合次数補完)とフェーズB(水素付加)に分割する

**フェーズA(結合次数補完)を先に、独立したタスクとして着手する。** フェーズBは幾何学的に大幅に難易度が高いため、フェーズA完了後に改めて計画する(本節では概要のみ記載)。

##### フェーズA: 結合トポロジー・結合次数の補完 (完了: 2026-09-22)

1. **テンプレートデータの抽出・同梱方針**: wwPDBの完全なCCD配布ファイル(`components.cif`)は数百MB規模で全実行時に読み込むのは非現実的なため、**標準アミノ酸20種・標準核酸(DNA/RNA各4種)・水(HOH)** の計29種を対象に、wwPDBの最新データ(`https://files.rcsb.org/ligands/view/{comp_id}.cif`)から抽出スクリプト(`scripts/build_ccd_bond_templates.py`)により抽出を実施した。データ形式はMessagePackバイナリ(`rust/crates/proteindf-bridge/src/data/ccd_bond_templates.msgpack`, 約9.2KB)として格納した。
   **同梱方式は`include_bytes!`によるバイナリ直接埋め込み**を採用し、wasm32ターゲットやPyO3 wheel配布でのポータビリティを確保した。初回参照時に`std::sync::OnceLock`で遅延デシリアライズし、静的なルックアップテーブル(`CcdTemplateDb`)として保持する。
2. **テンプレート構造体の設計**: `ccd_templates.rs`を新設し、`CcdBondTemplate { comp_id: String, atoms: Vec<String>, bonds: Vec<(String, String, usize)> }`および`CcdTemplateDb`(`global()`, `lookup()`)を実装した。
3. **結合情報への適用**: `AtomGroup::apply_ccd_bond_templates(&mut self, db: &CcdTemplateDb)`を実装した。各residueレベルグループについて、その`name`(残基名)でテンプレートDBを引き、原子名の対応が取れる結合ペアについて結合次数を設定する。
   **優先順位の保護**:
   (1) ファイル由来の明示的結合(CONECT等、すでに存在する結合)
   (2) CCDテンプレートによる結合(残基内の正準結合・結合次数)
   (3) `Bond::setup()`の共有結合半径ヒューリスティック(残基間ペプチド結合・非標準構造向け)
   CCDテンプレート適用時だけでなく、その後の`Bond::setup()`呼び出し時にも既存の結合ペア(`mol.get_bond_list()`)を事前収集して二重追加をスキップする仕様とし、重複登録や結合次数の不正な上書きを完全に防止した(§3.12参照)。
4. **検証**: 実PDBフィクスチャ(`1hls.pdb`)を用いたテスト(`tests/test_ccd_templates.rs`)において、`Bond::setup()`単独では次数1にしかならないGLU側鎖(CD=OE1)や主鎖カルボニル(C=O)、ARGグアニジノ基(CZ=NH2)が正しく二重結合(order 2)として設定されること、および既存結合の上書き防止、未知残基の安全なスキップを確認した。さらに、その後に`Bond::setup()`を実行しても同一原子ペアのレコードが重複登録されないことを検証した。

##### フェーズB: 水素付加(将来、フェーズA完了後に別途計画)

CCDの理想化座標(`pdbx_model_Cartn_*_ideal`)は単体コンポーネントの計算幾何であり、実際の(非理想的な)実験構造の重原子配置にそのまま重ね合わせることはできない。水素付加には、重原子の混成状態(sp3/sp2/sp、CCDのトポロジー情報から導出可能)に応じた幾何学的なH配置計算が必要になる(`hydrogen_bond.rs`のPhase 8実装にある主鎖疑似H座標計算——直前残基のC・現残基のN/CA座標からH位置を幾何学的に算出する手法——が同種のアプローチの前例になる。側鎖版はこれよりバリエーションが多く難易度が高い)。プロトネーション状態(pHによる荷電残基の水素数の違い等)の扱いも別途検討が必要。フェーズAの完了後、実装方針・スコープをあらためてこのドキュメントに追記してから着手する。

#### スコープ外(当面)

- CCD全件(`components.cif`全体)のデータベース化(非標準リガンド・稀少修飾残基まで含む網羅対応は、必要が生じた時点で個別追加する)。
- 水素付加(フェーズB、別途計画)。
- プロトネーション状態・互変異性体の推定。

### 3.10 `Bond::setup()`のVDWヒューリスティックを共有結合半径ベースに変更 (完了: 2026-09-22)

#### 背景

`Bond::setup()`(`bond.rs`、Python版`bond.py`から1:1移植)の現行判定式は以下の通りであった:

> `distance(p, q) <= vdw(p) + vdw(q) + 0.4Å`

この`vdw()`はBondiのファンデルワールス半径(`periodic_table.rs`の`VDW`テーブル)であり、**本来「非結合の接触距離」を表す値であって、結合検出用ではない**。実際に数値を確認すると、例えばC-C原子ペアのカットオフは`1.70 + 1.70 + 0.4 = 3.8Å`になる。実際のC-C共有結合長は約1.5Å程度であり、この閾値は水素結合(~2.7〜3.5Å)や単なるVDW接触(側鎖パッキング等)まで「結合あり」と誤検出しうる範囲まで踏み込んでいた。

OpenBabel・RDKit・ASEの`natural_cutoffs`・Jmol/PyMOL等、主要な構造化学ツールは、結合検出には**共有結合半径(covalent radius)の和+小さめの許容値**を用いるのが標準である。§3.9のCCDテンプレートが効くのは標準残基(名前でルックアップできるもの)に限られるため、テンプレートが無い残基・非標準構造・リガンド全般で使われ続ける`Bond::setup()`自体の精度を上げることは、CCDテンプレートと独立に価値がある。

#### 実施内容 (2026-09-22完了)

1. **共有結合半径テーブルの追加**: `periodic_table.rs`に`COVALENT_RADIUS`テーブル(HからCmまでの96元素)を新設した。値はpymatgen・ASE等で標準参照されている信頼性の高い文献値(Cordero et al., "Covalent radii revisited", *Dalton Trans.*, 2008, 2832–2838, DOI: 10.1039/B801115J)を採用した。また、`PeriodicTable::covalent_radius`および`Atom::covalent_radius`APIを追加した。
2. **許容値(トレランス)の決定**: OpenBabel(`OBAtom::ConnectsTo`)の標準慣習に基づき、加算マージン`COVALENT_BOND_TOLERANCE = 0.45` Åを定数定義した。熱振動や実験誤差を許容しつつ、非結合のVDW接触や水素結合(>2.6Å)を明確に排除する。
3. **`Bond::setup()`の判定式変更**: `distance(p, q) <= covalent_radius(p) + covalent_radius(q) + COVALENT_BOND_TOLERANCE`に変更。PR#35の`CellList`によるO(N)探索の動的セルサイズ計算も`2 * max_cov + COVALENT_BOND_TOLERANCE`に更新した。VDW半径テーブル・`PeriodicTable::vdw()`は将来の接触判定用途のため保持している。
4. **Python版からの意図的な改善の明記**: 本変更はPython版`bond.py`からの意図的な改善であることをdocコメント・コミットに明記した。
5. **検証**:
   - `tests/test_covalent_bond_detection.rs`: 実PDB(`1hls.pdb`)での主鎖ペプチド結合(C-N)・残基内結合・ジスルフィド結合(S-S)の検出維持を検証。
   - 旧VDW式では結合と誤判定されていた2.5ÅのC-Cペアや2.8Åの水素結合ペアが「結合なし」と正しく除外されることを検証。
   - 100万原子ベンチマーク(`test_benchmark_1m_atoms`)が6.36秒で正常動作することを確認。

#### スコープ外(当面)

- 結合次数の判定(引き続き§3.9のCCDテンプレート、またはテンプレートが無い場合は次数1のまま)。
- VDW半径テーブル自体の削除・置き換え(共有結合半径テーブルを追加するのみ)。

### 3.11 CCD結合テンプレートDBの実行時拡張(ユーザー提供の外部CCDデータ) (完了: 2026-09-22)

#### 背景

§3.9フェーズAで実装した`CcdTemplateDb`は、標準アミノ酸20種・標準核酸8種・水の計29残基に固定された組み込みデータベースであり(`include_bytes!`でバイナリに埋め込み)、これはクレートが**常に**持っているべき最小限のデフォルトとして妥当な設計である(サイズが小さく、wasm32/PyO3配布でポータビリティを損なわない)。

一方、wwPDBのCCD全件(`components.cif`、数百MB規模、非標準リガンド・修飾残基・補酵素等を含む数万コンポーネント)のような**任意選択・大容量のデータ**まで同じ方式でバイナリに焼き込むのは悪手である(ほとんどのユーザーが使わないデータで全員のバイナリを肥大化させる、wwPDB側の更新への追従に再ビルドが必須になる等)。この種のデータは、**ユーザー自身が管理する外部ファイルとして、実行時に明示的にロードする方式**が適切(RDKit・OpenBabel等、他の主要ツールもCCD全件を配布物には同梱していない)。

既存の`format/mmcif.rs`の`SimpleMmcif`は、CCD全件ファイル(複数`data_`ブロック)を含め、CCD形式全般をパースできる機能を持っていた。本タスクにより、`SimpleMmcif`でパースした任意のCCDデータブロックから`CcdBondTemplate`を生成し、`CcdTemplateDb`に動的に登録・合成できる仕組みを整備した。

#### 実施内容 (2026-09-22完了)

1. **`SimpleMmcif`データブロックから`CcdBondTemplate`を組み立てる変換関数の新設**:
   - `CcdBondTemplate::from_mmcif_block(block: &MmcifDataBlock, comp_id: &str) -> Result<CcdBondTemplate>`を`ccd_templates.rs`に実装。
   - `format/mmcif.rs`から結合次数変換関数`parse_chem_comp_bond_order`を抽出し、`format/mmcif.rs`と`ccd_templates.rs`で共通利用(重複実装を排除)。
   - `block.has_atom_site()`が`true`の場合はマクロ分子構造データと判定し、適切なエラーを返却。
   - また、`SimpleMmcif::get_data_block`および`get_atomgroup`において、`data_`プレフィックスの有無にかかわらず柔軟にブロックを検索できるよう改善。
2. **`CcdTemplateDb`への動的登録・マージ手段の追加**:
   - `CcdTemplateDb::new()`: 空のテンプレートDBを生成。
   - `CcdTemplateDb::insert(&mut self, template: CcdBondTemplate) -> Option<CcdBondTemplate>`: 1件追加(同名存在時は上書きし旧値を返却)。
   - `CcdTemplateDb::merge(&mut self, other: &CcdTemplateDb)`: 複数DBを合成。同名重複時は`other`のエントリが優先される後勝ち(last-write-wins)仕様を明記。
   - 組み込みの`CcdTemplateDb::global()`は不変(`&'static`)として保持し、ユーザーは`CcdTemplateDb::default()`(組み込み29種のクローン)または`new()`を起点に拡張する。
3. **利用パターンのドキュメント化**:
   - `ccd_templates.rs`のモジュールレベルdocコメントに、外部CCDファイルのロードからテンプレート変換・DB登録・`AtomGroup`への適用までの一連のコード実例を記載(doctestでコンパイル検証済み)。
4. **検証**:
   - `tests/data/ALA.cif`から抽出したテンプレートが、組み込みDBの`ALA`エントリ(13原子・12結合・C=O二重結合)と完全一致することを検証。
   - 架空の合成リガンドCIF(`LIG`)を動的変換・登録し、`AtomGroup::apply_ccd_bond_templates`によって二重結合・単結合が正しく付与されることを検証。
   - `1HLS.cif`等のマクロ分子構造ブロックが`from_mmcif_block`で適切に拒絶されることを検証。
   - `merge`における後勝ち優先順位を検証。

#### スコープ外(当面)

- CCD全件の自動ダウンロード・キャッシュ機構(あくまでユーザーが自分でファイルを用意する前提。ネットワーク取得をクレートに組み込むことはしない)。
- 組み込みデフォルトDB(29残基)自体の拡張(§3.9フェーズAのスコープ、変更しない)。

### 3.12 `Bond::setup()`における既存結合の重複登録防止 (完了: 2026-09-22)

#### 背景・課題

§3.9でCCDテンプレートDBによる結合次数補完（C=Oの二重結合付与など）を導入し、続いて残基間ペプチド結合や非標準残基の補完として`Bond::setup()`を組み合わせるワークフローが確立された。
しかし、`AtomGroup::add_bond`が無条件にレコードを追加する仕様であったため、CCDテンプレート適用後に`Bond::setup()`を呼ぶと、すでにCCDで登録されている結合ペア（例: GLUのC=O、次数2）に対してヒューリスティックでも同一ペアが検出され、次数1の`BondRecord`が重複して追加されてしまう問題が判明した。

#### 実施内容 (2026-09-22完了)

1. **`Bond::setup()`内での既存結合スキップ**:
   - `Bond::setup(&mut self, mol: &mut AtomGroup)` の開始時に `mol.update_paths()` を呼び出し、`mol.get_bond_list()` から既存の原子パスペアを正規化（`atom1_path <= atom2_path`）した集合（`HashSet<(String, String)>`）を構築。
   - `CellList` による近傍探索で共有結合判定を満たしたペアであっても、すでに集合に存在するペアはスキップし、未結合のペアのみを `mol.add_bond(..., 1)` で追加するよう改修した。
   - これにより、一部の結合が既に登録されている状態（CCD適用後やCONECT等の一部のみが存在する場合）で `Bond::setup()` を呼んでも、不足分のみを補うヒューリスティックとして安全に機能するようになった。
2. **検証**:
   - `tests/test_ccd_templates.rs` の `test_apply_ccd_bond_templates_1hls_real_pdb` において、単なる `.find()` ではなく、該当原子ペア（GLU 4 C=O）の結合レコードがちょうど1件であり、かつ次数2を保持していること、および全体として重複結合が0件であることを厳密にアサート。
   - 新規回帰テスト `test_bond_setup_after_ccd_does_not_duplicate_bonds` を追加し、CCD適用後に `Bond::setup()` を実行しても結合数が重複せず、CCDテンプレートの結合次数が維持されることを検証。

### 3.13 結合解決の統合エントリポイント(ファイル由来 > CCDテンプレート > ヒューリスティック) (完了: 2026-09-22)

#### 背景

§3.8で確立した優先順位(ファイル由来結合 > CCDテンプレート > `Bond::setup()`ヒューリスティック)は、これまで呼び出し側が正しい順序で複数のAPIを個別に呼ぶ必要があった(`apply_ccd_bond_templates()` → `Bond::setup()`)。この2段階呼び出しを1回の呼び出しに統合し、`proteindf-bridge`を使うプログラム側の実装をシンプルにするために統合エントリポイントを整備した。

`Bond::setup()`自体を「CCDテンプレートも内部で使う」ように変更することは行わない。理由:
- `Bond::setup()`は「純粋な共有結合半径ヒューリスティックのみ」を行う関数としてドキュメント化・テストされており、その意味を保つ。
- `bond.rs`は基礎プリミティブであり、高レベル機能(CCDテンプレートDB)に依存させないレイヤリングを維持する。

代わりに、両方を正しい順序で組み合わせる統合メソッド `AtomGroup::resolve_bonds` を追加した。

#### 実施内容 (2026-09-22完了)

1. **`AtomGroup::resolve_bonds` の実装**:
   - `AtomGroup::resolve_bonds(&mut self, db: &CcdTemplateDb) -> Result<()>` を `atom_group.rs` に追加。
   - 内部で `self.apply_ccd_bond_templates(db)` を実行した後、`Bond::setup(self)` を呼び出す。
   - 事前にファイル由来の明示的結合（CONECT、MOL2、PRMTOP等）が存在する場合も、既存結合の上書き防止および重複排除ロジック（§3.12）により、最優先で保護される。
2. **`db` 引数による拡張性の担保**:
   - 呼び出し側は組み込みデフォルトDB（`CcdTemplateDb::global()`）をそのまま渡すことも、§3.11 で拡張・マージした独自DBを渡すことも可能。
3. **検証**:
   - `tests/test_ccd_templates.rs` に `test_resolve_bonds_1hls_real_pdb` を追加。実PDB（`1hls.pdb`）に対して1回の呼び出しでCCD由来の二重結合（GLU 4 C=O、ARG 22 CZ=NH2等）とヒューリスティックによる残基間ペプチド結合の両方が重複なく得られることを検証。
   - `test_resolve_bonds_preserves_existing_file_bonds` を追加。事前に登録されたファイル由来結合（CONECT等）が上書き・重複されず保持されることを検証。
   - 既存の `Bond::setup()` 単体・`apply_ccd_bond_templates()` 単体のテストが全て引き続きパスすることを確認。

### 3.14 `AtomGroup::setup()`を賢いデフォルトの結合解決エントリポイントにする(計画中、未着手)

#### 背景

§3.13で`AtomGroup::resolve_bonds(db)`を追加したが、呼び出し側から見ると「結合を解決するなら`Bond::setup()`を呼べばよい」という直感に反し、低レベルの`Bond::setup()`（純ヒューリスティックのみ）と高レベルの`AtomGroup::resolve_bonds()`（CCD+ヒューリスティック）のどちらを呼ぶべきかが分かりにくい、というユーザーからの指摘があった。

「`Bond::setup()`自体を賢くする」案（`bond.rs`が`ccd_templates.rs`に依存する形になり、§3.13で明示的に避けたレイヤリング崩壊を招く）ではなく、**「`setup`という名前を、高レベルAPIである`AtomGroup`側の賢いデフォルトに割り当てる」**方針を採用する。`bond.rs`は`ccd_templates.rs`に依存しないレイヤリングを維持する。

さらに、「ファイルを読み込んだ時点で自動的に結合解決まで完了していてほしい」という要望があった。ただし「`AtomGroup`の内容が更新されるたびに暗黙的に再実行する」方式は、以下の理由により採用しない:
- **性能**: 原子を1個ずつ追加しながら構築するケースで、追加のたびにCCD照合・近傍探索が走ると計算量が爆発する(N原子でO(N²)相当)。
- **正確性**: CCDテンプレート照合は残基内の全原子が揃っていることを前提とする。構築途中の不完全な状態で結合解決を走らせると、誤った結合が確定し、§3.12の「既存結合は上書きしない」方針により後から自動修正されなくなる。
- **実装コスト**: `AtomGroup`は再帰的な木構造であり、「更新」をどこで検知するか（`set_atom`、`add_group`、子グループの更新の伝播...）が広範囲に及ぶ。

代わりに、**「ファイルローダーの読み込み完了」という自然な区切りに限定して暗黙実行する**。手動で`AtomGroup`をゼロから構築するケース（テストコード等）では、構築完了後に明示的に`setup()`を呼ぶ、という形を維持する。

#### 対象

1. **`bond.rs`のリネーム**: `Bond::setup()`（純粋な共有結合半径ヒューリスティックのみ）を`Bond::setup_heuristic()`にリネームする。ロジックは変更しない。`proteindf-bridge-py`（Pythonバインディング）内の呼び出し箇所も追従させる。
2. **`AtomGroup`側APIの整備**（`atom_group.rs`）:
   - `AtomGroup::resolve_bonds(&mut self, db: &CcdTemplateDb) -> Result<()>` を `AtomGroup::setup_with_db(&mut self, db: &CcdTemplateDb) -> Result<()>` にリネームする（ロジックは変更しない。§3.11の拡張DBを使う場合の明示的エントリポイントとして残す）。
   - 新規に `AtomGroup::setup(&mut self) -> Result<()>` を追加する。内部で `self.setup_with_db(CcdTemplateDb::global())` を呼ぶだけの薄いラッパーとする。これが「呼べば良きに計らってくれる」デフォルトの公開APIになる。
3. **各フォーマットローダーでの暗黙実行**: `get_atomgroup()`（またはパース処理の完了直前）で、`ag.get_bond_list().is_empty()` の場合のみ `ag.setup()?` を呼んでから返すようにする。対象:
   - `format/pdb.rs`の`get_atomgroup()`
   - `format/mmcif.rs`の`get_atomgroup()`
   - `format/amber_prmtop.rs`の`get_atomgroup()`
   - `format/gro.rs`の`get_atomgroup()`
   - `format/mol2.rs`は`get_atomgroup()`が`&AtomGroup`（参照）を返す設計のため、`parse_str()`内で`self.set_by_atomgroup(&ag)`を呼ぶ直前に同様のチェックを行う。
   - いずれも「ファイル由来の結合が既にある場合は何もしない」ため、MOL2/PRMTOP/PDB(CONECT)等、明示的結合情報を持つフォーマットでは実質的にno-opとなり、§3.8の優先順位方針を壊さない。
4. **ドキュメント更新**: `RUST_PORT_SPEC.md` §3.8の呼び出し側推奨パターンを、本セクションの内容に合わせて更新済み（先行して反映済み）。各ローダーファイル内の`Bond::setup()`を参照するdocコメント（`mol2.rs`・`pdb.rs`・`amber_prmtop.rs`）も新名称に更新する。

#### 完了の定義(想定)

1. `Bond::setup()`への参照が名称`Bond::setup_heuristic()`に統一され、`cargo build --workspace`（Pythonバインディング含む）が通ること。
2. `AtomGroup::setup()`（引数なし）を呼ぶだけで、CCDテンプレート＋ヒューリスティックによる結合解決が行われることを検証する回帰テストを追加すること（既存の`test_resolve_bonds_*`系テストを新API名に追従させる形でよい）。
3. 各ローダー（PDB・mmCIF・PRMTOP・GRO・MOL2）について、明示的結合情報を持たない入力に対して`get_atomgroup()`を呼んだだけで結合が自動解決されていることを検証する回帰テストを、フォーマットごとに最低1件追加すること。
4. 明示的結合情報を持つ入力（PDBのCONECT、MOL2のBONDセクション、PRMTOPのBONDS等）に対して`get_atomgroup()`を呼んでも、ファイル由来の結合が上書き・重複されないことを確認する既存テストが引き続きパスすること。
5. `cargo clippy` / `cargo fmt` を通すこと。

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
- **[高] 二次構造情報の`AtomGroup`への書き戻し (PR#29対応済み)**: 現状`calc_secondary_structure(chain: &AtomGroup) -> Vec<SecondaryStructure>`は結果を別のVecとして返すのみで、`AtomGroup`ツリー自体には反映されない
  （`AtomGroup`に汎用メタデータフィールドが無いため）。一方`Bond::setup()`は`mol.add_bond(...)`で
  結果を`AtomGroup`自体に書き戻す設計になっており、一貫していない。`bonds: Vec<BondRecord>`と
  同格の、生物学的に意味の明確な専用フィールド（例: 各residueレベルの`AtomGroup`が持つ
  `secondary_structure: Option<SsCode>`）を追加し、`calc_secondary_structure`と対になる
  `apply_secondary_structure(chain: &mut AtomGroup)`のような書き戻し関数を提供してほしい
  （汎用メタデータ袋ではなく、`bonds`と同じ「specific typed field」パターンを希望）。
  YUI側はこれが無い間、residueのpath文字列をキーとする一時的なサイドマップで代替する
  （フェーズ6e-ii、`atom_group.rs`/`selector.rs`と同様「bridge実装までの一時代替」と明記）。
  → `secondary_structure: Option<SsCode>`フィールド（`bonds`と同様private、`secondary_structure()`/
  `set_secondary_structure()`経由でのみアクセス）と、`apply_secondary_structure(chain: &mut AtomGroup)`
  を実装済み。`merge`/`BitAnd`/`BitOr`/`BitXor`/`Clone`全てで正しくハンドリングされることをテストで検証済み。
- **[中] パスベース`BondRecord`の効率的な解決 (PR#30対応済み)**: `BondRecord`の`atom1_path`/`atom2_path`が
  文字列パスのため、大規模構造でこれを原子への参照へ解決するコストを確認したい。
  パス文字列→原子への効率的なルックアップAPI（O(1)またはO(log n)）が既にあるか、
  なければ追加してほしい。
  → 計測の結果、`get_atom_by_path`は既に階層深さのみに依存するO(depth)（実質O(1)）であることを実証
  （75,000原子まで探索時間が変化しないことをベンチマークで確認）。ゼロアロケーション最適化も実施。
  利便性のため`AtomGroup::resolve_bond(&self, record: &BondRecord) -> Option<(&Atom, &Atom)>`を追加した。
- **[中] wasm32ターゲット向けのデフォルト設定 (PR#27対応済み)**: `cargo check --target wasm32-unknown-unknown`は
  `--no-default-features --features ruzstd`を指定すれば成功することを確認したが、これを消費側が
  毎回指定するのではなく、`Cargo.toml`側で`[target.'cfg(target_arch = "wasm32")'.dependencies]`を
  使い、wasm32ターゲットでは自動的に`ruzstd`が使われるよう構成してほしい（YUI自身の
  `core/Cargo.toml`が`zstd`に対して既に行っているのと同じパターン）。
  → `[target.'cfg(target_arch = "wasm32")'.dependencies]`でwasm32では`ruzstd`が自動選択されるよう構成し、
  ネイティブターゲットでは従来通り`zstd`/`ruzstd`をfeatureで明示選択できる状態を維持した
  （バージョン指定は`[workspace.dependencies]`に一元化し重複を排除）。
- **[低・将来] クレート配布方式 (PR#36対応完了、4.1節参照)**: 現時点では対応不要。YUI側は相対パス依存を前提とし、次の現実的な選択肢としてGitHub git依存の指針・トレードオフを4.1節に明記した。
- **[低・将来] Pythonバインディングの名前空間整理 (PR#36対応完了、4.2節参照)**: `proteindf-bridge-py`（`proteindf_bridge_rs`）とYUI独自の`core-py`（`yui`）の責務と使い分け判断基準を4.2節に明記した。
- 内部数値計算（`Vector`/`Matrix`）を自前実装のまま保つか、`nalgebra`等の既存クレートに置き換えるかの最終判断（1:1移植完了後に検討）
