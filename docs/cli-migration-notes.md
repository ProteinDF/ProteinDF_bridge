# docs/cli-migration-notes.md

# CLI移行対応表 (純Python版 → Rust版 `proteindf_bridge.rs`)

本ドキュメントは、`RUST_PORT_SPEC.md` §4.5 および `docs/tasks/TASK_cli-rust-migration.md` に基づき、構造を扱う26本のCLIスクリプト（`proteindf_bridge/cli/`）が使用している純Python版APIと、それに対応するRust版API（`proteindf_bridge.rs`）の有無、移行対象PR（PR#53〜55）、および不足しているAPIを整理した対応表である。

---

## 1. 対象外ツール（構造を扱わない汎用ツール6本）

以下の6本は構造データを扱わない汎用ユーティリティであり、`RUST_PORT_SPEC.md` §4.5の方針決定に基づき、純Python版のままとし移行対象外とする。

1. `db2txt.py`: SQLite3 DBのテキストダンプ
2. `doctest_runner.py`: doctest実行ランナー
3. `module_inspect.py`: Pythonモジュールのインスペクト
4. `mpac2txt.py`: MessagePackバイナリのYAML風テキスト表示
5. `mpac2yml.py`: MessagePack → YAML 変換
6. `yml2mpac.py`: YAML → MessagePack 変換

---

## 2. 移行対象26本のAPI対応表

凡例:
- **移行PR**: PR#53（フォーマット変換10本）、PR#54（構造の操作・表示14本）、PR#55（モデリング2本）
- **Rust API有無**: ○（利用可能）、△（一部利用可能/代替あり）、×（不足）

| コマンド | 移行PR | 純Python版で使用しているAPI | 対応するRust版API (`proteindf_bridge.rs`) | Rust API有無 | 備考・差異 |
| :--- | :--- | :--- | :--- | :---: | :--- |
| `pdb2brd.py` | PR#53 (完了) | `Pdb`, `save_msgpack`, `AtomGroup.get_raw_data` | `Pdb`, `save_atomgroup` | ○ | `save_atomgroup` で1ステップ保存。移行完了。基準ファイルと完全一致。 |
| `brd2pdb.py` | PR#53 (完了) | `load_msgpack`, `AtomGroup(dict)`, `Pdb.set_by_atomgroup` | `load_atomgroup`, `Pdb.set_by_atomgroup` | ○ | 移行完了。意図した違い: B-factorデフォルト値（純Python版1.00 → Rust版0.00）。 |
| `brd2gro.py` | PR#53 (完了) | `load_msgpack`, `AtomGroup(dict)`, `SimpleGro` | `load_atomgroup`, `SimpleGro` | ○ | 移行完了。意図した違い: ボックスベクトル行の書式（GRO規格準拠のf10.5固定幅および末尾改行）。 |
| `gro2brd.py` | PR#53 (完了) | `SimpleGro`, `save_msgpack`, `AtomGroup.get_raw_data` | `SimpleGro`, `save_atomgroup` | ○ | `SimpleGro.get_atomgroup()` 後に `save_atomgroup` で保存。移行完了。基準ファイルと完全一致。 |
| `brd2xyz.py` | PR#53 (完了) | `load_msgpack`, `AtomGroup(dict)`, `Xyz` | `load_atomgroup`, `Xyz` | ○ | `load_atomgroup` で読み込み、`Xyz.set_by_atomgroup` で出力。移行完了。基準ファイルと完全一致。 |
| `xyz2brd.py` | PR#53 (完了) | `Xyz`, `save_msgpack`, `AtomGroup.get_raw_data` | `Xyz`, `save_atomgroup` | ○ | `Xyz.get_atomgroup()` 後に `save_atomgroup` で保存。移行完了。基準ファイルと完全一致。 |
| `brd2txt.py` | PR#53 (完了) | `load_msgpack`, `AtomGroup(dict)`, `AtomGroup.save_csv` | `load_atomgroup`, `AtomGroup.save_csv` | ○ | `AtomGroup.save_csv` をPyO3バインディングに追加（PR#53）。基準ファイルと完全一致。 |
| `mmcif2txt.py` | PR#53 (完了) | `SimpleMmcif`, `SimpleMmcif.save_msgpack` | `SimpleMmcif.data`, `msgpack` (Python) | ○ | `SimpleMmcif.data` を追加し、スクリプト内で `msgpack` を直接使用。意図した違い: 引用符トークナイズ不具合修正。 |
| `mmcif2mol2.py` | PR#53 (完了) | `SimpleMmcif`, `SimpleMol2`, `SimpleMmcif.save_msgpack` | `SimpleMmcif`, `SimpleMol2`, `msgpack` (Python) | ○ | スクラッチ保存に `SimpleMmcif.data` + `msgpack`、読込に `from_data` を使用。基準ファイルと完全一致。 |
| `read_amber_prmtop.py` | PR#53 (完了) | `AmberPrmtop`, `save_msgpack` | `AmberPrmtop`, `msgpack` (Python) | ○ | `AmberPrmtop` で読込後、座標辞書を `msgpack` で書き出し。基準ファイルと完全一致。 |
| `brd-box.py` | PR#54 | `load_msgpack`, `AtomGroup(dict)`, `AtomGroup.weight`, `AVOGADRO_CONST` | `load_atomgroup`, `AtomGroup.weight` | × | 定数 `AVOGADRO_CONST`（アボガドロ定数）が `proteindf_bridge.rs` に未公開。 |
| `brd-density.py` | PR#54 | `load_msgpack`, `AtomGroup(dict)`, `AtomGroup.weight`, `AVOGADRO_CONST` | `load_atomgroup`, `AtomGroup.weight` | × | 定数 `AVOGADRO_CONST` が `proteindf_bridge.rs` に未公開。 |
| `brd-formula.py` | PR#54 | `load_msgpack`, `AtomGroup(dict)`, `AtomGroup.formula` | `load_atomgroup`, `AtomGroup.formula` | ○ | `AtomGroup.formula` はRust版にも実装・公開済み。 |
| `brd-show-res.py` | PR#54 | `load_msgpack`, `AtomGroup(dict)`, `Format.is_models`, `Format.is_protein` | `load_atomgroup`, `Format.is_models`, `Format.is_protein` | ○ | `Format.is_models`, `is_protein` はRust版にも公開済み。 |
| `brd-show-bonds.py` | PR#54 | `load_msgpack`, `AtomGroup(dict)`, `AtomGroup.get_bond_list` | `load_atomgroup`, `AtomGroup.get_bond_list` | ○ | `AtomGroup.get_bond_list` はRust版にも実装・公開済み。 |
| `brd-setup-bond.py` | PR#54 | `load_msgpack`, `AtomGroup(dict)`, `Bond.setup`, `AtomGroup.get_dict_data`, `save_msgpack` | `load_atomgroup`, `AtomGroup.setup()` (または `Bond.setup_heuristic`) | △ | 純Python版では未エクスポートの `bridge.Bond` 参照でAttributeErrorになるバグがある。Rust版では `AtomGroup.setup()` 呼出し後に `save_atomgroup` で保存可能。 |
| `brd-select.py` | PR#54 | `load_msgpack`, `AtomGroup`, `Select_Path_wildcard`, `AtomGroup.select`, `save_msgpack` | `load_atomgroup`, `AtomGroup`, `Select_Path_wildcard`, `AtomGroup.select`, `save_atomgroup` | ○ | `Select_Path_wildcard`, `select` はRust版に実装・公開済み。 |
| `brd-select-path.py` | PR#54 | `load_msgpack`, `AtomGroup`, `Select_Path_simple`, `AtomGroup.select`, `save_msgpack` | `load_atomgroup`, `AtomGroup`, `Select_Path_simple`, `AtomGroup.select`, `save_atomgroup` | ○ | `Select_Path_simple`, `select` はRust版に実装・公開済み。 |
| `brd-divide.py` | PR#54 | `load_atomgroup`, `AtomGroup`, `save_atomgroup` | `load_atomgroup`, `AtomGroup`, `save_atomgroup` | ○ | すべてRust版に実装・公開済み。 |
| `brd-divide-mainchain.py` | PR#54 | `load_atomgroup`, `AtomGroup`, `save_atomgroup` | `load_atomgroup`, `AtomGroup`, `save_atomgroup` | ○ | すべてRust版に実装・公開済み。 |
| `brd-renumber-resid.py` | PR#54 | `load_atomgroup`, `AtomGroup`, `save_atomgroup` | `load_atomgroup`, `AtomGroup`, `save_atomgroup` | ○ | すべてRust版に実装・公開済み。 |
| `brd-restructure.py` | PR#54 | `load_atomgroup`, `AtomGroup.restructure`, `save_atomgroup` | `load_atomgroup`, `save_atomgroup` | × | `AtomGroup.restructure` がRust版APIに未公開。 |
| `reorder.py` | PR#54 | `load_atomgroup`, `AtomGroup`, `save_msgpack` | `load_atomgroup`, `AtomGroup`, `save_atomgroup` | ○ | `save_msgpack(ag.get_raw_data(), ...)` は `save_atomgroup(ag, ...)` に置換可能。 |
| `crystallize.py` | PR#54 | `load_atomgroup`, `AtomGroup`, `Position`, `AtomGroup.shift_by`, `save_msgpack` | `load_atomgroup`, `AtomGroup`, `Position`, `AtomGroup.shift_by`, `save_atomgroup` | ○ | `Position`, `shift_by` はRust版にも実装・公開済み。 |
| `neutralize.py` | PR#55 | `load_atomgroup`, `Neutralize.neutralized`, `save_msgpack` | `load_atomgroup`, `Neutralize.neutralized`, `save_atomgroup` | ○ | `Neutralize` はRust版に実装・公開済み。 |
| `superposer.py` | PR#55 | `load_atomgroup`, `Superposer`, `Superposer_quaternion` | `load_atomgroup`, `Superposer`, `Superposer_quaternion` | ○ | `Superposer`, `Superposer_quaternion` はRust版に実装・公開済み。 |

---

## 3. 足りないRust版APIの一覧

移行（PR#53〜55）を進めるにあたり、`proteindf_bridge.rs`（Pythonバインディング `rust/crates/proteindf-bridge-py`）への追加が必要となるAPIは以下の通りである。

1. **`proteindf_bridge.rs.AVOGADRO_CONST`**
   - **用途**: `brd-box.py`, `brd-density.py`
   - **内容**: アボガドロ定数（`6.02214129e23` または純Python版 `common.py` の定義値 `6.0221415e23`）。モジュールレベル定数として公開。
2. **`AtomGroup.save_csv(path)`**
   - **用途**: `brd2txt.py`
   - **内容**: 原子リストをCSV形式（`path,symbol,name,x,y,z,charge`）でファイル出力するメソッド。
3. **`AtomGroup.restructure(reference_ag, range_tol)`**
   - **用途**: `brd_restructure.py`
   - **内容**: 参照構造の階層ツリーに合わせてターゲット構造の原子を再構築・整列するメソッド。
4. **汎用MessagePack I/Oヘルパー (`load_msgpack` / `save_msgpack`) または 代替設計**
   - **用途**: `read_amber_prmtop.py`, `mmcif2txt.py`, `mmcif2mol2.py`
   - **内容**:
     - `read_amber_prmtop.py` は原子座標・電荷等のPython辞書を `.brd` 形式のMessagePackとして書き出している。
     - `mmcif2txt.py` / `mmcif2mol2.py` は `SimpleMmcif.save_msgpack` でmmCIFブロックの辞書データをキャッシュしている。
     - Python側の標準 `msgpack` パッケージをスクリプト内で使用するか、`proteindf_bridge.rs` にMessagePack保存関数・メソッドを追加するかはPR#53の実装時に判断する。
5. **純Python版 `brd-setup-bond.py` のバグとPR#54移行時の正しさの基準**
   - **純Python版の状態**: `proteindf_bridge/__init__.py` で `Bond` がインポートされていないため、純Python版の実行時に `AttributeError: module 'proteindf_bridge' has no attribute 'Bond'` となる。このため、移行前の基準として使える正常出力が存在しない（本PR#52ではこの例外送出をゴールデンテストとして固定している）。
   - **PR#54での移行方針と正しさの基準**:
     - Rust版では `AtomGroup.setup()`（CCDテンプレートDB + 共有結合半径ヒューリスティックによる統合結合解決、§3.13〜3.15）または `Bond.setup_heuristic()`（純共有結合半径ヒューリスティックのみ、§3.10）が既に提供されている。
     - PR#54で `brd-setup-bond.py` をRust版に切り替える際は、`AtomGroup.setup()` を呼び出して結合を解決し、`save_atomgroup` で保存する設計とする。
     - **実在するテストに基づく正しさの検証基準値**:
       1. `rust/crates/proteindf-bridge/tests/test_ccd_templates.rs` の `test_atomgroup_setup_1hls_real_pdb` (710行目〜):
          実PDBフィクスチャ `1hls.pdb` に対し、`ag.setup()` 実行後にCCDテンプレートDB由来の二重結合が付与されること:
          - GLU 4 の主鎖 `C=O` 結合がちょうど1本で次数2であること:
            `assert_eq!(glu4_co_bonds.len(), 1, "GLU 4 C=O bond should appear exactly once");`
            `assert_eq!(glu4_co_bonds[0].order, 2, "GLU 4 C=O bond should have order 2 from CCD template");`
          - 鎖Bの ARG 22 に次数2の結合（`C=O` および `CZ=NH2`）が2本以上存在すること:
            `assert!(arg_double_bonds >= 2, "ARG 22 should have at least 2 double bonds (C=O and CZ=NH2)");`
          - 残基間ペプチド結合（次数1）がヒューリスティックにより検出され、全体として重複結合が存在しないこと:
            `assert!(seen.insert(key), "Duplicate bond found: {:?}", key);`
       2. `rust/crates/proteindf-bridge/tests/test_ccd_templates.rs` の `test_resolve_bonds_preserves_existing_file_bonds` (805行目〜):
          ファイル由来の事前登録結合（CONECTや明示的結合情報）が上書きされず保持され、重複しないこと:
          - 事前登録した次数1の `C-O` 結合がCCDの次数2で上書きされず維持されること:
            `assert_eq!(co_bonds.len(), 1, "C-O bond must be unique");`
            `assert_eq!(co_bonds[0].order, 1, "Pre-existing C-O bond order 1 must be strictly preserved");`
          - ALAの重原子結合4本が重複なく得られること:
            `assert_eq!(bonds.len(), 4, "Total bonds should be 4 without duplicates");`
       3. `rust/crates/proteindf-bridge/tests/test_covalent_bond_detection.rs` の `test_non_covalent_contact_excluded` (45行目〜):
          共有結合半径（Cordero et al. 2008）+ マージン（0.45Å）による幾何判定により、非結合の接触（2.5ÅのC-Cペアや2.8Åの水素結合ペア）が結合として誤検出されず除外されること:
          `assert_eq!(bonds.len(), 1, "Expected exactly 1 genuine covalent bond, but found {}", bonds.len());`
          `assert!(is_c3_c4, "The only detected bond should be C3-C4");`
       4. `rust/crates/proteindf-bridge/tests/test_covalent_bond_detection.rs` の `test_known_covalent_bonds_1hls` (91行目〜):
          `1hls.pdb` において残基間ペプチド結合（鎖A残基1-2間、鎖B残基1-2間）、残基内共有結合（GLU 4のCA-CB）、およびジスルフィド結合（鎖A CYS6 SG - CYS11 SG）が正しく検出されること:
          `assert!(peptide_a1_a2.is_some(), "Peptide bond between Chain A residue 1 C and residue 2 N should be detected");`
          `assert!(peptide_b1_b2.is_some(), "Peptide bond between Chain B residue 1 C and residue 2 N should be detected");`
          `assert!(glu4_ca_cb.is_some(), "CA-CB bond in GLU residue 4 should be detected");`
          `assert!(ssbond_a6_a11.is_some(), "Disulfide bond between Chain A CYS6 SG and CYS11 SG should be detected");`

---

## 4. PR#53（フォーマット変換10本）の移行結果と差分の分類・根拠

PR#53においてフォーマット変換10本（`pdb2brd.py`, `brd2pdb.py`, `brd2gro.py`, `gro2brd.py`, `brd2xyz.py`, `xyz2brd.py`, `brd2txt.py`, `mmcif2txt.py`, `mmcif2mol2.py`, `read_amber_prmtop.py`）を純Python版からRust版（`proteindf_bridge.rs`）へ移行した。

### 4.1 10本ごとの比較結果一覧

| コマンド | 移行結果 | 比較結果 | 違いの分類 | 根拠・内容 |
| :--- | :--- | :--- | :--- | :--- |
| `pdb2brd.py` | 移行完了 | **完全一致** | なし | `rs.save_atomgroup` による1ステップ保存。 |
| `brd2pdb.py` | 移行完了 | **違いあり** | **意図した違い** | ATOMレコードの温度因子（B-factor）未指定時デフォルト値。純Python版は読み込み時の既定値が `0.0`（`proteindf_bridge/biopdb.py` 216行目）なのに、`set_by_atomgroup` の経路でだけ `1.0`（同474行目の `setdefault("temp_factor", 1.0)`）になる不整合があり、Rust版は読み込み（`rust/crates/proteindf-bridge/src/format/pdb.rs` 266〜267行目）および `set_by_atomgroup`（同559〜563行目）の双方で既定値 `0.0` として一貫している。基準ファイルを更新。 |
| `brd2gro.py` | 移行完了 | **違いあり** | **意図した違い** | ボックスベクトル行（末尾行）の書式。純Python版はスペース区切りの浮動小数点（改行なし）を出力していたが、Rust版はGRO公式規格（GROMACS仕様）に準拠した固定幅 `{:10.5}` カラムおよび末尾改行を出力する。基準ファイルを更新。 |
| `gro2brd.py` | 移行完了 | **完全一致** | なし | 読み込み・`.brd` 出力ともに一致。 |
| `brd2xyz.py` | 移行完了 | **完全一致** | なし | 読み込み・`.xyz` 出力ともに一致。 |
| `xyz2brd.py` | 移行完了 | **完全一致** | なし | 読み込み・`.brd` 出力ともに一致。 |
| `brd2txt.py` | 移行完了 | **完全一致** | なし | 新設した `AtomGroup.save_csv` により列・書式（CRLF）・浮動小数点表記が純Python版と完全一致。 |
| `mmcif2txt.py` | 移行完了 | **違いあり** | **意図した違い (不具合修正)** | `_pdbx_chem_comp_identifier` 内のトークナイズ結果。純Python版パーサーは引用符付き空白文字列（例: `"(2S)-2-aminopropanoic acid"`）の分割に失敗しカラムが壊れていたが、Rust版の厳密なトークナイザにより正しくパースされるようになった。基準ファイルを更新。 |
| `mmcif2mol2.py` | 移行完了 | **完全一致** | なし | `.mol2` 出力およびスクラッチMessagePack保存・復元が一致。 |
| `read_amber_prmtop.py` | 移行完了 | **完全一致** | なし | PRMTOP読込、内部座標辞書変換、MessagePack出力すべて一致。 |

### 4.2 追加したAPI
1. **`AtomGroup.save_csv(path)`**:
   - `rust/crates/proteindf-bridge-py/src/atom_group.rs`
   - 一時ファイル出力後にアトミックリネームする安全設計。RFC 4180準拠のCRLF改行およびPythonの`repr(float)`互換の浮動小数点書式。
2. **`SimpleMmcif.data`**:
   - `rust/crates/proteindf-bridge-py/src/format/mmcif.rs`
   - mmCIFブロック辞書（`{block_name: [key_values, tables]}`）を返すプロパティ。
3. **`SimpleMmcif.from_data(dict)`**:
   - `rust/crates/proteindf-bridge-py/src/format/mmcif.rs`
   - mmCIFブロック辞書から `SimpleMmcif` を復元する静的メソッド。不正入力時は `BrInputError` を送出。

### 4.3 CLIオプションとAPIの扱いについての決定
- **`pdb2brd.py` の `-d/--debug` オプションと `Pdb.debug` プロパティの扱い**:
  純Python版の `Pdb` には `debug` 属性が存在したが、Rust版の `Pdb` コアにはデバッグフラグが存在しない。何もしないgetter/setterを公開APIに追加することは避け、`Pdb.debug` プロパティはバインディングに追加しない方針とした。`pdb2brd.py` 内の `pdb_obj.debug = debug` の代入行は削除したが、全PR共通ルール1（コマンド名・引数・オプションは変えない）に従い、CLIの `-d/--debug` 引数定義はそのまま維持している。
