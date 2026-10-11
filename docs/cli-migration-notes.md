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
| `pdb2brd.py` | PR#53 | `Pdb`, `save_msgpack`, `AtomGroup.get_raw_data` | `Pdb`, `save_atomgroup` (または `save_msgpack`) | △ | `save_msgpack(ag.get_raw_data(), path)` は Rust版では `save_atomgroup(ag, path)` で1ステップ実行可能。`save_msgpack` 自体は未公開。 |
| `brd2pdb.py` | PR#53 | `load_msgpack`, `AtomGroup(dict)`, `Pdb.set_by_atomgroup` | `load_atomgroup`, `Pdb.set_by_atomgroup` | ○ | `load_msgpack` + `AtomGroup` は `load_atomgroup` で置換可能。 |
| `brd2gro.py` | PR#53 | `load_msgpack`, `AtomGroup(dict)`, `SimpleGro` | `load_atomgroup`, `SimpleGro` | ○ | `load_atomgroup` で読み込み、`SimpleGro` で出力可能。 |
| `gro2brd.py` | PR#53 | `SimpleGro`, `save_msgpack`, `AtomGroup.get_raw_data` | `SimpleGro`, `save_atomgroup` | △ | `SimpleGro.get_atomgroup()` を取得後、`save_atomgroup` で保存可能。 |
| `brd2xyz.py` | PR#53 | `load_msgpack`, `AtomGroup(dict)`, `Xyz` | `load_atomgroup`, `Xyz` | ○ | `load_atomgroup` で読み込み、`Xyz.set_by_atomgroup` で出力可能。 |
| `xyz2brd.py` | PR#53 | `Xyz`, `save_msgpack`, `AtomGroup.get_raw_data` | `Xyz`, `save_atomgroup` | △ | `Xyz.get_atomgroup()` 後に `save_atomgroup` で保存可能。 |
| `brd2txt.py` | PR#53 | `load_msgpack`, `AtomGroup(dict)`, `AtomGroup.save_csv` | `load_atomgroup` | × | `AtomGroup.save_csv` がRust版APIに未公開。 |
| `mmcif2txt.py` | PR#53 | `SimpleMmcif`, `SimpleMmcif.save_msgpack` | `SimpleMmcif` | × | `SimpleMmcif.save_msgpack` がRust版APIに未公開。 |
| `mmcif2mol2.py` | PR#53 | `SimpleMmcif`, `SimpleMol2`, `SimpleMmcif.save_msgpack` | `SimpleMmcif`, `SimpleMol2` | △ | `save_msgpack`/`load_msgpack` スクラッチ機能以外は `SimpleMmcif`・`SimpleMol2` ともにRust版に揃っている。 |
| `read_amber_prmtop.py` | PR#53 | `AmberPrmtop`, `save_msgpack` | `AmberPrmtop` | × | ProteinDF内部座標辞書形式のMessagePack書き出しに `save_msgpack` を使用している。 |
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
5. **純Python版 `brd-setup-bond.py` のバグ解消と移行**
   - **純Python版の状態**: `proteindf_bridge/__init__.py` で `Bond` がインポートされていないため、純Python版の実行時に `AttributeError: module 'proteindf_bridge' has no attribute 'Bond'` となる。
   - **Rust版での移行方針**: Rust版では `AtomGroup.setup()` または `Bond.setup_heuristic()` が利用可能であり、PR#54にて正常に結合解決を行うCLIとして移行する。
