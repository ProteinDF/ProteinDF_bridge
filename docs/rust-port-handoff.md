# proteindf-bridge Rust移植 — antigravity向け作業指示

本ドキュメントは、`RUST_PORT_SPEC.md` に基づくRust移植作業を実装担当(antigravity)に委任するにあたっての、フェーズ単位の作業指示を記録する。Claude(このリポジトリでのレビュー担当)は各PRを本ドキュメント・`RUST_PORT_SPEC.md`・`SPEC.md` と突き合わせて仕様適合チェックを行う。

> **現在の運用ルールは本ファイル末尾の「運用ルールの変更(2026-09-15、2026.9.0リリース後): GitFlow運用への移行」を参照。** 以下の「運用ルール」節は`rust-port`統合ブランチ時代(Phase 1〜4)の記録として残しているが、現在は無効。

## 運用ルール(過去の記録、現在は無効)

- **正本**: 仕様は `RUST_PORT_SPEC.md`。移植元の正解データは `proteindf_bridge/*.py` と `tests/test_*.py`。
- **作業ブランチ**: 統合ブランチ `rust-port` を `docs/rust-port-spec`(本ドキュメントとRUST_PORT_SPEC.mdを含む)から作成し、そこにPRを積み上げる。Phaseが完了し、Claudeのレビューを通過した時点でまとめて `main` へマージする。Phase途中の未完成状態が `main` に混ざらないようにするため。
- **PRごとの機能ブランチ(MUST)**: `rust-port` へ直接コミットしない。`feature/phaseN-prM` のようなPR単位の機能ブランチを `rust-port` から切って作業する。
- **マージ前レビューゲート(MUST)**: 機能ブランチでの実装・テスト・`cargo clippy`/`cargo fmt` 確認が終わっても、**`rust-port` へは自分でマージしない**。ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。承認が出るまで次のPRの実装に着手しない。承認後にユーザーまたはClaudeが `rust-port` へマージする。
  - 背景: Phase 1(PR#1〜3)はこのルールに反して `rust-port` へ直接コミットされ、3PR分が未レビューのまま積み上がってしまった(2026-09-14に発覚)。同じ状態を繰り返さないための明文化。
- **PR粒度**: モジュール単位、または関連の強い2〜3モジュールをまとめた単位で小さく切る。
- **各PR説明に記載する項目**: (1) 対応するPythonモジュール、(2) 移植したテスト数/件数、(3) 意図的な差分があれば理由(例: Python側の型階層をRustのenumに統合した、Python側にない独自エラー種別を追加した、等)。「特になし」でもよいが、必ず言及すること。

### Cargo workspace構成(2026-09-14 決定)

```
rust/
├── Cargo.toml            # workspace定義
└── crates/
    └── proteindf-bridge/  # コアライブラリ(Phase 1の対象)
        ├── Cargo.toml
        └── src/
            ├── lib.rs
            ├── error.rs
            └── periodic_table.rs
```

単一クレート直下構成ではなく `crates/` 配下にコアクレートを置く構成を採用する。理由: §4のバインディング方針(C ABI用cdylib、PyO3用Pythonモジュール)は将来的にコアクレートとは別クレート(`crates/proteindf-bridge-capi/`、`crates/proteindf-bridge-py/` 等)として追加する想定のため、最初から `crates/` 構成にしておく。

### 命名変更(2026-09-15): `pdf-bridge` → `proteindf-bridge`

当初 `pdf-bridge` という名称を使っていたが、「PDF」(Portable Document Format)との混同を避けるため `proteindf-bridge` に変更した。Rustクレート名・ディレクトリ名(`rust/crates/pdf-bridge/` → `rust/crates/proteindf-bridge/`)・`Cargo.toml`のpackage名・全テストファイルの`use pdf_bridge::...`を`use proteindf_bridge::...`に、専用PRで対応すること(詳細は「Phase 4.5: プロジェクト名変更」参照、Phase 5着手前に完了させる)。Pythonバインディングのパッケージ名も`pdf_bridge`ではなく`proteindf_bridge_rs`とする。

## Phase 1: 基盤・データモデルの1:1移植(完了 2026-09-14)

`rust/` に Cargo workspace を新設し、`RUST_PORT_SPEC.md` §2の対応表のうち以下のみを移植する。

| Python | Rust | 推奨PR分割 |
| --- | --- | --- |
| `error.py` | `error.rs` | PR#1 |
| `periodictable.py` | `periodic_table.rs` | PR#1 |
| `vector.py` | `vector.rs` | PR#2 |
| `matrix.py`(`Matrix`, `SymmetricMatrix`) | `matrix.rs` | PR#2 |
| `position.py` | `position.rs` | PR#2 |
| `atom.py` | `atom.rs` | PR#3 |
| `bond.py` | `bond.rs` | PR#3 |
| `atomgroup.py` | `atom_group.rs` | PR#3 |

**スコープ外(今回は着手しない)**: フォーマットI/O(xyz/gro/mol2/mmcif/amber_prmtop/biopdb)、構造操作(select/aminoacid/ssbond/ionpair/neutralize/modeling/superposer系)、§3の新規機能(mmCIF堅牢化・二次構造推定・InteractionSet・QC結果I/O・CUBEパーサ)、§4の多言語バインディング。これらは次フェーズ以降で別途指示する。

## 完了の定義(Definition of Done)

1. 対応する `tests/test_*.py` を同粒度で `#[cfg(test)]` として移植し、全てpassすること。
2. 数値計算系(`vector`/`matrix`/`position`)はPython版と同一入力で同一出力になることをテストで担保すること(浮動小数点誤差の許容範囲はPRで明記)。
3. `cargo clippy` / `cargo fmt` を通すこと。
4. 内部実装の最適化(`nalgebra`採用など)は行わず、まず動作一致を優先すること(`RUST_PORT_SPEC.md` §1方針1)。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- `RUST_PORT_SPEC.md` §8の未決事項(内部数値計算を自前実装のままにするか等)を実装側の判断で勝手に確定させない。判断が必要な箇所は実装を止めてPR上で質問すること。
- Phase 1の範囲外(上記スコープ外の項目)には手を出さない。

## Phase 1 是正事項(2026-09-14 レビューで判明、Phase 2着手前に対応)

Phase 1(PR#1〜3)はブランチ運用ルール違反に加え、以下のテストカバレッジ不足・実バグが見つかった(コードレビューで実装を読んで確認済み)。`rust-port` への追加コミットで是正してから Phase 2 へ進むこと(この是正コミットのみ、後始末のため直接コミットを許可する)。**次の機能ブランチ運用に切り替えるのは、この是正が完了した後から。**

### テストカバレッジ不足

1. `tests/test_atomgroup.py` の `test_op_ior` / `test_op_ixor` / `test_ixor_operator` が未移植。`atom_group.rs` 側に `BitOrAssign`/`BitXorAssign`(`|=`, `^=`)の実装自体はあるため、対応するテストを追加すること。
2. `test_path_copy`(コピー構築時のpath再計算を検証するテスト)相当のテストがない。`AtomGroup` の `Clone` 実装がPython版のコピー構築子と同じpath結果になることを検証するテストを追加すること。

### 実バグ(要修正)

3. `Matrix::new`(`matrix.rs`)が `rows>0 && cols>0` をassertしており、0次元で`panic`する。`Bond::setup`が原子数0の`AtomGroup`に対してクラッシュする(Python版のnumpyベース`SymmetricMatrix(0)`はエラーにならない)。0次元を許容するよう修正すること。
4. `bond.rs`の`make_bond_matrix`が`self.atoms[p].vdw().unwrap_or(0.0)`としており、VDWテーブルに存在しない元素(atomic number >= 53: I, Xe, ランタノイド/アクチノイド、後半の遷移金属等)に対してエラーを握りつぶし0.0にフォールバックしている。該当元素を含む結合がサイレントに検出されなくなるため、エラーを適切に伝播させること(あるいはVDWテーブル自体を全元素分に拡充すること)。
5. `AtomGroup::merge`/`BitAnd`/`BitOr`/`BitXor`のいずれも`bonds: Vec<BondRecord>`フィールドをコピー・マージしていない。`bond.setup()`で検出した結合情報が、その後のマージや集合演算でサイレントに消える。全ての演算で`bonds`も適切にマージ/結合するよう修正すること。
6. `AtomGroup::merge`のサブグループ照合が`self.groups.get_mut(key)`というkey一致のみで、Python版`_merge_group`の`has_group(key)`(key一致 OR 対象グループの`name`一致)を再現していない。本来1つに統合されるべきグループが重複した兄弟グループとして分裂するケースがある。name一致のフォールバックも追加すること。
7. `Position::from_str`がPython版より厳格。Python版`Position(str)`はトークン不足時に末尾座標を0.0で補うが、Rust版はトークンが3未満だとエラーを返す。Python版と同じ「不足分は0.0補完」の挙動に合わせること。

### 中程度(可能なら是正、必須ではない)

8. `get_atom`の名前フォールバックが`a.name.trim() == key_or_name.trim()`とtrimする一方、Python版`get_atom`はtrimしない(trimするのは別メソッド`has_atomname`のみ)。空白パディングされた原子名での挙動がPython版と食い違う。
9. `BitXor`実装が`HashSet<String>`でキーを集めてから処理しているため、結果の子要素の順序が実行ごとに非決定的になる。このファイルの他の演算は`IndexMap`で順序を保証しているので、`BitXor`も同様にkeyの出現順を保つこと。
10. `SymmetricMatrix::eig()`のn=0分岐が形状不整合(空の`Vector`と1x1の`Matrix`のペアを返す)。#3の修正で0次元が到達可能になるため、あわせて修正すること。

### スコープ逸脱(2026-09-14 方針決定)

11. `SelectRange`/`Selector`が`atom_group.rs`に`pub`で実装されている件。方針:
    - `pub trait Selector`(`atom_group.rs:40`)は残してよい。`AtomGroup::select()`(Phase 1スコープの`atomgroup.py`が持つメソッド)のシグネチャに必要な最小限の抽象化のため。
    - `pub struct SelectRange`(`atom_group.rs:49`)は**`pub`から外し、`#[cfg(test)]`のテストモジュール内に移動**すること。これはPython版`select.py`の`Select_Range`クラスの機能そのものであり、Phase 1スコープ外(`select.py`→`selector.rs`は次フェーズ)。`test_select_range`を通すためだけに存在するので、本番APIとして公開する必要はない。
    - 次フェーズで`selector.rs`に正式な`Select_Range`(Python版と1:1)を実装する際、このテスト専用の暫定`SelectRange`は削除してよい(非公開なので削除は非破壊的)。

## Phase 1完了後の流れ

Phase 1は完了(全11件の是正事項を含め、Claudeレビュー通過。2026-09-14)。**ただし `main` へはまだマージしない** — Phase 2以降も完了するまで `rust-port` に積み上げていく方針とする(ユーザー判断、2026-09-14)。

## Phase 2: フォーマットI/Oの1:1移植(完了 2026-09-14)

PR#4(format/xyz/gro)・PR#5(mol2/amber_prmtop)・PR#6(biopdb)、全てClaudeレビュー通過・`rust-port`へマージ済み。PR#6は実PDBファイル3種(`1hls.pdb`/`2MGO.pdb`/`3i3zH.pdb`)をPython版で実行した結果とRust版のテスト期待値を直接突き合わせ、完全一致を確認済み(原子数・階層構造・SSBOND結合パスまで)。途中で発覚した`add_bond`/`get_bond_list`の設計差異(共通祖先ルーティング欠如)も`feature/phase2-atomgroup-bond-fix`で修正し、実データで検証済み。mmCIFは計画通りこのPhaseに含めていない(次項参照)。`main`へはまだマージしない(Phase 1と同じ方針)。

`RUST_PORT_SPEC.md` §2の対応表のうち以下を移植する。**mmCIFはこのPhaseに含めない**(下記「mmCIFを除外する理由」参照)。

| Python | Rust | 推奨PR分割 |
| --- | --- | --- |
| `format.py`(`Format`: `is_residue`/`is_chain`等の階層判定ヘルパー) | `format/mod.rs` | PR#4 |
| `xyz.py` | `format/xyz.rs` | PR#4 |
| `gro.py`(`SimpleGro`) | `format/gro.rs` | PR#4 |
| `mol2.py`(`SimpleMol2`) | `format/mol2.rs` | PR#5 |
| `amber_prmtop.py` | `format/amber_prmtop.rs` | PR#5 |
| `biopdb.py`(`Pdb`) | `format/pdb.rs` | PR#6(613行と大きいので単独PR、必要なら内部で複数コミットに分けてよい) |

`format.py`は`RUST_PORT_SPEC.md` §2の対応表に記載漏れがあるが、`xyz`/`gro`/`mol2`/`biopdb`等の実装が`Format.is_residue`/`is_chain`等の階層判定ヘルパーに依存しているため、Phase 2の対象に含める。

**スコープ外(今回は着手しない)**: mmCIF(下記参照)、構造操作(select/aminoacid/ssbond/ionpair/neutralize/modeling/superposer系)、§3の新規機能、§4の多言語バインディング。

### mmCIFを除外する理由

現行 `tests/test_mmcif.py` は合成データ1件が生のキー・バリュー辞書(`_data`)にパースできることしか検証しておらず、`AtomGroup`への変換・複数モデル・altloc・insertion code・100万原子規模の往復保証など、`RUST_PORT_SPEC.md` §3.1が要求する堅牢化の受け入れ基準を検証できるテストになっていない。「既存テストを1:1移植」するだけでは§3.1の目的を満たさないため、mmCIFは受け入れ基準とテストフィクスチャを先に整備してから独立したPhaseとして着手する。antigravityはこのPhaseでmmCIFに触れないこと。

### テストフィクスチャの扱い

`proteindf_bridge/data/`(`1hls.pdb`, `2MGO.pdb`, `3i3zH.pdb`, `ACE_ALA_NME.xyz`, `sample.gro`)を参照する既存Pythonテストがある。Rust側でも同じフィクスチャファイルを使うこと(`rust/crates/proteindf-bridge/tests/data/`等にコピーし、`env!("CARGO_MANIFEST_DIR")`基点の絶対パスで参照する。CWD依存にしないこと — Python版で過去にCWD依存のテストが壊れた実例があるため)。

### 完了の定義(Definition of Done)

1. 対応する `tests/test_*.py` を同粒度で `#[cfg(test)]` として移植し、全てpassすること。
2. `biopdb.py`(613行)はテストが3件と実装規模に対して薄いため、テスト移植に加えて、上記フィクスチャ3種のPDBファイルをPython版・Rust版の両方でパースし、原子数・座標・チェインID構成が一致することを確認するテストを追加すること(既存テストの1:1移植だけでは不十分と判断)。
3. `cargo clippy` / `cargo fmt` を通すこと。
4. Phase 1と同様、内部実装の最適化は行わず、まず動作一致を優先すること。

### やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- mmCIFには触れない(上記参照)。
- Phase 2の範囲外(構造操作・§3新規機能・§4バインディング)には手を出さない。
- **ブランチ運用ルール(上記「運用ルール」節のMUST項目)を厳守すること**: `feature/phase2-prM` ブランチで作業し、`rust-port`へは自分でマージせず、レビュー承認を待つ。Phase 1で一度違反があったため、Phase 2では特に徹底すること。

## PR#4(format + xyz + gro)レビュー結果(2026-09-14、要修正)

ブランチ運用ルールは今回遵守されていた(`feature/phase2-pr4`、`rust-port`への自己マージなし)。`cargo test`(56 passed)/`clippy`/`fmt`も報告通りクリーンで確認済み。一方、コードレビューで複数の実バグが見つかった。**`rust-port`へマージする前に、機能ブランチ上で以下を修正すること。**

### 実バグ(要修正)

1. **`AtomGroup::box()`がPhase 1で未移植**(`atomgroup.py`の`box()`メソッド自体がPhase 1で漏れていた)。このため`gro.rs`の`set_by_atomgroup`が`box_vectors`を常に`[0.0, 0.0, 0.0]`にハードコードしており(`gro.rs`の該当箇所)、GRO書き出し時に実際のボックスサイズが失われる。Python版は`atomgroup.box()`(bounding box)から計算する。`AtomGroup::box()`を`atom_group.rs`に追加してから`gro.rs`側を修正すること。
2. **`format/mod.rs`の`is_residue`/`is_chain`/`is_protein`の真偽条件がPython版と食い違う。** Python版は「直下に原子がない」ことと「サブグループがあるなら全て次階層の条件を満たす」ことだけを見ており、**サブグループが0個であること自体は失格条件ではない**(else節はloggerのみ)。`is_residue`も同様に「原子数」は真偽値に無関係(loggerのみ)。Rust版は`is_chain`/`is_protein`で`groups()==0`の場合に`false`を返しており、`is_residue`で`atoms()>0`を要求しており、いずれもPython版と異なる。空のresidue/chain/modelで判定が食い違う。`is_models`は正しく実装されているので、それを参考に3つとも修正すること。
3. **`xyz.rs`: 宣言原子数より実際の行数が少ない場合にエラーにならない。** `lines.take(num_of_atoms)`は単に反復を打ち切るだけで、行が足りなくても`Ok(())`を返し、`atoms.len() < num_of_atoms`のまま黙って抜ける。Python版は`fin.readline()`がEOFで空文字列を返し`words[0]`で`IndexError`になるため、少なくとも「エラーになる」点は一致させること。
4. **`gro.rs`: `residue_number`/`atom_number`のパース失敗時に`.unwrap_or(0)`で握りつぶしている。** Python版は`int(line[0:5])`が`ValueError`で例外を投げる。不正な列(例: GROMACSのオーバーフロー表示`*****`)が原子0番・残基0番として黙って読み込まれ、既存の原子と誤って同一グループに混入する。エラーを伝播させること。
5. **`gro.rs`: バイト単位スライス(`line[0..5]`等)が非ASCII文字を含む行で`panic`しうる。** Rustの`&str`インデックスはバイト境界でしか切れない。Python版は文字単位のスライス(`line[0:5]`)なのでpanicしない。`.chars()`ベースの切り出しに変更すること(残基名・原子名に非ASCII文字が入るケースは実際にありうる)。
6. **`gro.rs`の`get_text()`: 残基/原子番号のラップが`% 100000`(Python版は`% 10000`)。** フィールド幅(5桁)的にはRust版の方が自然だが、Python版と異なる出力になる。意図的差異として明記するか、`% 10000`に揃えるか判断すること(揃えることを推奨)。
7. **`gro.rs`: box vectorの読み込みパディングが3要素まで(Python版は6要素)。** `for i in range(len(box_vectors), 6)`相当にすること。

### PR#5のブロッカー(次PR着手前に対応必須)

8. **`AtomGroup::get_number_of_bonds()` / `get_bond_list()`がPhase 1で未移植。** これらは`atomgroup.py`にある(box()と同様の見落とし)。次PR#5で移植する`mol2.py`はこの2メソッドに依存しており(`mol2.py:54,104`)、既存`test_mol2.py`の`str(mol2)`呼び出しがこれを経由するため、**このままではPR#5が動かない**。PR#4の是正と合わせて`atom_group.rs`に追加すること。

### コード品質(必須ではないが推奨)

9. GRO固定カラム幅(5/5/5/5/8/8/8/8/8/8)や`% 100000`等のマジックナンバーが複数箇所に散らばっている。名前付き定数にまとめることを推奨。
10. nm↔Å変換で`Position`の`Mul<f64>`実装を使わず、x/y/zを手動展開している(2箇所)。`position * 10.0`のように書けるはず。
11. `get_text()`内のresidue_name/atom_name切り詰めロジックが重複している。
12. `xyz.rs`の座標パースが`Position::from_str`と同等のロジックを再実装している。再利用を検討。
13. `GroAtom.velocity`が生の`[f64; 3]`。`Position`型(既にAdd/Mul/length/dot等を持つ)の再利用を検討。

### 未確認の欠落メソッド(バックログ、都度対応)

Phase 1の`atomgroup.py`移植で他にも漏れているメソッドがある: `formula`(`get_formula`とは別)、`get_atom_keys`、`get_family`、`get_xyz`、`pickup_atoms`、`restructure`、`assign_charges`。`get_raw_data`/`set_by_dict_data`はbrd往復フォーマット用として後続フェーズで対応する想定なので今は保留でよい。上記以外は、今後のPRで依存が発生した時点で都度`atom_group.rs`に追加すること(今まとめて移植する必要はない)。

### `add_bond`/`get_bond_list`の設計差異(2026-09-14、修正完了)

**対応済み。** biopdb移植(PR#6)の事前検証でPython版と異なる結果を生む実シナリオ(SSBOND処理: `model.add_bond`後に`root.set_group`で再配置すると結合パスが古いまま取り残される)が実際に確認されたため、`get_common_path`/`get_family`(下方探索のみ、Rustの所有権木に親への逆参照がないための意図的な範囲限定)を追加し、`add_bond`が共通祖先グループへ相対パスでルーティングするよう修正した。実際のバグシナリオを再現する回帰テストで検証済み。

PR#4是正で`get_number_of_bonds`/`get_bond_list`を追加した際に判明。Python版`AtomGroup.add_bond`(`_add_bond_normalize`)は、2原子の**共通祖先グループ**(`get_family(common_path)`)を探し、そこに**相対パス**(`self.path`からの差分)で結合情報を格納する。`get_bond_list()`はこの相対パスに`self.path`を連結して絶対パスを復元する再帰処理になっている。

一方Rust版の`add_bond`(Phase 1由来)は共通祖先へのルーティングを行わず、呼び出された`self`にそのまま**絶対パス**(`atom.path`そのまま)で格納する。今回追加した`get_bond_list`は絶対パスかどうかを`starts_with('/')`で判定して連結をスキップする実装になっており、現状(結合がフラットな構造の最上位グループに追加されるケースのみ)ではPython版と同じ結果になるが、根本的な格納方式が異なる。

ネストした階層(model/chain/residue)の異なる枝にまたがる結合(例: ジスルフィド結合、PDBのCONECTレコードで表現される遠い残基間の結合)を扱うPR#6(biopdb)やPhase 3(ssbond.rs)着手前に、この差異がPython版と異なる結果を生まないか検証すること。必要であれば`add_bond`に共通祖先ルーティングを追加する。

## Phase 3: 構造操作の1:1移植(完了 2026-09-14)

PR#7(select)・PR#8(aminoacid)・PR#9(ssbond)・PR#10(ionpair)・PR#11(superposer/superposer_quaternion)、全てClaudeレビュー通過・`rust-port`へマージ済み(累計112テスト)。`modeling.py`/`neutralize.py`は計画通り除外(次項参照)。

**重要な発見(PR#11)**: `proteindf_bridge/superposer_quaternion.py`(現行Python版、`2026.8.0`)には実バグがある。`SymmetricMatrix.add()`が`get`/`set`と異なり下三角に正規化されておらず、`eig()`の`numpy.linalg.eigh(self._data, "L")`が下三角しか読まないため、四元数法の対称行列の非対角成分が実質無視される。実際にPythonを実行して確認(任意回転でKabsch法RMSD≈6e-16に対し四元数法RMSD≈0.549)。Rust版の`SymmetricMatrix::add`(Phase 1由来)は元々対称に書き込むため、`superposer_quaternion.rs`はこの影響を受けない(同一シナリオでRMSD≈3e-16、Kabsch法と一致することを確認済み)。§2で「実験的」と明記されたモジュールであるため、Python版の壊れた挙動を忠実に再現するのではなく正しい実装を優先した。**このバグはPython版`ProteinDF_bridge`本体にも存在するため、別途Python側での修正を検討する価値がある**(本Rust移植プロジェクトのスコープ外)。

## (旧)Phase 3: 構造操作の1:1移植(元のスコープ記述)

`RUST_PORT_SPEC.md` §2の対応表のうち以下を移植する。**`modeling.py`/`neutralize.py`はこのPhaseに含めない**(下記「除外する理由」参照)。

| Python | Rust | 推奨PR分割 |
| --- | --- | --- |
| `select.py`(`Select`/`Select_Symbol`/`Select_Name`/`Select_Path`系/`Select_Range`/`Select_Atom`/`Select_AtomGroup`) | `selector.rs` | PR#7 |
| `aminoacid.py` | `amino_acid.rs` | PR#8 |
| `ssbond.py`(`SSBond`) | `ssbond.rs` | PR#9 |
| `ionpair.py`(`IonPair`、`aminoacid.py`に依存) | `ion_pair.rs` | PR#10(PR#8の後) |
| `superposer.py` + `superposer_quaternion.py` | `superposer.rs` + `superposer_quaternion.rs` | PR#11 |

**スコープ外(今回は着手しない)**: `modeling.py`/`neutralize.py`(下記参照)、mmCIF(Phase 2から継続除外)、§3の新規機能、§4の多言語バインディング。

### `modeling.py`/`neutralize.py`を除外する理由

`modeling.py`(720行、ACE/NME末端キャッピング・中性化テンプレート等)には**専用テストファイルが存在せず、doctestも0件**——既存Pythonテストスイートでの検証が一切ない状態。`neutralize.py`は`from .modeling import Modeling`で直接依存しているため、`modeling.py`なしには移植できない。mmCIFと同じ理由(「既存テストの1:1移植」では正しさを担保できない)でPhase 3の対象外とし、受け入れ基準・テストフィクスチャ(例: 既知のアミノ酸構造に対する末端キャッピング結果の幾何学的検証)を別途整備してから独立したPhaseとして着手する。

### PR#7(select.py)の注意点

`atom_group.rs`には既にPhase 1由来の`pub trait Selector`と、`test_select_range`用の非公開(`#[cfg(test)]`内)の暫定`SelectRange`がある(`docs/rust-port-handoff.md`の「スコープ逸脱」節参照)。PR#7では`selector.rs`に正式な`Select_Range`(Python版と1:1)を実装し、`atom_group.rs`内の暫定実装は削除して`selector.rs`のものに置き換えること。`Selector`トレイト自体は`atom_group.rs`に残ったままでよい(`AtomGroup::select()`のシグネチャに必要なため)。

### PR#9(ssbond.py)の注意点

`tests/test_ssbond.py`の実質的なテストメソッドはコメントアウトされており(`# def test_check(self):`)、有効なのはdoctestのみ。このdoctestは`data/1hls.pdb`を使った実データ検証で、期待値`[('/model_1/A/6/', '/model_1/A/11/'), ('/model_1/A/7/', '/model_1/B/7/'), ('/model_1/A/20/', '/model_1/B/19/')]`(距離ベース、SG-SG距離 < 2.31Å)が明記されている。PR#6と同様、このdoctestをRust側でも実データ照合テストとして再現すること(空のテストメソッドの移植だけでは不十分)。

### 完了の定義(Definition of Done)

1. 対応する `tests/test_*.py`(および上記doctest)を同粒度で `#[cfg(test)]` として移植し、全てpassすること。
2. `superposer.rs`/`superposer_quaternion.rs`(Kabschアルゴリズム・四元数法)は数値計算系なので、Phase 1の`vector`/`matrix`/`position`と同様、Python版と同一入力で同一出力(RMSD・回転行列)になることをテストで担保すること。
3. `cargo clippy` / `cargo fmt` を通すこと。
4. 内部実装の最適化は行わず、まず動作一致を優先すること。

### やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- `modeling.py`/`neutralize.py`には触れない(上記参照)。
- Phase 3の範囲外(mmCIF・§3新規機能・§4バインディング)には手を出さない。
- **ブランチ運用ルール(MUST項目)を厳守**: `feature/phase3-prM` ブランチで作業し、`rust-port`へは自分でマージせず、レビュー承認を待つ。

## Phase 4: mmCIFサポート(完了 2026-09-15)

PR#12(CCD形式1:1移植)・PR#13(`_atom_site`形式新規実装)、ともにClaudeレビュー通過・`rust-port`へマージ済み(累計120テスト)。PR#13は`chain_id`/残基キーに`label_asym_id`/`label_seq_id`ではなく指示通り`auth_asym_id`/`auth_seq_id`を使うよう是正し(初版では水分子の残基が全て1つに潰れる実バグがあった)、`3I3Z.cif`で独立計算した正解値(chain A: 186原子/44残基、chain B: 274原子/64残基)と一致することを確認済み。1HLS.cif/2MGO.cifは既存のPDB形式パーサ(`pdb.rs`)の結果と原子単位で完全一致することも検証済み(座標含む782原子全件)。`_struct_conn`によるSSBOND検出も実装され、PDB SSBOND記録との整合を確認済み。

### 残された既知のギャップ

- insertion code(`pdbx_PDB_ins_code`)は値の保持のみ確認済みで、実データでの検証は未実施(適切なフィクスチャが見つかっていないため)。
- 100万原子規模での性能ベンチマークは未実施(別タスク)。

## (旧)Phase 4: mmCIFサポート(元のスコープ記述)

### 背景: スコープの再定義

`proteindf_bridge/mmcif.py`(`SimpleMmcif`)を精読した結果、既存の`get_atomgroup()`は**PDB Chemical Component Dictionary(CCD)形式**(`_chem_comp`/`_chem_comp_atom`/`_chem_comp_bond`カテゴリ、個々のリガンド/残基テンプレート定義1件を表す)にしか対応しておらず、通常のタンパク質全体構造で使われる**`_atom_site`カテゴリには一切対応していない**ことが判明した。実際の利用スクリプト(`scripts/mmcif2mol2.py`)も「複数の化学成分定義を個別のmol2ファイルに変換する」用途で、CCD形式の前提と一致する。

`RUST_PORT_SPEC.md` §3.1が目指す「PDBx/mmCIFをRust版の主力フォーマットとし、100万原子規模の構造を扱う」には、`_atom_site`ベースの全構造パーサが必要だが、これはPython版に存在しない**新規実装**である。そのため本Phaseは2つのPRに分割する。

### フィクスチャ(RCSB PDBから取得・検証済み、`proteindf_bridge/data/`と`rust/crates/proteindf-bridge/tests/data/`に配置済み)

| ファイル | 内容 | 検証済みの値 |
| --- | --- | --- |
| `ALA.cif` | CCD形式(アラニン単体定義) | 既存Python`SimpleMmcif`で実際に読み込み確認済み: 13原子、12結合(C=Oのみbond order 2、他は1)。座標は`pdbx_model_Cartn_*_ideal`列由来。 |
| `1HLS.cif` | `_atom_site`形式、20モデルNMRアンサンブル(インスリン) | model 1: 782原子、チェインA/B。既存`1hls.pdb`フィクスチャ(PR#6で検証済み)の値と完全一致。altloc/insertion codeなし。 |
| `2MGO.cif` | `_atom_site`形式、20モデルNMR(オキシトシン) | model 1: 134原子、チェインA、9残基。既存`2MGO.pdb`フィクスチャと完全一致。 |
| `3I3Z.cif` | `_atom_site`形式、X線構造・単一モデル・4チェイン | チェインA: 163原子/21残基、B: 258原子/30残基、C: 23原子/1残基(HETATM)、D: 36原子/1残基(HETATM、HOH)。**altloc "A"/"B"が20箇所ずつ存在**(`label_alt_id`列)。insertion codeは全て"?"(該当データなし、既知のギャップとして許容)。 |

### PR#12: 既存`SimpleMmcif`の1:1移植(CCD形式、`format/mmcif.rs`)

移植対象: `load`/`_get_line`(セミコロンブロック・引用符付き値・`#`コメントを含む汎用CIFトークナイザ)、`_load_data_block`/`_load_loop_block`(`loop_`テーブルの読み込み)、`get_atomgroup`(CCD形式: `_chem_comp.id`/`_chem_comp_atom.*`/`_chem_comp_bond.*`)。

**完了の定義**:
1. 既存`tests/test_mmcif.py`(薄い)を移植すること。
2. `ALA.cif`フィクスチャで実データ検証テストを追加すること。期待値: 13原子(座標は上記の通りPython版で実際に確認済み)、12結合(`C-O`のみbond order 2)。
3. `get_coordinates`相当のロジック(`pdbx_model_Cartn_*_ideal`を優先し、なければ`model_Cartn_*`にフォールバック)を1:1で再現すること。
4. `type_symbol`が`"D"`(重水素)の場合`"H"`に置き換える処理を再現すること。
5. `cargo clippy`/`cargo fmt`を通すこと。

### PR#13: `_atom_site`形式の新規実装(本来の「堅牢化」目標、`format/mmcif.rs`に追加)

**新規機能**なので、PR#12のような「Python版との1:1」ではなく、以下の受け入れ基準を満たすこと。階層構造・API設計は`biopdb.py`/`pdb.rs`(`Pdb::get_atomgroup`)と同じ規約(`model_<N>` → `<chain_id>` → `<res_seq>` → `<serial>_<name>`)に揃え、`select_model`/`select_altloc`パラメータも`Pdb::get_atomgroup`と同じ形にすること(フォーマットが違ってもダウンストリーム(ssbond.rs/ion_pair.rs等)が同じように扱えるようにするため)。

**カラムマッピングの方針**: `auth_asym_id`/`auth_seq_id`(レガシーPDB互換の著者番号付け)を`pdb.rs`の`chain_id`/`res_seq`相当として使うこと(`label_asym_id`/`label_seq_id`は内部番号付けでPDBの慣習と異なる場合があるため)。`pdbx_PDB_model_num`をモデル番号として使うこと。`label_alt_id`をaltlocとして使うこと(`.`または`?`は「altlocなし」として扱う)。

**完了の定義**:
1. `1HLS.cif`・`2MGO.cif`のmodel 1を解析した結果が、既存の`Pdb::get_atomgroup`(PDB形式、PR#6で検証済み)が`1hls.pdb`/`2MGO.pdb`から生成する結果と**完全一致**すること(原子数・チェイン/残基構造・最初の原子のsymbol/座標)。フォーマットが違っても同じ構造からは同じ結果が得られることを保証する強力な検証。
2. `3I3Z.cif`で、チェインごとの原子数・残基数(表の値)が一致すること。altlocのデフォルト選択(`"A"`または空欄)で、20件の`"B"`altloc原子が除外されることを検証すること。
3. カラム固定幅を使わないmmCIFの性質上、レガシーPDBの99,999原子上限は原理的に存在しないことをコードレビューで確認する(巨大原子数の実ベンチマークは本PRの必須要件ではなく、別タスクとする)。
4. `insertion code`(`pdbx_PDB_ins_code`)は値を保持・パースすること(専用フィクスチャがないため、実データでのテストは今回のスコープ外。既知のギャップとして`docs/rust-port-handoff.md`に残す)。
5. `cargo clippy`/`cargo fmt`を通すこと。

### スコープ外

- §3の他の新規機能(二次構造推定・InteractionSet・QC結果I/O・CUBEパーサ)、§4の多言語バインディング。
- insertion codeの実データ検証(専用フィクスチャ未確保、既知のギャップ)。
- 100万原子規模での性能ベンチマーク(別タスク)。

### やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- **ブランチ運用ルール(MUST項目)を厳守**: `feature/phase4-prM` ブランチで作業し、`rust-port`へは自分でマージせず、レビュー承認を待つ。PR#12を先に、PR#13をその後に。

## 運用ルールの変更(2026-09-15): `rust-port`統合ブランチの廃止

Phase 1〜4が`main`へマージされ(`superposer_quaternion.py`のPythonバグ修正も`main`へ直接マージ済み)、`rust-port`ブランチは`main`より遅れた状態になった。Phase 5以降は**`rust-port`を経由せず、`main`から直接機能ブランチを切り、レビュー承認後に`main`へ直接マージする**運用に切り替える。ブランチ命名・レビューゲート(MUST項目)等の他のルールは変更なし。

## 運用ルールの変更(2026-09-15、2026.9.0リリース後): GitFlow運用への移行

`2026.9.0`(Rust移植Phase 1〜5 + `superposer_quaternion.py`バグ修正 + SPDXヘッダー移行を含む)を`main`に直接コミット・タグ付けしてリリースした後、**今後はGitFlow的な運用に切り替える**(本プロジェクトが2024.03.0までは`develop`+`release/*`ブランチを使っていた慣習への回帰)。`main`から`develop`ブランチを作成済み。

**新しいブランチ運用(2026.9.0以降のPRから適用、MUST)**:
- `main`: 常にリリース済みの状態のみを反映する。直接コミットしない。バージョンタグ(例: `2026.9.0`、`v`プレフィックスなし)はここに打つ。
- `develop`: 通常の開発の統合ブランチ。**Phase 6以降の機能ブランチはここから切り、ここへマージする**(`main`ではない)。
- `feature/*`: `develop`から切り、レビュー承認後に`develop`へマージする(これまでの`feature/phaseN-prM`と同じ命名・レビューゲートのルールを維持)。
- `release/X.Y.Z`: リリース準備時に`develop`から切る。バージョン番号の更新等の最終調整をここで行い、`main`(タグ付け)と`develop`の両方にマージする。
- `hotfix/*`: リリース済みの`main`に対する緊急修正用。`main`から切り、`main`(タグ付け)と`develop`の両方にマージする。

**PRごとの機能ブランチ・マージ前レビューゲート等の既存MUST項目は維持。ブランチの起点と合流先が`main`→`develop`に変わるだけ。**

## Phase 4.5: プロジェクト名変更(`pdf-bridge` → `proteindf-bridge`、Phase 5着手前に必須)

「PDF」(Portable Document Format)との混同を避けるため、プロジェクト名を`pdf-bridge`から`proteindf-bridge`に変更する。Phase 1〜4で既にマージ済みの全ファイルに影響するため、専用PRとして先に対応すること。

**対象**:
1. ディレクトリ: `rust/crates/pdf-bridge/` → `rust/crates/proteindf-bridge/`(`git mv`)
2. `rust/Cargo.toml`: workspace memberパスを`crates/proteindf-bridge`に更新
3. `rust/crates/proteindf-bridge/Cargo.toml`: `package.name`を`proteindf-bridge`に変更
4. 全テストファイル(`tests/test_*.rs`、8ファイル)の`use pdf_bridge::...`を`use proteindf_bridge::...`に変更(Cargoは`-`を`_`に自動変換するため、Rustコード内の参照はこの形になる)
5. `RUST_PORT_SPEC.md`・`docs/rust-port-handoff.md`は既にClaude側で`proteindf-bridge`表記に更新済み。追加の対応は不要。

**完了の定義**: `cargo build --workspace`・`cargo test --workspace`(120件)・`cargo clippy`・`cargo fmt`が全て通ること。振る舞いの変更は一切ないので、既存テストが全てそのままpassすることを確認するだけでよい。

**ブランチ**: `fix/rename-to-proteindf-bridge`を`main`から切り、レビュー承認後`main`へマージ。この対応が完了してからPhase 5(PR#14)に着手すること。

## Phase 5: Pythonバインディング(PyO3)(完了 2026-09-15)

PR#14(基盤・データモデル)・PR#15(フォーマットI/O)・PR#16(構造操作)、全てClaudeレビュー通過・`main`へマージ済み。`proteindf_bridge_rs`パッケージとして`import`可能。pytestベースの検証(`tests/test_rs_phase1.py`〜`test_rs_phase3.py`、累計42件)は全て既存の純Python版`proteindf_bridge`との直接比較になっており、実データフィクスチャ(1hls.pdb/2MGO.pdb/1HLS.cif等)による相互検証も実施済み。

レビューで見つかった主な実バグ(いずれも修正済み): `cargo test --workspace`のリンクエラー(`extension-module`featureの扱い)、`AtomGroup.__getitem__`の非多態性(グループしか返さない)、`AtomGroup.select()`のカスタムPythonセレクタで例外が握りつぶされる問題。

### 背景・目的

`RUST_PORT_SPEC.md` §4の多言語バインディング方針のうち、Pythonバインディングに着手する。既存`ProteinDF_bridge`ユーザーが最小コストで移行できるよう、既存Python API(クラス名・メソッド名)を可能な限り踏襲する。`selector.rs`で既に確立されているパターン(`Select_Symbol`等のPython互換エイリアス)を踏襲すること。

### Cargo workspace構成

```
rust/
├── Cargo.toml
└── crates/
    ├── proteindf-bridge/     # コアライブラリ(既存)
    └── proteindf-bridge-py/  # 新規: PyO3バインディング
        ├── Cargo.toml     # crate-type = ["cdylib"], pyo3依存
        ├── pyproject.toml # maturin設定
        └── src/lib.rs
```

Pythonパッケージ名は`proteindf_bridge_rs`とし、既存の純Python版`proteindf_bridge`パッケージと共存できるようにすること(名前が衝突しないため、移行期間中に両方インストールして比較検証できる)。

### エラー変換方針

`BridgeError`を、既存Python例外階層(`BrError`基底、`BrInputError`、`BrValueError`)に対応するPython例外クラスとして`pyo3::create_exception!`で定義し、変換すること。`BridgeError::General`→`BrError`、`InputError`→`BrInputError`、`ValueError`→`BrValueError`。`PeriodicTable`関連の独自エラー種別(`AtomicNumberNotFound`等、Python版に対応物がない)は`BrValueError`にマッピングすること(Python版periodictable.pyは例外を素通しするだけで独自メッセージを持たないため、最も意味的に近いものとして扱う)。

### 完了の定義(Definition of Done、全PR共通)

1. `maturin develop`でビルドでき、Pythonから`import proteindf_bridge_rs`できること。
2. pytestベースのテストを追加し、**既存`proteindf_bridge`の同等クラスと同じ操作をして結果を比較する**(例: 同じ分子構造を新旧両方のAPIで構築し、原子数・座標が一致することを確認)。単にバインディングが動くことだけでなく、既存Python版との挙動一致を検証すること(これまでのRust移植と同じ検証方針)。
3. クラス名・メソッド名は既存Python API(`SPEC.md`参照)に合わせること。
4. `cargo clippy`/`cargo fmt`に加え、Python側のテストも実行して報告すること。

### PR#14: 基盤・データモデルのバインディング(Phase 1相当)

対象: `error`(例外階層)、`PeriodicTable`、`Vector`、`Matrix`/`SymmetricMatrix`、`Position`、`Atom`、`Bond`、`AtomGroup`。

### PR#15: フォーマットI/Oのバインディング(Phase 2・4相当)

対象: `Xyz`、`SimpleGro`、`SimpleMol2`、`AmberPrmtop`、`Pdb`、`SimpleMmcif`、`Format`。PR#14完了後に着手。

### PR#16: 構造操作のバインディング(Phase 3相当)

対象: `Selector`系(`Select_*`)、`AminoAcid`、`SSBond`、`IonPair`、`Superposer`、`SuperposerQuaternion`。PR#14完了後に着手(PR#15と並行可)。

### スコープ外

- C/C++バインディング(cbindgen、§4の別項目)。
- `modeling.py`/`neutralize.py`、§3新機能、§3.4-3.6のQC結果I/O。

### やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- **ブランチ運用ルール(MUST項目)を厳守**: `feature/phase5-prM` ブランチを`main`から切って作業し、`main`へは自分でマージせず、レビュー承認を待つ。PR#14から着手すること。

## ブランチ運用(2026-09-15以降、GitFlow)

`CONTRIBUTING.md`を参照。**Phase 6以降の機能ブランチは`main`ではなく`develop`から切り、`develop`へマージすること。** PRごとの機能ブランチ・マージ前レビューゲート等の既存MUST項目は変更なし。

## Phase 6: `modeling.py`/`neutralize.py`(完了 2026-09-15)

PR#17(`brd.rs`)・PR#18(`modeling.rs`)・PR#19(`neutralize.rs`)、全てClaudeレビュー通過・`develop`へマージ済み(累計140テスト)。`get_ACE`/`get_NME`は実フィクスチャでPython版と座標が完全一致することを、`neutralize`は実PDBフィクスチャ(1hls.pdb、10件のイオン追加)でPython版と完全一致することを、それぞれ実際にPythonを再実行して確認済み。`_exempt_list`のデッドコード挙動(Python版で実質機能していない)も忠実に再現されている。

**これで`RUST_PORT_SPEC.md` §2の1:1移植対応表(全モジュール)が完了した。**

### 背景: 隠れた前提条件の発見

`modeling.py`を精読した結果、`Modeling.__init__`が**コンストラクタの時点で無条件に**4つの`.brd`(MessagePack)参照構造ファイル(`proteindf_bridge/data/ACE_ALA_NME_{trans1,trans2,cis1,cis2}.brd`)を読み込むことが判明した。これは`get_ACE`/`get_NME`だけでなく、`neutralize.py`が使う`neutralize_Nterm`等のメソッドを使うだけでも発生する(`Neutralize.__init__`が内部で`Modeling()`を生成するため)。

つまり`functions.py`のMessagePack I/O(`load_msgpack`/`save_msgpack`)と、`AtomGroup`/`Atom`の`get_raw_data`/辞書コンストラクタ(`set_by_dict_data`/`set_by_raw_data`)が前提条件だが、これは`RUST_PORT_SPEC.md` §2の対応表に`functions.py` → `brd.rs`として元々計画されていたにもかかわらず、Phase 1〜5のどこでも着手されていなかった。

さらに`RUST_PORT_SPEC.md` §1方針3は、Rust版のネイティブ`.brd`形式がYUI側の`MessagePack + zstd`ヘッダー設計(`[Magic: "YUI\0"(4B)] + [Version(1B)] + [Compression Flag(1B)] + [Payload]`)と相互運用できることを求めている。しかし実際に確認したところ、**既存のPython版`.brd`ファイル(このPhaseで使う参照構造フィクスチャ含む)はこのヘッダーを持たないプレーンなMessagePack**だった。したがって本Phaseは、mmCIFの時と同様に「1:1移植部分」と「新規実装部分」に分かれる。3つのPRに分割する。

### フィクスチャ(検証済み、`proteindf_bridge/data/`に既存)

| ファイル | 内容 | 検証済みの値 |
| --- | --- | --- |
| `ACE_ALA_NME_trans1.brd` / `_trans2.brd` / `_cis1.brd` / `_cis2.brd` | ACE-ALA-NMEの4種コンフォーマー参照構造(`Modeling.__init__`が読み込む本体) | 各22原子、3グループ(residue "1"=ACE, "2"=ALA, "3"=NME)。Python版で実際にロード確認済み。 |
| `ACE_ALA_NME.brd` | 同上、未分類版 | 22原子、3グループ。 |
| `NML.brd` / `NML_trans.brd` | 他の参照構造(`modeling.py`では未使用と思われるが存在確認) | 各12原子、0グループ(フラット構造)。 |

**参考(Python版で実際に確認したget_ACEの出力例)**: `ACE_ALA_NME_trans1.brd`のresidue "2"(ALA)をresに、residue "3"(NME)をnext_aaに渡すと、`get_ACE`は6原子のAtomGroupを返す(CH3, 3×H, C, O相当)。Rust版でも同じ入力で同じ出力(座標を含む)になることを確認すること。

### PR#17: `brd.rs`(MessagePack往復フォーマット、`functions.py`の一部 + 新規YUIヘッダー)

**1:1移植部分**:
- `functions.py`の`load_msgpack`/`save_msgpack`(プレーンMessagePack読み書き、ヘッダーなし)。
- `Atom::get_raw_data()` / 辞書からの構築(Python版`Atom.set_by_raw_data`相当)。スキーマ: `{"Z": atomic_number(int), "name": str, "Q": charge(float), "xyz": [x,y,z], "force": [x,y,z]}`。
- `Position::get_raw_data()`: `[x, y, z]`。
- `AtomGroup::get_raw_data()` / 辞書からの構築(Python版`AtomGroup.set_by_dict_data`相当)。スキーマ: `{"name": str, "groups": {key: <AtomGroup再帰>, ...}, "atoms": {key: <Atom>, ...}, "bonds": [...]}`(`groups`/`atoms`は空なら省略)。
- 未知のキーはPython版と同様、エラーにせずログ出力のみで無視すること(実際に検証したところ、既存`.brd`データの読み込み時に`AtomGroup::set_by_dict_data(): unknown key: Q=None`という警告がPython版でも出るが、クラッシュはしない。この寛容な挙動を1:1で再現すること)。

**新規実装部分(RUST_PORT_SPEC.md §1.3)**:
- YUI互換ヘッダー付き形式: `[Magic: "YUI\0"(4B)] + [Version(1B)] + [Compression Flag(1B)] + [Payload(MessagePack, オプションでzstd圧縮)]`の読み書き。既存のプレーンMessagePack形式とは別のAPI(例: `save_brd_yui`/`load_brd_yui`)として実装し、既存`.brd`ファイルの読み込みパス(1:1移植部分)を壊さないこと。zstd圧縮には適切なクレート(`zstd`クレート等)を使用すること。

**完了の定義**:
1. 上記フィクスチャ全て(`ACE_ALA_NME_{trans1,trans2,cis1,cis2}.brd`、`ACE_ALA_NME.brd`、`NML.brd`、`NML_trans.brd`)をRust版で読み込み、Python版(`load_msgpack`+`AtomGroup(data)`)と原子数・グループ数・パスリストが一致することをテストで検証すること。
2. `AtomGroup`/`Atom`の`get_raw_data`→再構築のラウンドトリップ(構造が保持されること)をテストすること。
3. YUI互換ヘッダー形式は、書き込み→読み込みの自己整合性(ヘッダーのマジックバイト・バージョン・圧縮フラグが正しく解釈されること)をテストすること(他言語実装との相互運用テストは本PRの範囲外、将来のC/C++バインディングやYUI側との統合時に別途検証)。
4. `cargo clippy`/`cargo fmt`を通すこと。

### PR#18: `modeling.py`本体(`modeling.rs`、PR#17完了後)

対象: `get_ACE`/`get_NME`/`get_ACE_simple`/`get_NME_simple`/`_match_ACE`/`_match_NME`/`_match_residues`(ACE/NME末端キャッピング)、`add_methyl`/`get_NH3`/`arbitary_rotate_matrix`/`select_residues`/`get_last_index`(幾何ヘルパー)、`neutralize_Nterm`/`neutralize_Cterm`/`neutralize_GLU`/`neutralize_ASP`/`neutralize_LYS`/`neutralize_ARG`/`neutralize_FAD`/`_get_neutralize_pos_{NH3,NH2,COO,POO}_type`(中性化イオン位置計算)。

**完了の定義**:
1. `get_ACE`/`get_NME`は、上記フィクスチャを使い、Python版と同じ入力(residue "2"をres、residue "3"またはNoneをnext_aaに)で同じ出力(原子数・座標・RMSD)になることをテストすること。4種のコンフォーマー全てで検証すること。
2. `neutralize_*`系は、期待する原子(例: `neutralize_Nterm`ならN/H1/H2/HXTまたはH3を持つ残基)を合成データで用意し、Python版と同じイオン位置(座標)になることをテストすること。
3. `neutralize_FAD`のPython版は`OP1`/`O1P`のどちらの命名でも対応する分岐があるが、どちらも存在しない場合`raise`(引数なしの再送出、実質クラッシュ)する。この「該当なしならエラー」という挙動を1:1で再現すること(サイレントなフォールバックにしないこと — これまでのレビューで繰り返し指摘した`unwrap_or`パターンと同じ考え方)。
4. `cargo clippy`/`cargo fmt`を通すこと。

### PR#19: `neutralize.py`(`neutralize.rs`、PR#18完了後)

対象: `Neutralize`クラス。`modeling.rs`(PR#18)と既存の`ion_pair.rs`(Phase 3で完了済み)に依存する。

**注意**: Python版`_neutralize`メソッド内の`exempt_list = []  # self._exempt_list()`は、`_exempt_list()`の呼び出しがコメントアウトされており、**`_exempt_list`/`_divide_path`は実質デッドコードで、除外リストの仕組みは常に空リストとして動作する**(実際には機能していない)。これはPython版の実際の挙動なので、1:1移植方針に従いそのまま(常に空の除外リストとして)再現すること。「動いていない機能を直す」ことはスコープ外。

**完了の定義**:
1. 荷電残基(GLU/ASP/LYS/ARG、N末端/C末端)を含む合成データまたは既存PDBフィクスチャで、Python版と同じ数・位置のイオンが追加されることをテストすること。
2. `cargo clippy`/`cargo fmt`を通すこと。

### スコープ外

- §3の新規機能(二次構造推定・InteractionSet)、§3.4-3.6のQC結果I/O、§4 C/C++バインディング。

### やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- **ブランチ運用ルール(MUST項目)を厳守**: `feature/phase6-prM`ブランチを**`develop`から**切って作業し(GitFlow運用、上記参照)、`develop`へは自分でマージせず、レビュー承認を待つ。PR#17→PR#18→PR#19の順に着手すること。

## Phase 7: バックボーン二面角(φ/ψ)計算 — Ramachandranプロット対応(完了 2026-09-15)

PR#20(`dihedral_angle` + `ramachandran.rs`)・PR#21(Pythonバインディング)、全てClaudeレビュー通過・`develop`へマージ済み。基準値表(1hls.pdb、chain A、残基4/5/10のφ/ψ)はRust版・Pythonバインディング双方で誤差1e-3度以内の一致を確認済み。PR#20レビュー時に指摘した`IndexMap`挿入順依存の問題(`calc_phi_psi`が`sort_nicely`による明示ソートをせず、挿入順に依存していた)は、挿入順を意図的に崩したチェインでの回帰テストとともに修正済み。PR#21のPythonバインディングは例外を投げる経路がなく(`extract_position`の`?`伝播のみ)、`unwrap`/`unwrap_or`によるサイレントなフォールバックも無いことをコードレビューで確認した。pytest側でも幾何学的サニティ・実PDB基準値・欠損主鎖の安全スキップ・挿入順スクランブルの4ケースを検証し、既存の回帰テスト(46件)を含め全てパス。

**これでRamachandranプロット対応(新規機能)が完了した。**

### 背景

ユーザーからRamachandranプロット(タンパク質主鎖のφ/ψ二面角の散布図によるコンフォメーション評価)のサポート可否を問われた。調査の結果、**既存Python版・Rust版のどちらにも二面角計算のコードが一切存在しない**ことを確認した(`grep -rn "dihedral\|phi\|psi" proteindf_bridge/*.py`で該当なし)。これは`RUST_PORT_SPEC.md` §3.2(二次構造推定)とは別の、独立した新規機能である(§3.2のDSSP相当は主鎖水素結合ベースの判定で、φ/ψ角ベースの判定とは別のアプローチ)。プロット描画そのもの(可視化)は本ライブラリのスコープ外とする(構造I/Oライブラリであり、可視化はYUI側またはPython側の別スクリプトの役割)。

**Python版の実装が存在しないため、他のPhaseと異なり「Python版との1:1比較」という受け入れ基準が使えない。** 代わりに以下の2段階で検証する: (1) 幾何学的に正解が既知の合成データでの符号・大きさの検証、(2) 実PDBフィクスチャに対して独立に計算した基準値との比較。

### 二面角の定義(IUPAC標準、数値的に安定な`atan2`方式を使用)

4点 `p1, p2, p3, p4` に対して:

```
b1 = p2 - p1
b2 = p3 - p2
b3 = p4 - p3
n1 = b1 × b2
n2 = b2 × b3
m1 = n1 × (b2 / |b2|)
x = n1 · n2
y = m1 · n2
angle = atan2(y, x)  (ラジアン、度に変換する場合は180/πを乗じる。範囲: -180°〜+180°)
```

φ(phi) = dihedral(前残基のC, 現残基のN, 現残基のCA, 現残基のC)
ψ(psi) = dihedral(現残基のN, 現残基のCA, 現残基のC, 次残基のN)

最初の残基はφが定義できず(前の残基のCがない)、最後の残基はψが定義できない(次の残基のNがない)ので、`Option<f64>`で表現すること。

### 検証済みの基準値(独立計算、`proteindf_bridge/data/1hls.pdb`、model_1、chain A)

上記の`atan2`方式の二面角計算をnumpyで独立実装し、実際に計算して得た値:

| 残基番号 | φ(度) | ψ(度) |
| --- | --- | --- |
| 4 | 70.5993 | 3.3945 |
| 5 | 121.6311 | 26.5618 |
| 10 | 84.5507 | -99.6314 |

Rust版でも同じ入力(同じPDBファイル、同じ残基)から同じ値(誤差1e-3度程度まで)が得られることを確認すること。

### PR#20: `dihedral_angle`関数 + Ramachandran計算(`ramachandran.rs`、新規)

**対象**:
1. `Position`(または新規`geometry.rs`)に汎用の二面角計算関数を追加: `dihedral_angle(p1: &Position, p2: &Position, p3: &Position, p4: &Position) -> f64`(上記の`atan2`方式、度単位で返す)。
2. `ramachandran.rs`(新規): チェイン(`AtomGroup`)を受け取り、連続する残基(整数キー、`get_group_list`等で自然順ソート済みのキーを使う)ごとにφ/ψを計算する関数。例: `calc_phi_psi(chain: &AtomGroup) -> Vec<RamachandranAngle>`、`RamachandranAngle { residue_key: String, residue_name: String, phi: Option<f64>, psi: Option<f64> }`。N/CA/C原子が欠けている残基は安全にスキップすること(エラーにしない。可視化目的のデータ収集なので、部分的に欠損した構造でも可能な範囲で結果を返す方が実用的)。

**完了の定義**:
1. 合成データによる幾何学的サニティテスト: 平面上に4点を配置し、角度が解析的に既知の値(0°, 90°, 180°, -90°等)になるケースを複数用意し、`dihedral_angle`の符号・大きさが正しいことを検証すること。
2. 上記の基準値表(1hls.pdb、chain A、残基4/5/10)と一致することを検証する実データテストを追加すること。
3. 最初の残基でφが`None`、最後の残基でψが`None`になることを検証すること。
4. `cargo clippy`/`cargo fmt`を通すこと。

### PR#21: Pythonバインディング(`proteindf-bridge-py`、PR#20完了後)

`ramachandran.rs`をPyO3で公開する。Phase 5の既存パターン(クラス名・メソッド名をなるべく分かりやすく、エラーは`BrError`系にマッピング)を踏襲すること。pytestで、上記基準値表と同じ検証(1hls.pdbの残基4/5/10のφ/ψ)を行うこと。

### スコープ外

- プロット描画(matplotlib等での可視化)。
- §3.2のDSSP相当(水素結合ベースの二次構造推定)は別Phaseとする。
- C/C++バインディング。

### やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない(この機能はPython版に存在しないため、新規Pythonコードの追加も不要)。
- **ブランチ運用ルール(MUST項目)を厳守**: `feature/phase7-prM`ブランチを`develop`から切って作業し、`develop`へは自分でマージせず、レビュー承認を待つ。PR#20を先に、PR#21をその後に。

## Phase 8: 主鎖水素結合検出 + 二次構造推定(DSSP相当、§3.2、今回のスコープ、2026-09-15 受け入れ基準確定)

### 背景

`RUST_PORT_SPEC.md` §3.2は「Kabsch-Sanderの静電モデルによる主鎖 N-H...O=C 水素結合パターンからヘリックス/シートを推定する、オリジナルDSSPと同様のアルゴリズムを`secondary_structure.rs`に新規実装する」ことを求めている。同§3.3は側鎖水素結合・CH-π等の`InteractionSet`検出を求めており、「この主鎖水素結合検出ロジックは3.3節の水素結合検出と共有する」と明記されている。そのため本Phaseでは**主鎖水素結合検出を独立モジュール(`hydrogen_bond.rs`)として先に実装し(PR#22)、それを土台に二次構造推定(`secondary_structure.rs`、PR#23)を実装する**。側鎖水素結合判定・CH-π判定は§3.3として別Phaseに残す(スコープ外、下記参照)。

Phase 7と同様、**Python版に対応する実装が一切存在しない新規機能**のため、「Python版との1:1比較」という受け入れ基準は使えない。代わりに、Kabsch & Sander (1983)のオリジナルDSSPアルゴリズムを1983年の原論文の定義通りに再実装したオープンソース参照実装([PyDSSP](https://github.com/ShintaroMinami/PyDSSP)、MITライセンス、PyTorch/NumPy実装)を独立に動かし、実PDBフィクスチャ(`1hls.pdb`)に対する基準値(残基ごとのH-bondエネルギー・二次構造ラベル)を算出した。**この参照実装のコードそのものを移植するのではなく、以下に記述するアルゴリズム定義(Kabsch-Sander 1983の一次情報に基づく)に従って独立に実装し、算出された基準値と一致することをテストで検証すること。**

### アルゴリズム定義

#### 1. 主鎖アミドH原子の疑似座標(PDBに水素原子がない場合)

残基 `i`(N, CA, C, O座標を持つ)について、直前の残基 `i-1` のC原子座標が必要:

```
vec_cn  = normalize(N(i) - C(i-1))
vec_can = normalize(N(i) - CA(i))
vec_nh  = normalize(vec_cn + vec_can)
H(i)    = N(i) + 1.01 * vec_nh
```

チェインの最初の残基(直前のCがない)は疑似H座標を計算できないため、その残基はH-bondのドナーになれない(後述のHBond判定で常にfalseとする)。

#### 2. Kabsch-Sander静電エネルギー

残基 `d`(ドナー側、N-H)と残基 `a`(アクセプター側、C=O)の間のH-bondエネルギー(kcal/mol):

```
E(d, a) = q1*q2 * (1/r(O_a,N_d) + 1/r(C_a,H_d) - 1/r(O_a,H_d) - 1/r(C_a,N_d)) * 332
```

`q1*q2 = 0.084`(部分電荷 q1=0.42, q2=0.20 の積。332はkcal/mol変換定数)。`r(X,Y)`は原子間距離(Å)。

**HBond(d, a)** は次の条件を全て満たすときtrue:
1. `E(d, a) < -0.5`(kcal/mol)
2. `d != a` かつ `d`と`a`が隣接残基でない(`|d - a| <= 2` の場合は自明な骨格內H-bondとして除外。DSSPの標準的な扱い)
3. 残基`d`が疑似H座標を計算可能(直前残基が存在する)

#### 3. ターン(turn)の定義

`n ∈ {3, 4, 5}` について、残基`i`に対する「n-turn」は `HBond(i+n, i)` が成立すること(残基`i+n`のN-Hが残基`i`のC=Oに水素結合)。

#### 4. ヘリックス判定

2つの連続するn-turn(位置`i`と`i+1`両方でn-turnが成立)がある場合、その2つのターンが覆う`n`残基(`i`から`i+n-1`、実装により`i+1`から`i+n`の場合もあるので基準値表と突き合わせて調整すること)を「n-ヘリックスのコア」としてマークする。**4-turn(α-ヘリックス相当)を最優先とし、既に4-ヘリックスとしてマークされた残基は3-ヘリックス・5-ヘリックスの判定から除外する**(DSSP標準の優先順位)。最終的に3-ヘリックス・4-ヘリックス・5-ヘリックスのいずれかに属する残基を一律「H(ヘリックス)」として出力する(**本Phaseでは3-10/α/πの区別はせず統合した1状態として扱う**。スコープ簡略化、下記参照)。

#### 5. ブリッジ(β-シート)判定

残基`i`と`j`(`i`, `j`は互いに疎な位置)について:

**平行ブリッジ**: `[HBond(j+1, i) かつ HBond(i+2, j+1)]` または `[HBond(i+1, j) かつ HBond(j+2, i+1)]`

**逆平行ブリッジ**: `[HBond(j+1, i+1) かつ HBond(i+1, j+1)]`(相互H-bond) または `[HBond(j+2, i) かつ HBond(i+2, j)]`

いずれかが成立する`(i,j)`ペアを「ラダー」とする。残基`i`が何らかの`j`とラダーを形成する場合、その残基を「E(ストランド/シート)」として出力する。

#### 6. 最終分類

各残基は「H(ヘリックス)」「E(ストランド)」「-(ループ/コイル)」のいずれか1つに分類される(優先順位: ヘリックスとストランドが両方成立することは通常ないが、もし競合したらヘリックス優先でよい)。

### スコープの簡略化(意図的な設計判断)

- 本来のDSSPは8状態(H, G, I, E, B, T, S, -)に分類するが、本Phaseでは**H(ヘリックス統合)/E(ストランド)/-(ループ)の3状態**に簡略化する。理由: `RUST_PORT_SPEC.md` §3.2はYUI側のリボン(カートゥーン)表示での利用を目的としており、そこで必要なのはヘリックス/シート/コイルの区別のみで、3-10ヘリックスとαヘリックスの区別等は現時点で要求されていない。より細かい分類が必要になった場合は将来のPhaseで拡張する。
- プロリンなど、本来アミドHを持たない特殊残基の扱い(DSSPはプロリンをH-bondドナーから除外する)は、本Phaseでは特別扱いしない(全残基を一律ドナー候補として扱う簡易実装)。理由: 検証用フィクスチャにプロリンがドナー側になる問題ケースが含まれておらず、正しい特殊扱いの実装・検証には追加のケーススタディが必要なため、後続Phaseの課題とする。

### 検証済みの基準値(独立計算、`proteindf_bridge/data/1hls.pdb`)

上記アルゴリズムをNumPyで独立実装した参照コードを実際に動かして得た値。

**chain A(21残基)**:

| 残基 | 残基名 | 判定 |
| --- | --- | --- |
| 1 | GLY | - |
| 2 | ILE | - |
| 3 | VAL | H |
| 4 | GLU | H |
| 5 | GLN | H |
| 6 | CYS | H |
| 7 | CYS | - |
| 8 | THR | - |
| 9 | SER | - |
| 10 | ILE | - |
| 11 | CYS | - |
| 12 | SER | - |
| 13 | LEU | - |
| 14 | TYR | - |
| 15 | GLN | - |
| 16 | LEU | - |
| 17 | GLU | H |
| 18 | ASN | H |
| 19 | TYR | H |
| 20 | CYS | - |
| 21 | ASN | - |

**chain B(30残基)**:

| 残基 | 残基名 | 判定 |
| --- | --- | --- |
| 1 | PHE | - |
| 2 | VAL | - |
| 3 | ASN | - |
| 4 | GLN | - |
| 5 | HIS | - |
| 6 | LEU | - |
| 7 | CYS | - |
| 8 | GLY | - |
| 9 | SER | H |
| 10 | HIS | H |
| 11 | LEU | H |
| 12 | VAL | H |
| 13 | GLU | H |
| 14 | ALA | H |
| 15 | LEU | H |
| 16 | HIS | H |
| 17 | LEU | H |
| 18 | VAL | H |
| 19 | CYS | H |
| 20 | GLY | - |
| 21 | GLU | - |
| 22 | ARG | - |
| 23 | GLY | - |
| 24 | PHE | - |
| 25 | PHE | - |
| 26 | TYR | - |
| 27 | THR | - |
| 28 | PRO | - |
| 29 | LYS | - |
| 30 | THR | - |

このPhaseのフィクスチャに`E(ストランド)`判定の残基が1件も現れない(インスリンA/B鎖はヘリックスのみの小さな構造のため)。ブリッジ/ストランド判定ロジックのテストは、上記の実データに加えて、**幾何学的に構築した合成データ(逆平行β-シートを模した2本の短いストランドを人工的に配置し、教科書的なジグザグ骨格水素結合パターンを持たせたもの)でE判定が出ることを別途検証すること**(合成データの設計はPR#23の完了の定義に含む)。

**参考: 検証済みの主鎖H-bondエネルギー(chain A、抜粋、kcal/mol)**:

| ドナー残基 | アクセプター残基 | E |
| --- | --- | --- |
| 3 | 2 | -3.8074 |
| 4 | 3 | -3.7912 |
| 5 | 4 | -3.7696 |
| 6 | 5 | -3.7979 |
| 7 | 6 | -3.7884 |
| 8 | 7 | -3.7836 |

(残基2→1、3→2等の隣接残基間の弱いH-bondは上記の条件2(`|d-a|<=2`除外)によりH-bond判定から除外される。実際に`E(2,1) = -3.8103`のように閾値-0.5を超える強いエネルギーが出るが、これは自明な骨格内相互作用でありターン検出の対象外とすること。)

### PR#22: `hydrogen_bond.rs`(主鎖Kabsch-Sander水素結合検出、新規、§3.2と§3.3で共有)

**対象**:
1. 疑似H原子座標計算関数(上記1節)。
2. `HydrogenBond { donor_residue_key: String, acceptor_residue_key: String, energy: f64 }`構造体。
3. `calc_backbone_hbonds(chain: &AtomGroup) -> Vec<HydrogenBond>`: チェイン内の全残基ペアについてKabsch-Sanderエネルギーを計算し、上記2節のHBond条件(エネルギー閾値 + 隣接除外 + ドナー可否)を満たすペアのみ返す関数。残基は`sort_nicely`でソート済みの順序で扱うこと(Phase 7と同じ理由でIndexMap挿入順に依存しないこと)。
4. N/CA/C/Oのいずれかが欠けている残基は安全にスキップすること(Ramachandranと同じ設計方針。エラーにしない)。

**完了の定義**:
1. 上記の基準値表(chain A、ドナー→アクセプターのエネルギー抜粋6件)と誤差1e-3 kcal/mol以内で一致することをテストすること。
2. 隣接残基間(`|d-a|<=2`)の強いエネルギーがHBond判定から除外されることを検証すること(例: 残基2→1のエネルギーは-0.5未満でも`calc_backbone_hbonds`の結果に含まれないこと)。
3. 最初の残基がドナーになれない(疑似H座標が計算不能)ことを検証すること。
4. `cargo clippy`/`cargo fmt`を通すこと。

### PR#23: `secondary_structure.rs`(二次構造推定、PR#22完了後)

**対象**: 上記3〜6節のターン/ヘリックス/ブリッジ判定ロジック。`SsCode { Helix, Strand, Loop }`、`SecondaryStructure { residue_key: String, residue_name: String, code: SsCode }`、`calc_secondary_structure(chain: &AtomGroup) -> Vec<SecondaryStructure>`(内部で`hydrogen_bond.rs`の`calc_backbone_hbonds`を使う)。

**完了の定義**:
1. 上記の基準値表(chain A全21残基、chain B全30残基)と完全一致することをテストすること。
2. ブリッジ/ストランド判定について、逆平行β-シートを模した合成データ(2本の短いストランド、ジグザグ骨格水素結合)でE判定が出ることを検証する専用テストを追加すること(実データにE判定の例がないため)。
3. `cargo clippy`/`cargo fmt`を通すこと。

**Pythonバインディングは本Phaseでは対象外**(PR#22/PR#23完了後、必要であれば別PRとして追加検討する)。

### スコープ外

- §3.3の側鎖水素結合判定・CH-π相互作用検出(`hydrogen_bond.rs`の側鎖拡張・`ch_pi.rs`)は別Phaseとする。
- 8状態DSSP分類(3-10/α/πヘリックスの区別、B/T/S等)への拡張。
- §3.4-3.6のQC結果I/O、§4 C/C++バインディング。
- Pythonバインディング(PyO3)。

### やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない(この機能はPython版に存在しないため)。
- **ブランチ運用ルール(MUST項目)を厳守**: `feature/phase8-prM`ブランチを`develop`から切って作業し、`develop`へは自分でマージせず、レビュー承認を待つ。PR#22を先に、PR#23をその後に。
