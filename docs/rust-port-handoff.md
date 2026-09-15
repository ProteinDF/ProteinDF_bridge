# proteindf-bridge Rust移植 — antigravity向け作業指示

本ドキュメントは、`RUST_PORT_SPEC.md` に基づくRust移植作業を実装担当(antigravity)に委任するにあたっての、フェーズ単位の作業指示を記録する。Claude(このリポジトリでのレビュー担当)は各PRを本ドキュメント・`RUST_PORT_SPEC.md`・`SPEC.md` と突き合わせて仕様適合チェックを行う。

## 運用ルール

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
