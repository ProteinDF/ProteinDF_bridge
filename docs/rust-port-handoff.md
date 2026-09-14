# pdf-bridge Rust移植 — antigravity向け作業指示(Phase 1)

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
    └── pdf-bridge/        # コアライブラリ(Phase 1の対象)
        ├── Cargo.toml
        └── src/
            ├── lib.rs
            ├── error.rs
            └── periodic_table.rs
```

単一クレート直下構成ではなく `crates/` 配下にコアクレートを置く構成を採用する。理由: §4のバインディング方針(C ABI用cdylib、PyO3用Pythonモジュール)は将来的にコアクレートとは別クレート(`crates/pdf-bridge-capi/`、`crates/pdf-bridge-py/` 等)として追加する想定のため、最初から `crates/` 構成にしておく。

## Phase 1: 基盤・データモデルの1:1移植(今回のスコープ)

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

Phase 1の全PRがマージされ `rust-port` ブランチ上でビルド・テストが通った時点で、Claudeが `RUST_PORT_SPEC.md` §2対応表との突き合わせレビューを行う。問題なければ `rust-port` → `main` へのマージを提案し、Phase 2(フォーマットI/O)の指示を別途作成する。
