# TASK: mmCIF CCDパースの結合次数・座標欠損の是正

> ユーザーから「wwPDB Chemical Component Dictionary (CCD, https://www.wwpdb.org/data/ccd) のパースが現在の実装で可能か」という調査依頼があり、Claudeが`format/mmcif.rs`を調査した結果、構造的には対応済み(複数`data_`ブロックを含むCCD配布ファイルの読み込み・`AtomGroup`構築は動作する)だが、実データで問題になる2つのギャップが見つかった。本タスクはこれを修正する。

## ブランチ運用(MUST)

- `develop` から `fix/mmcif-ccd-bond-order` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。
- 依存関係: なし。Phase 10とは無関係の独立した不具合修正。

## 背景・調査結果

`SimpleMmcif::get_atomgroup`(CCD形式、`_chem_comp`/`_chem_comp_atom`/`_chem_comp_bond`カテゴリ)は、複数コンポーネントの`data_`ブロックを含むファイル(wwPDBの完全なCCD配布ファイル`components.cif`相当)を正しくパースできることを、合成データ(ALA・架空の芳香族コンポーネントBNZの2ブロック)で実際に検証済み(`get_molecule_names()`が両ブロックを返し、`get_atomgroup()`がそれぞれ正しく`AtomGroup`を構築する)。

しかし、以下2点のギャップが見つかった:

### 1. 芳香族結合(`_chem_comp_bond.value_order = "AROM"`)が結合次数0として記録される(要修正、優先度高)

`mmcif.rs`(432〜459行目付近)の結合次数変換ロジックは以下の通り:
```rust
let bond_order = match bond_order_str {
    "SING" => 1,
    "DOUB" => 2,
    "TRIP" => 3,
    _ => 0,
};
```
`AROM`(芳香族)はこのmatchのいずれにも該当せず、デフォルトの`_ => 0`に落ちる。実際に合成データ(`BNZ C1 C2 AROM`)で検証したところ、`BondRecord { atom1_path: "C1", atom2_path: "C2", order: 0 }`が生成されることを確認した。

**再現手順**: `_chem_comp_bond.value_order`が`AROM`の行を含むCCDデータを`get_atomgroup()`でパースすると、当該結合の`order`が0になる。`order: 0`はこのコードベースの他の箇所(`bond.rs`の`bondmat`等)で「結合なし」を意味する値として扱われており、紛らわしい・下流処理で無視されるリスクがある。

wwPDB CCDの多くのリガンド(ベンゼン環を含む化合物、核酸塩基、多くの補酵素等)は芳香族結合を持つため、実データでは高頻度に発生する問題。

**修正方針**: `mol2.rs`の`@<TRIPOS>BOND`パーサに既に前例がある(`"ar" => 1`、SYBYL芳香族結合タイプをbond order 1として登録)。**同じ方針でCCDの`AROM`もbond order 1として登録すること**(整数の結合次数しか表現できない現状の`BondRecord`スキーマでの整合性を優先。将来的に`is_aromatic`フラグ等を追加する設計変更は本タスクのスコープ外)。あわせて`QUAD`(四重結合、稀だが仕様上存在する)も`4`として登録すること。それ以外の未知の`value_order`値が来た場合にフォールバックする`_ => 0`は、意図的なフォールバックであることをdocコメントで明記すること(サイレントに見えないようにする)。

### 2. 座標が`ideal`/`model`両方とも欠損している場合、エラーにならず座標(0,0,0)に静かにフォールバックする(要修正、優先度中)

`mmcif.rs`の`get_coordinate`(659〜675行目付近)は、`pdbx_model_Cartn_{axis}_ideal`と`model_Cartn_{axis}`の両方が存在しない、またはパース不能(`"?"`等)な場合に`None`を返す。呼び出し元(`extract_atoms_and_name`、640〜645行目付近)は`if let (Some(x), Some(y), Some(z)) = (x, y, z) { atom.xyz = ...; }`という構造になっており、**いずれか欠けていれば`atom.xyz`は`Atom::new()`のデフォルト値(原点)のまま、エラーも警告もなく処理が続行される**。

このプロジェクトでは`unwrap_or`等によるサイレントなエラー握りつぶしが繰り返し問題になってきた経緯があり(`docs/rust-port-handoff.md`の教訓カタログ参照)、同じパターンを踏襲すべきではない。

**修正方針**: 座標が両方とも取得できない場合は、その原子をサイレントに原点へフォールバックさせず、**明示的なエラーを返す**ように修正すること(`Result`を返す設計に変更するか、該当原子をスキップしてログ出力する等、既存のエラーハンドリング方針に合わせること。判断に迷う場合はコード変更前にユーザー経由でClaudeに相談すること)。

## 完了の定義(Definition of Done)

1. 芳香族結合(`AROM`)を含む合成CCDデータで、結合次数が1として正しく登録されることを検証する回帰テストを追加すること。四重結合(`QUAD`)についても同様に4になることを検証すること。
2. 座標が両方とも欠損しているCCDデータ(合成)で、期待通りエラーになる(またはスキップされてログが残る)ことを検証する回帰テストを追加すること。
3. 既存のCCDテスト(`ALA.cif`等、既存フィクスチャ)が引き続き全てパスすること(回帰確認)。
4. 複数`data_`ブロックを含む合成CCDデータ(例: 2〜3コンポーネント)を新規フィクスチャとして追加し、`get_molecule_names()`・各ブロックの`get_atomgroup()`が正しく動作することを検証するテストを追加すること(既存のALA.cif単体テストは単一ブロックのみなので、複数ブロックのケースが未検証)。
5. `cargo clippy` / `cargo fmt` を通すこと。

## スコープ外

- `BondRecord`に`is_aromatic`のような専用フィールドを追加する設計変更(将来必要になれば別タスク)。
- 実際のwwPDB完全CCD配布ファイル(`components.cif`、数百MB規模)でのパフォーマンス検証(今回は正しさの検証のみが目的)。
- `_chem_comp.type`/`.formula`等、`AtomGroup`構築に不要なメタデータカテゴリの追加パース。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない(CCDのマルチブロック対応はPython版に対応物がない可能性が高いため、着手前に`proteindf_bridge/mmcif.py`を確認し、対応物があれば1:1方針を優先すること)。
- 座標欠損時のエラーハンドリング方針(エラーを返す vs スキップ)を独断で決めず、判断に迷えばClaudeに相談すること。
