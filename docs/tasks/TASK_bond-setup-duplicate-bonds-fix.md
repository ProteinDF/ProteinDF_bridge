# TASK: `Bond::setup()`が既存の結合を重複登録するバグの修正

> ユーザーから「CCDテンプレートとヒューリスティックを両方使う場合の挙動」について質問があり、Claudeが検証した結果、`apply_ccd_bond_templates()`の後に`Bond::setup()`を呼ぶと、同じ原子ペアに対して重複した`BondRecord`が作られることが判明した。本タスクはこれを修正する。

## ブランチ運用(MUST)

- `develop` から `fix/bond-setup-duplicate-bonds` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。
- 依存関係: なし。§3.9(CCDテンプレートDB)・§3.10(共有結合半径)いずれもマージ済みの状態を前提にする。

## 背景・再現手順

`apply_ccd_bond_templates()`→`Bond::setup()`の順で呼ぶと、既にCCDで登録済みの結合に対して`Bond::setup()`が同じ原子ペアをヒューリスティックでも検出し、別のレコードとして追加してしまう(実際に検証済み):

```
CCD適用後:        bonds = [BondRecord { C, O, order: 2 }]
Bond::setup()後:  bonds = [BondRecord { C, O, order: 2 }, BondRecord { C, O, order: 1 }]  ← 重複
```

**原因**: `AtomGroup::add_bond`(内部の`add_bond_direct`、`atom_group.rs`)は無条件に`self.bonds.push(...)`するだけで、同一原子ペアの既存レコードをチェックしていない。`apply_ccd_bond_templates`側には「既存結合を上書きしない」ための`existing_bonds`チェックが実装されている(`atom_group.rs`の`apply_ccd_bond_templates_recursive`参照)が、`Bond::setup()`側(`bond.rs`)には同様のチェックが無い。

**なぜフェーズAのテストで検出されなかったか**: `tests/test_ccd_templates.rs`の`test_apply_ccd_bond_templates_1hls_real_pdb`は、まさにこの「CCD適用後にBond::setup()を呼ぶ」パターンをテストしているが、アサーションが`final_bonds.iter().find(...)`で最初に見つかったレコードの`order`のみを検証しており、`final_bonds.len()`(総件数)や重複の有無は検証していなかった。`get_bond_list()`が挿入順を保持するため、CCD由来のレコード(先に追加された)が`.find()`で最初にヒットし、後から追加された重複レコードの存在に気づけなかった。

**§3.8との関係**: §3.8の推奨利用パターン(`if ag.get_bond_list().is_empty() { Bond::setup(&mut ag)?; }`)は、`Bond::setup()`を「結合情報が一切無い場合のみ」呼ぶことを前提にしており、この前提を守っている限りは重複は起きない。しかし§3.9で確立した「CCDテンプレート→(補完として)Bond::setup()」という組み合わせ方は、この前提(空である場合のみ呼ぶ)の外側にあるパターンであり、現状無防備である。

## 対象

- `Bond::setup()`(`bond.rs`)を、**呼び出し時点で既に登録されている結合(原子パスのペア)をスキップする**ように修正する。`apply_ccd_bond_templates_recursive`(`atom_group.rs`)が使っている「`self.get_bond_list()`から既存の原子パスペアの集合を事前収集し、それに含まれるペアは追加しない」というパターンを踏襲すること(重複実装せず、可能なら共通化を検討してもよい)。
- これにより、`Bond::setup()`は「結合情報が一切無い場合」だけでなく、「一部の結合が既に登録されている状態(CCDテンプレート適用後、または将来的にファイル由来結合の一部だけがある場合等)」で呼んでも安全になる(不足分のみを補うヒューリスティックとして機能する)。

## 完了の定義(Definition of Done)

1. 今回発見した再現手順(CCD適用後にBond::setup()を呼ぶ)を回帰テストとして追加し、`ag.bonds().len()`が重複なく、既存のCCD結合次数(order 2)がそのまま維持されることを検証すること。
2. 既存の`tests/test_ccd_templates.rs`の`test_apply_ccd_bond_templates_1hls_real_pdb`のアサーションを、単なる`.find()`ではなく総結合数の重複が無いこと(または該当原子ペアのレコード数が1件であること)まで検証するように強化すること。
3. 既存の`test_bond_setup`系・`test_covalent_bond_detection`系のテストが引き続き全てパスすることを確認すること(回帰確認)。
4. `cargo clippy` / `cargo fmt` を通すこと。
5. `RUST_PORT_SPEC.md` §3.9または適切な箇所に、この修正内容を追記すること。

## スコープ外

- CCDテンプレート自体・共有結合半径テーブル自体の変更。
- `add_bond`/`add_bond_direct`の汎用的な重複排除(全呼び出し元に影響する変更)は、影響範囲が広いため本タスクでは`Bond::setup()`側での対応に留める。`add_bond`自体の重複排除が必要だと判断した場合は、実装前にユーザー経由でClaudeに相談すること。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- `add_bond`/`add_bond_direct`自体の挙動(重複チェックの有無)を、影響範囲を確認せずに変更しない(他の呼び出し元、特に`apply_ccd_bond_templates`・各フォーマットローダーのCONECT/BONDS処理への影響を先に確認すること)。
