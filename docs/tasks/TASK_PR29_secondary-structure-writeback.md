# TASK_PR29: 二次構造情報の`AtomGroup`への書き戻し

> 本タスクは `docs/rust-port-handoff.md` の「Phase 10」節(RUST_PORT_SPEC.md §9対応)から抽出したものです。全体の背景・優先順位・他PRとの関係は同ドキュメントを参照してください。

## ブランチ運用(MUST)

- `develop` から `feature/phase10-pr29` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。承認が出るまで次のタスクの実装に着手しない。
- 依存関係: なし(PR#28とは独立)。着手推奨順序ではPR#28の次。

## 背景

2026-09-19に `RUST_PORT_SPEC.md` §9へ追記された新規「高」優先度項目。`calc_secondary_structure(chain: &AtomGroup) -> Vec<SecondaryStructure>`(Phase 8 PR#23で実装済み、`secondary_structure.rs`)は結果を別のVecとして返すのみで、`AtomGroup` ツリー自体には反映されない(`AtomGroup` に汎用メタデータフィールドが無いため)。一方 `Bond::setup()` は `mol.add_bond(...)` で結果を `AtomGroup` 自体に書き戻す設計になっており、一貫していない。YUI側は現状、residueのpath文字列をキーとする一時的なサイドマップで代替している(bridge側の対応までの暫定措置)。

## 対象

1. `bonds: Vec<BondRecord>` と同格の、residueレベルの `AtomGroup` が持つ専用フィールド `secondary_structure: Option<SsCode>` を追加する(汎用メタデータ袋ではなく、`bonds` と同じ「specific typed field」パターンをYUI側は希望している)。
2. `calc_secondary_structure` と対になる `apply_secondary_structure(chain: &mut AtomGroup)`(算出結果を対応するresidueグループの `secondary_structure` フィールドに書き戻す関数)を追加する。

## 完了の定義(Definition of Done)

1. `apply_secondary_structure` を呼んだ後、chain内の各residueグループの `secondary_structure` フィールドが、`calc_secondary_structure` が返すVecの対応するエントリと一致することをテストすること(Phase 8の基準値データ、`1hls.pdb` のchain A/Bで検証)。
2. `AtomGroup` の `Clone`・マージ演算(`merge`/`BitAnd`/`BitOr`/`BitXor`)が新フィールドを正しく扱う(消えない・上書きロジックが妥当)ことを確認すること。**Phase 1是正事項5で `bonds` フィールドのマージ漏れが実際にバグとして見つかった前例があるため、同じ轍を踏まないこと。**
3. `cargo clippy` / `cargo fmt` を通すこと。

## スコープ外

- 8状態DSSP分類への拡張(Phase 8のスコープ外のまま)。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない(この機能はPython版に存在しないため)。
- **`AtomGroup` のマージ・集合演算全て(`merge`/`BitAnd`/`BitOr`/`BitXor`/`Clone`)を洗い出してから着手すること。** 新フィールド追加のたびに一部の演算だけ対応漏れするパターンがPhase 1是正事項5で実際に起きている。
- Phase 10の他タスク(schema検証・ファイル由来結合・`Bond::setup()`のスケーラビリティ等)には手を出さない。

## レビュー結果(2026-09-19、ブロッカーなし・要修正2件)

`feature/phase10-pr29`(コミット`9a22d25`)をClaudeがレビューした。ビルド・`cargo clippy -D warnings`・`cargo fmt --check`・`cargo test --workspace`はいずれも問題なくパス。`apply_secondary_structure`は実データ(`1hls.pdb`のchain A/B)で`calc_secondary_structure`の結果と一致することを確認済み。Phase 1是正事項5で起きた「bondsフィールドのマージ漏れ」の再発は無く、`merge`/`BitAnd`/`BitXor`いずれにも`secondary_structure`のハンドリングがある(`BitOr`は`clone()+merge()`委譲のため自動的にカバーされている)。**マージをブロックする実バグはない**が、以下2件のMinor〜Moderateな指摘がある。

1. **設計一貫性(`atom_group.rs`、131行目付近)**: `secondary_structure`フィールドが`pub`になっており、同格のはずの`bonds`フィールド(private、`bonds()`/`set_bonds()`経由でのみアクセス)と一貫していない。タスク要件は「`bonds`と同格の専用フィールド」だったが、このPR自身が追加した`secondary_structure()`/`set_secondary_structure()`ゲッター・セッターを迂回できてしまう。実際、追加されたテスト自身も`res1.secondary_structure = Some(SsCode::Helix)`のように直接フィールドを書き換えており、セッターを使っていない。
   - **提案**: フィールドをprivateにし、`secondary_structure()`/`set_secondary_structure()`経由でのみアクセスする形に統一する。テスト側の直接代入もセッター呼び出しに直すこと。

2. **集合演算の仕様(`atom_group.rs`、944行目・1069〜1075行目)**: `BitAnd`は`self.secondary_structure.or(rhs.secondary_structure)`、`BitXor`は明示的なmatch式だが、**両方に値があり、かつ値が食い違う場合**(例: 片方がHelix、もう片方がStrand)、どちらも無条件に`self`側の値を採用する。実際に検証すると`Helix ∩ Strand → Some(Helix)`となり、`bonds`/`atoms`の交差判定(両方に存在するものだけ残す)とは異なる「食い違いを無視して片方を優先する」挙動になっている。追加されたテスト(`test_secondary_structure_merge_and_set_operations`)は「両方が同じ値」のケースしか検証しておらず、この食い違いケースは未テスト。
   - **提案**: 意図的な単純化であれば、その旨をdocコメントに明記する。そうでなければ`bonds`同様「一致する場合のみ残す」(`BitAnd`なら食い違い時は`None`)形に揃え、食い違いケースの回帰テストを追加する。

### 完了の定義(対応する場合)

1. 上記1・2への対応方針を決め、実装する(またはdocコメントで意図的な設計として明記する)。
2. 食い違いケース(`BitAnd`/`BitXor`で両者が異なる`Some`値を持つ場合)の回帰テストを追加する。
3. `cargo clippy` / `cargo fmt` / `cargo test`を通すこと。
4. 対応後、同じ`feature/phase10-pr29`ブランチに追加コミットし、再度ユーザー経由でClaudeにレビュー依頼すること。

## レビュー指摘対応結果 (2026-09-19)

指摘1・2への対応を完了しました。

1. **フィールド private 化とゲッター/セッター統一**:
   - `AtomGroup` の `secondary_structure` フィールドを private に変更し、`bonds` と同様に `secondary_structure(&self)` / `set_secondary_structure(&mut self, ...)` 経由でのみアクセスする設計に統一。
   - `apply_secondary_structure` およびテストコード内の直接フィールド代入をすべて `set_secondary_structure` / `secondary_structure()` 呼び出しに置換。
2. **集合演算（`BitAnd` / `BitXor`）での食い違いケースの厳密化**:
   - `BitAnd`（交差）: 両者が同一の `Some` 値を持つ場合のみその値を残し、値が食い違う場合や片方にしかない場合は `None` となるよう修正。
   - `BitXor`（対称差）: 片方のみが `Some` を持つ場合のみその値を残し、両者が値を持つ場合（一致・不一致問わず）は相殺されて `None` となるよう修正。
3. **回帰テスト追加**:
   - `test_secondary_structure_merge_and_set_operations` に、異なる二次構造値（`Helix` と `Strand`）を持つグループ同士の `BitAnd` (-> `None`)、`BitOr` (-> `Some(Strand)`)、`BitXor` (-> `None`) のテストを追加。

