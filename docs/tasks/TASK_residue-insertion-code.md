# TASK: `get_atomgroup()`系ローダーが挿入コード（insertion code）を無視しているバグの修正

> `yui`リポジトリのフェーズ6e-vii実データ検証の過程で発見された。詳細な根本原因分析は`RUST_PORT_SPEC.md` 9章「結(YUI)からの要求リスト」の該当項目（2026-09-23付）を参照。本タスクはその修正案を実装する。

## ブランチ運用(MUST)

- `develop` から `fix/residue-insertion-code` を作成して作業する。
- **`develop` へは自分でマージしない。** 実装・テスト・`cargo clippy`/`cargo fmt`確認が終わったら、ブランチ名と完了内容をユーザー経由でClaudeに報告し、レビュー承認(または修正依頼)を待つこと。
- 依存関係: なし。

## 背景・再現手順

`src/format/pdb.rs`の`Pdb::get_atomgroup()`と`src/format/mmcif.rs`の対応する読み込み処理の両方で、残基キー（`res_key`）の構築が挿入コード（insertion code）を無視している:

- `pdb.rs`（418行目付近）: `let res_key = format!("{}", item.res_seq);` — `PdbRecord.i_code`（`i_code: String`、パース自体は正しく行われている）を全く参照していない。
- `mmcif.rs`（566行目付近）: `let res_key = format!("{res_seq}");` — `pdbx_pdb_ins_code`（`_atom_site.pdbx_PDB_ins_code`、こちらもパース自体は正しく`"."`/`"?"`をフィルタして空文字既定値にしている）を全く参照していない。

同一`res_seq`（または`auth_seq_id`/`label_seq_id`）だが異なる挿入コードを持つ複数の残基（例: 抗体のCDRループ等で一般的な`"52A"`/`"52B"`のような番号付け）を含む実ファイルを読み込むと、`chain.has_group(&res_key)`が既に真になっているため後続の挿入コード違いの残基がサイレントにマージ・上書きされ、データが破損する。

## 対象

### 1. 主経路の`res_key`構築の修正

```rust
// src/format/pdb.rs 内、get_atomgroup()
// 変更前: let res_key = format!("{}", item.res_seq);
let res_key = format!("{}{}", item.res_seq, item.i_code.trim());
// i_code の既定値は半角スペース " " なので trim() が必須（通常ケースでは res_key は
// 従来通り素の数値文字列のままになり、後方互換性が保たれる）。
```

```rust
// src/format/mmcif.rs 内、対応する読み込み処理
// 変更前: let res_key = format!("{res_seq}");
let res_key = format!("{res_seq}{}", item.pdbx_pdb_ins_code);
// pdbx_pdb_ins_code は取得時点で "."/"?" が既にフィルタ済み・空文字既定なので trim 不要。
```

### 2. SSBOND/`_struct_conn`結合解決部分の追従修正（見落としやすい）

主経路の`res_key`だけ直しても、SSBOND/CONECT由来の結合の両端原子が属する残基を`chain.get_group(&res_key)`で引き直す処理が別途あり、**そちらも全く同じフォーマットに揃えないと、挿入コード付き残基に対するSSBOND/CONECT結合が「該当残基が見つからない」形でサイレントに解決失敗する**（新しいキー形式と一致しなくなるため）。ただし、この2箇所は現状、挿入コード自体をパースすらしていないため、まず以下の拡張が必要:

- **`pdb.rs`の`SsBondRecord`構造体**（55行目付近）: `icode1: String`, `icode2: String`フィールドを追加し、SSBOND行パース箇所（175行目付近、`seq_num1`/`seq_num2`をパースしている箇所）で対応するカラムから挿入コードもパースする。wwPDB PDB Format v3.3のSSBONDレコード仕様（カラム22がicode1、カラム36がicode2、いずれも1文字）を確認の上、既存の0-indexed `slice_chars(start, end)`呼び出しパターン（`seq_num1 = slice_chars(&chars, 17, 21)`は1-indexedカラム18-21に対応）に合わせて実装すること。
  - 修正後、455行目付近の`let res_key1 = format!("{}", ssbond.seq_num1);` / `let res_key2 = format!("{}", ssbond.seq_num2);` を、主経路と同じ`"{}{}"` + `trim()`形式に更新する。
- **`mmcif.rs`の`StructConnRecord`構造体**（206行目付近）: `ptnr1_ins_code: String`, `ptnr2_ins_code: String`フィールドを追加し、`from_row()`（219行目付近）で`_struct_conn.pdbx_ptnr1_PDB_ins_code` / `_struct_conn.pdbx_ptnr2_PDB_ins_code`（実際のmmCIFカテゴリ定義でのタグ名を要確認）から、主経路の`pdbx_pdb_ins_code`と同じ`"."`/`"?"`フィルタ・空文字既定ロジックでパースする。
  - 修正後、606行目付近の`let res_key1 = format!("{seq1}");` / `let res_key2 = format!("{seq2}");` を、主経路と同じ形式に更新する。

修正時は、ファイル内で`res_seq`/`auth_seq_id`/`label_seq_id`から`res_key`相当の文字列を組み立てている箇所を全て洗い出し、共通のヘルパー関数（例: `fn build_residue_key(seq: i32, ins_code: &str) -> String`）に切り出すことを推奨する（`pdb.rs`・`mmcif.rs`それぞれに1つずつでよい。両ファイル間で共有する必要はない）。

## 完了の定義(Definition of Done)

1. 同一`res_seq`だが異なる挿入コードを持つ複数残基（例: `"52"`, `"52A"`, `"52B"`）を含む合成PDBフィクスチャ・合成mmCIFフィクスチャをそれぞれ作成し、`get_atomgroup()`で読み込んだ結果、3つの残基が別々のグループとして正しく存在すること（マージ・上書きされていないこと）を検証する回帰テストを追加すること。
2. 挿入コード付き残基に対してSSBOND（PDB）/`_struct_conn`（mmCIF）由来のジスルフィド結合が正しく解決されること（該当残基が見つからずサイレントに結合が欠落しないこと）を検証する回帰テストを追加すること。
3. 挿入コードを持たない通常のファイル（既存の`tests/data/1hls.pdb`等）に対する既存テストが全て引き続きパスすること（`res_key`のフォーマットが空の挿入コードに対して後方互換であることの回帰確認）。
4. `cargo clippy --workspace --all-targets -- -D warnings` / `cargo fmt --check` を通すこと。
5. `RUST_PORT_SPEC.md` 9章の該当項目に、実施内容・完了ステータスを追記すること。

## スコープ外

- `yui`側の`core::atom_group::parse_residue_key()`の変更（yui側は既に挿入コード付き残基キー`"52A"` → `(52, Some('A'))`を扱える設計になっているため、本タスクでの対応は不要）。
- 挿入コード以外の残基キー構築ロジックの変更。

## やってはいけないこと

- 既存Pythonコード(`proteindf_bridge/`)は変更しない。
- 挿入コードを持たない既存ファイルに対する`res_key`のフォーマット（素の数値文字列）を変更しない（後方互換性を壊さないこと）。
