# TASK: mmCIF読み込み時に解決できなかった`_struct_conn`結合を黙って捨てない

> PR#43のレビュー(2026-10-03)で判明した。`SimpleMmcif::get_structure_atomgroup_for_block`は、`_struct_conn`(`disulf`・`covale`)の相手原子がモデル内に見つからない場合、その結合を黙って捨てている。`disulf`では元からあった挙動だが、PR#43で`covale`にも広がった。「黙って不完全な結果を作らない」方針(`docs/rust-port-handoff.md`教訓13)に合わないため、別タスクとして起票した(ユーザー承認、2026-10-03)。

## 見つからなくなる典型的な原因

- altLocの選択(`select_altloc`、既定は`A`)で、相手の原子が除外された。
- `select_model`で特定のモデルだけを読んだ(この場合は結合の相手は見つかるはずなので、原因にならないことを確認する)。
- 鎖ID・残基番号・挿入コードの解決方法が`_atom_site`と`_struct_conn`で食い違っている(auth/labelの扱いの違いなど)。
- ファイル自体の不整合。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaude(`/code-review`)が担当する。
- **PR#44(`TASK_mmcif-writer.md`)の完了後に着手する。**
- `develop`から専用のブランチを切って作業する。**`develop`へは自分でマージしない。**

## 対象

1. 解決できなかった`_struct_conn`の行を、呼び出し側が確認できる形で返す。エラーにするか、警告として報告するかは、次の点を考慮して設計し、完了報告で理由を説明する。
   - altLocで原子を除外した結果、結合が解決できなくなるのは実データでは普通に起こりうる(エラーにすると、そうした実構造が読めなくなる)。
   - 既存の`get_structure_atomgroup`のシグネチャ(戻り値`Result<AtomGroup>`)を壊さないこと。例えば、詳細を返す別のメソッドを追加し、既存のメソッドはそれを呼ぶ形にする。
2. どの行(`_struct_conn.id`)が、なぜ(どちらの相手の、鎖・残基・原子のどれが見つからなかったか)解決できなかったのかがわかるようにする。

## 完了の定義

1. 実データまたは合成データで、相手原子がaltLocの選択によって除外されるケースを作り、その結合が報告されることをテストで確認する。
2. 既存の実データ(1HLS・2FB4・1WCT等)では、解決できない結合が0件であることを確認する。
3. `cargo test --workspace`、`cargo clippy --workspace --all-targets -- -D warnings`、`cargo fmt --all -- --check`が通る。
