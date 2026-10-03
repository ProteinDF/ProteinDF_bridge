# TASK: Pythonバインディング比較テスト(phase1)の結合本数テストの更新

> `tests/test_rs_phase1.py`の`test_atom_group_hierarchy_and_bonds`が、developで失敗し続けている(`AssertionError: 2 != 3`、`len(rs_bonds) == len(py_bonds)`の比較)。pyo3の警告対応のレビュー(2026-10-03)で判明した。ユーザー確認のうえ、別タスクとして起票した。

## 原因(Claude確認済み、2026-10-03)

テストは N(0,0,0)・CA(1.45,0,0)・C(2.0,1.4,0) の3原子で結合を判定している。N-C間は約2.44Åで、化学的には結合していない。

- Python版`Bond.setup`(ファンデルワールス半径によるヒューリスティック)は N-CA・CA-C・N-C の3本を結合と判定する。
- Rust版`Bond::setup`は§3.10(2026-09-22)で共有結合半径による判定に変わり、化学的に正しい N-CA・CA-C の2本を返す。

つまり**Rust版が正しく、Python版との違いは§3.10で意図して作ったもの**であり、直すべきはテストである。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaude(`/code-review`)が担当する。
- `develop`から`fix/phase1-bond-count-test`(名前はagyの判断でよい)を切って作業する。
- **`develop`へは自分でマージしない。** 完了したら、ユーザー経由でClaudeに報告し、レビュー承認を待つこと。

## 対象

1. `test_atom_group_hierarchy_and_bonds`の結合本数の比較を、Python版との一致ではなく、Rust版が化学的に正しい2本(N-CA、CA-C)を返すことの検証に置き換える。どの原子の組が結合しているかまで確認する。
2. Python版と結果が異なる理由(§3.10の意図した違い)をテストのコメントに書く。
3. 同じテスト内の、結合本数以外のPython版との比較(原子数・グループ数・重心など)はそのまま残す。
4. 他の`tests/test_rs_phase*.py`に同じ種類の比較(結合本数のPython版との一致)がないか確認し、あれば同様に扱う。なければ、ないことを報告する。

## 完了の定義

1. `tests/test_rs_phase{1,2,3,7}.py`がすべて成功する。実行方法(拡張モジュールのビルド方法、Python環境)と結果を完了報告に書く。参考: Claudeは`cargo build -p proteindf-bridge-py --features extension-module`で作った`.so`を`proteindf_bridge_rs.abi3.so`としてコピーし、`PYTHONPATH`に置いたうえで`uv run --no-project --with numpy --with pyyaml --with msgpack python -m unittest tests/test_rs_phaseN.py`で実行した。
2. Rust側のコードは変更しない(テストの更新のみ)。

## レビュー結果(1回目、2026-10-03、収束・マージ済み)

`fix/phase1-bond-count-test`(`9111556`)をレビューした。変更は`tests/test_rs_phase1.py`のみで、Rust側の変更はない。結合本数の比較を、N-CA・CA-Cの2本を原子の組まで含めて確認する形に置き換え、Python版との違いの理由をコメントに書いている。コメント中の閾値(Python版`vdw_p + vdw_q + 0.4`、Rust版`COVALENT_BOND_TOLERANCE = 0.45`)は実装と一致することを確認した。他のphaseのテストに同じ種類の比較はなかった。`tests/test_rs_phase{1,2,3,7}.py`(10・16・16・4件)がすべて成功することをClaudeが確認した。ユーザー承認のうえ、2026-10-03にdevelopへマージした。**本タスクは完了。**
