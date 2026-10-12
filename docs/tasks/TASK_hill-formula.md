# TASK: 化学式をHill方式にする (PR#57)

> Rust版の化学式(`AtomGroup::get_formula`)を、原子番号順からHill方式に変える。ユーザー確認済み(2026-10-12): `get_formula`そのものをHill方式にする(新しい関数は作らない)、個数が1の元素は数字を省略する(`CH4`)。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaudeが担当する。
- `develop`から`feature/hill-formula`を切って作業する。**`develop`へは自分でマージしない。** タグやReleaseは作らない。
- 完了報告は、起動時に指示されたパス(リポジトリの外)に、決まった見出しの短いファイルとして書く。コマンドの出力は要約・再構成せず、実際の出力行をそのまま貼る。
- テスト名・ファイル名・行番号は、必ず実際のファイルを開いて確認してから書く。確かめていない推測を事実として書かない。

## 仕様(Hill方式)

1. 炭素を含む場合: `C`、続いて`H`(含む場合)、続いて残りの元素を元素記号のアルファベット順。
2. 炭素を含まない場合: `H`も含めて、すべての元素を元素記号のアルファベット順。
3. 個数が1の元素は数字を書かない(`CH4`、`NaCl`、`H2O`)。
4. 元素記号が不明な原子(`X`)は、これまでどおり最後に`X{個数}`(個数1なら`X`)を付ける。
5. 原子がない場合は空文字列。
6. アルファベット順は元素記号の文字列の順(`Cl` < `Na` < `O` など、大文字小文字は元素記号の表記どおりで比べる)。

例: 1hlsは`C254H379N67O76S6`、水は`H2O`、塩化ナトリウムは`ClNa`、メタンは`CH4`。

## 対象

1. `rust/crates/proteindf-bridge/src/atom_group.rs`の`get_formula`を上の仕様に変える。ドキュメントコメントにHill方式であることと規則を書く。
2. Rustのテスト: 既存の`test_get_formula`(期待値`"H2C1"`)を新しい仕様に直し、上の例(炭素あり、炭素なし、個数1の省略、`X`あり、空)を網羅するテストを追加する。
3. Pythonバインディング: `AtomGroup.get_formula`・`AtomGroup.formula`はコアをそのまま呼ぶので、変更は不要のはず。Pythonのテストで、`proteindf_bridge.rs`の`get_formula()`・`formula()`がHill方式を返すことを確かめる。
4. 純Python版との一致テスト: `tests/test_rs_phase1.py`の`self.assertEqual(rs_root.get_formula(), py_root.get_formula())`は、純Python版(原子番号順、変更しない)と一致しなくなる。このアサートを、Rust版がHill方式で正しい値を返すことの確認に置き換え、純Python版とは意図して違うことをコメントに書く。純Python版のパッケージ・テスト(`tests/test_atomgroup.py`)は変更しない。
5. 基準ファイル: `tests/cli_golden/brd_formula.stdout.txt`を`C254H379N67O76S6: MW=...`に更新する(化学式以外の部分が変わらないことを確かめる)。`docs/cli-migration-notes.md`の`brd-formula.py`の記録を「Hill方式に変更(本タスク)」に更新する。

## 守ること

- 純Python版のパッケージは変更しない。
- `RUST_PORT_SPEC.md`・`docs/tasks/`の過去の記録は書き換えない(本ファイルへの追記もしない)。
- 意味のまとまりごとにコミットする。

## 完了の定義

1. 対象1〜5に対応している。
2. 次のテストが通る: 新しい一時環境にパッケージを`pip install`して`python -P -m unittest discover -s tests`、`cargo nextest run --workspace --status-level fail --final-status-level fail`、`cargo test --doc`、`cargo clippy --workspace --all-targets -- -D warnings`、`cargo fmt --all -- --check`。
3. `brd-formula.py`を1hlsの`.brd`で実行した実際の出力を報告に貼る。

## PR#57 完了(2026-10-12)

`42f553f`・`fdf18a3`・`0051080`をレビューし、1回で収束した。`get_formula`の実装が仕様1〜6どおりであること、Rustのテストが仕様の例を網羅すること、`brd-formula.py`の1hlsの出力が`C254H379N67O76S6: MW=5779.662000000002`で基準ファイルの違いが化学式だけであること、純Python版のパッケージと`tests/test_atomgroup.py`が変更されていないことをClaudeが確認した。Pythonテスト一式237件OK(新しい環境)、`cargo nextest` `320 tests run: 320 passed, 2 skipped`、doctest・clippy・fmt OK。`--no-ff`でdevelopにマージした(`52c4414`)。リリースの時期はユーザーの指示待ち。

## リリース 2026.10.4(2026-10-12)

PR#57(Hill方式)を含めてリリースした。Releaseワークフローが成功し、wheel 3つとsdistが添付された。Claudeが添付ファイルをすべて展開して確認した。含まれるパスはCIのランナー(`/home/runner`・`/Users/runner`)のみ、メールアドレスはnoreplyと、SBOM(`dist-info/sboms/*.cyclonedx.json`)に載る依存クレートの公開作者情報のみで、ローカルのパス・個人のメールアドレスは含まれていない。`--find-links`でインストールし、`brd-formula.py`が1hlsで`C254H379N67O76S6: MW=5779.662000000002`を出すことを確認した。
