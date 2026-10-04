# TASK: 2026.10.2リリース前の個人情報の整理

> 2026.10.2(wheel付きの最初のGitHub Release)の前に、公開される配布物に個人情報や秘密情報が含まれていないかをClaudeが調査した(2026-10-05)。ローカルのパス・トークン・秘密鍵・ローカル設定ファイルは見つからなかったが、次の2点についてユーザーが判断した。

## 調査で見つかったこと

1. **Gaussianの出力ファイル`proteindf_bridge/data/ACE-ALA-NME/*.out`(4ファイル、計約3.2MB)**: 計算機のホスト名(`VEDA`)、ユーザー名(`HIRANO`)、`/opt/g16`・`/tmp/scratch`などのパス、Gaussian社の著作権表示を含む。2022年からgitリポジトリにあるが、以前の配布物(setuptoolsで`data/*.brd`のみを同梱)には入っていなかった。PR#49でmaturinに切り替えた際、パッケージのディレクトリ全体が同梱されるようになり、wheelとsdistに入るようになった。コード・テストからの参照はない。入力ファイル`*.gjf`にはパスやユーザー名は含まれない。
2. **作者のメールアドレス`hiracchi@gmail.com`**: `pyproject.toml`の作者情報と、`docs/locale/ja/LC_MESSAGES/{index,installation,usage}.po`の`Last-Translator`にある。

**ユーザー判断(2026-10-05)**: 1はリポジトリからも削除する(過去の履歴には残る)。2はGitHubのnoreplyアドレス`hiracchi@users.noreply.github.com`(gitのコミット・タグで使っているもの)に変える。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaude(`/code-review`)が担当する。
- `develop`から`fix/release-privacy-cleanup`を切って作業する。**`develop`へは自分でマージしない。** タグやReleaseは作らない。

## 対象

1. `proteindf_bridge/data/ACE-ALA-NME/`の`AAN_cis1.out`・`AAN_cis2.out`・`AAN_trans1.out`・`AAN_trans2.out`を`git rm`で削除する。`*.gjf`と`README.md`は残す。`README.md`に、Gaussianの出力ファイルは個人情報(計算機名・ユーザー名・パス)を含むためリポジトリに置いていないこと、`*.gjf`から再計算できることを書き足す。
2. `pyproject.toml`の作者のメールアドレスと、上記3つの`.po`ファイルの`Last-Translator`のメールアドレスを`hiracchi@users.noreply.github.com`に変える。
3. 配布物に必要のないファイルが今後また入り込まないよう、`pyproject.toml`の`[tool.maturin]`の設定で、パッケージに同梱するデータを確認する。少なくとも、実行時に必要なデータ(`Modeling`が読む`*.brd`など)がwheelに含まれ続けることを確認する。テスト用の構造ファイル(`1HLS.cif`など、PDBの公開データ)を同梱から外すかどうかは、テストや実行時の利用状況を調べて報告するだけにとどめ、変更はしない。

## 完了の定義

1. このブランチから作ったsdistとwheelの中身の一覧を確認し、`*.out`が含まれないこと、`hiracchi@gmail.com`が含まれないこと(`PKG-INFO`・`METADATA`を含む)を、実際のコマンドと出力で示す。
2. `git grep hiracchi@gmail.com`の結果が空であることを示す。
3. 新しい一時環境に`pip install`して、`python -P -m unittest discover -s tests`が通る。`cargo test --workspace`・`cargo clippy --workspace --all-targets -- -D warnings`・`cargo fmt --all -- --check`が通る。
