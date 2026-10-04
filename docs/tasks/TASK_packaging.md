# TASK: 1つのパッケージへの統合とwheel配布 (PR#49〜50)

> 純Python版`proteindf_bridge`とRust版のPythonバインディング`proteindf_bridge_rs`を1つのパッケージに統合し、ビルド済みのwheelをGitHub Releaseで配る。方針はユーザー確認済み(2026-10-04)で、設計の根拠は`RUST_PORT_SPEC.md` §4.4にある。**着手前に§4.4を必ず精読すること。** 本ファイルは実行チェックリストである。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaude(`/code-review`)が担当する。
- 着手順はPR#49 → PR#50。PR#50はPR#49がdevelopにマージされてから始める。
- `develop`から`feature/packaging-prNN`(名前はagyの判断でよい)を切って作業する。
- **`develop`へは自分でマージしない。タグを打たない・プッシュしない。** 完了したら、ユーザー経由でClaudeに報告し、レビュー承認を待つこと。
- **完了報告では、コマンドの出力を要約・再構成せず、実際の出力行をそのまま貼ること。**
- 既存の`.venv`は変更しない。検証用の環境は`uv venv`などで一時的に作り、作業後に削除する。

## PR#49: パッケージの統合

### 対象

1. **ルートの`pyproject.toml`をmaturinに切り替える。**
   - `[project]`: 名前は`proteindf_bridge`。バージョン(`2026.10.1`など、暦ベース)・説明・作者・ライセンス(GPLv3)・依存パッケージ(`msgpack`・`pyyaml`・`numpy`など、現在の`setup.cfg`の`install_requires`)・`docs`の追加指定(extras)を、`setup.cfg`から移す。`requires-python = ">=3.9"`とし、分類(classifiers)のPythonバージョンも3.9以降に直す。
   - `[tool.maturin]`: Rustのクレートは`rust/crates/proteindf-bridge-py/Cargo.toml`を指す。拡張モジュールは`proteindf_bridge.rs`としてパッケージ内に置く(Rust側の`#[pymodule]`の名前もこれに合わせる)。`pyo3/extension-module`を有効にする。`proteindf_bridge/data/*.brd`がwheelに含まれること。
   - **バージョン番号の置き場所**: `pyproject.toml`と`proteindf_bridge/_version.py`の両方に書く場合は、2つが一致していることを確認するテストを追加する(リリース時に片方だけ更新する事故を防ぐ)。Rustクレートのバージョン(0.1.0)はこれとは別のままでよい。
2. **スクリプト**: 現在`setup.cfg`の`scripts=`でインストールしている全スクリプト(コメントアウトされているものを除く)が、これまでと同じ名前でインストールされ、実行できること。スクリプトを関数(エントリーポイント)に書き換えるのではなく、maturinのデータディレクトリ(`<data>/scripts/`)の仕組みでそのままインストールする方法を第一候補にする。ディレクトリの移動が必要なら`git mv`で行う。
3. **古いビルド設定の整理**: `setup.cfg`・`setup.py`(と、`rust/crates/proteindf-bridge-py/pyproject.toml`)は、二重の定義にならないよう削除する。残す必要があるものがあれば、理由を報告する。
4. **テスト**: `tests/test_rs_*.py`の`import proteindf_bridge_rs`を`proteindf_bridge.rs`に書き換える。テストの実行方法を、パッケージを一時的な環境に`pip install -e .`(またはwheelをビルドしてインストール)してから`python -m unittest discover -s tests`を実行する形に変え、その手順を`CONTRIBUTING.md`(または適切な開発者向け文書)に書く。
5. **ドキュメント**: `README.md`・`docs/installation.md`のインストール手順を更新する(ソースからのインストールにはRustツールチェーンが必要になること、wheelでの導入はPR#50で用意すること)。`docs/`のAPIリファレンスなどで`proteindf_bridge_rs`に触れている箇所があれば直す。
6. **ドキュメントのCI**(`.github/workflows/docs.yml`の`pip install -e .[docs]`)がmaturinへの移行後も動くこと。GitHubのUbuntuランナーにはRustが入っているが、必要なら`dtolnay/rust-toolchain`などで明示する。

### 完了の定義

1. `maturin build --release`(または`pip wheel .`)でwheelが作れ、それを**新しい一時環境**にインストールして次が通る。
   - `python -m unittest discover -s tests`(既存の純Python版のテストと`tests/test_rs_*.py`のすべて)
   - `import proteindf_bridge`、`import proteindf_bridge.rs`、`proteindf_bridge.__version__`がパッケージのバージョンと一致すること
   - 代表的なスクリプト(`pdb2brd.py`・`brd2pdb.py`・`neutralize.py`など)が`--help`などで起動できること。インストールされたスクリプトの一覧を報告する
   - wheelに`proteindf_bridge/data/*.brd`が含まれていること(`unzip -l`などの出力を報告する)
2. `pip install -e .`(編集可能インストール)でも同じテストが通る。
3. `cargo test --workspace`、`cargo clippy --workspace --all-targets -- -D warnings`、`cargo fmt --all -- --check`が通る。
4. ドキュメントのビルド(`sphinx-build`、英語・日本語)がローカルで通る。

## PR#50: wheelのビルドとGitHub Releaseへの自動公開

### 対象

1. **リリース用のワークフロー**(例: `.github/workflows/release.yml`)を追加する。
   - きっかけ: バージョンタグ(既存のタグの形式。例: `2026.10.2`)のプッシュ。加えて、Releaseを作らずにビルドとテストだけを行う手動実行(`workflow_dispatch`)を用意する。
   - wheel: Linux x86_64・Linux aarch64(manylinux2014)、macOS(Apple Silicon)。maturin公式の`PyO3/maturin-action`を使う。ソース配布物(sdist)も作る。
   - テスト: ビルドしたwheelを各環境(少なくともLinux x86_64とmacOS。aarch64はエミュレーションで可能なら)でインストールし、`python -m unittest discover -s tests`を実行する。
   - 公開: すべて成功したら、タグ名のGitHub Releaseを作り、wheelとsdistを添付する。Releaseの本文には、注釈付きタグのメッセージ(変更内容の要約)を使う。
   - 使うGitHub Actionsはバージョンを固定する。必要な権限(`contents: write`など)は最小限にする。
2. **通常のCI**(例: `.github/workflows/ci.yml`): `develop`・`main`へのプッシュとプルリクエストで、Linux x86_64のwheelをビルドしてテストを実行する。`cargo test`・`clippy`・`fmt`も実行する。
3. **インストール手順の文書化**: `README.md`・`docs/installation.md`に、GitHub Releaseのwheelからインストールする方法を書く。`pip install <wheelのURL>`に加えて、`pip install proteindf_bridge --find-links <Releaseのページ>`のように、pipに環境に合うwheelを選ばせる方法が使えるかを実際に確かめ、使えるならそれを第一の方法として書く。

### 完了の定義

1. 手動実行(`workflow_dispatch`)でワークフローを実際に動かし、全環境のwheelのビルドとテストが成功することを確認する(実行のURLと結果を報告する)。**この段階ではタグを打たず、Releaseも作らない。**
2. 手動実行で作られたwheel(Actionsの成果物)を手元(Linux x86_64)の新しい一時環境にインストールし、テストとスクリプトの起動が通ることを確認する。
3. タグをきっかけにしたReleaseの作成は、PR#50のマージ後、次のリリースのときに、Claudeがユーザーの承認を得て最初に実行し、確認する。

## 全PR完了後

- `RUST_PORT_SPEC.md` §4.4に「実施内容・検証」と既知の限界を追記する。
- `docs/rust-port-handoff.md`に記録する。
- PyPIのアカウントが用意できたら、PyPIへの公開(Trusted Publishing)を別タスクとして追加する。
- スクリプトをRust版を使う形に書き換える計画を、別に立てる。

## PR#49 レビュー結果(1回目、2026-10-04、要修正)

`feature/packaging-pr49`(`b6f066a`〜`8d343ae`の5コミット)をレビューした。Claudeが手元で`maturin build`・`maturin sdist`を実行し、それぞれを新しい一時環境にインストールして確認した。

- **wheel**: 問題なし。旧`setup.cfg`の32本のスクリプトがそのままインストールされ(一覧が旧`setup.cfg`と完全に一致することを確認)、`pdb2brd.py --help`が起動する。`python -P -m unittest discover -s tests`で197件すべて成功。`proteindf_bridge.__version__`は`2026.10.1`、`import proteindf_bridge.rs`も通る。`proteindf_bridge/data/*.brd`もwheelに含まれる。
- 移動・削除したファイル(`scripts/` → `proteindf_bridge.data/scripts/`、`setup.cfg`・`setup.py`・バインディング側の`pyproject.toml`の削除)と、バージョン一致のテスト(`tests/test_version.py`)は妥当。

### 修正依頼

1. **【不具合】ソース配布物(sdist)にスクリプトのディレクトリ`proteindf_bridge.data/`が含まれず、sdistからのインストールが失敗する。** Claudeが`maturin sdist`で作ったsdistを新しい環境に`uv pip install`したところ、`Caused by: No such data directory .../proteindf_bridge.data`でビルドが失敗した。§4.4は「wheelのない環境ではソースからビルドできる」ことを求めており、PR#50でもsdistをReleaseに置く。`[tool.maturin]`の`include`にsdist向けの指定(例: `{ path = "proteindf_bridge.data/**/*", format = "sdist" }`)を加えるなどで直す。**sdistを新しい一時環境にインストールし、スクリプトがインストールされ、テストが通ることを確認する**(完了の定義に追加)。
2. **【文書】`README.md`・`docs/installation.md`の「`rustc` 1.80+」に根拠がない。** Cargo.tomlに`rust-version`(最小対応バージョン)は宣言されていない。実際に確かめた最小バージョンがあるならCargo.tomlに`rust-version`として宣言して文書と一致させ、なければ具体的な数字を書かない(「最新の安定版のRust」など)。
3. **【文書】利用者向けの文書(`README.md`・`docs/installation.md`)に「PR#50」という内部の作業番号が書かれている。** 利用者には意味がないので、「今後GitHub Releaseで提供予定」のような書き方にする。
4. **【文書】`docs/usage.md`の冒頭が、スクリプトが`scripts/`の下にインストールされると書いている。** 実際には`pip install`で実行ファイルとして(環境の`bin/`に)インストールされる。現状に合わせて直す。

### 完了の定義(修正後)

1. 上記1〜4に対応し、同じブランチに追加コミットする。
2. PR#49の完了の定義1〜4に加え、sdistからのインストールとテストが通る。

## PR#49 レビュー結果(2回目、2026-10-04、収束・マージ済み)

修正コミット`bf09909`を確認した。sdistに`proteindf_bridge.data/scripts/`の32本が含まれるようになり、Claudeがwheelとsdistのそれぞれを新しい一時環境にインストールして、どちらでもスクリプト32本のインストール、`python -P -m unittest discover -s tests`(197件)、`neutralize.py --help`の起動を確認した。文書の指摘2〜4(根拠のないRustのバージョン、内部の作業番号、`docs/usage.md`の古い説明)も直った。`cargo test --workspace`(318 passed)、clippy、fmtも確認した。

**参考**: `maturin sdist`の最初に出る`error: manifest path \`Cargo.toml\` does not exist`は、maturinがカレントディレクトリの`Cargo.toml`を最初に探すときのメッセージで、その後`pyproject.toml`の`manifest-path`のクレートを読んで正常に終了する(終了コード0)。害はない。

**注意**: マージ後は、ソースからの`pip install .`にRustツールチェーンが必要になる。ビルド済みwheelの配布はPR#50で用意する。

ユーザー承認のうえ、2026-10-04にdevelopへマージした(`6af3a4d`)。**PR#49は完了。** 次はPR#50。
