# TASK: CIをローカルで回す (PR#63)

> ユーザーの指示(2026-10-12): 「CIは基本的にローカルで回してください」。GitHub Actionsの`ci.yml`は手動実行(`workflow_dispatch`)だけを残す(ユーザー確認済み)。リリース用の`release.yml`は変えない。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaudeが担当する。
- `develop`から`feature/local-ci`を切って作業する。**`develop`へは自分でマージしない。** タグやReleaseは作らない。
- 完了報告は、起動時に指示されたパス(リポジトリの外)に、決まった見出しの短いファイルとして書く。コマンドの出力は要約・再構成せず、実際の出力行をそのまま貼る。
- テスト名・ファイル名・行番号は、必ず実際のファイルを開いて確認してから書く。確かめていない推測を事実として書かない。

## 対象

1. **ローカルCIのスクリプト** `scripts/ci-local.sh`(bash、`set -euo pipefail`)。リポジトリのルートで実行する(どのディレクトリから呼ばれてもルートに移動する)。引数でジョブを選べ、引数なしなら既定のジョブをすべて実行する。ジョブ:
   - `rust`: `cargo fmt --all -- --check`、`cargo clippy --workspace --all-targets -- -D warnings`、`cargo nextest run --workspace --status-level fail --final-status-level fail`、`cargo test --doc --workspace`(いずれも`--manifest-path rust/Cargo.toml`)。
   - `hdf5`: 機能`proteindf-bridge/hdf5`を有効にしたclippyと`cargo nextest`。
   - `wasm`: `cargo check --target wasm32-unknown-unknown -p proteindf-bridge --no-default-features --features ruzstd`(既存の手順と同じ。ターゲットが入っていなければ、入れ方を示して失敗する)。
   - `wheel`: リポジトリの外の一時ディレクトリに新しい仮想環境を作り、`maturin build --release`でwheelを作ってインストールし、`python -P -m unittest discover -s tests`と、`ci.yml`と同じスクリプトの`--help`を実行する。
   - `sdist`: 同様に一時ディレクトリでsdistを作ってインストールし、同じテストを実行する。
   - `docs`: 一時ディレクトリの仮想環境に`.[docs]`を入れ、英語・日本語のHTMLを一時ディレクトリにビルドする(警告の扱いは`ci.yml`と同じ)。
   - `qc-data`(既定では実行しない): 環境変数`PROTEINDF_BRIDGE_QC_DATA`が設定されているときに、機能`hdf5`で`#[ignore]`のテストを実行する。このテストはPR#61で入る予定なので、本PRの時点で該当するテストがなくても失敗しないようにし、その旨をスクリプトのコメントに書く。
   - Pythonは環境変数`PYTHON`(既定`python3`)で選べる。`maturin`・`cargo-nextest`など必要なコマンドがないときは、何が足りないかを出して失敗する。仮想環境のツールは一時環境に`pip install`してよい。
   - 一時ディレクトリは終了時(失敗時も)に消す。リポジトリの中には、cargoの`target/`以外のファイルを作らない(`dist/`なども作らない)。既存の`.venv`などには触れない。
   - 最後に、ジョブごとの成否と所要時間の一覧を出す。どれかが失敗したら0以外で終わる(失敗したジョブで止めるか、残りを続けて最後にまとめて失敗にするかは、どちらかに決めて`--help`に書く)。
2. **`ci.yml`**: `on:`を`workflow_dispatch`だけにし、ジョブの中身は変えない。ファイルの先頭に、通常の確認は`scripts/ci-local.sh`で行うことをコメントで書く。
3. **文書**: `CONTRIBUTING.md`の「Running Tests」を、`scripts/ci-local.sh`を使う手順に書き換える(ジョブの一覧、前提となるツール、`qc-data`の使い方)。`ci.yml`が手動実行だけになったことも書く。

## 守ること

- `release.yml`・`docs.yml`は変更しない。
- 純Python版のパッケージ・Rustのコードは変更しない。
- `RUST_PORT_SPEC.md`・`docs/tasks/`の過去の記録は書き換えない(本ファイルへの追記もしない)。
- 意味のまとまりごとにコミットする。

## 完了の定義

1. 対象1〜3に対応している。
2. `scripts/ci-local.sh`(引数なし)がすべて成功し、最後の一覧を実際の出力で示す。`scripts/ci-local.sh rust`のようにジョブを選んで実行できることも示す。
3. 実行の前後で`git status --porcelain`が変わらない(リポジトリの中にファイルが増えない)ことを、実際の出力で示す。
4. `shellcheck`が使えれば、`shellcheck scripts/ci-local.sh`の結果を示す(使えなければその旨を書く)。

## PR#63 レビュー結果(1回目、2026-10-12、要修正)

`feature/local-ci`(`42651b6`・`6225afa`・`f5578aa`)をレビューした。ジョブの構成、`--fail-fast`、前提ツールの確認、`ci.yml`を`workflow_dispatch`だけにした変更、`CONTRIBUTING.md`の書き換えは指示どおり。

### 修正依頼

1. **【失敗を見逃す】`scripts/ci-local.sh`の313〜324行目は、`set +e`にしてからジョブの関数を呼んでいるため、関数の中でも`errexit`が無効になり、最後のコマンド以外の失敗が無視される。** 例えば`rust`ジョブで`cargo fmt --check`や`clippy`が失敗しても、最後の`cargo test --doc`が成功すれば`SUCCESS`になる(Claudeが`bash -c 'set -euo pipefail; f(){ false; echo continued; }; set +e; f; echo rc=$?'`で`rc=0`になることを確認した)。各ジョブをサブシェルで`set -e`を有効にして実行する(`( set -e; run_job_x )`の終了コードを取る)などで直す。**ジョブの途中のコマンドが失敗したときに`FAILED`になることを確かめるテスト**を用意する(例: 環境変数で途中のコマンドをわざと失敗させる仕組みは公開しない。代わりに、`cargo fmt --check`が失敗するようにリポジトリの外の一時コピーで書式を崩して実行する、など。方法は任せるが、リポジトリに失敗用の仕組みを足さない)。
2. **【リポジトリの中にファイルを作っている】`git status --porcelain`は変わらないが、`.gitignore`で無視されたファイルがリポジトリの中にできている。** `docs`ジョブの`pip install -e .[docs]`(編集可能インストール)で`proteindf_bridge/rs.abi3.so`が、Sphinxの日本語ビルドで`docs/locale/ja/LC_MESSAGES/**/*.mo`ができる(Claudeが`git status --porcelain --ignored`で確認した)。TASKの「リポジトリの中には、cargoの`target/`以外のファイルを作らない」に反する。`docs`ジョブは、編集可能インストールをやめ、作業ツリーの必要なファイルをリポジトリの外の一時ディレクトリにコピーして、そこでインストールとビルドを行うなどで直す。完了の定義3の確認は`git status --porcelain --ignored`の実行前後の比較で行う(`rust/target/`は除いてよい)。
3. **【後片付け】worktreeに、作業中に作られた`.venv/`・`proteindf_bridge/rs.abi3.so`・`__pycache__/`・`docs/locale/ja/LC_MESSAGES/**/*.mo`が残っている。** 消してから完了の定義3の確認を行う。
4. **【コメントの誤り】`run_job_qc_data`のコメント(294〜296行目)は「対象のテストが0件のこともある」としているが、既に`#[ignore]`のテストが2件あり、`--run-ignored only`で実行される(Claudeの確認: `2 tests run: 2 passed, 325 skipped`)。実際の状態に合わせて直す。`--help`の`PROTEINDF_BRIDGE_QC_DATA`の説明「Path or flag」も、作り直した`.h5`のディレクトリのパスであることに直す。
5. **【未確認の記述】報告の「気になった点」の、sdistで出る`error: manifest path 'Cargo.toml' does not exist`について、どのコマンドのどの段階で出ているのかを実際の出力行で示し、原因を調べて書く。** 害がないと判断する場合は、その根拠(どの段階で、なぜ最終的に正しく動くのか)を示す。
6. 一時ディレクトリは`mktemp -d -t proteindf-bridge-ci-XXXXXX`のように`TMPDIR`に従う形にする。

### 完了の定義(修正後)

1. 上記1〜6に対応し、同じブランチに追加コミットする。
2. 完了の定義2〜3(PR#63)を、修正後のスクリプトで改めて示す。修正依頼1のテストで`FAILED`になる実際の出力を示す。

## PR#63 完了(2026-10-12)

修正コミット(`c18be9f`)で修正依頼1〜6に対応済み。各ジョブを`set -e`を有効にしたサブシェルで実行し、途中のコマンドの失敗で`FAILED`・終了コード1になること(リポジトリの外のコピーで書式を崩して確認)、全ジョブの実行の前後で`git status --porcelain --ignored`にファイルが増えないことをClaudeが確認した。worktreeにできる`.venv/`は、ユーザーのmiseの設定(`python.uv_venv_auto`)によるもので、スクリプトとは関係ない。sdistの最初に出る``error: manifest path `Cargo.toml` does not exist``は、その直後に`pyproject.toml`の`manifest-path`でビルドが進み、インストールとテストが通るため害はない。`--no-ff`でdevelopにマージした(`6531714`)。マージ後のdevelopで`scripts/ci-local.sh all qc-data`が全ジョブ成功した。
