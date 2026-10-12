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
