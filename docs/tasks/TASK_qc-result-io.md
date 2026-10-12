# TASK: 計算結果の読み込み (PR#58〜) — RUST_PORT_SPEC.md §3.4〜3.6

> ProteinDFの計算結果を読み込む機能をコアのクレート(`rust/crates/proteindf-bridge`)に追加する。対象は、Gaussian CUBEファイル(§3.6)、ProteinDFのHDF5(§3.4)、軌道ごとのMulliken寄与(§3.5)。

## 方針(ユーザー確認済み、2026-10-12)

1. **進める順**: CUBEの読み込み(PR#58) → HDF5の形式の調査(PR#59) → HDF5の読み込み(PR#60) → 軌道ごとのMulliken寄与(PR#61)。PR#60・61の詳細は、PR#59の調査結果を見てから本ファイルに追記する。
2. **HDF5はRustのオプション機能にする**: コアのクレートに機能フラグ`hdf5`を設け、既定では無効にする(wasm32のビルドにも影響させない)。依存は`hdf5-metno`。Pythonのwheelには入れない(Pythonからは当面ProteinDF_pytoolsの`h5py`版を使う)。必要になったらwheelへの同梱を別途検討する。
3. **テストデータ**: ProteinDF_test(GitHubで公開、GPLv3、本リポジトリと同じライセンス)の`serial_dev/*/reference.h5`のうち小さなものを2〜3個、`rust/crates/proteindf-bridge/tests/data/`にコピーしてよい。コピーする前に、ファイル内にローカルのパス・ユーザー名・ホスト名などの個人データがないことを確かめる。
4. **CUBE**: ProteinDFが出力した実データは手元にない。テストでは、仕様どおりの小さなCUBEをテストの中で作って使う。
5. 参照する既存実装: ProteinDF_pytools(`/home/toshi/work/dev/pdf-dev/ProteinDF_pytools`、`proteindf_tools/pdfparam_hdf5.py`・`matrix.py`・`vector.py`・`basisset.py`・`basis2.py`・`orbinfo.py`・`poputils.py`)。読むだけで変更しない。

## 役割分担・ブランチ運用(MUST)

- 実装はagy、レビューはClaudeが担当する。各PRは前のPRがレビューで承認され、developにマージされてから始める。
- `develop`から`feature/<名前>`を切って作業する。**`develop`へは自分でマージしない。** タグやReleaseは作らない。
- 完了報告は、起動時に指示されたパス(リポジトリの外)に、決まった見出しの短いファイルとして書く。コマンドの出力は要約・再構成せず、実際の出力行をそのまま貼る。
- テスト名・ファイル名・行番号は、必ず実際のファイルを開いて確認してから書く。確かめていない推測を事実として書かない。
- 既存の設計方針に従う: エラーは`BridgeError`(Pythonでは`BrError`の系統)、ファイルの書き出しを伴う場合は一時ファイル+名前の変更、Pythonバインディングは`RUST_PORT_SPEC.md` §4.3「共通の設計方針」に従う。

## 全PR共通のテスト

新しい一時環境にパッケージを`pip install`して`python -P -m unittest discover -s tests`、`cargo nextest run --workspace --status-level fail --final-status-level fail`、`cargo test --doc`、`cargo clippy --workspace --all-targets -- -D warnings`、`cargo fmt --all -- --check`。PR#60以降は、機能フラグ`hdf5`を有効にしたテスト・clippyも行う。

## PR#58: Gaussian CUBEファイルの読み込み

ブランチ: `feature/cube-reader`

### 対象

1. `src/format/cube.rs`を新設し、Gaussian CUBE形式の読み込みを実装する。
   - 2行のコメント行、原子数と原点、3軸それぞれの格子点数と刻みベクトル、原子の行(原子番号・電荷・座標)、格子点の値。
   - 原子数が負の場合は、原子の行の後に軌道番号の行がある(分子軌道のCUBE)。軌道の数と番号を読み、格子点の値は軌道ごとに持つ。
   - 格子点数の符号で単位が決まる(正ならbohr、負ならÅ)。読み込んだ値は単位の情報とともに保持し、勝手に換算しない。
   - 値の並びは、1つ目の軸が最も外側、3つ目の軸が最も内側(1行あたりの値の数は決め打ちしない)。
   - 形式の誤り(行が足りない、値の数が合わない、数値でないなど)は、行番号を含む`BridgeError`にする。
2. 原子を`AtomGroup`として取り出す関数(座標はÅに換算する。換算係数はコアに既にある定数があればそれを使い、なければ定数として定義して根拠を書く)と、格子点の座標・値を参照する関数を用意する。
3. Pythonバインディングにも公開する(外部の依存がないため)。クラス名・メソッド名は既存の形式(`SimpleGro`・`Xyz`など)にならう。
4. テスト(Rust・Python): 小さなCUBE(例: 2×3×4の格子、原子2個)をテスト内で作り、ヘッダー・原子・各格子点の値・単位を確かめる。分子軌道のCUBE(原子数が負)、Å単位のCUBE、形式の誤り(各種)も確かめる。

### 完了の定義

1. 対象1〜4に対応している。
2. 全PR共通のテストが通る。
3. `RUST_PORT_SPEC.md` §3.6に実施内容を追記する(完了日、主な関数、単位の扱い)。

## PR#59: HDF5の形式の調査(文書のみ)

ブランチ: `feature/qc-hdf5-survey`

### 対象

1. テストデータの候補(`serial_dev/H2O_HF`・`serial_dev/N2_UHF`など)の`reference.h5`の中身(グループ・データセット・属性の木構造、型、形)を調べ、`docs/qc-result-format.md`に書く。調べるための道具(一時環境の`h5py`、または`hdf5-metno`を使った使い捨てのプログラム)は、リポジトリに残さない。
2. ProteinDF_pytoolsの`pdfparam_hdf5.py`が読み書きする項目と、1の実際の中身を対応づける。基底関数系、分子(原子と座標)、軌道エネルギー、占有数、MO係数行列、重なり行列(S)、Mulliken電荷などが、HDF5の中にあるのか、別ファイル(ProteinDFの`.mat`・`.vtr`のバイナリ)にあるのかを明らかにする。別ファイルの場合は、その形式(`matrix.py`・`vector.py`の読み込み処理)も書く。
3. §3.5の軌道ごとのMulliken寄与の計算に必要なデータ(MO係数・S行列・基底関数と原子の対応)が、テストデータからそろうかを判断する。そろわない場合は、何が足りないかを書く。
4. 1〜3を踏まえて、PR#60・61の範囲の案を書く(Claudeとユーザーが判断する)。
5. テストデータとして使うファイルを`rust/crates/proteindf-bridge/tests/data/`にコピーする(方針3。個人データの有無の確認方法と結果を報告する)。

### 完了の定義

1. 対象1〜5に対応している。コードの変更はない(テストデータの追加と文書のみ)。

## PR#58 レビュー結果(1回目、2026-10-12、要修正)

`feature/cube-reader`(`19790cb`・`63bdb8b`・`1a2aad3`)をレビューした。CUBEの読み込みは仕様どおり(単位を格子点数の符号で判定し換算しない、原子数が負なら軌道番号の行を読み、格子点の値は軌道が最も内側に並ぶ形で軌道ごとに分ける、1行あたりの値の数を決め打ちしない、形式の誤りは行番号付きの`BridgeError`)。テストも分子軌道の並び(`1.0 2.0 1.1 2.1 ...`)を区別して確かめている。Pythonテスト一式242件、`cargo nextest`(`325 tests run: 325 passed, 2 skipped`)、clippy・fmtをClaudeが新しい環境で確認した。

### 修正依頼

1. **【方針】純Python版のパッケージを変更している。** `proteindf_bridge/cube.py`を新設し、`proteindf_bridge/__init__.py`に`from .cube import Cube`を足している。これにより、純Python版の`import proteindf_bridge`が初めてRustの拡張モジュールを読み込むことになる(これまで純Python版のモジュールは`proteindf_bridge.rs`をimportしていない)。Rust版のクラスは`proteindf_bridge.rs`から使う方針(`RUST_PORT_SPEC.md` §4.4)なので、`proteindf_bridge/cube.py`を削除し、`__init__.py`の変更を元に戻す。Pythonのテストも`proteindf_bridge.rs.Cube`だけを使うようにする。
2. **【API】同じ機能に2つの名前を付けている。** 新しいクラスなので互換のための別名は不要。既存のクラス(`SimpleGro`の`num_of_atoms`・`get_atomgroup`)にそろえ、Pythonの`Cube`からは`num_atoms`と`get_atom_group`を削除する。コアの`Cube`も、`get_atomgroup`と`get_atom_group`のように同じ機能の別名があれば1つにする。

### 完了の定義(修正後)

1. 上記1・2に対応し、同じブランチに追加コミットする。
2. 全PR共通のテストが通る。`git diff develop -- proteindf_bridge/`の出力が空であることを示す。

## PR#58 完了(2026-10-12)

修正コミット(`7e7de3a`・`fbbc0d9`)で修正依頼1・2に対応済み。`git diff develop -- proteindf_bridge/`が空であること、Pythonのテストが`proteindf_bridge.rs.Cube`だけを使うこと、別名(`num_atoms`・`get_atom_group`)が削除されたことをClaudeが確認した。Pythonテスト一式242件OK(新しい環境)、`cargo nextest` `325 tests run: 325 passed, 2 skipped`、clippy・fmt OK。`--no-ff`でdevelopにマージした(`64d88f0`)。
