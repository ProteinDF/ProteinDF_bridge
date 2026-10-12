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

### 実施内容 (2026-10-12)

1. **HDF5形式の調査**: `ProteinDF_test/serial_dev` の `reference.h5` を調査し、ルート属性・グループ・データセットの木構造、型、形を `docs/qc-result-format.md` にまとめた。
2. **ProteinDF_pytoolsとの対応づけ**: `pdfparam_hdf5.py`, `matrix.py`, `vector.py`, `orbinfo.py`, `basisset.py` の実装と対応づけ、各データの所在（HDF5内 vs 外部バイナリファイル）およびバイナリ行列・ベクトルのフォーマット仕様を記載した。
3. **Mulliken寄与計算の充足性判断**: MO係数行列 $C$ とAO-原子対応はそろうが、重なり行列 $S$ は `reference.h5` に含まれていない（不足）ことを突き止め、その根拠と理由を明記した。
4. **PR#60・61の範囲の案**: HDF5読み込み機能（PR#60）および軌道別Mulliken寄与計算・S行列入力（PR#61）の設計案を提示した。
5. **テストデータのコピー**: `H2O_HF_reference.h5`, `Ne_HF_reference.h5`, `N2_UHF_reference.h5` を `rust/crates/proteindf-bridge/tests/data/` に追加。コピー前に個人データがないことを確認済み。

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

## PR#59 レビュー結果(1回目、2026-10-12、要修正・文書のみ)

`feature/qc-hdf5-survey`(`daec449`・`9761241`・`c2c6498`)をレビューした。`H2O_HF_reference.h5`の木構造(ルート属性、`C/rks_49`(24×24)、`energy_level/rks_N`、`occ/rks`、`pop/mulliken/rks_49`、`molecule/xyzc`(固定長文字列、座標はbohr)、`basisset/{name,cgto,pgto}`)をClaudeが`h5py`で確認し、`docs/qc-result-format.md`の記述と一致した。テストデータ3ファイルに個人データがないことも確認した。

### 修正依頼

1. **【結論の誤り】「重なり行列SはHDF5内には存在しない」は、テストデータ(`reference.h5`)についてのみ正しい。** ProteinDF_pytoolsの`scripts/pdf-archive-h5.py`(ProteinDFの計算結果をHDF5にまとめる標準のスクリプト)は、`archive_s`(129〜133行目付近、`logger.debug("save S")`→`entry.save_s_matrix(S)`)で常にS行列をデータセット`s`として書き出す。Claudeが、このスクリプトで作られた実際のファイル(Ne原子のB3LYP計算の`pdfresults.h5`)を確認したところ、`s`(14×14)・`h`・`h2`・`C/rks_1`・`P/rks_1`・`F/rks_1`があった(`pop/mulliken`はなかった)。`docs/qc-result-format.md`と報告を、次のように事実どおりに直す。
   - `reference.h5`(ProteinDF_testの比較用の参照ファイル)と、`pdf-archive-h5.py`が作る`pdfresults.h5`は、同じ形式の木構造だが、含まれるデータセットが違う。後者には`s`・`h`・`h2`・`P`・`F`があり、前者にはない。各データセットがどちらに含まれるかを表にする(確認していないものは「未確認」)。`pdf-archive-h5.py`の該当行を実物で確かめて書く。
   - 軌道ごとのMulliken寄与に必要なデータは、`pdf-archive-h5.py`が作るHDF5ならそろう(S・C・基底関数系・分子)。
2. **【PR#61の案の修正】** テストは合成のS行列ではなく、実際の計算結果で行う。テストデータ(`s`を含む多原子分子のHDF5)はユーザーが用意する(ユーザー確認済み、2026-10-12)。PR#61の案のうち、外部ファイル`Spq.mat`の読み込みは、HDF5の`s`で足りるため範囲から外す(必要になったら別途)。

### 完了の定義(修正後)

1. 上記1・2に対応し、同じブランチに追加コミットする。コードの変更はない。

## PR#59 完了(2026-10-12)

修正コミット`da4e2df`で修正依頼1・2に対応済み(文書のみ)。`reference.h5`と`pdf-archive-h5.py`が作る`pdfresults.h5`の違いの表と、`pdf-archive-h5.py`の引用行(`archive_s` 128〜133行目ほか)をClaudeが実物で確認した。`pdf-archive-h5.py`は`pop/mulliken`を書き出さない。`--no-ff`でdevelopにマージした(`9db1415`)。

## PR#60: ProteinDFのHDF5の読み込み(計画、ユーザー確認済み 2026-10-12)

ブランチ: `feature/qc-hdf5-reader`。形式は`docs/qc-result-format.md`による。

### 対象

1. **機能フラグ**: コアのクレートに機能`hdf5`を設ける(依存は`hdf5-metno`)。既定では無効。wasm32とPythonのwheelには入れない。`hdf5-metno`のHDF5を、ソースからの静的ビルド(`static`機能)にするか、システムのライブラリを使うかは、ローカル(システムにlibhdf5はない。cmakeとCコンパイラはある)とCIの両方でビルドできる方を選び、理由を報告する。
2. **読み込む内容**(`src/qc_result/`、型の名前は`docs/qc-result-format.md` §4の案を基にしてよい):
   - ルート属性(手法`method`、交換相関汎関数、原子数・AO数・MO数、反復回数、収束したか)。
   - 分子(`molecule/xyzc`。座標はbohrのまま保持し、`AtomGroup`への変換ではÅに換算する。換算係数は`format/cube.rs`の`BOHR_TO_ANGSTROM`を使う)。
   - 基底関数系(`basisset/{name,cgto,pgto}`)と、AOと原子の対応(ProteinDF_pytoolsの`orbinfo.py`と同じ規則。AOの並び順・角運動量の扱いを実物で確かめてドキュメントコメントに書く)。
   - 軌道エネルギー・占有数・MO係数。反復ごとのデータは既定で最終反復(`iterations`)を使い、反復を指定して読むこともできる。RKS(`rks`)とUKS(`uks-alpha`・`uks-beta`など、実際の名前を確かめる)の両方に対応する。
   - ファイルにあれば読むもの: `s`・`P`・`F`・`h`・`pop/mulliken`。ないものは`Option`で表し、エラーにしない。
   - 数値が固定長文字列として格納されている部分(`xyzc`・`cgto`・`pgto`)は数値に変換する。変換できない値は`BridgeError`にする。
3. **テスト**: `tests/data/`の3ファイル(`H2O_HF_reference.h5`・`Ne_HF_reference.h5`・`N2_UHF_reference.h5`)で、上の各項目を確かめる。期待値は、リポジトリの外の一時環境で`h5py`やProteinDF_pytoolsを使って読んだ値と比べて決め、どう決めたかをテストのコメントに書く。`s`などがないファイルで`None`になることも確かめる。
4. **CI**: `.github/workflows/ci.yml`に、機能`hdf5`を有効にしたテストとclippyを加える。
5. `RUST_PORT_SPEC.md` §3.4に実施内容を追記する。

### 完了の定義

1. 対象1〜5に対応している。
2. 全PR共通のテストに加え、`cargo nextest run --workspace --features proteindf-bridge/hdf5 --status-level fail --final-status-level fail`と、同じ機能でのclippyが通る。機能`hdf5`なしのビルド(既定)と、`cargo check --target wasm32-unknown-unknown`(既存の手順どおり)が、これまでどおり通る。
3. CIの変更は、ローカルで同じコマンドが通ることを示す(CIでの結果はマージ後にClaudeが確認する)。

## PR#61: 軌道ごとのMulliken寄与(予定)

ユーザーが用意する、`s`を含む多原子分子の`pdfresults.h5`を受け取ってから始める。答え合わせは、(1)各軌道の寄与の原子についての和が1、(2)占有軌道の寄与を占有数で重み付けした和が、ファイル内の`P`と`s`から計算した原子ごとのMulliken gross populationと一致、(3)`pop/mulliken`があればそれとも一致、で行う。

## PR#60 レビュー結果(1回目、2026-10-12、要修正)

`feature/qc-hdf5-reader`(`0aaa323`〜`1b8dc14`の5コミット)をレビューした。機能`hdf5`は既定で無効、`hdf5-metno`の`static`機能でHDF5をソースからビルドする選択(ローカルにlibhdf5がなく、CIでも追加パッケージが不要)は妥当。AOの展開順(`s`、`px py pz`、`dxy dyz dzx dxx-yy dzz`)はProteinDF_pytoolsの`basisset.py` 156行目と、実行種別の名前(`rks`、`uks_alpha`/`uks_beta`、`roks_close`/`roks_open`)は`pdfparam_object.py` 62〜67行目と一致することをClaudeが確認した。`cargo nextest`(機能`hdf5`あり: `328 tests run: 328 passed, 2 skipped`)、機能ありのclippy、wasm32のcheck、Pythonテスト一式242件をClaudeが確認した。

### 修正依頼

1. **【重複する型】`src/qc_result/types.rs`に独自の`SymmetricMatrix`(全要素を持つ行優先の`Vec<f64>`)を定義しているが、コアには既に`crate::matrix::SymmetricMatrix`(`src/matrix.rs` 516行目)と`crate::matrix::Matrix`がある。** 同じ名前の別の型が公開APIに2つあると、利用者が取り違える。S・P・F・h・h2は`crate::matrix::SymmetricMatrix`、MO係数は`crate::matrix::Matrix`(または`Matrix`を内部に持つ型)で表す。HDF5のデータセットが全要素(N×N)で格納されていることは、読み込み時に対称性を確かめる(許容誤差を決めて、外れたら`BridgeError`)か、下三角だけを使うかを決めて、ドキュメントコメントに書く。
2. **【テスト不足】UKS・ROKSの軌道データ(MO係数・軌道エネルギー・占有数・Mulliken電荷)と、`s`・`P`・`F`・`h`・`h2`の読み込みが、一度もテストで通っていない。** テストデータ3ファイルのうち、`N2_UHF_reference.h5`には軌道データがなく、どのファイルにも`s`などがない。テストの中で、`hdf5-metno`を使って一時ディレクトリに小さなHDF5を書き、読み込みを確かめる。書くデータの木構造(データセットの名前、`C/uks_alpha_{反復}`などのパス、`s`などの属性`type`・`row`・`col`)は、ProteinDF_pytoolsの`pdfparam_hdf5.py`の`_get_h5path_*`(643〜690行目付近)と`_set_matrix`の書き方に合わせ、テストのコメントに該当行を書く。UKS(α・βで値が違うこと)、ROKS、`s`ありで対称でない値を入れた場合(修正依頼1で決めた扱い)を確かめる。実際の計算結果での確認は、PR#61でユーザーが用意するファイルで行う。

### 完了の定義(修正後)

1. 上記1・2に対応し、同じブランチに追加コミットする。
2. 完了の定義2(PR#60)のテストが通る。

## PR#60 完了(2026-10-12)

修正コミット(`4fba5c3`・`2907e4b`)で修正依頼1・2に対応済み。独自の`SymmetricMatrix`が削除され、S・P・F・h・h2はコアの`matrix::SymmetricMatrix`、MO係数はコアの`Matrix`を持つ型で表すこと、全要素で格納された行列を読み込み時に対称性検査(許容誤差`1e-8`)すること、合成HDF5によるUKS・ROKS・`s`等・非対称`s`のテスト(`test_synthetic_uks_hdf5`・`test_synthetic_roks_hdf5`・`test_synthetic_asymmetric_s_matrix_error`)をClaudeが確認した。`cargo nextest`は機能`hdf5`ありで`331 tests run: 331 passed, 2 skipped`、なしで`325 tests run: 325 passed, 2 skipped`。clippy(機能あり)・fmt・wasm32のcheck OK。`hdf5-metno`の依存`proc-macro-error2`について将来のRustでの非互換の警告が出る(現状のビルドには影響なし)。`--no-ff`でdevelopにマージした(`4d0333f`)。

## 計画の変更(ユーザー確認済み、2026-10-12)

PR#61の検証データは、ユーザーが用意する結果ファイルではなく、ProteinDF_testの入力ファイル(`fl_Userinput`)をコピーし、ProteinDFを実行して作り直せるようにする。結果ファイル(`.h5`)はコミットしない。Claudeが手元で試した結果は次のとおり。

- ProteinDF(`ProteinDF` `bf581a8`のビルド)で、ProteinDF_test(`8150fc0`)の`serial_dev/CH3COOX_HF.mmap`と`parallel/GLY.omp`の入力を実行すると、それぞれ約18秒・約46秒で収束し、最終の全エネルギーはProteinDF_testの`pdfresults.std.h5`と1e-11の差で一致した(CH3COOX_HF `-227.32196883082707`、GLY `-282.86416358371594`、反復回数88・49)。
- `pdf-archive-h5.py`(ProteinDF_pytools)の出力には`s`・`h`・`h2`・`P`・`F`がある。**ただし、`type`がSYの行列はすべて、下三角だけに値が入り、上三角は0で格納されている**(ProteinDF_pytoolsの`pdfparam_hdf5.py`の`_set_matrix`が`matrix.data`をそのまま書くため)。下三角を折り返した`s`は、(CCᵀ)⁻¹と1e-14で一致し、これと`C`・占有数から計算したMulliken電荷は、ProteinDF_testの`pdfresults.std.h5`の`pop/mulliken`と小数点以下8桁で一致した(ダミー原子`X`は基底関数を持たないので、電荷は`xyzc`の電荷の列の値1.0)。
- PR#60の読み込みは、全要素の対称性を確かめて外れたらエラーにするため、`pdf-archive-h5.py`の実際の出力のSY行列をすべて拒否する(PR#60のテストデータ`reference.h5`にはSYの行列がなかったため、気づかなかった)。
- 力の計算(`force`)を含む入力(`serial_dev/H2O_HF`)では、`pdf-archive-h5.py`が`pdfparam_object.py`の`set_by_raw_data`(1355行目)で`KeyError: 0`になる。ProteinDF_pytools側の問題で、本タスクの対象外。

これに伴い、PR#61を「SY行列の読み込みの修正と検証データの作り直しの仕組み」、PR#62を「軌道ごとのMulliken寄与」とする。

## PR#61: SY行列の読み込みの修正と、検証データを作り直す仕組み(計画、ユーザー確認済み 2026-10-12)

ブランチ: `feature/qc-regen-testdata`。

### 対象

1. **SY行列の読み込みの修正**(`src/qc_result/pdfparam_hdf5.rs`の`read_symmetric_matrix`と、グループ内の行列を読む部分): `type`がSYのデータセットは、下三角(行≥列)の値から対称行列を作る。上三角(行<列)は、(a)すべて0なら下三角だけの格納として受け入れ、(b)0でない値があれば、対応する下三角の値と許容誤差`MATRIX_SYMMETRY_TOLERANCE`以内で一致することを確かめ、外れたら`BridgeError`にする。この規則と、`pdf-archive-h5.py`が下三角だけを書くこと(根拠の行番号は実際のファイルで確かめる)をドキュメントコメントと`docs/qc-result-format.md`に書く。GEの行列(`C`)の読み方は変えない。
2. **合成HDF5のテスト**(`tests/test_qc_hdf5.rs`): 下三角だけで格納した`s`・`P`・`F`を読めて、折り返した値になることを確かめるテストを追加する。既存の`test_synthetic_asymmetric_s_matrix_error`は、上三角に下三角と合わない0でない値がある場合のエラーとして残す(必要なら名前とコメントを直す)。
3. **入力ファイルのコピー**: ProteinDF_test(`8150fc0`)の`serial_dev/CH3COOX_HF.mmap/fl_Userinput`と`parallel/GLY.omp/fl_Userinput`を、`rust/crates/proteindf-bridge/tests/qc_inputs/CH3COOX_HF/fl_Userinput`・`.../GLY/fl_Userinput`に、内容を変えずにコピーする。同じディレクトリの`README.md`に、出典(リポジトリのURL・コミット・元のパス)、ライセンス(GPLv3)、作り直し方、ProteinDF_testの`pdfresults.std.h5`の最終の全エネルギーと反復回数(上記)を書く。
4. **作り直すスクリプト**: `rust/crates/proteindf-bridge/tests/qc_inputs/regenerate.sh`(bash)。出力先ディレクトリを引数で受け取り、入力ごとに一時ディレクトリへ`fl_Userinput`をコピーして`pdf-setup.sh`・`PDF.x`・`pdf-archive-h5.py`を実行し、`<出力先>/<名前>.h5`に置く。`PDF_HOME`が未設定のとき、または必要なコマンド(`PDF.x`・`pdf-setup.sh`・`pdf-archive-h5.py`)が見つからないときは、理由を出して0以外で終わる。SCFが収束しなかったとき(`pdf.log`に`SCF is well converged`がない)もエラーにする。ローカルのパスを書き込まない。リポジトリ内に出力する場合の置き場として`rust/crates/proteindf-bridge/tests/qc_generated/`を`.gitignore`に加える。
5. **実データのテスト**(`tests/test_qc_hdf5.rs`、機能`hdf5`): `#[ignore]`を付け、環境変数`PROTEINDF_BRIDGE_QC_DATA`(作り直した`.h5`のディレクトリ)から読む。`--run-ignored`で実行したときに環境変数がない・ファイルがない場合は、黙って飛ばさずに失敗させる。2つのファイルについて、次を確かめる。
   - 読み込みが成功し、`s`・`P`・`F`・`h`があり、原子数・AO数・MO数・反復回数が期待どおり(CH3COOX_HF: 原子8(ダミー原子`X`を含む)・AO42・MO42・反復88、GLY: 原子10・AO95・MO95・反復49。値は作り直したファイルで確かめる)。
   - 最終の全エネルギーが上記の値と1e-6以内で一致する。
   - CᵀSCが単位行列と1e-8以内で一致する。
   - 最終反復の`P`が、`C`と占有数から作った密度行列と1e-8以内で一致し、tr(PS)が電子数(CH3COOX_HF 32、GLY 40)と1e-8以内で一致する。
   - 期待値をどう決めたか(どのツールで、どのファイルから)をテストのコメントに書く。
6. **CI**: 実データのテストは`#[ignore]`のため実行されない。合成HDF5のテストは、既存の機能`hdf5`のステップで実行される。CIの設定の変更は不要のはず(変更した場合は理由を報告する)。

### 完了の定義

1. 対象1〜6に対応している。
2. 全PR共通のテストと、機能`hdf5`ありの`cargo nextest`・clippy、wasm32のcheckが通る。
3. `regenerate.sh`でリポジトリの外に2つの`.h5`を作り、`PROTEINDF_BRIDGE_QC_DATA`を設定して`cargo nextest run --workspace --features proteindf-bridge/hdf5 --run-ignored ignored-only --status-level fail --final-status-level fail`(または同等)が通ることを、実際の出力で示す。環境変数なしで同じコマンドを実行すると失敗することも示す。

## PR#62: 軌道ごとのMulliken寄与(予定)

PR#61のマージ後に詳細を書く。PR#61で作り直した`.h5`の`s`を使い、(1)各軌道の寄与の原子についての和が1、(2)占有数で重み付けした和を`xyzc`の電荷の列から引いた値が、ProteinDF_testの`pdfresults.std.h5`の`pop/mulliken`と一致、を確かめる。

## PR#61 レビュー結果(1回目、2026-10-12、要修正)

`feature/qc-regen-testdata`(`79c32a1`・`7bac411`・`36e3de0`)をレビューした。SY行列の読み方(下三角から作り、上三角がすべて0なら下三角だけの格納として受け入れ、0でない値があれば許容誤差で照合)は計画どおり。根拠の行番号(`pdfparam_hdf5.py` 605行目の`data=matrix.data`、純Python版`matrix.py`の`SymmetricMatrix.set`)をClaudeが確認した。Claudeが`regenerate.sh`でリポジトリの外に2つの`.h5`を作り直し、実データのテスト(`4 tests run: 4 passed, 332 skipped`)、機能`hdf5`ありの`332 tests run: 332 passed, 4 skipped`、なしの`325 tests run: 325 passed, 2 skipped`、clippy(機能あり)・fmt・wasm32のcheckが通ることを確認した。追加したファイルにローカルのパスは含まれていない。

### 修正依頼

1. **【どの反復の行列か分からない】`QcResult`の`density_matrices`・`fock_matrices`は`Vec<SymmetricMatrix>`で、実行種別と反復の情報を持たない。** 実際のファイルには`P/rks_1`〜`P/rks_88`のように全反復の行列があり、実データのテストは並び順(`read_group_matrices`の並べ替え)に頼って`.last()`を最終反復として使っている。MO係数(`MoCoefficients`、`types.rs` 446行目)と同じく、実行種別・反復・行列を持つ型にし、`get_mo_coefficients`と同じ形の取り出し関数(反復を省略すると最終反復)を用意する。実データのテストは、その関数で最終反復の`P`を取り出す形に直す。合成HDF5のテストで、UKSのα・βと反復の取り違えがないことを確かめる。
2. **【読めないデータを黙って捨てる】`read_group_matrices`は、名前を解釈できないメンバー、正方でないデータセット、読み込みに失敗したデータセットを、エラーにせず飛ばしている。** 読み込みの失敗と正方でない行列は`BridgeError`にする。名前を解釈できないメンバーを飛ばすなら、その理由(ProteinDF_pytoolsが書く名前の規則)をドキュメントコメントに書く。
3. **【基準値の出典の取り違え】`tests/qc_inputs/README.md`の表「ProteinDF_test基準値 (pdfresults.std.h5)」と、実データのテストのコメント「Source: ... pdfresults.std.h5」「Expected final energy: -227.32196883082707」は、作り直したファイルの値(`bf581a8`のビルド)を、ProteinDF_testの`pdfresults.std.h5`の値として書いている。** `pdfresults.std.h5`の最終の全エネルギーは、CH3COOX_HFが`-227.32196883081747`、GLYが`-282.86416358371673`(Claudeが`h5py`で確認)。作り直した値は、CH3COOX_HFが`-227.32196883082707`、GLYが`-282.86416358371594`。両方を区別して書き、テストの期待値は`pdfresults.std.h5`の値にする(許容誤差1e-6はそのまま)。この取り違えは、本ファイルの「計画の変更」の書き方(Claude)があいまいだったことにもよる。
4. **【ダミー原子の扱いのテスト】`OrbitalInfo::build`で、基底関数の定義がない原子(ダミー原子`X`)をAO数0として飛ばすように変えたが、これを確かめるのはCIで動かない実データのテストだけ。** 合成HDF5で、基底関数のない原子を途中に含む分子を作り、AOと原子の対応がずれないこと(`X`の後ろの原子のAOが正しい原子に対応すること)を確かめるテストを追加する。
5. **【期待値の決め方の記述】実データのテストのコメント「Evaluated using h5py / numpy in isolated scratchpad environment.」は、何をどう確かめたのかが分からない。** どのファイルのどのデータセットを読んで、どの値を確かめたかを書く。
6. **【後片付け】worktreeに`.venv/`・`__pycache__/`が残っている。** 消す。

### 完了の定義(修正後)

1. 上記1〜6に対応し、同じブランチに追加コミットする。
2. 完了の定義2・3(PR#61)を改めて示す。

## PR#61 完了(2026-10-12)

修正コミット(`7bc65a7`)で修正依頼1〜6に対応済み。`P`・`F`が実行種別と反復を持つ`DensityMatrix`・`FockMatrix`になり、`get_density_matrix`・`get_fock_matrix`で取り出せること、`read_group_matrices`が読み込みの失敗と正方でない行列をエラーにすること、基準値をProteinDF_testの値と作り直した値に分けて書き、テストの期待値をProteinDF_testの値にしたこと、ダミー原子の合成テスト(`test_synthetic_dummy_atom_orbital_mapping`)をClaudeが確認した。Claudeが作り直した`.h5`で実データのテスト`4 tests run: 4 passed, 333 skipped`、機能`hdf5`ありで`333 tests run: 333 passed, 4 skipped`、なしで`325 tests run: 325 passed, 2 skipped`、Pythonテスト一式242件OK(新しい環境)。`--no-ff`でdevelopにマージした(`61c9903`)。マージ後のdevelopで`scripts/ci-local.sh all qc-data`が全ジョブ成功した。
