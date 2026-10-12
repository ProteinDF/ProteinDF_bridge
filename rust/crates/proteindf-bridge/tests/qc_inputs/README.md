# ProteinDF QC Test Inputs

本ディレクトリには、ProteinDF計算結果の検証テストデータ（`.h5`）を再生成するための入力ファイル（`fl_Userinput`）および生成スクリプトが格納されています。

## 出典

- リポジトリ: [ProteinDF_test](https://github.com/ProteinDF/ProteinDF_test.git)
- コミット: `8150fc0`
- 元のパス:
  - `serial_dev/CH3COOX_HF.mmap/fl_Userinput` -> `CH3COOX_HF/fl_Userinput`
  - `parallel/GLY.omp/fl_Userinput` -> `GLY/fl_Userinput`

## ライセンス

GPL-3.0-or-later (GPLv3)
ProteinDF_testリポジトリおよび本プロジェクト（`proteindf-bridge`）と同一ライセンスです。

## ProteinDF_test基準値 (pdfresults.std.h5)

ProteinDF本体（コミット `bf581a8` ビルド）にて各入力を実行し、ProteinDF_testの `pdfresults.std.h5` と比較した基準値：

| 体系 | 手法 | 原子数 | AO数 | MO数 | 最終全エネルギー (Hartree) | 収束反復回数 |
| --- | --- | --- | --- | --- | --- | --- |
| `CH3COOX_HF` | RKS (HF) | 8 (ダミー原子X含む) | 42 | 42 | `-227.32196883082707` | 88 |
| `GLY` | RKS (HF) | 10 | 95 | 95 | `-282.86416358371594` | 49 |

## 作り直し方 (データの再生成手順)

`regenerate.sh` スクリプトを実行して、指定した出力先ディレクトリに各体系の `.h5` ファイル（`CH3COOX_HF.h5`, `GLY.h5`）を生成します。
生成された `.h5` ファイルはリポジトリにはコミットせず、環境変数 `PROTEINDF_BRIDGE_QC_DATA` を介してテストから参照します。

### 前提環境
1. `PDF_HOME` に ProteinDF のインストールディレクトリが設定されていること。
   - `bin` が `PATH` に、`lib` が `LD_LIBRARY_PATH` に含まれていること（`PDF.x`, `pdf-setup.sh` が実行可能）。
2. `pdf-archive-h5.py` および `h5py` がインストールされた Python 環境の `bin` が `PATH` に含まれていること。

### 実行例

```bash
# 出力先ディレクトリを指定して実行
./rust/crates/proteindf-bridge/tests/qc_inputs/regenerate.sh /path/to/qc_data_dir

# テストの実行
PROTEINDF_BRIDGE_QC_DATA=/path/to/qc_data_dir cargo nextest run --workspace --features proteindf-bridge/hdf5 --run-ignored ignored-only
```
