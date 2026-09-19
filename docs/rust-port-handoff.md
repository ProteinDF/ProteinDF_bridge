# proteindf-bridge Rust移植 — アーキテクチャ決定・フェーズ実績・教訓記録

> **ドキュメントの位置づけ**: 本ドキュメントは、Phase 1〜10（基盤1:1移植から「結 (YUI)」統合要求対応まで、全36 PR）が2026-09-20にすべて完了したことに伴い、従来の「作業指示書・引き継ぎ手順書」から、今後の開発・保守・YUI統合において参照されるべき**「アーキテクチャ方針・各フェーズ実績・教訓記録（Lessons Learned & Architecture Decision Record）」**として集約・再編したものである。

---

## 1. 共通運用ルール・リポジトリ構造

### 1.1 GitFlow 開発運用ルール (MUST)
Phase 1初期にブランチ運用違反による未レビューの積み上がりが起きた反省を踏まえ、マージ前レビューゲートが厳格にルール化された。`2026.9.0` リリース後は GitFlow 運用へと回帰している。

- **`main`**: 常にリリース済みの安定状態のみを反映（直接コミット禁止）。バージョンタグ（例: `2026.9.0`）はここに打つ。
- **`develop`**: 日常開発の統合ブランチ。機能ブランチはここから切り、ここへ合流する。
- **機能ブランチ（`feature/*`, `fix/*`, `chore/*`）**: `develop` から切り、実装・テスト・`cargo clippy`/`fmt` 通過後、**自己マージせず** Claude のレビュー承認を得てから合流する。

### 1.2 Cargo Workspace 構成と先読み設計
単一クレート直下ではなく、最初から `crates/` 配下にコアクレートを置く構成を採用した。
```text
rust/
├── Cargo.toml            # workspace定義
└── crates/
    ├── proteindf-bridge/     # コアライブラリ
    └── proteindf-bridge-py/  # PyO3 Pythonバインディング (Phase 5新設)
```
- **設計判断**: 将来のバインディング（PyO3、C ABI等）が別クレートとして追加されることを見越し、Phase 1開始時からマルチクレート構成を選択。これにより Phase 5（`proteindf-bridge-py`）追加時に既存の構成を一切破壊することなくスムーズに拡張できた。

### 1.3 命名規約
- **コアクレート**: `proteindf-bridge`（「PDF: Portable Document Format」との混同を避けるため、Phase 4.5で `pdf-bridge` から改名）。
- **Pythonバインディングパッケージ**: `proteindf_bridge_rs`（既存の純Python版 `proteindf_bridge` と共存インストール・並行比較できるように命名）。

---

## 2. 開発で得られた重要教訓・ハマりどころ集 (Lessons Learned & Pitfalls)

将来の実装や YUI 統合で同じ不具合や落とし穴を繰り返さないための、重要知見カタログ。

### A. データモデル・メモリ管理・所有権
1. **`AtomGroup` へのフィールド追加時は全マージ・集合演算を漏れなく更新する（Phase 1是正事項5、Phase 10 PR#29）**:
   - `AtomGroup` に新フィールドを追加した際、`merge`/`BitAnd`/`BitOr`/`BitXor`/`Clone` でコピー・合成処理が漏れると、マージや演算時にその情報がサイレントに消失する（Phase 1で `bonds` が消失していた実バグ）。Phase 10で `secondary_structure` フィールドを追加した際も、この教訓を引用して全演算の保持を検証した。
2. **階層構造を跨ぐ結合（SSBOND等）は共通祖先グループに相対パスで格納する（Phase 2 PR#4〜6）**:
   - ネストした階層（model/chain/residue）の異なる枝にまたがる結合を呼び出し元グループに絶対パスで格納すると、親グループの再配置時にパスが古いまま取り残される。共通祖先（`get_common_path`）を探索し、そこからの相対パスで格納して `get_bond_list()` で復元するルーティングが必須。
3. **0次元（空の系）の正常許容（Phase 1是正事項3, 10）**:
   - `Matrix::new` が `rows > 0 && cols > 0` を assert していたため、原子数0の系に対する `Bond::setup` でクラッシュした。空の構造や0次元行列は正常系として許容し、固有値計算（`eig`）の n=0 分岐の形状整合も保つ必要がある。
4. **サブグループ照合における name フォールバック（Phase 1是正事項6）**:
   - `AtomGroup::merge` で key 一致のみ見ていると、同一とみなすべきグループが兄弟グループとして分裂する。key 一致 OR `name` 一致のフォールバックが必要。
5. **コレクションの順序決定性（Phase 1是正事項9、Phase 7 PR#20、Phase 8 PR#22）**:
   - `HashSet` の利用や `IndexMap` の挿入順への暗黙的な依存は、残基走査や集合演算の順序を実行ごとに非決定的にする。残基走査では `sort_nicely` による自然順ソートを明示的に適用すること。
6. **大容量データでのメモリ保護（Phase 10 PR#35）**:
   - 100万原子規模では全原子ペアの $O(N^2)$ 密行列（`SymmetricMatrix`）は約 4 TB のメモリを消費し即座に OOM となる。公開定数 `MAX_DENSE_MATRIX_ATOMS = 2000` を設け、2,000原子超では密行列（`distmat`/`bondmat`）の構築をスキップして `None` を設定するフォールバックが不可欠。

### B. 幾何計算・アルゴリズム
7. **Python版 `superposer_quaternion.py` の非対称バグ発見とRust版での是正（Phase 3 PR#11）**:
   - Python版では `SymmetricMatrix.add()` が下三角に正規化されておらず、`eigh` が下三角のみを読むため四元数法の非対角成分が実質無視されていた（Kabsch法 RMSD≈6e-16 に対し四元数法 RMSD≈0.549 と乖離）。Rust版では実験的モジュールのバグを模倣せず、正しい対称書き込みを優先して Kabsch 法と完全一致（RMSD≈3e-16）させた。
8. **二面角計算の数値安定性と定義（Phase 7 PR#20）**:
   - IUPAC 標準に従い数値的に安定な `atan2` 方式を採用。最初の残基の φ、最後の残基の ψ は定義不能なため `Option<f64>`（`None`）とし、欠損残基は可視化目的のため安全にスキップする設計とした。
9. **Kabsch-Sander 水素結合・DSSP の判定バグ是正（Phase 8 PR#22, PR#23）**:
   - 欠損残基フィルタ後の配列インデックスで隣接残基判定（`|d-a| <= 2`）を行うと、鎖内に欠損残基が連続した場合に非隣接ペアが誤除外される。必ず元の残基番号（`orig_idx`）ベースで距離判定を行うこと。
   - 逆平行ブリッジの定義式におけるインデックスのずれを手計算で再導出して修正。また、可視化リボン用途に合わせ DSSP を 3状態（H/E/-）に意図的に簡略化した。
10. **CH-π 幾何判定の落とし穴（Phase 9 PR#25）**:
    - 環法線ベクトルは最小二乗平面フィット（SVD）で求めるが、法線の向き（符号）は数学的に不定であるため、角度判定では必ず内積の絶対値（`|dot|`）を評価すること。
    - `is_carbon_atom` が原子名 "C" 始まりだけで判定していたため、塩化物イオン "CL" を誤認していた。必ず `atomic_number() == 6` で判定すること。
11. **セルリスト近傍探索のスケーラビリティ（Phase 10 PR#35）**:
    - 動的セルサイズ（`(2 * max_vdw + 0.4).max(3.0)`）を用いた half-neighborhood（13方向）走査により、結合探索を $O(N)$ 化。コールバック走査 API ではインデックス順序規約 `i < j` を保証。

### C. フォーマット・パーサ・I/O
12. **mmCIF の実態と著者番号（auth）採用理由（Phase 4 PR#12, PR#13）**:
    - Python版 `SimpleMmcif` は単一リガンド（CCD形式）のみ対応で、全体構造（`_atom_site`）に未対応だった。
    - `label_asym_id`/`label_seq_id` をチェイン/残基キーに使うと、水分子等の残基がすべて 1 つに潰れてしまう。PDB互換の著者番号 `auth_asym_id`/`auth_seq_id` を採用することが不可欠。
13. **文字列処理・パーサのエラー握りつぶし防止（Phase 1是正事項4、Phase 2 PR#4）**:
    - `unwrap_or(0.0)` や `unwrap_or(0)` によるサイレントなフォールバックは、未定義元素の結合見落としや、GROのオーバーフロー（`*****`）が原子0番・残基0番として誤混入する原因となる。エラーは適切に伝播させること。
    - GRO 固定長パースでのバイト単位スライス（`line[0..5]`）は非ASCII文字で panic するため、`.chars()` ベースで処理すること。
14. **`.brd` (MessagePack) の互換性と YUI ヘッダー（Phase 6 PR#17, PR#18）**:
    - 既存 Python 版 `.brd` はヘッダーを持たないプレーン MessagePack であり、`modeling.py` が無条件にロードしていた。既存プレーン形式と YUI 互換ヘッダー（Magic+Version+zstd）形式を別 API として両立させた。
15. **明示的結合情報とヒューリスティック結合の優先順位ポリシー（Phase 10 PR#34）**:
    - MOL2, PRMTOP, PDB CONECT 由来の明示的な結合情報が存在する場合はそれを優先し、ヒューリスティックな `Bond::setup()` による上書き・重複定義を避ける利用規約を `RUST_PORT_SPEC.md` §3.8 に確立。
16. **テストフィクスチャの CWD 非依存パス解決（Phase 2）**:
    - フィクスチャの相対パス参照は実行時カレントディレクトリに依存して壊れやすいため、`env!("CARGO_MANIFEST_DIR")` を基点とする絶対パス解決に統一。

### D. 依存クレート選定・配布方針
17. **メンテナンス終了クレートの回避（Phase 9 PR#26）**:
    - YAML パーサの選定にあたり、メンテナンスが終了した `serde_yaml` を回避し、保守されているフォーク版 `serde_yaml_ng` を採用。
18. **wasm32 ターゲット依存の自動切り替え（Phase 10 PR#27）**:
    - wasm32 環境では C バインディングの `zstd` がビルドできないため、`Cargo.toml` の `target.'cfg(target_arch = "wasm32")'.dependencies` で純 Rust 実装の `ruzstd` が自動選択されるよう構成。
19. **クレート配布における GitHub git 依存の指針（Phase 10 PR#36）**:
    - 相対パス依存（ローカル開発）と crates.io 公開の中間の現実的選択肢として、Cargo の GitHub git 依存（`tag`, `branch` 指定）の指針とトレードオフ（SemVer 解決不可、プライベートリポジトリ時の認証要件等）を明文化。

### E. 検証方法論 (Methodology)
20. **「新規機能・テスト不足モジュールは独立検証で担保する」横断方針**:
    - テストが一切存在しないモジュール（`modeling.py`/`neutralize.py`）やテストが極めて薄いモジュール（`mmcif.py`）は、拙速な移植を避け、独立フェーズとして受け入れ基準とフィクスチャを整備してから着手した。
    - Python版に対応コードのない新規機能（Ramachandran、DSSP、側鎖H-bond、CH-π、明示結合I/O）では、「Python版との1:1比較」が使えないため、**「(1) 一次アルゴリズム・文献値に基づく独立手計算基準値（実データ照合）」＋「(2) 境界値・幾何学的サニティを保証する合成データテスト」** の2段階検証を徹底。これにより Python 版の潜在バグに引きずられず、高い信頼性を確保した。

---

## 3. 各フェーズの実績と設計記録 (Phase 1〜10)

全10フェーズ、計36 PRの完了実績サマリー。

### Phase 1: 基盤・データモデルの1:1移植 (PR#1〜3, 完了 2026-09-14)
- **移植モジュール**: `error.rs`, `periodic_table.rs`, `vector.rs`, `matrix.rs` (`Matrix`, `SymmetricMatrix`), `position.rs`, `atom.rs`, `bond.rs`, `atom_group.rs`。
- **成果**: 基礎データ構造と幾何数値計算の移植。発覚した11件の是正事項（0次元行列対応、`bonds`マージ漏れ修正、未定義元素VDWエラー伝播、`Position::from_str`補完等）を完了。

### Phase 2: フォーマットI/Oの1:1移植 (PR#4〜6, 完了 2026-09-14)
- **移植モジュール**: `format/mod.rs`, `format/xyz.rs`, `format/gro.rs`, `format/mol2.rs`, `format/amber_prmtop.rs`, `format/pdb.rs`。
- **成果**: 実PDB構造（`1hls.pdb`, `2MGO.pdb`, `3i3zH.pdb`）でPython版と完全一致を検証。階層結合ルーティング（共通祖先探索・相対パス格納）を確立。

### Phase 3: 構造操作の1:1移植 (PR#7〜11, 完了 2026-09-14)
- **移植モジュール**: `selector.rs`, `amino_acid.rs`, `ssbond.rs`, `ion_pair.rs`, `superposer.rs`, `superposer_quaternion.rs`。
- **成果**: 重ね合わせ演算（Kabsch法・四元数法）の実装。Python版四元数法の実バグ（対称行列加算不整合）を特定し是正。

### Phase 4: mmCIFサポート (PR#12〜13, 完了 2026-09-15)
- **移植・新規モジュール**: `format/mmcif.rs`。
- **成果**: 既存CCD形式（`ALA.cif`）の1:1移植（PR#12）に加え、`_atom_site` カテゴリの全体構造パーサを新規実装（PR#13）。著者番号（`auth_asym_id`/`auth_seq_id`）を採用し、PDBパーサ結果と原子単位で完全一致を達成。

### Phase 4.5: プロジェクト名変更 (完了 2026-09-15)
- **成果**: `pdf-bridge` から `proteindf-bridge` への名称変更（Cargo.toml, ディレクトリ, テストコード一括更新）。

### Phase 5: Pythonバインディング (PyO3) (PR#14〜16, 完了 2026-09-15)
- **新規クレート**: `rust/crates/proteindf-bridge-py`（パッケージ名: `proteindf_bridge_rs`）。
- **成果**: 基盤データモデル、フォーマットI/O、構造操作の Python バインディングを提供。既存 Python 版との pytest 直接比較（42件）で一致を検証。

### Phase 6: `modeling.py`/`neutralize.py` (PR#17〜19, 完了 2026-09-15)
- **新規モジュール**: `brd.rs`, `modeling.rs`, `neutralize.rs`。
- **成果**: MessagePack 往復 I/O と YUI ヘッダー形式（zstd圧縮対応）の実装。ACE/NME 末端キャッピングおよび中性化イオン配置アルゴリズムを完全再現。

### Phase 7: バックボーン二面角計算 — Ramachandranプロット対応 (PR#20〜21, 完了 2026-09-15)
- **新規モジュール**: `ramachandran.rs`, Pythonバインディング。
- **成果**: IUPAC 標準 `atan2` による主鎖 φ/ψ 二面角計算関数 `calc_phi_psi` を新規実装。実データ基準値（`1hls.pdb` chain A）と誤差 1e-3 度以内で一致。

### Phase 8: 主鎖水素結合検出 + 二次構造推定 (DSSP) (PR#22〜23, 完了 2026-09-15)
- **新規モジュール**: `hydrogen_bond.rs`, `secondary_structure.rs`。
- **成果**: Kabsch-Sander 静電エネルギーモデルによる水素結合検出および DSSP アルゴリズム（H/E/- の 3状態分類）を新規実装。PyDSSP 参照値と完全一致。

### Phase 9: 側鎖水素結合 + CH-π相互作用検出 (PR#24〜26, 完了 2026-09-15)
- **新規モジュール**: `hydrogen_bond.rs`（側鎖拡張）, `ch_pi.rs`, `interaction_set.rs`。
- **成果**: 側鎖水素結合、CH-π 相互作用（SVD平面フィット法線ベクトル評価）、および相互作用を集約・MessagePack/YAML で往復可能な `InteractionSet` を実装。

### Phase 10: 「結 (YUI)」統合要求への対応 (PR#27〜36, 完了 2026-09-20)
- **PR#27**: wasm32 ターゲットでの `ruzstd` 自動選択設定。
- **PR#28**: protein schema 階層規約（`/model_N/chain_id/res_key/atom_key`）の形式化、位置判定 `is_*_level()` と構造判定の責務分離、違反検出ヘルパー `validate_schema()`。
- **PR#29**: 二次構造情報の `AtomGroup` への書き戻し（`secondary_structure: Option<SsCode>` フィールド追加、`apply_secondary_structure`）。
- **PR#30**: パスベース `BondRecord` 解決コストの検証（$O(\text{深さ})$ ルックアップの実証）と `resolve_bond` ヘルパー追加。
- **PR#31〜33**: ファイル由来明示結合のパース（MOL2, PRMTOP `BONDS_*`, PDB `CONECT`）。
- **PR#34**: 明示的結合情報とヒューリスティック判定の優先順位ポリシー確立（`RUST_PORT_SPEC.md` §3.8）。
- **PR#35**: `Bond::setup()` のスケーラビリティ改善（空間セルリスト探索導入、100万原子ベンチマーク達成）。
- **PR#36**: クレート配布方式（相対パス・GitHub git依存）および Python バインディング役割分担指針の策定。

---

## 4. スケーラビリティ・ベンチマーク実績

Phase 10 において実証された、大規模構造処理に関する性能指標。

### 100万原子規模での `Bond::setup()` ベンチマーク (`test_benchmark_1m_atoms`)
- **条件**: 1,000,000 原子（100×100×100 グリッド、格子間隔 2.0 Å、炭素原子）
- **原子生成時間**: 308.94 ms
- **`Bond::setup()` 実行時間**: **3.39 秒**（リリースビルド）
- **検出結合数**: 12,731,796 本
- **メモリ消費量（概算）**:
  - 従来の $O(N^2)$ 密行列（約 4 TB 必要）を回避。
  - 座標配列（24 MB）とセルリストのグリッドを含め、数十 MB 程度で完走。
- **安全保護**: `MAX_DENSE_MATRIX_ATOMS = 2000` を超える原子数では `distmat`/`bondmat` を `None` とし OOM を防止。

### パスベース原子解決性能 (`test_get_atom_by_path_depth_scaling`)
- `AtomGroup::get_atom_by_path` の計算量は原子総数 $N$ に依存せず、パスの深さ $D$（通常 4〜5）と各ノードの `IndexMap` ルックアップ $O(1)$ の積 $O(D)$ で動作することを実証済み。

---

## 5. スコープ外・将来課題（未決事項）

本移植プロジェクト（Phase 1〜10）では着手せず、今後の独立タスクまたは必要時に検討すべき事項。

- **量子化学計算（QC）結果 I/O・CUBE パーサ（`RUST_PORT_SPEC.md` §3.4〜3.6）**:
  - `OrbitalContribution`, HDF5 結果ファイルの I/O, Gaussian/ProteinDF CUBE 形式のパース。
- **C/C++ バインディング（`RUST_PORT_SPEC.md` §4）**:
  - `cdylib` + `cbindgen` による C ABI ヘッダー生成。
- **内部数値計算ライブラリの検討（`RUST_PORT_SPEC.md` §8）**:
  - 現在の自前実装（`Vector`, `Matrix`）を `nalgebra` 等の外部クレートへ置き換えるかの是非（現状で十分な性能が出ており、必要性が生じた時点で検討）。
- **近傍探索最適化の他モジュール展開**:
  - `secondary_structure.rs`, `hydrogen_bond.rs`, `ch_pi.rs` へのセルリスト適用（現時点では対象構造規模において十分高速であるため未適用）。
