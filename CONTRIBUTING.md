# Contributing

## Branching model (GitFlow)

This project follows a GitFlow-style branching model (2026-09-15 onward, reverting to the
convention used through `2024.03.0` after a period of committing directly to `main`).

- **`main`** — always reflects the latest released state. Nothing is committed here directly;
  it only receives merges from `release/*` and `hotfix/*` branches, each merge tagged with the
  version number (e.g. `2026.9.0`, no `v` prefix).
- **`develop`** — the integration branch for ongoing development. Day-to-day work targets this
  branch, not `main`.
- **`feature/*`** — branch off `develop` for a unit of work, merge back into `develop` once
  reviewed.
- **`release/X.Y.Z`** — branched from `develop` when preparing a release. Final adjustments
  (version bump, changelog, etc.) happen here, then it merges into both `main` (tagged) and
  `develop`.
- **`hotfix/*`** — branched from `main` for urgent fixes to an already-released version. Merges
  into both `main` (tagged) and `develop`.

Versioning follows `YYYY.M.PATCH` (calendar versioning): the first release in a given
year/month is `PATCH=0` (e.g. `2026.9.0`), subsequent releases in the same month increment
`PATCH` (e.g. `2026.9.1`). Keep `proteindf_bridge/_version.py` and the `version` field in
`pyproject.toml` in sync (checked by `tests/test_version.py`).

## Running Tests & CI

Continuous Integration is primarily run locally via `scripts/ci-local.sh`. GitHub Actions CI (`.github/workflows/ci.yml`) is retained for manual verification (`workflow_dispatch`) only.

### Prerequisites

- **Python 3.9+** with `venv` module (specified via `$PYTHON`, defaults to `python3`)
- **Rust toolchain** (`cargo`, `rustc`, `rustup`)
- **cargo-nextest**: `cargo install cargo-nextest`
- **wasm32-unknown-unknown target** (for `wasm` job): `rustup target add wasm32-unknown-unknown`
- **maturin** (used for `wheel` and `sdist` jobs; automatically installed into the isolated temporary virtualenv if not in `$PATH`)

### Running Local CI

Run all default CI jobs (`rust`, `hdf5`, `wasm`, `wheel`, `sdist`, `docs`):
```bash
./scripts/ci-local.sh
```

You can also run specific jobs by name:
```bash
./scripts/ci-local.sh rust hdf5
```

Available jobs:
- `rust`: Rust formatting check (`cargo fmt`), clippy (`-D warnings`), unit tests (`cargo nextest`), and doc tests (`cargo test --doc`).
- `hdf5`: Clippy and tests with `proteindf-bridge/hdf5` feature enabled (`cargo nextest`).
- `wasm`: Compilation check for `wasm32-unknown-unknown` target.
- `wheel`: Builds a release wheel in an isolated temporary virtual environment, installs it, runs `unittest`, and verifies CLI entry points (`--help`).
- `sdist`: Builds a source distribution in an isolated temporary virtual environment, installs it, runs `unittest`, and verifies CLI entry points (`--help`).
- `docs`: Builds English and Japanese Sphinx documentation in an isolated temporary virtual environment.
- `qc-data`: Runs `#[ignore]` tests with the `hdf5` feature. Not run by default; requires setting `PROTEINDF_BRIDGE_QC_DATA` to the directory containing regenerated `.h5` files:
  ```bash
  PROTEINDF_BRIDGE_QC_DATA=/path/to/qc_data ./scripts/ci-local.sh qc-data
  ```

Options:
- `--fail-fast`: Stop immediately upon the first failed job (by default, all requested jobs run and a summary table is printed at the end).
- `--help`: Display usage and detailed information.

See `docs/rust-port-handoff.md` for the Rust port (`rust/`) subproject's specific PR/review
workflow, which follows this same branching model.
