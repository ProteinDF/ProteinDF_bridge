#!/usr/bin/env bash
# scripts/ci-local.sh: Run CI checks locally
set -euo pipefail

# Always run from the repository root (CI_LOCAL_ROOT checks another worktree with
# this copy of the script; devtool/review.sh uses it)
REPO_ROOT="$(cd "${CI_LOCAL_ROOT:-$(dirname "${BASH_SOURCE[0]}")/..}" && pwd)"
cd "$REPO_ROOT"

# Prevent mise from auto-creating .venv in the repository
export MISE_PYTHON_UV_VENV_AUTO=false
export MISE_PYTHON_VENV_AUTO=false

# Prevent Python bytecode caching (__pycache__) in repository
export PYTHONDONTWRITEBYTECODE=1

# Configuration
PYTHON="${PYTHON:-python3}"
FAIL_FAST=0

# Temporary directory outside the repository (cleaned up on exit)
TMP_ROOT="$(mktemp -d -t proteindf-bridge-ci-XXXXXX)"
cleanup() {
    rm -rf "$TMP_ROOT"
}
trap cleanup EXIT

# Available and default jobs
ALL_DEFAULT_JOBS=("rust" "hdf5" "wasm" "lint" "wheel" "sdist" "docs")
REQUESTED_JOBS=()

show_help() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] [JOB...]

Run ProteinDF_bridge CI checks locally.
If no jobs are specified, all default jobs are executed:
  ${ALL_DEFAULT_JOBS[*]}

Options:
  -h, --help       Show this help message and exit
  --fail-fast      Stop immediately on the first failed job.
                   By default, all requested jobs are executed and a summary
                   is displayed at the end.

Available Jobs:
  rust             Run cargo fmt, cargo clippy, cargo nextest, and cargo test --doc
  hdf5             Run cargo clippy and cargo nextest with proteindf-bridge/hdf5 feature
  wasm             Run cargo check for wasm32-unknown-unknown target
  lint             Run ruff check (F821, PLE0704) on proteindf_bridge in isolated venv
  wheel            Build release wheel in isolated venv and run unit tests + CLI checks
  sdist            Build sdist package in isolated venv and run unit tests + CLI checks
  docs             Build English and Japanese Sphinx documentation in isolated venv
  qc-data          Run #[ignore] tests with hdf5 feature (requires PROTEINDF_BRIDGE_QC_DATA)
                   (Note: Not run by default; must be explicitly specified)
  all              Run all default jobs (${ALL_DEFAULT_JOBS[*]})

Environment Variables:
  PYTHON                      Python interpreter to use (default: python3)
  CI_LOCAL_ROOT               Repository root to check (default: the one containing this script)
  PROTEINDF_BRIDGE_QC_DATA    Path to directory containing regenerated QC .h5 files (required for qc-data job)

Prerequisites:
  - cargo, rustc (Rust toolchain)
  - cargo-nextest (cargo install cargo-nextest)
  - wasm32-unknown-unknown target for wasm job (rustup target add wasm32-unknown-unknown)
  - Python 3.9+ with venv support
  - maturin (auto-installed in isolated venv if not found in PATH)

Failure Handling:
  By default, all requested jobs are executed even if an earlier job fails.
  At the end of execution, a summary table showing the result (SUCCESS/FAILED)
  and duration for each job is displayed. If any job failed, the script exits
  with status code 1.
EOF
}

# Parse options and arguments
while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help)
            show_help
            exit 0
            ;;
        --fail-fast)
            FAIL_FAST=1
            shift
            ;;
        all)
            REQUESTED_JOBS+=("${ALL_DEFAULT_JOBS[@]}")
            shift
            ;;
        rust|hdf5|wasm|lint|wheel|sdist|docs|qc-data)
            REQUESTED_JOBS+=("$1")
            shift
            ;;
        *)
            echo "Error: Unknown argument '$1'" >&2
            echo "Run '$(basename "$0") --help' for usage." >&2
            exit 1
            ;;
    esac
done

if [ ${#REQUESTED_JOBS[@]} -eq 0 ]; then
    REQUESTED_JOBS=("${ALL_DEFAULT_JOBS[@]}")
fi

# Tool requirement check
check_tool() {
    local cmd="$1"
    local hint="$2"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Error: Required command '$cmd' not found in PATH." >&2
        echo "Hint: $hint" >&2
        return 1
    fi
    return 0
}

# Verify prerequisites for requested jobs
verify_prerequisites() {
    local has_rust=0
    local has_hdf5=0
    local has_wasm=0
    local has_py=0
    local has_qc=0
    local missing=0

    for job in "${REQUESTED_JOBS[@]}"; do
        case "$job" in
            rust) has_rust=1 ;;
            hdf5) has_hdf5=1 ;;
            wasm) has_wasm=1 ;;
            lint|wheel|sdist|docs) has_py=1 ;;
            qc-data) has_qc=1 ;;
        esac
    done

    if [ $has_rust -eq 1 ] || [ $has_hdf5 -eq 1 ] || [ $has_wasm -eq 1 ] || [ $has_qc -eq 1 ]; then
        check_tool "cargo" "Please install Rust toolchain (https://rustup.rs/)" || missing=1
    fi

    if [ $has_rust -eq 1 ] || [ $has_hdf5 -eq 1 ] || [ $has_qc -eq 1 ]; then
        check_tool "cargo-nextest" "Install with: cargo install cargo-nextest" || missing=1
    fi

    if [ $has_wasm -eq 1 ]; then
        check_tool "rustup" "rustup is required to check wasm32-unknown-unknown target." || missing=1
        if command -v rustup >/dev/null 2>&1; then
            if ! rustup target list --installed | grep -q "^wasm32-unknown-unknown$"; then
                echo "Error: Target 'wasm32-unknown-unknown' is not installed." >&2
                echo "Hint: Install it with: rustup target add wasm32-unknown-unknown" >&2
                missing=1
            fi
        fi
    fi

    if [ $has_py -eq 1 ]; then
        check_tool "$PYTHON" "Python interpreter '$PYTHON' is required." || missing=1
        if command -v "$PYTHON" >/dev/null 2>&1; then
            if ! "$PYTHON" -m venv -h >/dev/null 2>&1; then
                echo "Error: '$PYTHON' does not have 'venv' module available." >&2
                echo "Hint: Install python3-venv package or use a Python with venv support." >&2
                missing=1
            fi
        fi
    fi

    if [ $has_qc -eq 1 ]; then
        if [ -z "${PROTEINDF_BRIDGE_QC_DATA:-}" ]; then
            echo "Error: Environment variable 'PROTEINDF_BRIDGE_QC_DATA' is not set." >&2
            echo "Hint: Set PROTEINDF_BRIDGE_QC_DATA to run qc-data tests." >&2
            missing=1
        fi
    fi

    if [ $missing -ne 0 ]; then
        exit 1
    fi
}

verify_prerequisites

# Job implementations
run_job_rust() {
    echo "==> Running cargo fmt check..."
    cargo fmt --manifest-path rust/Cargo.toml --all -- --check

    echo "==> Running cargo clippy..."
    cargo clippy --workspace --all-targets --manifest-path rust/Cargo.toml -- -D warnings

    echo "==> Running cargo nextest..."
    cargo nextest run --workspace --manifest-path rust/Cargo.toml --status-level fail --final-status-level fail

    echo "==> Running cargo doc tests..."
    cargo test --doc --workspace --manifest-path rust/Cargo.toml
}

run_job_hdf5() {
    echo "==> Running cargo clippy (with hdf5 feature)..."
    cargo clippy --workspace --all-targets --manifest-path rust/Cargo.toml --features proteindf-bridge/hdf5 -- -D warnings

    echo "==> Running cargo nextest (with hdf5 feature)..."
    cargo nextest run --workspace --manifest-path rust/Cargo.toml --features proteindf-bridge/hdf5 --status-level fail --final-status-level fail
}

run_job_wasm() {
    echo "==> Checking wasm32-unknown-unknown target build..."
    cargo check --manifest-path rust/Cargo.toml --target wasm32-unknown-unknown -p proteindf-bridge --no-default-features --features ruzstd
}

# Helper to create temporary virtualenv and prepare maturin
setup_python_env() {
    local env_dir="$1"
    local need_maturin="${2:-0}"

    mkdir -p "$env_dir"
    "$PYTHON" -m venv "$env_dir/venv"
    local venv_pip="$env_dir/venv/bin/pip"

    if [ "$need_maturin" -eq 1 ]; then
        if command -v maturin >/dev/null 2>&1; then
            MATURIN_BIN="maturin"
        else
            echo "==> Installing maturin in temporary virtualenv..."
            "$venv_pip" install maturin
            MATURIN_BIN="$env_dir/venv/bin/maturin"
        fi
    fi
}

run_job_lint() {
    local job_dir="$TMP_ROOT/lint"
    setup_python_env "$job_dir" 0

    local venv_pip="$job_dir/venv/bin/pip"
    local venv_ruff="$job_dir/venv/bin/ruff"

    echo "==> Installing ruff in isolated venv..."
    "$venv_pip" install --quiet "ruff==0.17.0"

    echo "==> Running ruff check (F821, PLE0704) on proteindf_bridge..."
    "$venv_ruff" check --select F821,PLE0704 proteindf_bridge
}

run_job_wheel() {
    local job_dir="$TMP_ROOT/wheel"
    setup_python_env "$job_dir" 1

    local venv_pip="$job_dir/venv/bin/pip"
    local venv_python="$job_dir/venv/bin/python"
    local dist_dir="$job_dir/dist"
    mkdir -p "$dist_dir"

    echo "==> Building wheel with maturin..."
    "$MATURIN_BIN" build --release --out "$dist_dir"

    echo "==> Installing wheel in isolated venv..."
    "$venv_pip" install "$dist_dir"/*.whl

    echo "==> Running Python unit tests..."
    "$venv_python" -B -P -m unittest discover -s tests

    echo "==> Verifying CLI scripts (--help)..."
    "$job_dir/venv/bin/pdb2brd.py" --help >/dev/null
    "$job_dir/venv/bin/brd2pdb.py" --help >/dev/null
    "$job_dir/venv/bin/neutralize.py" --help >/dev/null
}

run_job_sdist() {
    local job_dir="$TMP_ROOT/sdist"
    setup_python_env "$job_dir" 1

    local venv_pip="$job_dir/venv/bin/pip"
    local venv_python="$job_dir/venv/bin/python"
    local dist_dir="$job_dir/dist"
    mkdir -p "$dist_dir"

    echo "==> Building sdist with maturin..."
    "$MATURIN_BIN" sdist --out "$dist_dir"

    echo "==> Installing sdist in isolated venv..."
    "$venv_pip" install "$dist_dir"/*.tar.gz

    echo "==> Running Python unit tests..."
    "$venv_python" -B -P -m unittest discover -s tests

    echo "==> Verifying CLI scripts (--help)..."
    "$job_dir/venv/bin/pdb2brd.py" --help >/dev/null
    "$job_dir/venv/bin/brd2pdb.py" --help >/dev/null
    "$job_dir/venv/bin/neutralize.py" --help >/dev/null
}

run_job_docs() {
    local job_dir="$TMP_ROOT/docs"
    setup_python_env "$job_dir" 0

    local venv_pip="$job_dir/venv/bin/pip"
    local venv_sphinx="$job_dir/venv/bin/sphinx-build"

    local job_src="$job_dir/src"
    mkdir -p "$job_src"
    echo "==> Copying working tree files to temporary directory for isolated build..."
    cp -R docs pyproject.toml README.md proteindf_bridge "$job_src/"
    if [ -f LICENSE ]; then
        cp LICENSE "$job_src/"
    fi
    mkdir -p "$job_src/rust"
    cp -R rust/Cargo.toml rust/Cargo.lock rust/crates "$job_src/rust/"

    echo "==> Installing package with docs dependencies in isolated venv..."
    "$venv_pip" install "$job_src"[docs]

    local out_en="$job_dir/out/en"
    local out_ja="$job_dir/out/ja"
    mkdir -p "$out_en" "$out_ja"

    echo "==> Building English Sphinx HTML..."
    "$venv_sphinx" -b html -D language=en "$job_src/docs" "$out_en"

    echo "==> Building Japanese Sphinx HTML..."
    "$venv_sphinx" -b html -D language=ja "$job_src/docs" "$out_ja"
}

run_job_qc_data() {
    echo "==> Running #[ignore] tests with proteindf-bridge/hdf5 feature..."
    # Note: Runs existing #[ignore] benchmark/heavy tests (and future QC tests
    # requiring PROTEINDF_BRIDGE_QC_DATA planned for PR#61).
    cargo nextest run --workspace --manifest-path rust/Cargo.toml --features proteindf-bridge/hdf5 --run-ignored only --status-level fail --final-status-level fail
}

# Results tracking
JOB_RESULTS=()
JOB_DURATIONS=()
TOTAL_FAILED=0
SUITE_START=$SECONDS

for job in "${REQUESTED_JOBS[@]}"; do
    echo "============================================================"
    echo "Starting job: $job"
    echo "============================================================"
    job_start=$SECONDS
    job_status="SUCCESS"

    set +e
    (
        set -euo pipefail
        case "$job" in
            rust) run_job_rust ;;
            hdf5) run_job_hdf5 ;;
            wasm) run_job_wasm ;;
            lint) run_job_lint ;;
            wheel) run_job_wheel ;;
            sdist) run_job_sdist ;;
            docs) run_job_docs ;;
            qc-data) run_job_qc_data ;;
        esac
    )
    exit_code=$?
    set -e

    job_duration=$((SECONDS - job_start))

    if [ $exit_code -ne 0 ]; then
        job_status="FAILED"
        TOTAL_FAILED=$((TOTAL_FAILED + 1))
        echo "==> Job '$job' FAILED (exit code: $exit_code)" >&2
    else
        echo "==> Job '$job' SUCCESS (${job_duration}s)"
    fi

    JOB_RESULTS+=("$job_status")
    JOB_DURATIONS+=("$job_duration")

    if [ $exit_code -ne 0 ] && [ $FAIL_FAST -eq 1 ]; then
        echo "==> --fail-fast specified, aborting remaining jobs." >&2
        break
    fi
done

TOTAL_DURATION=$((SECONDS - SUITE_START))

# Summary output
echo ""
echo "============================================================"
echo "CI Local Summary"
echo "============================================================"
printf "%-12s %-10s %s\n" "Job" "Status" "Duration"
printf "%-12s %-10s %s\n" "---" "------" "--------"

for i in "${!JOB_RESULTS[@]}"; do
    job_name="${REQUESTED_JOBS[$i]}"
    status="${JOB_RESULTS[$i]}"
    duration="${JOB_DURATIONS[$i]}s"
    printf "%-12s %-10s %6s\n" "$job_name" "$status" "$duration"
done

echo "============================================================"
if [ $TOTAL_FAILED -eq 0 ]; then
    echo "Result: ALL PASSED (total: ${TOTAL_DURATION}s)"
    exit 0
else
    echo "Result: $TOTAL_FAILED JOB(S) FAILED (total: ${TOTAL_DURATION}s)" >&2
    exit 1
fi
