#!/usr/bin/env bash
# Compares the speed of two releases of proteindf_bridge from Python, for example the
# last pure-Python AtomGroup release (2026.10.13) against a Rust-backed one.
# Downloads each release's wheel for this platform from GitHub (gh), installs it in a
# temporary venv, runs devtool/bench_switch.py on each structure file, and prints a
# table of median seconds per operation and the ratio new/old (< 1 means faster).
#
# Usage: devtool/bench-switch.sh [-r N] [-o op1,op2] <old> <new> [file...]
#   -r N   repeat count per operation (default 5; bond_setup takes tens of seconds)
#   -o     comma-separated subset of operations (default: all)
#   file   PDB files (default: 2FB4.pdb and 2MGO.pdb from the test data)
# Raw JSON results are kept in the printed temporary directory; the venvs are removed.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
repeat=5
ops=""
while getopts "r:o:h" opt; do
    case "$opt" in
        r) repeat=$OPTARG ;;
        o) ops=$OPTARG ;;
        *) sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
    esac
done
shift $((OPTIND - 1))
(($# >= 2)) || { sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
old=$1 new=$2
shift 2
files=()
for f in "$@"; do files+=("$(cd "$(dirname "$f")" && pwd)/$(basename "$f")"); done
if ((${#files[@]} == 0)); then
    files=("$root/rust/crates/proteindf-bridge/tests/data/2FB4.pdb"
           "$root/rust/crates/proteindf-bridge/tests/data/2MGO.pdb")
fi

case "$(uname -s)-$(uname -m)" in
    Linux-x86_64) pattern="*manylinux*x86_64.whl" ;;
    Linux-aarch64) pattern="*manylinux*aarch64.whl" ;;
    Darwin-arm64) pattern="*macosx*arm64.whl" ;;
    *) echo "no release wheel for $(uname -s)-$(uname -m)" >&2; exit 1 ;;
esac

work=$(mktemp -d -t bench-switch-XXXXXX)
echo "==> work dir: $work" >&2
trap 'rm -rf "$work/$old" "$work/$new"' EXIT
for v in "$old" "$new"; do
    mkdir -p "$work/$v/dl"
    gh release download "$v" -R ProteinDF/ProteinDF_bridge -p "$pattern" -D "$work/$v/dl" >&2
    python3 -m venv "$work/$v/venv"
    "$work/$v/venv/bin/pip" install -q "$work/$v/dl/"*.whl >&2
done

# Copy the runner out of the repository and run it from the work dir, so the installed
# package is imported and not the working tree.
cp "$here/bench_switch.py" "$work/bench_switch.py"
for f in "${files[@]}"; do
    for v in "$old" "$new"; do
        echo "==> $v $(basename "$f")" >&2
        (cd "$work" && "$work/$v/venv/bin/python" -I bench_switch.py "$f" \
            --repeat "$repeat" --ops "$ops") > "$work/$v-$(basename "$f").json"
    done
done

python3 - "$work" "$old" "$new" "${files[@]}" <<'EOF'
import json, os, sys
work, old, new, files = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
for f in files:
    name = os.path.basename(f)
    a = json.load(open(os.path.join(work, f"{old}-{name}.json")))
    b = json.load(open(os.path.join(work, f"{new}-{name}.json")))
    print(f"\n{name} ({a['atoms']} atoms, median of {a['repeat']})")
    print(f"{'operation':<18}{old:>12}{new:>12}{'new/old':>9}")
    for op in list(a["seconds"]) + [k for k in b["seconds"] if k not in a["seconds"]]:
        ta = a["seconds"].get(op, {}).get("median")
        tb = b["seconds"].get(op, {}).get("median")
        cell = lambda t: f"{t:>11.4f}s" if t is not None else f"{'-':>12}"
        ratio = f"{tb / ta:>9.2f}" if ta and tb is not None else f"{'-':>9}"
        print(f"{op:<18}{cell(ta)}{cell(tb)}{ratio}")
EOF
