#!/usr/bin/env bash
# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

set -euo pipefail

if [ $# -lt 1 ] || [ -z "$1" ]; then
    echo "Error: Output directory argument required." >&2
    echo "Usage: $0 <output_directory>" >&2
    exit 1
fi

OUTPUT_DIR="$1"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

# Verify PDF_HOME is set
if [ -z "${PDF_HOME:-}" ]; then
    echo "Error: PDF_HOME environment variable is not set." >&2
    echo "Please set PDF_HOME to the ProteinDF installation directory." >&2
    exit 1
fi

# Verify required commands exist in PATH
for cmd in PDF.x pdf-setup.sh pdf-archive-h5.py; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Error: Required command '$cmd' not found in PATH." >&2
        exit 1
    fi
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGETS=("CH3COOX_HF" "GLY")

for target in "${TARGETS[@]}"; do
    input_file="$SCRIPT_DIR/$target/fl_Userinput"
    if [ ! -f "$input_file" ]; then
        echo "Error: Input file '$input_file' not found." >&2
        exit 1
    fi

    echo "=== Regenerating $target ==="
    work_dir="$(mktemp -d -t "pdf_regen_${target}_XXXXXX")"
    
    # Ensure temporary working directory is cleaned up
    cleanup() {
        rm -rf "$work_dir"
    }
    trap cleanup EXIT INT TERM

    cp "$input_file" "$work_dir/fl_Userinput"

    (
        cd "$work_dir"
        pdf-setup.sh
        PDF.x > pdf.log 2>&1 || {
            echo "Error: PDF.x execution failed for $target." >&2
            cat pdf.log >&2
            exit 1
        }

        if ! grep -q "SCF is well converged" pdf.log; then
            echo "Error: SCF did not converge for $target ('SCF is well converged' not found in pdf.log)." >&2
            cat pdf.log >&2
            exit 1
        fi

        pdf-archive-h5.py -f -o "$OUTPUT_DIR/$target.h5"
    )

    rm -rf "$work_dir"
    trap - EXIT INT TERM

    echo "Successfully generated: $OUTPUT_DIR/$target.h5"
done

echo "All targets regenerated successfully in $OUTPUT_DIR"
