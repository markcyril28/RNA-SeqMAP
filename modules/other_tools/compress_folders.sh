#!/usr/bin/env bash
# compress_folders.sh - Compress specified folders using 7z LZMA2 at maximum compression with 64 threads
set -euo pipefail

SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)"
# Resolve project root (this script lives in modules/other_tools/)
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Use available cores (nproc), fallback to 64 for systems without nproc
THREADS=$(nproc 2>/dev/null || echo 64)
printf -v TIMESTAMP '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || TIMESTAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p HPC
OUTPUT="HPC/HeatSeq_archive_${TIMESTAMP}.7z"

FOLDERS=(
    "1_SRRs/3_FastQC_v2"
    "1_SRRs/C_FastQC"
    "2_ALIGNMENT_RESULTs"
    "3_POST_PROC"
    "logs"
    "z_archive"
)

# Verify all folders exist
for folder in "${FOLDERS[@]}"; do
    if [[ ! -d "$folder" ]]; then
        echo "ERROR: Directory not found: $folder"
        exit 1
    fi
done

echo "=== Compressing ${#FOLDERS[@]} folders into $OUTPUT ==="
echo "    Compression: LZMA2, ultra (mx=9), threads=$THREADS"
echo "    Folders: ${FOLDERS[*]}"
echo ""

# Show total size before compression — single du pass instead of two O(tree) traversals
mapfile -t _du_lines < <(du -shc "${FOLDERS[@]}" 2>/dev/null)
# Print per-folder lines (all but last "total" line)
printf '%s\n' "${_du_lines[@]::${#_du_lines[@]}-1}"
echo ""
TOTAL="${_du_lines[-1]%%$'\t'*}"
echo "Total size: $TOTAL"
echo ""

time 7z a -t7z -m0=lzma2 -mx=9 -mfb=273 -md=64m -ms=on -mmt="$THREADS" \
    "$OUTPUT" "${FOLDERS[@]}"

echo ""
echo "=== Done ==="
ls -lh "$OUTPUT"