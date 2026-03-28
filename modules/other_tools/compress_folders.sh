#!/usr/bin/env bash
# compress_folders.sh - Compress specified folders using 7z LZMA2 at maximum compression with 64 threads
set -euo pipefail

SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)" || { echo "[ERROR] compress_folders.sh: Failed to resolve script directory" >&2; exit 1; }
# Resolve project root (this script lives in modules/other_tools/)
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)" || { echo "[ERROR] compress_folders.sh: Failed to resolve PROJECT_ROOT" >&2; exit 1; }
cd "$PROJECT_ROOT" || { echo "[ERROR] compress_folders.sh: Cannot cd to PROJECT_ROOT: $PROJECT_ROOT" >&2; exit 1; }

# Respect pre-set THREADS; check HPC scheduler vars before nproc
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-${PBS_NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 64)}}}"
printf -v TIMESTAMP '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || TIMESTAMP=$(date +%Y%m%d_%H%M%S)
mkdir -p "$PROJECT_ROOT/HPC"
OUTPUT="$PROJECT_ROOT/HPC/HeatSeq_archive_${TIMESTAMP}.7z"

FOLDERS=(
    "${SRR_OUTPUT_ROOT:-$PROJECT_ROOT/1_SRRs}/C_FastQC"
    "${ALIGNMENT_RESULTS_ROOT:-$PROJECT_ROOT/2_ALIGNMENT_RESULTs}"
    "${POST_PROCESSING_ROOT:-$PROJECT_ROOT/3_POST_PROC}"
    "${CONCORDANCE_OUTPUT_ROOT:-$PROJECT_ROOT/4_CONCORDANCE_ANALYSIS}"
    "$PROJECT_ROOT/logs"
    "$PROJECT_ROOT/z_archive"
)

# Filter to existing folders, warn about missing ones
_valid_folders=()
for folder in "${FOLDERS[@]}"; do
    if [[ -d "$folder" ]]; then
        _valid_folders+=("$folder")
    else
        echo "[WARN] compress_folders.sh: Skipping missing directory: $folder" >&2
    fi
done
if [[ ${#_valid_folders[@]} -eq 0 ]]; then
    echo "[ERROR] compress_folders.sh: No directories found to compress" >&2
    exit 1
fi
FOLDERS=("${_valid_folders[@]}")

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