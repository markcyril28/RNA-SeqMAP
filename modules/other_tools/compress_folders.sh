#!/usr/bin/env bash
# compress_folders.sh - Compress specified folders using 7z LZMA2 at maximum compression with 64 threads
set -euo pipefail

SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)"
cd "$SCRIPT_DIR"

THREADS=64
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT="HeatSeq_archive_${TIMESTAMP}.7z"

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

# Show total size before compression
du -sh "${FOLDERS[@]}" 2>/dev/null
echo ""
TOTAL=$(du -sc "${FOLDERS[@]}" 2>/dev/null | tail -1 | cut -f1)
echo "Total size: $(numfmt --to=iec --from-unit=1024 "$TOTAL" 2>/dev/null || echo "${TOTAL}K")"
echo ""

time 7z a -t7z -m0=lzma2 -mx=9 -mfb=273 -md=64m -ms=on -mmt="$THREADS" \
    "$OUTPUT" "${FOLDERS[@]}"

echo ""
echo "=== Done ==="
ls -lh "$OUTPUT"