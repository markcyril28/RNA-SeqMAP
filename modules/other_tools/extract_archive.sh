#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)"
ARCHIVE_DIR="${SCRIPT_DIR}/HPC"

# Source logging utilities for consistent output
source "${SCRIPT_DIR}/modules/logging/logging_utils.sh" 2>/dev/null || {
    log_info()  { echo "[INFO]  $*"; }
    log_warn()  { echo "[WARN]  $*"; }
    log_error() { echo "[ERROR] $*" >&2; }
    log_step()  { echo ""; echo "==> $*"; }
}

# Find the latest .7z archive (avoid hardcoding timestamped filenames)
if [[ -n "${1:-}" && -f "$1" ]]; then
	ARCHIVE="$1"
elif [[ -d "$ARCHIVE_DIR" ]]; then
	ARCHIVE="$(find "$ARCHIVE_DIR" -maxdepth 1 -name 'HeatSeq_archive_*.7z' -type f | sort -r | head -n1)"
else
	ARCHIVE=""
fi

if [[ -z "$ARCHIVE" || ! -f "$ARCHIVE" ]]; then
	log_error "No .7z archive found in ${ARCHIVE_DIR}/"
	log_error "Usage: $0 [path/to/archive.7z]"
	exit 1
fi

THREADS="${THREADS:-12}"
log_info "Extracting: $ARCHIVE (threads=$THREADS)"
if ! 7z x "$ARCHIVE" -o"${SCRIPT_DIR}" -aoa -mmt="${THREADS}"; then
	log_error "Archive extraction failed"
	exit 1
fi
log_info "Extraction completed successfully"
