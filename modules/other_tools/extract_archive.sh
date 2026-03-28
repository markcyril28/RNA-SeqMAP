#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)" || { echo "[ERROR] extract_archive.sh: Failed to resolve script directory" >&2; exit 1; }
# Resolve project root (this script lives in modules/other_tools/)
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)" || { echo "[ERROR] extract_archive.sh: Failed to resolve PROJECT_ROOT" >&2; exit 1; }
ARCHIVE_DIR="${PROJECT_ROOT}/HPC"

# Source logging utilities for consistent output
source "${PROJECT_ROOT}/modules/logging/logging_utils.sh" 2>/dev/null || {
    log_info()  { echo "[INFO]  $*"; }
    log_warn()  { echo "[WARN]  $*"; }
    log_error() { echo "[ERROR] $*" >&2; }
    log_step()  { echo ""; echo "==> $*"; }
}

# Find the latest .7z archive (avoid hardcoding timestamped filenames)
if [[ -n "${1:-}" && -f "$1" ]]; then
	ARCHIVE="$1"
elif [[ -d "$ARCHIVE_DIR" ]]; then
	# Bash glob + loop replaces find|sort|head pipeline (3 subshells → 0)
	# O(N) single pass over glob results; lexicographic sort by shell is sufficient
	ARCHIVE=""
	for _f in "$ARCHIVE_DIR"/HeatSeq_archive_*.7z; do
		[[ -f "$_f" ]] && [[ "$_f" > "${ARCHIVE:-}" ]] && ARCHIVE="$_f"
	done
else
	ARCHIVE=""
fi

if [[ -z "$ARCHIVE" || ! -f "$ARCHIVE" ]]; then
	log_error "No .7z archive found in ${ARCHIVE_DIR}/"
	log_error "Usage: $0 [path/to/archive.7z]"
	exit 1
fi

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-${PBS_NCPUS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 12)}}}"
log_info "Extracting: $ARCHIVE (threads=$THREADS)"
if ! 7z x "$ARCHIVE" -o"${PROJECT_ROOT}" -aoa -mmt="${THREADS}"; then
	log_error "Archive extraction failed"
	exit 1
fi
log_info "Extraction completed successfully"
