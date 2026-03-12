#!/usr/bin/env bash
set -euo pipefail

ARCHIVE="$(dirname "$(realpath "$0")")/HPC/HeatSeq_archive_20260312_063948.7z"

THREADS=12
7z x "$ARCHIVE" -o"$(dirname "$ARCHIVE")" -aoa -mmt="${THREADS}"
