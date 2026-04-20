#!/bin/bash
#===============================================================================
# CONCORDANCE VIEWER — ORCHESTRATOR
#===============================================================================
# Generates (and optionally serves) an interactive HTML viewer for all
# concordance analysis results stored under 4_CONCORDANCE_ANALYSIS/.
#
# The viewer provides:
#   - Figure gallery:  card-based layout grouped by concordance mode
#   - Report viewer:   rendered Markdown concordance reports
#   - Data tables:     sortable CSV tables with correlation heatmap coloring
#   - Sidebar filters: concordance mode, analysis type, gene group
#   - Full-screen modal with keyboard navigation
#   - Zoom controls (buttons + Ctrl+scroll)
#
# Usage (from project root):
#   bash modules_gea/other_tools/run_concordance_viewer.sh                    # Generate viewer
#   bash modules_gea/other_tools/run_concordance_viewer.sh --serve            # Generate + serve
#   bash modules_gea/other_tools/run_concordance_viewer.sh --serve --port 9090
#   bash modules_gea/other_tools/run_concordance_viewer.sh --manifest-only    # Only rebuild manifest
#   bash modules_gea/other_tools/run_concordance_viewer.sh --open             # Generate + open
#   bash modules_gea/other_tools/run_concordance_viewer.sh --output path.html # Custom output path
#   bash modules_gea/other_tools/run_concordance_viewer.sh --dry-run          # Preview without changes
#
# Output:
#   4_CONCORDANCE_ANALYSIS/concordance_viewer.html
#   4_CONCORDANCE_ANALYSIS/concordance_manifest.json
#===============================================================================

set -o pipefail

# Resolve script directory: honour pre-set BASE_DIR from orchestrators (Nextflow/Snakemake)
if [[ -z "${BASE_DIR:-}" ]]; then
    SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
    [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
    SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)" || { echo "[ERROR] run_concordance_viewer.sh: Failed to resolve script directory" >&2; exit 1; }
    # Navigate up two levels: modules_gea/other_tools/ -> project root
    BASE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)" || { echo "[ERROR] run_concordance_viewer.sh: Failed to resolve project root from $SCRIPT_DIR" >&2; exit 1; }
else
    if [[ "$BASE_DIR" != /* ]]; then
        BASE_DIR="$(cd "$BASE_DIR" 2>/dev/null && pwd)" || { echo "[ERROR] BASE_DIR is set but invalid: $BASE_DIR" >&2; exit 1; }
    fi
    SCRIPT_DIR="$BASE_DIR"
fi

#===============================================================================
# CONFIGURATION
#===============================================================================

CONCORDANCE_DIR="$BASE_DIR/4_CONCORDANCE_ANALYSIS"
VIEWER_SCRIPT="$BASE_DIR/modules_gea/c_post_processing/utilities/generate_concordance_viewer.py"
OUTPUT_HTML=""          # Empty = default
SERVE_PORT="${SERVE_PORT:-8081}"
SERVE_HOST="${SERVE_HOST:-0.0.0.0}"

DO_SERVE=false
DO_OPEN=false
MANIFEST_ONLY=false
DRY_RUN=false

#===============================================================================
# ARGUMENT PARSING
#===============================================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --serve)          DO_SERVE=true ;;
        --open)           DO_OPEN=true ;;
        --manifest-only)  MANIFEST_ONLY=true ;;
        --dry-run)        DRY_RUN=true ;;
        --port)           shift; [[ $# -gt 0 ]] || { echo "[ERROR] --port requires a value" >&2; exit 1; }; SERVE_PORT="$1" ;;
        --output)         shift; [[ $# -gt 0 ]] || { echo "[ERROR] --output requires a value" >&2; exit 1; }; OUTPUT_HTML="$1" ;;
        --concordance-dir) shift; [[ $# -gt 0 ]] || { echo "[ERROR] --concordance-dir requires a value" >&2; exit 1; }; CONCORDANCE_DIR="$1" ;;
        --help|-h)
            sed -n '2,/^[^#]/{ /^#/!q; s/^# \?//p }' "$0"
            exit 0 ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            echo "Run 'bash $0 --help' for usage." >&2
            exit 1 ;;
    esac
    shift
done

# Resolve user-supplied paths to absolute (defensive: works from any cwd)
if [[ "$CONCORDANCE_DIR" != /* ]]; then
    _orig_concordance_dir="$CONCORDANCE_DIR"
    CONCORDANCE_DIR="$(cd "$CONCORDANCE_DIR" 2>/dev/null && pwd)" || { echo "[ERROR] --concordance-dir path is invalid: $_orig_concordance_dir" >&2; exit 1; }
fi

#===============================================================================
# LOGGING
#===============================================================================

source "$BASE_DIR/modules_gea/logging/logging_utils.sh" 2>/dev/null || {
    log_info()  { echo "[INFO]  $*"; }
    log_warn()  { echo "[WARN]  $*"; }
    log_error() { echo "[ERROR] $*" >&2; }
    log_step()  { echo ""; echo "==> $*"; echo ""; }
}

LOG_DIR="$CONCORDANCE_DIR/logs/log_files"
mkdir -p "$LOG_DIR" 2>/dev/null || true
export LOG_DIR

setup_logging "FALSE" 2>/dev/null || true

#===============================================================================
# CONDA ENVIRONMENT
#===============================================================================

# Skip conda activation when orchestrator manages the environment (Nextflow/Snakemake)
if [[ -z "${WF_MANAGED_ENV:-}" ]]; then
    if [[ "${CONDA_DEFAULT_ENV:-}" != "gea" ]] || ! command -v python3 &>/dev/null; then
        eval "$(conda shell.bash hook 2>/dev/null)" 2>/dev/null || true
        conda activate gea 2>/dev/null || log_warn "conda env 'gea' not found, using current env"
    fi
fi

#===============================================================================
# VALIDATE
#===============================================================================

log_step "Concordance Viewer — Orchestrator"

if [[ ! -d "$CONCORDANCE_DIR" ]]; then
    log_error "4_CONCORDANCE_ANALYSIS directory not found: $CONCORDANCE_DIR"
    log_error "Run the concordance pipeline first: bash run_concordance_analysis.sh"
    exit 1
fi

if [[ ! -f "$VIEWER_SCRIPT" ]]; then
    log_error "Viewer generator script not found: $VIEWER_SCRIPT"
    exit 1
fi

if ! command -v python3 &>/dev/null; then
    log_error "python3 is required but not found in PATH"
    exit 1
fi

# Cache browser command once — avoids repeated command -v PATH lookups in open/serve blocks
_BROWSER_CMD=""
if   command -v xdg-open     &>/dev/null; then _BROWSER_CMD="xdg-open"
elif command -v open         &>/dev/null; then _BROWSER_CMD="open"
elif command -v wslview      &>/dev/null; then _BROWSER_CMD="wslview"
elif command -v explorer.exe &>/dev/null; then _BROWSER_CMD="explorer.exe"
fi

# Count available outputs — single find pass for files (3 traversals → 1)
read -r _png_count _csv_count _md_count < <(
    find "$CONCORDANCE_DIR" \( -name "*.png" -o -path "*/tables/*.csv" -o -name "concordance_report.md" \) -print0 2>/dev/null \
    | awk -v RS='\0' '
        /\.png$/                       { png++ }
        /\/tables\/[^\/]+\.csv$/       { csv++ }
        /concordance_report\.md$/      { md++  }
        END { print png+0, csv+0, md+0 }
    ')
# Bash glob for top-level dirs — avoids find+wc fork pair (O(1) stat vs 2 forks)
_mode_count=0
for _d in "$CONCORDANCE_DIR"/*/; do
    [[ -d "$_d" ]] || continue
    case "${_d%/}" in */logs|*/_viewer_cache) continue ;; esac
    (( _mode_count++ ))
done

log_info "Concordance dir   : $CONCORDANCE_DIR"
log_info "Concordance modes : $_mode_count"
log_info "Figures found     : $_png_count PNG"
log_info "Tables found      : $_csv_count CSV"
log_info "Reports found     : $_md_count Markdown"

if [[ "$_png_count" -eq 0 && "$_csv_count" -eq 0 && "$_md_count" -eq 0 ]]; then
    log_warn "No concordance outputs found — viewer will be empty."
    log_warn "Run the concordance pipeline first: bash run_concordance_analysis.sh"
fi

#===============================================================================
# BUILD GENERATOR ARGS
#===============================================================================

_gen_args=("$CONCORDANCE_DIR")
[[ -n "$OUTPUT_HTML" ]] && _gen_args+=("--output" "$OUTPUT_HTML")
[[ "$MANIFEST_ONLY" == true ]] && _gen_args+=("--manifest-only")

_html_out="${OUTPUT_HTML:-$CONCORDANCE_DIR/concordance_viewer.html}"

#===============================================================================
# GENERATE VIEWER
#===============================================================================

log_step "Generating Concordance Viewer"

if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY RUN] Would execute: python3 $VIEWER_SCRIPT ${_gen_args[*]}"
else
    _gen_log="$LOG_DIR/concordance_viewer_gen.log"
    python3 "$VIEWER_SCRIPT" "${_gen_args[@]}" >"$_gen_log" 2>&1 && {
        log_info "Manifest : $CONCORDANCE_DIR/concordance_manifest.json"
        [[ "$MANIFEST_ONLY" == false ]] && log_info "Viewer   : $_html_out"
        log_info "Log      : $_gen_log"
    } || {
        log_error "Viewer generation failed — see $_gen_log"
        cat "$_gen_log" >&2
        exit 1
    }
fi

#===============================================================================
# OPEN IN BROWSER (optional)
#===============================================================================

if [[ "$DO_OPEN" == true && "$DRY_RUN" == false && -f "$_html_out" ]]; then
    log_info "Opening viewer in browser..."
    if [[ -n "$_BROWSER_CMD" ]]; then
        if [[ "$_BROWSER_CMD" == "explorer.exe" ]]; then
            explorer.exe "$(wslpath -w "$_html_out" 2>/dev/null || echo "$_html_out")" &
        else
            "$_BROWSER_CMD" "$_html_out" &
        fi
    else
        log_warn "Cannot auto-open browser — open manually: $_html_out"
    fi
fi

#===============================================================================
# SERVE (optional)
#===============================================================================

if [[ "$DO_SERVE" == true ]]; then
    if [[ ! -f "$_html_out" && "$DRY_RUN" == false ]]; then
        log_error "Viewer HTML not found at $_html_out — generation may have failed"
        exit 1
    fi

    # Serve from the directory containing the HTML file (handles custom --output paths)
    _serve_dir="$(cd "$(dirname "$_html_out")" && pwd)"
    _serve_url="http://localhost:${SERVE_PORT}/${_html_out##*/}"

    log_step "Serving Concordance Viewer"
    log_info "Directory : $_serve_dir"
    log_info "URL       : $_serve_url"
    log_info "Stop      : Ctrl+C"
    log_info ""

    # Auto-open after small delay (reuses cached _BROWSER_CMD — no repeated PATH lookups)
    if [[ "$DO_OPEN" == false && "$DRY_RUN" == false && -n "$_BROWSER_CMD" ]]; then
        (sleep 1; "$_BROWSER_CMD" "$_serve_url") &>/dev/null &
    fi

    if [[ "$DRY_RUN" == true ]]; then
        log_info "[DRY RUN] Would execute: python3 -m http.server $SERVE_PORT --bind $SERVE_HOST --directory $_serve_dir"
    else
        python3 -m http.server "$SERVE_PORT" --bind "$SERVE_HOST" --directory "$_serve_dir"
    fi
fi

#===============================================================================
# SUMMARY
#===============================================================================

log_step "Done"
if [[ "$MANIFEST_ONLY" == true ]]; then
    log_info "Manifest updated: $CONCORDANCE_DIR/concordance_manifest.json"
else
    log_info "Open the viewer:"
    log_info "  File browser : $_html_out"
    log_info "  Local server : bash modules_gea/other_tools/run_concordance_viewer.sh --serve"
fi
