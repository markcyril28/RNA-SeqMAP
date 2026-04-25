#!/bin/bash
#===============================================================================
# HTML RESULTS VIEWER — ORCHESTRATOR
#===============================================================================
# Generates (and optionally serves) an interactive HTML viewer for all heatmap
# results stored under II_RESULTS/3_POST_PROC/.
#
# The viewer provides:
#   - Grid layout: rows = alignment methods (M1–M5), columns = references
#   - Sidebar toggles: count type, normalization, processing level, gene group,
#     row orientation, sort order, analysis type
#   - Per-image download, full-screen modal, keyboard navigation
#   - viewer_manifest.json — pre-cached index of all figures
#
# Usage (from project root):
#   bash modules_gea/other_tools/run_post_processing_html_viewer.sh                    # Generate viewer (default)
#   bash modules_gea/other_tools/run_post_processing_html_viewer.sh --serve            # Generate + serve on localhost
#   bash modules_gea/other_tools/run_post_processing_html_viewer.sh --serve --port 9090
#   bash modules_gea/other_tools/run_post_processing_html_viewer.sh --manifest-only    # Only rebuild manifest JSON
#   bash modules_gea/other_tools/run_post_processing_html_viewer.sh --open             # Generate + open in browser
#   bash modules_gea/other_tools/run_post_processing_html_viewer.sh --output path.html # Custom output path
#   bash modules_gea/other_tools/run_post_processing_html_viewer.sh --dry-run          # Show what would be done
#
# Output:
#   II_RESULTS/3_POST_PROC/alignment_results_viewer.html
#   II_RESULTS/3_POST_PROC/viewer_manifest.json
#===============================================================================

set -o pipefail

# Resolve script directory: honour pre-set BASE_DIR from orchestrators (Nextflow/Snakemake)
if [[ -z "${BASE_DIR:-}" ]]; then
    SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
    [[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
    SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)" || { echo "[ERROR] run_post_processing_html_viewer.sh: Failed to resolve script directory" >&2; exit 1; }
    # Navigate up two levels: modules_gea/other_tools/ -> project root
    BASE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)" || { echo "[ERROR] run_post_processing_html_viewer.sh: Failed to resolve project root from $SCRIPT_DIR" >&2; exit 1; }
else
    if [[ "$BASE_DIR" != /* ]]; then
        BASE_DIR="$(cd "$BASE_DIR" 2>/dev/null && pwd)" || { echo "[ERROR] BASE_DIR is set but invalid: $BASE_DIR" >&2; exit 1; }
    fi
    SCRIPT_DIR="$BASE_DIR"
fi

#===============================================================================
# CONFIGURATION
#===============================================================================

POST_PROC_DIR="$BASE_DIR/II_RESULTS/3_POST_PROC/${CURRENT_GENE_GROUP:-_active}"
VIEWER_SCRIPT="$BASE_DIR/modules_gea/c_post_processing/utilities/generate_html_viewer.py"
OUTPUT_HTML=""          # Empty = default (II_RESULTS/3_POST_PROC/alignment_results_viewer.html)
SERVE_PORT="${SERVE_PORT:-8080}"
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
        --post-proc-dir)  shift; [[ $# -gt 0 ]] || { echo "[ERROR] --post-proc-dir requires a value" >&2; exit 1; }; POST_PROC_DIR="$1" ;;
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
if [[ "$POST_PROC_DIR" != /* ]]; then
    _orig_post_proc_dir="$POST_PROC_DIR"
    POST_PROC_DIR="$(cd "$POST_PROC_DIR" 2>/dev/null && pwd)" || { echo "[ERROR] --post-proc-dir path is invalid: $_orig_post_proc_dir" >&2; exit 1; }
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

LOG_DIR="$POST_PROC_DIR/logs/log_files"
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

log_step "HTML Results Viewer — Orchestrator"

if [[ ! -d "$POST_PROC_DIR" ]]; then
    log_error "II_RESULTS/3_POST_PROC directory not found: $POST_PROC_DIR"
    log_error "Run the post-processing pipeline first: bash run_post_processing.sh"
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

# Count available figures — POSIX-compatible: one filename per line, counted by wc -l.
# (Replaces GNU-only find -printf which is unavailable on macOS/BSD.)
_png_count=0
_png_count=$(find "$POST_PROC_DIR" -name "*.png" 2>/dev/null | wc -l)
log_info "Post-processing dir : $POST_PROC_DIR"
log_info "Figures found       : $_png_count PNG files"
if [[ "$_png_count" -eq 0 ]]; then
    log_warn "No PNG figures found — viewer will be empty."
    log_warn "Run the pipeline first: bash run_post_processing.sh"
fi

#===============================================================================
# BUILD GENERATOR ARGS
#===============================================================================

_gen_args=("$POST_PROC_DIR")
[[ -n "$OUTPUT_HTML" ]] && _gen_args+=("--output" "$OUTPUT_HTML")
[[ "$MANIFEST_ONLY" == true ]] && _gen_args+=("--manifest-only")

# Resolve output path for reporting
_html_out="${OUTPUT_HTML:-$POST_PROC_DIR/alignment_results_viewer.html}"

#===============================================================================
# GENERATE VIEWER
#===============================================================================

log_step "Generating Viewer"

if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY RUN] Would execute: python3 $VIEWER_SCRIPT ${_gen_args[*]}"
else
    _gen_log="$LOG_DIR/html_viewer_gen.log"
    python3 "$VIEWER_SCRIPT" "${_gen_args[@]}" >"$_gen_log" 2>&1 && {
        log_info "Manifest : $POST_PROC_DIR/viewer_manifest.json"
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

    log_step "Serving Results Viewer"
    log_info "Directory : $_serve_dir"
    log_info "URL       : $_serve_url"
    log_info "Stop      : Ctrl+C"
    log_info ""

    # Auto-open after small delay so the server is ready (reuses cached _BROWSER_CMD)
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
    log_info "Manifest updated: $POST_PROC_DIR/viewer_manifest.json"
else
    log_info "Open the viewer:"
    log_info "  File browser : $_html_out"
    log_info "  Local server : bash modules_gea/other_tools/run_post_processing_html_viewer.sh --serve"
fi
