#!/bin/bash
#===============================================================================
# SETUP: HTML RESULTS VIEWER — MINIMAL DEPENDENCIES
#===============================================================================
# Installs the minimal conda/mamba environment required to run the HTML viewer
# generator independently of the full 'gea' pipeline environment.
#
# The viewer generator (generate_html_viewer.py) uses only Python 3.11 stdlib,
# but this script also installs an optional lightweight HTTP server so the viewer
# can be served locally (avoids browser file:// CORS restrictions on images).
#
# Usage:
#   bash setup_html_viewer_env.sh              # Create env and install
#   bash setup_html_viewer_env.sh --serve      # Also start HTTP server after setup
#   bash setup_html_viewer_env.sh --generate   # Generate viewer then open browser
#   bash setup_html_viewer_env.sh --dry-run    # Show what would be done
#
# Environment: html_viewer  (lightweight, <200 MB)
#===============================================================================

set -euo pipefail

#===============================================================================
# CONFIGURATION
#===============================================================================

ENV_NAME="gea"
PYTHON_VERSION="3.11"
SERVE_PORT="${SERVE_PORT:-8080}"
DRY_RUN=false
DO_SERVE=false
DO_GENERATE=false

# Parameter expansion avoids nested $() subshell fork
SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)"

POST_PROC_DIR="$SCRIPT_DIR/3_POST_PROC"
VIEWER_SCRIPT="$SCRIPT_DIR/modules/c_post_processing/utilities/generate_html_viewer.py"
VIEWER_HTML="$POST_PROC_DIR/alignment_results_viewer.html"

#===============================================================================
# ARGUMENT PARSING
#===============================================================================

for arg in "$@"; do
    case "$arg" in
        --serve)    DO_SERVE=true ;;
        --generate) DO_GENERATE=true ;;
        --dry-run)  DRY_RUN=true ;;
        --help|-h)
            # Single awk pass replaces grep|head|sed 3-process pipe
            awk 'NR>20{exit} /^#/{sub(/^# ?/,""); print} !/^#/{exit}' "$0"
            exit 0
            ;;
        *) echo "WARNING: Unknown argument: $arg"; echo "Valid: --serve, --generate, --dry-run"; exit 1 ;;
    esac
done

#===============================================================================
# UTILITIES
#===============================================================================

log_info()  { echo "[INFO]  $*"; }
log_warn()  { echo "[WARN]  $*"; }
log_error() { echo "[ERROR] $*" >&2; }

run_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        echo "[DRY RUN] $*"
    else
        "$@"
    fi
}

check_command() { command -v "$1" &>/dev/null; }

#===============================================================================
# DETECT PACKAGE MANAGER
#===============================================================================

if check_command mamba; then
    PKG_MGR="mamba"
    log_info "Using mamba"
else
    PKG_MGR="conda"
    log_warn "mamba not found, using conda (slower)"
fi

#===============================================================================
# CONDA INIT
#===============================================================================

if [[ -z "${CONDA_EXE:-}" ]]; then
    for _conda_sh in \
        "$HOME/miniconda3/etc/profile.d/conda.sh" \
        "$HOME/anaconda3/etc/profile.d/conda.sh" \
        "/opt/conda/etc/profile.d/conda.sh"; do
        [[ -f "$_conda_sh" ]] && { source "$_conda_sh"; break; }
    done
fi

[[ -z "${CONDA_EXE:-}" ]] && { log_error "conda not found"; exit 1; }
eval "$(conda shell.bash hook)"

#===============================================================================
# PACKAGES
#===============================================================================

# Minimal packages: Python + optional HTTP server utilities
VIEWER_PACKAGES=(
    "python=${PYTHON_VERSION}"
    "uvicorn=0.34.0"   # lightweight ASGI server for local file serving
    "jq=1.7.1"         # optional: inspect/filter viewer_manifest.json
)

#===============================================================================
# CREATE OR UPDATE ENVIRONMENT
#===============================================================================

log_info "========================================"
log_info "HTML Viewer env: $ENV_NAME"
log_info "========================================"

if ${PKG_MGR} env list | grep -q "^${ENV_NAME} "; then
    log_info "Environment '${ENV_NAME}' already exists — updating..."
    run_cmd "${PKG_MGR}" install -n "${ENV_NAME}" -c conda-forge -y "${VIEWER_PACKAGES[@]}"
else
    log_info "Creating environment '${ENV_NAME}'..."
    run_cmd "${PKG_MGR}" create -n "${ENV_NAME}" -c conda-forge -y "${VIEWER_PACKAGES[@]}"
fi

#===============================================================================
# GENERATE VIEWER (optional)
#===============================================================================

if [[ "$DO_GENERATE" == true ]]; then
    log_info "Generating HTML viewer..."
    if [[ ! -d "$POST_PROC_DIR" ]]; then
        log_warn "3_POST_PROC/ not found at $POST_PROC_DIR — skipping generation"
    elif [[ ! -f "$VIEWER_SCRIPT" ]]; then
        log_error "Viewer script not found: $VIEWER_SCRIPT"
        exit 1
    else
        run_cmd conda run -n "${ENV_NAME}" python3 "$VIEWER_SCRIPT" "$POST_PROC_DIR"
        log_info "Viewer generated: $VIEWER_HTML"
    fi
fi

#===============================================================================
# SERVE (optional)
#===============================================================================

if [[ "$DO_SERVE" == true ]]; then
    if [[ ! -f "$VIEWER_HTML" ]]; then
        log_info "Viewer HTML not found — generating first..."
        conda run -n "${ENV_NAME}" python3 "$VIEWER_SCRIPT" "$POST_PROC_DIR" \
            || { log_error "Generation failed"; exit 1; }
    fi
    log_info "Serving results viewer at http://localhost:${SERVE_PORT}"
    log_info "Open:  http://localhost:${SERVE_PORT}/alignment_results_viewer.html"
    log_info "Stop:  Ctrl+C"
    run_cmd conda run -n "${ENV_NAME}" \
        python3 -m http.server "${SERVE_PORT}" --directory "$POST_PROC_DIR"
fi

#===============================================================================
# COMPLETION
#===============================================================================

log_info "========================================"
log_info "Setup complete: $ENV_NAME"
log_info "========================================"
log_info ""
log_info "Generate viewer:  python3 $VIEWER_SCRIPT $POST_PROC_DIR"
log_info "Serve locally:    bash setup_html_viewer_env.sh --serve"
log_info "Open viewer:      $VIEWER_HTML"
log_info ""
log_info "Or activate env and run manually:"
log_info "  conda activate $ENV_NAME"
log_info "  python3 $VIEWER_SCRIPT $POST_PROC_DIR"
