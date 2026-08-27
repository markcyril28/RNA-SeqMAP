#!/bin/bash

#===============================================================================
# COMBINED CONDA/MAMBA ENVIRONMENT SETUP SCRIPT
#===============================================================================
# Creates/updates conda environment 'gea' with all dependencies for:
#   - RNA-seq alignment & preprocessing       (GEA pipeline)
#   - Post-processing & statistical analysis  (R / Bioconductor)
#
# Environment name : gea
# Lockfile         : setup_conda_gea.yml  (generated — do not edit by hand)
#
# Quick start:
#   1. Ensure Conda/Miniconda is installed
#   2. Run     : bash setup_conda_gea.sh
#   3. Activate: conda activate gea
#
# Run 'bash setup_conda_gea.sh --help' for the full flag list.
#===============================================================================

set -euo pipefail

#===============================================================================
# CONFIGURATION
#===============================================================================

ENV_NAME="gea"
PYTHON_VERSION="3.11"

# Parameter expansion avoids nested $(dirname) subshell fork
SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)" || { echo "[ERROR] setup_conda_gea.sh: Failed to resolve script directory" >&2; exit 1; }

# Single source of truth for generated/consumed paths
LOCKFILE="$SCRIPT_DIR/setup_conda_gea.yml"
R_INSTALLER="$SCRIPT_DIR/modules_gea/c_post_processing/utilities/install_R_packages.R"

# Array form so every call site expands identically without word-splitting
CHANNELS=(-c conda-forge -c bioconda)

UPDATE_MODE=false
ENV_RESTART_MODE=false
DRY_RUN=false           # If true, only show what would be done without installing
SKIP_UPDATE_CHECK=false # If true, skip slow update availability check

#===============================================================================
# USAGE
#===============================================================================

usage() {
    cat << EOF
Usage: bash setup_conda_gea.sh [OPTIONS]

Creates or updates the '${ENV_NAME}' conda environment for the GEA pipeline.

Options:
  --update              Update an existing environment
  --restart             Remove and recreate the environment
  --dry-run             Show what would be done without installing
  --skip-update-check   Skip the slow update availability check
  -h, --help            Show this help and exit

With no flags, the environment is created if missing, and packages are
(re)installed into it if it already exists.
EOF
}

#===============================================================================
# ARGUMENT PARSING
#===============================================================================

for arg in "$@"; do
    case "$arg" in
        --update)             UPDATE_MODE=true ;;
        --restart)            ENV_RESTART_MODE=true ;;
        --dry-run)            DRY_RUN=true ;;
        --skip-update-check)  SKIP_UPDATE_CHECK=true ;;
        -h|--help)            usage; exit 0 ;;
        *)  echo "ERROR: Unknown argument: $arg" >&2; usage >&2; exit 1 ;;
    esac
done

#===============================================================================
# PACKAGE LISTS
#===============================================================================
# These arrays are the source of truth for the environment's DIRECT dependencies.
# setup_conda_gea.yml is generated from the resolved environment — edit here, not there.

# Preprocessing & data acquisition tools  (GEA pipeline)
# IMPORTANT: Versions are pinned to tested versions for reproducibility.
# Update only after validation on a test dataset.
PREPROCESSING_TOOLS=(
    "aria2=1.37.0"
    "parallel-fastq-dump=0.6.7"
    "sra-tools=3.2.1"
    "entrez-direct=24.0"
    "kingfisher=0.4.1"
    "trim-galore=0.6.11"
    "trimmomatic=0.40"
    "cutadapt=5.2"
    "fastqc=0.12.1"
    "multiqc=1.33"
    "parallel=20260122"
    "wget=1.25.0"      # ENA FTP fallback downloader
    "curl=8.18.0"      # ENA portal API queries
    "dos2unix=7.5.4"   # Line-ending normalization
    "pigz=2.8"         # Multi-threaded gzip; preprocessing silently falls back to gzip if absent
)

# Core alignment / quantification tools
ALIGNMENT_TOOLS=(
    "hisat2=2.2.2"
    "stringtie=3.0.3"
    "samtools=1.22.1"
    "salmon=1.10.3"
    "bowtie2=2.5.5"
    "rsem=1.3.3"
    "star=2.7.11b"
    "trinity=2.15.2"
    "gffread=0.12.7"               # Transcript FASTA generation from genome+GTF
    "rseqc=5.0.4"                  # infer_experiment.py for strandness auto-detection
    "ucsc-gtftogenepred=482"       # GTF -> genePred conversion for BED12
    "ucsc-genepredtobed=482"       # genePred -> BED12 conversion for infer_experiment.py
)

# R base and essentials
R_BASE=(
    "r-base=4.3.3"
    "r-essentials=4.3"
    "r-biocmanager=1.30.26"
)

# Bioconductor packages (pinned to Bioconductor 3.18 / R 4.3)
BIOCONDUCTOR=(
    "bioconductor-deseq2=1.42.0"
    "bioconductor-complexheatmap=2.18.0"
    "bioconductor-tximport=1.30.0"
    "bioconductor-tximeta=1.20.1"
    "bioconductor-annotationdbi=1.64.1"
    "bioconductor-ballgown=2.34.0"
    "bioconductor-clusterprofiler=4.10.0"
    "bioconductor-enrichplot=1.22.0"
    "bioconductor-dose=3.28.1"
    "bioconductor-fgsea=1.28.0"
)

# CRAN / WGCNA packages
CRAN_PACKAGES=(
    "r-wgcna=1.73"
    "r-dynamictreecut=1.63_1"
    "r-fastcluster=1.3.0"
    "r-tidyverse=2.0.0"
    "r-dplyr=1.1.4"
    "r-tibble=3.3.0"
    "r-readr=2.1.5"
    "r-ggplot2=3.5.2"
    "r-ggrepel=0.9.6"
    "r-rcolorbrewer=1.1_3"
    "r-circlize=0.4.16"
    "r-pheatmap=1.0.13"
    "r-rtsne=0.17"
    "r-umap=0.2.10.0"
    "r-factoextra=1.0.7"
    "r-igraph=2.1.4"
    "r-reshape2=1.4.4"
    "r-getopt=1.20.4"
    "r-visnetwork=2.1.4"
    "r-networkd3=0.4.1"
    "r-htmlwidgets=1.6.4"
    "r-heatmaply=1.6.0"
    "r-corrplot=0.95"
    "r-dendextend=1.19.1"
    "r-gridextra=2.3"
    "r-scales=1.4.0"
)

# Combined package list — what actually gets installed
ALL_PACKAGES=(
    "${PREPROCESSING_TOOLS[@]}"
    "${ALIGNMENT_TOOLS[@]}"
    "${R_BASE[@]}"
    "${BIOCONDUCTOR[@]}"
    "${CRAN_PACKAGES[@]}"
)

#===============================================================================
# VERIFICATION LISTS
#===============================================================================

# Executables expected on PATH once the environment is active
VERIFY_CMDS=(
    R Rscript python3
    samtools salmon hisat2 STAR bowtie2 stringtie rsem-calculate-expression
    fastqc trim_galore trimmomatic
    prefetch fasterq-dump pigz
    gffread infer_experiment.py gtfToGenePred genePredToBed dos2unix
)

# R packages spot-checked after install.
# The authoritative, complete R package lists live in install_R_packages.R.
VERIFY_R_PKGS=(
    DESeq2 ComplexHeatmap WGCNA tximport clusterProfiler
    ggplot2 pheatmap igraph corrplot dendextend gridExtra scales
)

#===============================================================================
# UTILITY FUNCTIONS
#===============================================================================

log_info()  { echo "[INFO] $*"; }
log_warn()  { echo "[WARN] $*"; }
log_error() { echo "[ERROR] $*" >&2; }

# Banner used for the major section headings printed at runtime
log_section() {
    log_info "========================================"
    log_info "$*"
    log_info "========================================"
}

check_command() { command -v "$1" &> /dev/null; }

# Run a command, or print it in dry-run mode.
# NOTE: redirections are applied by the caller, so this cannot wrap a command
#       whose output is redirected into a file (see the lockfile export below).
run_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        echo "[DRY RUN] Would execute: $*"
    else
        "$@"
    fi
}

#===============================================================================
# DETECT PACKAGE MANAGER  (prefer mamba, fall back to conda)
#===============================================================================

if check_command mamba; then
    PKG_MGR="mamba"
    log_info "Using mamba for faster installation"
else
    PKG_MGR="conda"
    log_warn "Mamba not found, using conda. Install mamba for faster setup:"
    log_warn "  conda install -c conda-forge mamba"
fi

#===============================================================================
# CONDA INITIALIZATION
#===============================================================================

log_info "Initializing conda..."

if [[ -z "${CONDA_EXE:-}" ]]; then
    if   [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    elif [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
        source "/opt/conda/etc/profile.d/conda.sh"
    else
        log_error "Cannot find conda installation"
        exit 1
    fi
fi

if [[ "${WF_MANAGED_ENV:-}" != "true" ]]; then
    eval "$(conda shell.bash hook)"
fi

#===============================================================================
# HELPERS: ENVIRONMENT / PACKAGE INSPECTION
#===============================================================================

# True if the conda environment already exists.
# Matches on the name column rather than a line anchor: conda prints
# "gea  /path" while mamba 2.x prints an indented "  gea  Active  /path" table,
# so "^${ENV_NAME} " silently never matched under mamba.
# Kept as a single awk rather than `awk | grep -q`: grep -q exits on first match,
# which can SIGPIPE awk (exit 141); under `set -o pipefail` that made an existing
# environment look missing, and the resulting `create` then failed on "prefix
# already exists".
env_exists() {
    conda env list 2>/dev/null | awk -v want="${ENV_NAME}" '$1 == want { found = 1 } END { exit !found }'
}

# Returns a space-separated list of package names missing from the environment
check_packages_installed() {
    # Build associative array for O(1) lookup (avoids echo|grep fork per package)
    local -A _installed_set=()
    local _line
    while IFS= read -r _line; do
        [[ -n "$_line" ]] && _installed_set["${_line%%=*}"]=1
    done < <("${PKG_MGR}" list -n "${ENV_NAME}" --export 2>/dev/null)

    local missing=()
    for pkg in "${ALL_PACKAGES[@]}"; do
        local pkg_name="${pkg%%[><=]*}"   # strip version constraint for comparison
        if [[ -z "${_installed_set[$pkg_name]+x}" ]]; then
            missing+=("$pkg_name")
        fi
    done
    echo "${missing[*]:-}"
}

#===============================================================================
# ENVIRONMENT CREATION / UPDATE
#===============================================================================

log_section "Setting up environment: $ENV_NAME"
log_info "Total packages: ${#ALL_PACKAGES[@]}"

if env_exists; then
    log_info "Environment '${ENV_NAME}' exists."

    if [[ "$ENV_RESTART_MODE" == true ]]; then
        log_info "Removing existing environment '${ENV_NAME}'..."
        run_cmd "${PKG_MGR}" env remove -n "${ENV_NAME}" -y
        log_info "Recreating environment '${ENV_NAME}'..."
        run_cmd "${PKG_MGR}" create -n "${ENV_NAME}" "${CHANNELS[@]}" python="${PYTHON_VERSION}" -y

    elif [[ "$UPDATE_MODE" == true ]]; then
        MISSING_PKGS=$(check_packages_installed)
        if [[ -n "$MISSING_PKGS" ]]; then
            log_info "Missing packages detected: $MISSING_PKGS"
            log_info "Installing/updating all packages..."
        else
            if [[ "$SKIP_UPDATE_CHECK" == true ]]; then
                log_info "All packages installed. Skipping update check (--skip-update-check)."
                exit 0
            fi
            log_info "All packages installed. Running update..."
            run_cmd "${PKG_MGR}" update -n "${ENV_NAME}" "${CHANNELS[@]}" --all -y
            log_info "Update complete."
            exit 0
        fi
    fi
else
    log_info "Creating new environment '${ENV_NAME}'..."
    run_cmd "${PKG_MGR}" create -n "${ENV_NAME}" "${CHANNELS[@]}" python="${PYTHON_VERSION}" -y
fi

#===============================================================================
# INSTALL PACKAGES
#===============================================================================

log_info "Installing packages into '${ENV_NAME}'..."
if ! run_cmd "${PKG_MGR}" install -n "${ENV_NAME}" "${CHANNELS[@]}" -y "${ALL_PACKAGES[@]}"; then
    # Only worth retrying when the first attempt was mamba — re-running the
    # identical conda command would just repeat the same failed solve.
    if [[ "${PKG_MGR}" == "conda" ]]; then
        log_error "conda installation failed"
        exit 1
    fi
    log_warn "${PKG_MGR} installation failed, falling back to conda..."
    run_cmd conda install -n "${ENV_NAME}" "${CHANNELS[@]}" -y "${ALL_PACKAGES[@]}"
fi

#===============================================================================
# CONFIGURE SRA TOOLS
#===============================================================================

log_info "Configuring SRA tools..."
run_cmd conda run -n "${ENV_NAME}" vdb-config --prefetch-to-cwd || \
    log_warn "vdb-config failed — SRA prefetch-to-cwd not set. Pipeline will still work but prefetch may use default cache location."

#===============================================================================
# ACTIVATE AND VERIFY
#===============================================================================

if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY-RUN] Would activate environment '${ENV_NAME}' and verify tooling"
else
    log_info "Activating environment..."
    if ! conda activate "${ENV_NAME}" 2>/dev/null; then
        # Not fatal: the environment is already built. Failing hard here would
        # skip the R-package installer and the lockfile export below.
        log_warn "Could not activate '${ENV_NAME}' in this shell — skipping verification."
        log_warn "Verify manually with: conda activate ${ENV_NAME}"
    else
        log_info "Verifying key executables..."
        for cmd in "${VERIFY_CMDS[@]}"; do
            if check_command "$cmd"; then
                log_info "  ✓ $cmd"
            else
                log_warn "  ✗ $cmd not found"
            fi
        done

        log_info "Checking R packages..."
        if check_command Rscript; then
            Rscript -e '
for (pkg in commandArgs(trailingOnly = TRUE)) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(paste0("  ✓ ", pkg, "\n"))
    } else {
        cat(paste0("  ✗ ", pkg, " NOT FOUND\n"))
    }
}
' "${VERIFY_R_PKGS[@]}"
        else
            log_warn "Rscript not available — skipping R package check"
        fi
    fi
fi

#===============================================================================
# INSTALL ADDITIONAL R PACKAGES
#===============================================================================
# install_R_packages.R holds the authoritative CRAN + Bioconductor lists and
# already installs only what is missing, so it is the only fallback needed here.

if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY-RUN] Would install missing R packages via $R_INSTALLER"
elif [[ ! -f "$R_INSTALLER" ]]; then
    log_warn "Not found: $R_INSTALLER — skipping R package fallback install"
elif ! check_command Rscript; then
    log_warn "Rscript not on PATH — skipping R package fallback install"
    log_warn "Run manually after activating: Rscript --no-save \"$R_INSTALLER\""
else
    log_info "Installing any missing R packages via install_R_packages.R..."
    Rscript --no-save "$R_INSTALLER" || log_warn "Some R packages may have failed"
fi

#===============================================================================
# EXPORT LOCKFILE
#===============================================================================
# Cannot use run_cmd here: the redirection happens in this shell, so a dry run
# would still truncate the lockfile.

if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY-RUN] Would export environment lockfile to: $LOCKFILE"
else
    log_info "Exporting environment lockfile..."
    # Sibling temp file + atomic move: a failed export must not destroy the
    # previous, known-good lockfile by truncating it to the header.
    LOCKFILE_TMP="${LOCKFILE}.tmp"
    {
        echo "#==============================================================================="
        echo "# CONDA ENVIRONMENT LOCKFILE — '${ENV_NAME}'"
        echo "#==============================================================================="
        echo "# GENERATED FILE — do not edit by hand; setup_conda_gea.sh overwrites it."
        echo "#"
        echo "# Direct, pinned dependencies live in setup_conda_gea.sh (PACKAGE LISTS)."
        echo "# Everything below is the full resolved set, transitive deps included,"
        echo "# produced by 'conda env export --no-builds'. The machine-specific"
        echo "# 'prefix:' line is stripped so the file stays portable."
        echo "#"
        echo "# Regenerate  : bash setup_conda_gea.sh"
        echo "# Recreate env: conda env create -f setup_conda_gea.yml"
        echo "#==============================================================================="
    } > "$LOCKFILE_TMP"

    if conda env export -n "$ENV_NAME" --no-builds | grep -v '^prefix:' >> "$LOCKFILE_TMP"; then
        mv -f "$LOCKFILE_TMP" "$LOCKFILE"
        log_info "Lockfile saved to: $LOCKFILE"
    else
        rm -f "$LOCKFILE_TMP"
        log_warn "Failed to export environment lockfile — kept existing: $LOCKFILE"
    fi
fi

#===============================================================================
# COMPLETION
#===============================================================================

log_section "Setup complete!"
log_info ""
log_info "Environment : $ENV_NAME"
log_info "Lockfile    : $LOCKFILE"
log_info "Activate    : conda activate $ENV_NAME"
log_info "Deactivate  : conda deactivate"
log_info ""
log_info "Run post-processing with:"
log_info "  cd \"$SCRIPT_DIR\" && bash run_post_processing.sh"
log_info ""
log_info "HTML results viewer is auto-generated by run_post_processing.sh."
log_info "For standalone viewer setup, see: modules_gea/other_tools/setup_html_viewer_env.sh"
log_info "  bash modules_gea/other_tools/setup_html_viewer_env.sh --generate   # generate viewer only"
log_info "  bash modules_gea/other_tools/setup_html_viewer_env.sh --serve      # serve on localhost:8080"
log_info ""
