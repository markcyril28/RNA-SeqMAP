#!/bin/bash

#===============================================================================
# COMBINED CONDA/MAMBA ENVIRONMENT SETUP SCRIPT
#===============================================================================
# Creates/updates conda environment 'gea' with all dependencies for:
#   - RNA-seq alignment & preprocessing  (GEA pipeline)
#   - Post-processing & statistical analysis (R / Bioconductor)
#
# Environment Name: gea
#
# Usage:
#   ./setup_conda_gea.sh                      # Standard setup / install
#   ./setup_conda_gea.sh --update             # Update existing env
#   ./setup_conda_gea.sh --restart            # Remove and recreate env
#   ./setup_conda_gea.sh --dry-run            # Show what would be done
#   ./setup_conda_gea.sh --skip-update-check  # Skip update availability check
#
# Instructions:
#   1. Ensure Conda/Miniconda is installed
#   2. Run: bash setup_conda_gea.sh
#   3. Activate: conda activate gea
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
CHANNELS="-c conda-forge -c bioconda"

UPDATE_MODE=false
ENV_RESTART_MODE=false
DRY_RUN=false           # If true, only show what would be done without installing
SKIP_UPDATE_CHECK=false # If true, skip slow update availability check

#===============================================================================
# ARGUMENT PARSING
#===============================================================================

for arg in "$@"; do
    case "$arg" in
        --update)             UPDATE_MODE=true ;;
        --restart)            ENV_RESTART_MODE=true ;;
        --dry-run)            DRY_RUN=true ;;
        --skip-update-check)  SKIP_UPDATE_CHECK=true ;;
        *)  echo "WARNING: Unknown argument: $arg"; echo "Valid: --update, --restart, --dry-run, --skip-update-check"; exit 1 ;;
    esac
done

#===============================================================================
# PACKAGE LISTS
#===============================================================================

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
    "dos2unix"         # Line-ending normalization (version not pinned — not yet installed)
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
    "gffread=0.12.7"              # Transcript FASTA generation from genome+GTF
    "rseqc"                        # infer_experiment.py for strandness auto-detection (version not pinned — not yet installed)
    "ucsc-gtftogenepred"           # GTF -> genePred conversion for BED12 (version not pinned — not yet installed)
    "ucsc-genepredtobed"           # genePred -> BED12 conversion for infer_experiment.py (version not pinned — not yet installed)
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

# Combined package list
ALL_PACKAGES=(
    "${PREPROCESSING_TOOLS[@]}"
    "${ALIGNMENT_TOOLS[@]}"
    "${R_BASE[@]}"
    "${BIOCONDUCTOR[@]}"
    "${CRAN_PACKAGES[@]}"
)

#===============================================================================
# UTILITY FUNCTIONS
#===============================================================================

log_info()  { echo "[INFO] $*"; }
log_warn()  { echo "[WARN] $*"; }
log_error() { echo "[ERROR] $*" >&2; }

check_command() { command -v "$1" &> /dev/null; }

# Run a command, or print it in dry-run mode
run_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        echo "[DRY RUN] Would execute: $*"
    else
        "$@"
    fi
}

#===============================================================================
# DETECT PACKAGE MANAGER  (prefer mamba > micromamba > conda)
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
# HELPER: CHECK INSTALLED PACKAGES
#===============================================================================

# Returns a space-separated list of package names missing from the environment
check_packages_installed() {
    # Build associative array for O(1) lookup (avoids echo|grep fork per package)
    local -A _installed_set=()
    local _line
    while IFS= read -r _line; do
        [[ -n "$_line" ]] && _installed_set["${_line%%=*}"]=1
    done < <(${PKG_MGR} list -n "${ENV_NAME}" --export 2>/dev/null)

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

log_info "========================================"
log_info "Setting up environment: $ENV_NAME"
log_info "========================================"
log_info "Total packages: ${#ALL_PACKAGES[@]}"

if ${PKG_MGR} env list | grep -q "^${ENV_NAME} "; then
    log_info "Environment '${ENV_NAME}' exists."

    if [[ "$ENV_RESTART_MODE" == true ]]; then
        log_info "Removing existing environment '${ENV_NAME}'..."
        run_cmd "${PKG_MGR}" env remove -n "${ENV_NAME}" -y
        log_info "Recreating environment '${ENV_NAME}'..."
        run_cmd "${PKG_MGR}" create -n "${ENV_NAME}" python="${PYTHON_VERSION}" -y

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
            if [[ "$DRY_RUN" == true ]]; then
                echo "[DRY RUN] Would execute: ${PKG_MGR} update -n ${ENV_NAME} ${CHANNELS} --all -y"
            else
                "${PKG_MGR}" update -n "${ENV_NAME}" ${CHANNELS} --all -y
            fi
            log_info "Update complete."
            exit 0
        fi
    fi
else
    log_info "Creating new environment '${ENV_NAME}'..."
    run_cmd "${PKG_MGR}" create -n "${ENV_NAME}" python="${PYTHON_VERSION}" -y
fi

#===============================================================================
# INSTALL PACKAGES
#===============================================================================

log_info "Installing packages into '${ENV_NAME}'..."
if [[ "$DRY_RUN" == true ]]; then
    echo "[DRY RUN] Would execute: ${PKG_MGR} install -n ${ENV_NAME} ${CHANNELS} ${ALL_PACKAGES[*]} -y"
else
    "${PKG_MGR}" install -n "${ENV_NAME}" ${CHANNELS} -y "${ALL_PACKAGES[@]}" || {
        log_warn "${PKG_MGR} installation failed, falling back to conda..."
        conda install -n "${ENV_NAME}" ${CHANNELS} -y "${ALL_PACKAGES[@]}"
    }
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

log_info "Activating environment..."
if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY-RUN] Would activate environment '${ENV_NAME}'"
else
    conda activate "${ENV_NAME}"
fi

log_info "Verifying key executables..."
VERIFY_CMDS=("R" "Rscript" "python3" "samtools" "salmon" "hisat2" "STAR" "fastqc" "trim_galore" "prefetch" "fasterq-dump" "trimmomatic" "stringtie" "bowtie2" "rsem-calculate-expression" "gffread" "infer_experiment.py" "gtfToGenePred" "genePredToBed" "dos2unix")
for cmd in "${VERIFY_CMDS[@]}"; do
    if check_command "$cmd"; then
        log_info "  ✓ $cmd"
    else
        log_warn "  ✗ $cmd not found"
    fi
done

log_info "Checking R packages..."
Rscript -e '
pkgs <- c("DESeq2", "ComplexHeatmap", "WGCNA", "tximport",
          "clusterProfiler", "ggplot2", "pheatmap", "igraph",
          "corrplot", "dendextend", "gridExtra", "scales")
for (pkg in pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(paste0("  \u2713 ", pkg, "\n"))
    } else {
        cat(paste0("  \u2717 ", pkg, " NOT FOUND\n"))
    }
}
'

#===============================================================================
# INSTALL ADDITIONAL R PACKAGES (fallback via BiocManager)
#===============================================================================

if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY-RUN] Would install missing R packages via BiocManager"
else
log_info "Installing any missing R packages via BiocManager..."

if [[ -f "$SCRIPT_DIR/modules/c_post_processing/utilities/install_R_packages.R" ]]; then
    log_info "Running install_R_packages.R..."
    Rscript --no-save \
        "$SCRIPT_DIR/modules/c_post_processing/utilities/install_R_packages.R" \
        || log_warn "Some R packages may have failed"
fi

# Fallback: install critical packages directly if still missing
Rscript -e '
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

pkgs_to_check <- c(
    # Bioconductor packages (critical for pipeline)
    "tximport", "tximeta", "DESeq2", "ComplexHeatmap",
    "AnnotationDbi", "ballgown", "clusterProfiler",
    "enrichplot", "DOSE", "fgsea",
    # CRAN packages
    "WGCNA", "dynamicTreeCut", "fastcluster",
    "Rtsne", "umap", "factoextra", "ggrepel",
    "pheatmap", "igraph", "reshape2",
    "corrplot", "dendextend", "gridExtra", "scales"
)
for (pkg in pkgs_to_check) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(paste0("Installing ", pkg, "...\n"))
        tryCatch(
            BiocManager::install(pkg, ask = FALSE, update = FALSE),
            error = function(e)
                install.packages(pkg, repos = "https://cloud.r-project.org")
        )
    }
}
'
fi  # end DRY_RUN gate for R package installation

#===============================================================================
# COMPLETION
#===============================================================================

if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY-RUN] Would export environment lockfile"
else
    log_info "Exporting environment lockfile..."
    conda env export -n "$ENV_NAME" --no-builds > "$SCRIPT_DIR/environment.yml"
    log_info "Lockfile saved to: $SCRIPT_DIR/environment.yml"
fi

log_info "========================================"
log_info "Setup complete!"
log_info "========================================"
log_info ""
log_info "Environment : $ENV_NAME"
log_info "Activate    : conda activate $ENV_NAME"
log_info "Deactivate  : conda deactivate"
log_info ""
log_info "Run post-processing with:"
log_info "  cd $SCRIPT_DIR && bash run_post_processing.sh"
log_info ""
log_info "HTML results viewer is auto-generated by run_post_processing.sh."
log_info "For standalone viewer setup, see: modules/other_tools/setup_html_viewer_env.sh"
log_info "  bash modules/other_tools/setup_html_viewer_env.sh --generate   # generate viewer only"
log_info "  bash modules/other_tools/setup_html_viewer_env.sh --serve      # serve on localhost:8080"
log_info ""
