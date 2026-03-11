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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHANNELS="-c conda-forge -c bioconda"

UPDATE_MODE=true
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
    esac
done

#===============================================================================
# PACKAGE LISTS
#===============================================================================

# Preprocessing & data acquisition tools  (GEA pipeline)
PREPROCESSING_TOOLS=(
    "aria2"
    "parallel-fastq-dump"
    "sra-tools>=3.0"
    "entrez-direct"
    "kingfisher"
    "trim-galore"
    "trimmomatic"
    "cutadapt>=4.1"
    "fastqc"
    "multiqc"
    "parallel"
    "wget"    # ENA FTP fallback downloader (download_srrs_wget)
    "curl"    # ENA portal API queries (download_srrs_wget)
)

# Core alignment / quantification tools
ALIGNMENT_TOOLS=(
    "hisat2"
    "stringtie"
    "samtools"
    "salmon"
    "bowtie2"
    "rsem"
    "star"
    "trinity"
    "gffread"  # Transcript FASTA generation from genome+GTF (required to build inputs/fasta/ reference files for M3/M4/M5)
)

# R base and essentials
R_BASE=(
    "r-base=4.3"
    "r-essentials"
    "r-biocmanager"
)

# Bioconductor packages
BIOCONDUCTOR=(
    "bioconductor-deseq2"
    "bioconductor-complexheatmap"
    "bioconductor-tximport"
    "bioconductor-tximeta"
    "bioconductor-annotationdbi"
    "bioconductor-ballgown"
    "bioconductor-clusterprofiler"
    "bioconductor-enrichplot"
    "bioconductor-dose"
    "bioconductor-fgsea"
)

# CRAN / WGCNA packages
CRAN_PACKAGES=(
    "r-wgcna"
    "r-dynamictreecut"
    "r-fastcluster"
    "r-tidyverse"
    "r-dplyr"
    "r-tibble"
    "r-readr"
    "r-ggplot2"
    "r-ggrepel"
    "r-rcolorbrewer"
    "r-circlize"
    "r-pheatmap"
    "r-rtsne"
    "r-umap"
    "r-factoextra"
    "r-igraph"
    "r-reshape2"
    "r-getopt"
    "r-visnetwork"
    "r-networkd3"
    "r-htmlwidgets"
    "r-heatmaply"
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

eval "$(conda shell.bash hook)"

#===============================================================================
# HELPER: CHECK INSTALLED PACKAGES
#===============================================================================

# Returns a space-separated list of package names missing from the environment
check_packages_installed() {
    local installed_pkgs
    installed_pkgs=$(${PKG_MGR} list -n "${ENV_NAME}" --export 2>/dev/null \
        | cut -d'=' -f1 | sort -u)

    local missing=()
    for pkg in "${ALL_PACKAGES[@]}"; do
        local pkg_name="${pkg%%[><=]*}"   # strip version constraint for comparison
        if ! echo "$installed_pkgs" | grep -qxF "$pkg_name"; then
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
        run_cmd ${PKG_MGR} env remove -n "${ENV_NAME}" -y
        log_info "Recreating environment '${ENV_NAME}'..."
        run_cmd ${PKG_MGR} create -n "${ENV_NAME}" python="${PYTHON_VERSION}" -y

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
                ${PKG_MGR} update -n "${ENV_NAME}" ${CHANNELS} --all -y
            fi
            log_info "Update complete."
            exit 0
        fi
    fi
else
    log_info "Creating new environment '${ENV_NAME}'..."
    run_cmd ${PKG_MGR} create -n "${ENV_NAME}" python="${PYTHON_VERSION}" -y
fi

#===============================================================================
# INSTALL PACKAGES
#===============================================================================

log_info "Installing packages into '${ENV_NAME}'..."
if [[ "$DRY_RUN" == true ]]; then
    echo "[DRY RUN] Would execute: ${PKG_MGR} install -n ${ENV_NAME} ${CHANNELS} ${ALL_PACKAGES[@]} -y"
else
    ${PKG_MGR} install -n "${ENV_NAME}" ${CHANNELS} -y "${ALL_PACKAGES[@]}" || {
        log_warn "${PKG_MGR} installation failed, falling back to conda..."
        conda install -n "${ENV_NAME}" ${CHANNELS} -y "${ALL_PACKAGES[@]}"
    }
fi

#===============================================================================
# CONFIGURE SRA TOOLS
#===============================================================================

log_info "Configuring SRA tools..."
run_cmd ${PKG_MGR} run -n "${ENV_NAME}" vdb-config --prefetch-to-cwd

#===============================================================================
# ACTIVATE AND VERIFY
#===============================================================================

log_info "Activating environment..."
conda activate "${ENV_NAME}"

log_info "Verifying key executables..."
VERIFY_CMDS=("R" "Rscript" "samtools" "salmon" "hisat2" "STAR" "fastqc" "trim_galore" "prefetch" "fasterq-dump" "trimmomatic" "stringtie" "bowtie2" "rsem-calculate-expression" "gffread")
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
          "clusterProfiler", "ggplot2", "pheatmap", "igraph")
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
    "pheatmap", "igraph", "reshape2"
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

#===============================================================================
# COMPLETION
#===============================================================================

log_info "========================================"
log_info "Setup complete!"
log_info "========================================"
log_info ""
log_info "Environment : $ENV_NAME"
log_info "Activate    : conda activate $ENV_NAME"
log_info "Deactivate  : conda deactivate"
log_info ""
log_info "Run post-processing with:"
log_info "  cd $SCRIPT_DIR && bash run_all_post_processing.sh"
log_info ""
