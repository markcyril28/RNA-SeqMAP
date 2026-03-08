#!/bin/bash

#===============================================================================
# COMBINED CONDA/MAMBA ENVIRONMENT SETUP SCRIPT
#===============================================================================
# Creates conda environment 'gea' with all dependencies for RNA-seq analysis
# Uses mamba if available, falls back to conda
#
# Environment Name: gea
#
# Usage:
#   ./setup_conda_post_proc.sh           # Standard setup
#   ENABLE_GPU=TRUE ./setup_conda_post_proc.sh  # With R torch GPU support
#
# GPU Acceleration:
#   R torch uses system CUDA (e.g., /usr/local/cuda) for GPU operations.
#   After setup, run: ./prep_GPU_post-proc.sh --install
#
# Instructions:
#   1. Ensure Conda/Miniconda is installed
#   2. Run: bash setup_conda_post_proc.sh
#   3. Activate: conda activate gea
#   4. For GPU: ./prep_GPU_post-proc.sh --install
#===============================================================================

set -euo pipefail

#===============================================================================
# CONFIGURATION
#===============================================================================

ENV_NAME="gea"
PYTHON_VERSION="3.11"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#===============================================================================
# PACKAGE LISTS
#===============================================================================

# Core alignment/quantification tools
ALIGNMENT_TOOLS=(
    "hisat2"
    "stringtie"
    "samtools"
    "salmon"
    "bowtie2"
    "rsem"
    "star"
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
    "bioconductor-org.mm.eg.db"
)

# WGCNA and dependencies (CRAN package, not Bioconductor)
WGCNA_PACKAGES=(
    "r-wgcna"
)

# CRAN packages
CRAN_PACKAGES=(
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
    "r-dynamictreecut"
    "r-fastcluster"
    "r-visnetwork"
    "r-networkd3"
    "r-htmlwidgets"
    "r-heatmaply"
)

# GPU support packages (optional, enable with ENABLE_GPU=TRUE)
# R torch downloads its own libtorch binaries - these are optional dependencies
# For full GPU setup, run: ./prep_GPU_post-proc.sh --install
GPU_PACKAGES=(
    # No conda packages needed - R torch uses system CUDA
    # The prep_GPU_post-proc.sh script handles R torch GPU installation
)

#===============================================================================
# UTILITY FUNCTIONS
#===============================================================================

log_info() { echo "[INFO] $*"; }
log_warn() { echo "[WARN] $*"; }
log_error() { echo "[ERROR] $*" >&2; }

check_command() {
    command -v "$1" &> /dev/null
}

#===============================================================================
# DETECT PACKAGE MANAGER
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

# Find conda installation
if [[ -z "${CONDA_EXE:-}" ]]; then
    if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
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
# ENVIRONMENT CREATION/UPDATE
#===============================================================================

ALL_PACKAGES=(
    "${ALIGNMENT_TOOLS[@]}"
    "${R_BASE[@]}"
    "${BIOCONDUCTOR[@]}"
    "${WGCNA_PACKAGES[@]}"
    "${CRAN_PACKAGES[@]}"
)

# Add GPU packages if requested (currently empty - R torch uses system CUDA)
if [[ "${ENABLE_GPU:-FALSE}" == "TRUE" ]]; then
    # GPU_PACKAGES is empty - R torch installs its own libtorch
    # Just log that GPU mode is enabled
    GPU_CHANNELS=""
    log_info "GPU mode enabled - R torch will use system CUDA"
    log_info "Run ./prep_GPU_post-proc.sh --install after setup for R torch GPU"
else
    GPU_CHANNELS=""
fi

log_info "========================================"
log_info "Setting up environment: $ENV_NAME"
log_info "========================================"
log_info "Total packages: ${#ALL_PACKAGES[@]}"

if conda env list | grep -q "^${ENV_NAME}\s"; then
    log_info "Environment '$ENV_NAME' exists. Updating..."
    $PKG_MGR install -n "$ENV_NAME" -c conda-forge -c bioconda $GPU_CHANNELS --yes "${ALL_PACKAGES[@]}"
else
    log_info "Creating new environment '$ENV_NAME'..."
    $PKG_MGR create -n "$ENV_NAME" python=$PYTHON_VERSION --yes
    
    log_info "Installing packages..."
    $PKG_MGR install -n "$ENV_NAME" -c conda-forge -c bioconda $GPU_CHANNELS --yes "${ALL_PACKAGES[@]}"
fi

#===============================================================================
# ACTIVATE AND VERIFY
#===============================================================================

log_info "Activating environment..."
conda activate "$ENV_NAME"

log_info "Verifying installation..."

# Check key packages
VERIFY_PACKAGES=("R" "Rscript" "samtools" "salmon")
for pkg in "${VERIFY_PACKAGES[@]}"; do
    if check_command "$pkg"; then
        log_info "  ✓ $pkg"
    else
        log_warn "  ✗ $pkg not found"
    fi
done

# Check R packages
log_info "Checking R packages..."
Rscript -e '
pkgs <- c("DESeq2", "ComplexHeatmap", "WGCNA", "tximport", "clusterProfiler")
for (pkg in pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(paste0("  ✓ ", pkg, "\n"))
    } else {
        cat(paste0("  ✗ ", pkg, " NOT FOUND\n"))
    }
}
'

#===============================================================================
# INSTALL ADDITIONAL R PACKAGES (if needed)
#===============================================================================

log_info "Installing any missing R packages via BiocManager..."

# Run external R package installer if available
if [[ -f "$SCRIPT_DIR/modules/utilities/install_R_packages.R" ]]; then
    log_info "Running install_R_packages.R..."
    if [[ "${ENABLE_GPU:-FALSE}" == "TRUE" ]]; then
        Rscript --no-save "$SCRIPT_DIR/modules/utilities/install_R_packages.R" --with-gpu || log_warn "Some R packages may have failed"
    else
        Rscript --no-save "$SCRIPT_DIR/modules/utilities/install_R_packages.R" || log_warn "Some R packages may have failed"
    fi
fi

# Fallback: install critical packages directly
Rscript -e '
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Packages that might need manual installation
pkgs_to_check <- c(
    "WGCNA",
    "dynamicTreeCut", 
    "fastcluster",
    "Rtsne",
    "umap",
    "factoextra",
    "ggrepel",
    "pheatmap",
    "igraph",
    "reshape2"
)

for (pkg in pkgs_to_check) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(paste0("Installing ", pkg, "...\n"))
        tryCatch({
            BiocManager::install(pkg, ask = FALSE, update = FALSE)
        }, error = function(e) {
            install.packages(pkg, repos = "https://cloud.r-project.org")
        })
    }
}
'

#===============================================================================
# COMPLETION
#===============================================================================

# Verify GPU setup if enabled
if [[ "${ENABLE_GPU:-FALSE}" == "TRUE" ]]; then
    log_info "Checking GPU availability..."
    if command -v nvidia-smi &>/dev/null; then
        log_info "  ✓ nvidia-smi available"
        nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null || true
    else
        log_warn "  ✗ nvidia-smi not found - GPU may not be available"
    fi
    
    log_info ""
    log_info "To install R torch with GPU support, run:"
    log_info "  ./prep_GPU_post-proc.sh --install"
fi

log_info "========================================"
log_info "Setup complete!"
log_info "========================================"
log_info ""
log_info "Environment: $ENV_NAME"
log_info "Activate with:"
log_info "  conda activate $ENV_NAME"
log_info ""
log_info "Run post-processing with:"
log_info "  cd $SCRIPT_DIR && bash run_all_post_processing.sh"
log_info ""
log_info "Deactivate with:"
log_info "  conda deactivate"
log_info ""
if [[ "${ENABLE_GPU:-FALSE}" == "TRUE" ]]; then
    log_info "GPU mode enabled."
    log_info "To complete R torch GPU setup, run:"
    log_info "  ./prep_GPU_post-proc.sh --install"
    log_info ""
    log_info "R torch GPU acceleration enables:"
    log_info "  - Correlation matrix computation (WGCNA)"
    log_info "  - PCA/dimensionality reduction"
    log_info "  - Large matrix operations"
else
    log_info "To enable GPU acceleration, run:"
    log_info "  ENABLE_GPU=TRUE ./setup_conda_post_proc.sh"
    log_info "  ./prep_GPU_post-proc.sh --install"
fi
log_info ""
