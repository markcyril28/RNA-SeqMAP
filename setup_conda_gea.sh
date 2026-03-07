#!/bin/bash

# Conda Environment Setup for GEA Pipeline
set -euo pipefail

# ===================== SETTINGS =====================
#ENV_NAME="GEA_ENV"
ENV_NAME="gea"
UPDATE_MODE=true
ENV_RESTART_MODE=false
DRY_RUN=false          # If true, only show what would be done without installing
SKIP_UPDATE_CHECK=false # If true, skip slow update availability check

# Channels
CHANNELS="-c conda-forge -c bioconda"

# Define required packages (use arrays for reliability)
# ====================================================
PACKAGES=(
    aria2
    parallel-fastq-dump
    sra-tools
    entrez-direct
    kingfisher
    hisat2
    stringtie
    samtools
    bowtie2
    rsem
    salmon
    star
    trim-galore
    trimmomatic
    cutadapt
    fastqc
    multiqc
    parallel
    r-wgcna
    r-dynamictreecut
    r-fastcluster
)

# Version constraints (only for packages that need them)
declare -A VERSION_CONSTRAINTS=(
    ["sra-tools"]=">=3.0"
    ["cutadapt"]=">=4.1"
)

# Parse arguments
for arg in "$@"; do
    case "$arg" in
        --update) UPDATE_MODE=true ;;
        --restart) ENV_RESTART_MODE=true ;;
        --dry-run) DRY_RUN=true ;;
        --skip-update-check) SKIP_UPDATE_CHECK=true ;;
    esac
done

# Helper function for running commands (respects DRY_RUN)
run_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        echo "[DRY RUN] Would execute: $*"
    else
        "$@"
    fi
}

# Check which package manager is available (prefer mamba > micromamba > conda)
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "Using mamba for faster installation"
elif command -v micromamba &> /dev/null; then
    CONDA_CMD="micromamba"
    echo "Using micromamba for faster installation"
else
    CONDA_CMD="conda"
    echo "Mamba/micromamba not found, using conda"
fi

# Function to check if all packages are installed (optimized: single list call)
check_packages_installed() {
    local installed_pkgs
    installed_pkgs=$(${CONDA_CMD} list -n "${ENV_NAME}" --export 2>/dev/null | cut -d'=' -f1 | sort -u)
    
    local missing=()
    for pkg in "${PACKAGES[@]}"; do
        if ! echo "$installed_pkgs" | grep -qx "$pkg"; then
            missing+=("$pkg")
        fi
    done
    echo "${missing[*]:-}"
}

# Function to build package spec with versions
build_package_spec() {
    local specs=()
    for pkg in "${PACKAGES[@]}"; do
        if [[ -v VERSION_CONSTRAINTS[$pkg] ]]; then
            specs+=("\"${pkg}${VERSION_CONSTRAINTS[$pkg]}\"")
        else
            specs+=("$pkg")
        fi
    done
    echo "${specs[*]}"
}

# Create or activate environment
if ${CONDA_CMD} env list | grep -q "^${ENV_NAME} "; then
    echo "Environment '${ENV_NAME}' exists"
    if [[ "$ENV_RESTART_MODE" == true ]]; then
        echo "Removing existing environment '${ENV_NAME}'..."
        run_cmd ${CONDA_CMD} env remove -n "${ENV_NAME}" -y
        echo "Recreating environment '${ENV_NAME}'..."
        run_cmd ${CONDA_CMD} create -n "${ENV_NAME}" -y
    elif [[ "$UPDATE_MODE" == true ]]; then
        # Check if all packages are installed first
        MISSING_PKGS=$(check_packages_installed)
        if [[ -n "$MISSING_PKGS" ]]; then
            echo "Missing packages detected: $MISSING_PKGS"
            echo "Installing missing packages..."
        else
            if [[ "$SKIP_UPDATE_CHECK" == true ]]; then
                echo "All packages installed. Skipping update check (--skip-update-check)."
                exit 0
            fi
            # All packages installed, run update (mamba/conda will skip if nothing to update)
            echo "All packages installed. Running update..."
            if [[ "$DRY_RUN" == true ]]; then
                echo "[DRY RUN] Would execute: ${CONDA_CMD} update -n ${ENV_NAME} ${CHANNELS} --all -y"
            else
                ${CONDA_CMD} update -n "${ENV_NAME}" ${CHANNELS} --all -y
            fi
            echo "Update complete"
            exit 0
        fi
    fi
else
    echo "Creating environment '${ENV_NAME}'..."
    run_cmd ${CONDA_CMD} create -n "${ENV_NAME}" -y
fi

PACKAGE_SPEC=$(build_package_spec)

echo "Installing packages into '${ENV_NAME}'..."
if [[ "$DRY_RUN" == true ]]; then
    echo "[DRY RUN] Would execute: ${CONDA_CMD} install -n ${ENV_NAME} ${CHANNELS} ${PACKAGE_SPEC} -y"
else
    eval "${CONDA_CMD} install -n ${ENV_NAME} ${CHANNELS} ${PACKAGE_SPEC} -y" || {
        echo "${CONDA_CMD} installation failed, falling back to conda..."
        eval "conda install -n ${ENV_NAME} ${CHANNELS} ${PACKAGE_SPEC} -y"
    }
fi

# Configure SRA tools
echo "Configuring SRA tools..."
run_cmd ${CONDA_CMD} run -n "${ENV_NAME}" vdb-config --prefetch-to-cwd

echo "Installation complete. Activate with: conda activate ${ENV_NAME}"
