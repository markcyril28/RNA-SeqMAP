#!/bin/bash

# ===============================================
# GENERIC CONFIGURATION - SHARED ACROSS METHODS
# ===============================================
# This file can be sourced by method-specific configs
# to provide common default values
# ===============================================

# Base directories (will be overridden by method-specific configs)
BASE_DIR="${BASE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
SCRIPT_DIR="${SCRIPT_DIR:-$BASE_DIR}"

# Master reference name (override via environment or method config)
MASTER_REFERENCE="${MASTER_REFERENCE:-All_Smel_Genes}"

# Common directory structure patterns
GENE_GROUPS_DIR="${GENE_GROUPS_DIR:-$BASE_DIR/0_INPUTs/gene_groups}"
ANALYSIS_MODULES_DIR="${ANALYSIS_MODULES_DIR:-$BASE_DIR/modules/c_post_processing/analysis_modules}"
UTILITIES_DIR="${UTILITIES_DIR:-$BASE_DIR/modules/c_post_processing/utilities}"
PREPROCESSING_DIR="${PREPROCESSING_DIR:-$BASE_DIR/modules/c_post_processing/preprocessing}"

# Output directory naming convention
FIGURES_OUTPUT_BASE="${FIGURES_OUTPUT_BASE:-Figure_Outputs}"
HEATMAP_SUBDIR="${HEATMAP_SUBDIR:-I_Basic_Heatmap}"
CV_HEATMAP_SUBDIR="${CV_HEATMAP_SUBDIR:-II_Heatmap_with_CV}"
BAR_GRAPH_SUBDIR="${BAR_GRAPH_SUBDIR:-III_Bar_Graphs}"
WGCNA_SUBDIR="${WGCNA_SUBDIR:-III_Coexpression_WGCNA}"
DE_SUBDIR="${DE_SUBDIR:-V_Differential_Expression}"
GSEA_SUBDIR="${GSEA_SUBDIR:-VI_Gene_Set_Enrichment}"
PCA_SUBDIR="${PCA_SUBDIR:-VII_PCA}"
CLUSTER_SUBDIR="${CLUSTER_SUBDIR:-VIII_Sample_Clustering}"
SPECIFICITY_SUBDIR="${SPECIFICITY_SUBDIR:-IX_Tissue_Specificity}"

# Processing parameters
THREADS="${THREADS:-12}"
ENABLE_GPU="${ENABLE_GPU:-FALSE}"

# GPU/CUDA configuration
# If conda environment has CUDA configured, source the activation script
if [[ -n "${CONDA_PREFIX:-}" && -f "$CONDA_PREFIX/etc/conda/activate.d/cuda_env.sh" ]]; then
    source "$CONDA_PREFIX/etc/conda/activate.d/cuda_env.sh"
fi

# Export CUDA-related variables for R torch compatibility
export TORCH_CUDA_VERSION="${TORCH_CUDA_VERSION:-}"
export CUDA_HOME="${CUDA_HOME:-$CONDA_PREFIX}"
export CUDA_PATH="${CUDA_PATH:-$CONDA_PREFIX}"

# Export for R scripts
export MASTER_REFERENCE THREADS ENABLE_GPU
export ANALYSIS_MODULES_DIR GENE_GROUPS_DIR UTILITIES_DIR
