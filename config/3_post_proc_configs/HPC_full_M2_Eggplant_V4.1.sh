#!/bin/bash

# ==============================================================================
# HPC FULL CONFIG: M2 — HISAT2 De Novo (Eggplant_V4.1_transcripts.function)
# ==============================================================================
# Purpose: Production post-processing for M2 (HISAT2 de novo transcript
#          assembly via StringTie).  Uses the full eggplant tissue atlas.
#
# Method:
#   M2 — HISAT2 De Novo  →  uses Eggplant_V4.1_transcripts.function
#
# Preprocessing:  stringtie_matrix_builder.sh runs automatically.
#                 Enable "Stringtie_Matrix" only when matrices are NOT yet built.
# DE support:     NO  — StringTie normalized expression (not integer counts).
# NOT applicable: Tximport_STAR, Tximport_Salmon, Tximport_RSEM, Matrix_Creation
#
# Usage:
#   Uncomment this config in PIPELINE_CONFIGS inside run_all_post_processing.sh
# ==============================================================================

# ==============================================================================
# MASTER REFERENCE
# ==============================================================================
#
# M2 uses transcript-level references.

MASTER_REFERENCES=(
    "Eggplant_V4.1_transcripts.function"
)

# ==============================================================================
# GENE GROUPS
# ==============================================================================

GENE_GROUPS=(
    "SmelDMPs_with_1_18s_rRNA"
    #"SmelDMPs"
    #"SmelDMPs_with_2_18s_rRNA"
    #"SmelDMPs_with_SmelCyclo"
    #"SmelGRF-GIF"              # Needs stringtie_matrix_builder.sh first
    #"SmelGRF-GIF_with_1_18s_rRNA"
    #"SmelGRFs"
    #"SmelGIF"
    #"Selected_GRF_GIF_Genes_vAll_GIF"
    "Selected_GRF_GIF_Genes_vTwo_GIF"
)

# ==============================================================================
# SRR DATASETS
# ==============================================================================

SRR_DATASETS=(
    #"PRJNA328564"             # Main Dataset — Eggplant tissue atlas (full)
    "PRJNA328564_selected"
    #"PRJNA865018"            # SET_1: Good Dataset for SmelDMP GEA
    #"PRJNA941250"            # SET_2: Good Dataset for SmelDMP GEA
    #"PRJNA865018_and_PRJNA941250"
    #"SAMN28540077"           # Chinese Dataset 1 — replicability
    #"SAMN28540068"           # Chinese Dataset 2 — replicability
    #"OTHER_SRR_LIST"
)

# ==============================================================================
# ALIGNMENT METHODS
# ==============================================================================

METHODS=(
    "M2_HISAT2_DeNovo"
)

# ==============================================================================
# ANALYSES
# ==============================================================================
#
# Enable "Stringtie_Matrix" only on the first run (when matrices are not yet built).
# Tximport_* and Matrix_Creation are NOT applicable to M2.
# Differential_Expression is NOT supported (StringTie TPM, not integer counts).
# Tissue_Specificity is applicable to transcript-level expression values.

ANALYSES=(
    # ---- M2 preprocessing (enable on first run only) ----
    #"Stringtie_Matrix"          # M2 only — skip if matrices already built

    # ---- Visualisation ----
    "Basic_Heatmap"
    "Heatmap_with_CV"
)
