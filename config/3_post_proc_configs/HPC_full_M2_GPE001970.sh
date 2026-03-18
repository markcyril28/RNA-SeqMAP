#!/bin/bash

# ==============================================================================
# HPC FULL CONFIG: M2 — HISAT2 De Novo (GPE001970_transcripts)
# ==============================================================================
# Purpose: Production post-processing for M2 (HISAT2 de novo transcript
#          assembly via StringTie).  Uses the full eggplant tissue atlas.
#
# Method:
#   M2 — HISAT2 De Novo  →  uses GPE001970_transcripts
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
    "GPE001970_transcripts"
)

# ==============================================================================
# GENE GROUPS
# ==============================================================================

GENE_GROUPS=(
    #"SmelDMPs_v5"
    "SmelDMPs_v5_with_18s_and_HAP2"
    #"SmelGRF-GIF_with_Control"
    "Selected_SmelGRF-GIF_with_two_GIF"
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
# Stringtie_Matrix preprocessing runs automatically via run_method_preprocessing.
# Tximport_* and Matrix_Creation are NOT applicable to M2.
# Differential_Expression is NOT supported (StringTie TPM, not integer counts).
# Tissue_Specificity is applicable to transcript-level expression values.

ANALYSES=(
    # ---- M2 preprocessing (runs automatically via run_method_preprocessing) ----
    #"Stringtie_Matrix"          # Legacy name — preprocessing handles this automatically

    # ---- Visualisation ----
    "Basic_Heatmap"
    "Heatmap_with_CV"
)
