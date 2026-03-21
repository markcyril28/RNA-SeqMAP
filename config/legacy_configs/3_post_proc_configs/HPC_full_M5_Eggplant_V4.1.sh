#!/bin/bash

# ==============================================================================
# HPC FULL CONFIG: M5 — RSEM + Bowtie2 (Eggplant_V4.1_transcripts.function)
# ==============================================================================
# Purpose: Production post-processing for M5 (RSEM probabilistic quantification
#          via Bowtie2).  Uses the full eggplant tissue atlas.
#
# Method:
#   M5 — RSEM + Bowtie2  →  uses Eggplant_V4.1_transcripts.function
#
# Preprocessing:  Tximport_RSEM → Matrix_Creation (Tximport is a no-op when
#                 Matrix_Creation is active; keep both for flexibility).
# DE support:     Possible via tximport counts, but not enabled by default.
# NOT applicable: Stringtie_Matrix, Tximport_STAR, Tximport_Salmon
#
# Usage:
#   Uncomment this config in PIPELINE_CONFIGS inside run_post_processing.sh
# ==============================================================================

# ==============================================================================
# MASTER REFERENCE
# ==============================================================================
#
# M5 uses transcript-level references.

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
    #"SmelGRF-GIF"
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
    "M5_RSEM_Bowtie2"
)

# ==============================================================================
# ANALYSES
# ==============================================================================
#
# Tximport_RSEM is a no-op when Matrix_Creation is also active (keep for
# standalone use).  Stringtie_Matrix and Tximport_STAR/Salmon are NOT applicable.
# Tissue_Specificity is applicable to transcript-level quantification.

ANALYSES=(
    # ---- M5 preprocessing ----
    "Tximport_RSEM"             # standalone; no-op when Matrix_Creation is active
    "Matrix_Creation"           # builds count matrix from RSEM quant output

    # ---- Visualisation ----
    "Basic_Heatmap"
    "Heatmap_with_CV"
)
