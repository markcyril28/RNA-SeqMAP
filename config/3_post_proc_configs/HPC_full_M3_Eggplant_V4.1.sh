#!/bin/bash

# ==============================================================================
# HPC FULL CONFIG: M3 — STAR Splice-Aware Alignment (Eggplant_V4.1)
# ==============================================================================
# Purpose: Production post-processing for M3 (STAR genome alignment).
#          Uses the full eggplant tissue atlas.
#
# Method:
#   M3 — STAR Splice-Aware Alignment  →  uses Eggplant_V4.1
#
# Preprocessing:  Tximport_STAR → Matrix_Creation (Tximport is a no-op when
#                 Matrix_Creation is active; keep both for flexibility).
# DE support:     YES — tximport counts feed into DESeq2.
# NOT applicable: Stringtie_Matrix, Tximport_Salmon, Tximport_RSEM
#
# Usage:
#   Uncomment this config in PIPELINE_CONFIGS inside run_all_post_processing.sh
# ==============================================================================

# ==============================================================================
# MASTER REFERENCE
# ==============================================================================
#
# M3 uses genome-level references.
#
# Reference → Source files:
#   Eggplant_V4.1  →  inputs/fasta/reference_genomes/Eggplant_V4.1.fa
#                      inputs/gtf/reference/Eggplant_V4.1_function_IPR_final_stringtie.gtf

MASTER_REFERENCES=(
    "Eggplant_V4.1"
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
    "Selected_GRF_GIF_Genes_vTwo_GIF.csv"
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
    "M3_STAR_Align"
)

# ==============================================================================
# ANALYSES
# ==============================================================================
#
# Tximport_STAR is a no-op when Matrix_Creation is also active (keep for
# standalone use).  Stringtie_Matrix is NOT applicable to M3.
# Differential_Expression IS supported (tximport counts → DESeq2).

ANALYSES=(
    # ---- M3 preprocessing ----
    "Tximport_STAR"             # standalone; no-op when Matrix_Creation is active
    "Matrix_Creation"           # builds count matrix from STAR quant output

    # ---- Visualisation ----
    "Basic_Heatmap"
    "Heatmap_with_CV"
    #"BarGraph"

    # ---- Statistics ----
    #"Differential_Expression"   # Supported for M3 (tximport counts)
    #"PCA_Dimensionality_Reduction"
    #"Sample_Correlation_Clustering"
    #"Tissue_Specificity"       # Not applicable — genome-based counts
    #"Coexpression_using_WGCNA"
    #"Gene_Set_Enrichment"
)
