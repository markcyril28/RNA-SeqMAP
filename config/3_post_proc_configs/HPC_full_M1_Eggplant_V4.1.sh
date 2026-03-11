#!/bin/bash

# ==============================================================================
# HPC FULL CONFIG: M1 — HISAT2 Reference-Guided (Eggplant_V4.1)
# ==============================================================================
# Purpose: Production post-processing for M1 (HISAT2 ref-guided genome
#          alignment).  Uses the full eggplant tissue atlas.
#
# Method:
#   M1 — HISAT2 Reference-Guided  →  uses Eggplant_V4.1
#
# Preprocessing:  prepde_matrix_linker.sh runs automatically (integer counts).
# DE support:     YES — prepDE.py integer counts feed directly into DESeq2.
# NOT applicable: Tximport_*, Matrix_Creation, Stringtie_Matrix
#
# Usage:
#   Uncomment this config in PIPELINE_CONFIGS inside run_all_post_processing.sh
# ==============================================================================

# ==============================================================================
# SYSTEM RESOURCES
# ==============================================================================

THREADS=96
ENABLE_GPU="FALSE"
ENABLE_GNU_PARALLEL="TRUE"
DESIRED_CPU_PER_JOB=2
AVAILABLE_RAM_GB=64
GPU_VRAM_GB=8

# ==============================================================================
# LOGGING AND OUTPUT
# ==============================================================================

CLEAR_LOGS="FALSE"
CLEAR_OUTPUT_FOLDER="FALSE"

# ==============================================================================
# MASTER REFERENCE
# ==============================================================================
#
# M1 uses genome-level references.
#
# Reference → Source files:
#   Eggplant_V4.1  →  inputs/fasta/reference_genomes/Eggplant_V4.1.fa
#                      inputs/gtf/reference/Eggplant_V4.1_function_IPR_final.gtf

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
    #"SmelGRF-GIFs"              # Needs stringtie_matrix_builder.sh first
    #"SmelGRF-GIFs_with_1_18s_rRNA"
    #"SmelGRFs"
    #"SmelGIFs"
    #"Selected_GRF_GIF_Genes_vAll_GIFs"
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
    "M1_HISAT2_RefGuided"
)

# ==============================================================================
# ANALYSES
# ==============================================================================
#
# M1 uses prepde_matrix_linker.sh (runs automatically) for integer count matrices.
# Tximport_* and Matrix_Creation are NOT applicable to M1.
# Differential_Expression IS supported (prepDE integer counts → DESeq2).

ANALYSES=(
    # ---- M1 preprocessing (automatic — do not add Tximport_* or Matrix_Creation) ----

    # ---- Visualisation ----
    "Basic_Heatmap"
    "Heatmap_with_CV"
    #"BarGraph"

    # ---- Statistics ----
    #"Differential_Expression"   # Supported for M1 (prepDE integer counts)
    #"PCA_Dimensionality_Reduction"
    #"Sample_Correlation_Clustering"
    #"Tissue_Specificity"       # Not applicable — genome-based counts
    #"Coexpression_using_WGCNA"
    #"Gene_Set_Enrichment"
)
