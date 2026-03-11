#!/bin/bash

# ==============================================================================
# HPC FULL CONFIG: GENOME-BASED METHODS (M1 + M3) — FULL DATASET
# ==============================================================================
# Purpose: Production post-processing run for genome-alignment methods.
#          Uses the full eggplant tissue atlas and Eggplant V4.1 genome reference.
#          (Matches alignment config: config/2_alignment/HPC_full_ref_guided.sh)
#
# Methods:
#   M1 — HISAT2 Reference-Guided  →  uses Eggplant_V4.1
#   M3 — STAR Splice-Aware Alignment  →  uses Eggplant_V4.1
#
# Usage:
#   bash run_all_post_processing.sh --config config/3_post_proc_configs/HPC_full_ref_guided.sh
# ==============================================================================

# ==============================================================================
# SYSTEM RESOURCES
# ==============================================================================

THREADS=64
ENABLE_GPU="false"
ENABLE_GNU_PARALLEL="TRUE"
DESIRED_CPU_PER_JOB=4
AVAILABLE_RAM_GB=32
GPU_VRAM_GB=8

# ==============================================================================
# LOGGING AND OUTPUT
# ==============================================================================

CLEAR_LOGS="FALSE"
CLEAR_OUTPUT_FOLDER="TRUE"

# ==============================================================================
# MASTER REFERENCE (uncomment ONE)
# ==============================================================================

MASTER_REFERENCES=(
    #"All_Smel_Genes"
    #"Eggplant_V4.1_transcripts.function"
    "Eggplant_V4.1"               # M1_HISAT2_RefGuided / M3_STAR_Align  (Eggplant_V4.1.fa) — run alignment first
    #"GPE001970_genome"            # M1_HISAT2_RefGuided / M3_STAR_Align  (GPE001970_genome.fa)  — test only
    #"GPE001970_transcripts"       # M2_HISAT2_DeNovo only
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
    "PRJNA328564"             # Main Dataset — Eggplant tissue atlas (full)
    #"PRJNA328564_selected"
    #"PRJNA865018"            # SET_1: Good Dataset for SmelDMP GEA
    #"PRJNA941250"            # SET_2: Good Dataset for SmelDMP GEA
    #"PRJNA865018_and_PRJNA941250"
    #"SAMN28540077"           # Chinese Dataset 1 — replicability
    #"SAMN28540068"           # Chinese Dataset 2 — replicability
    #"OTHER_SRR_LIST"
)

# ==============================================================================
# ALIGNMENT METHODS (comment/uncomment to enable/disable)
# ==============================================================================

METHODS=(
    "M1_HISAT2_RefGuided"
    #"M2_HISAT2_DeNovo"
    "M3_STAR_Align"
    #"M4_Salmon_Saf"
    #"M5_RSEM_Bowtie2"
)

# ==============================================================================
# ANALYSES (comment/uncomment to enable/disable)
# ==============================================================================

ANALYSES=(
    #"Stringtie_Matrix"         # M2 only — NOT for M1

    "Tximport_STAR"             # M3 standalone (no-op if Matrix_Creation active)
    #"Tximport_Salmon"          # M4 only
    #"Tximport_RSEM"            # M5 only
    "Matrix_Creation"           # M3/M4/M5 only — NOT for HISAT2

    "Basic_Heatmap"
    "Heatmap_with_CV"
    #"BarGraph"
    "Differential_Expression"   # M1 supported via prepDE.py integer counts
    #"Gene_Set_Enrichment"
    "PCA_Dimensionality_Reduction"
    "Sample_Correlation_Clustering"
    #"Tissue_Specificity"
    #"Coexpression_using_WGCNA"
)
