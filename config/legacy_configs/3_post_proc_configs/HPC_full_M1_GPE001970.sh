#!/bin/bash

# ==============================================================================
# HPC FULL CONFIG: M1 — HISAT2 Reference-Guided (GPE001970_genome)
# ==============================================================================
# Purpose: Production post-processing for M1 (HISAT2 ref-guided genome
#          alignment).  Uses the full eggplant tissue atlas.
#
# Method:
#   M1 — HISAT2 Reference-Guided  →  uses GPE001970_genome
#
# Preprocessing:  prepde_matrix_linker.sh runs automatically (integer counts).
# DE support:     YES — prepDE.py integer counts feed directly into DESeq2.
# NOT applicable: Tximport_*, Matrix_Creation, Stringtie_Matrix
#
# Usage:
#   Uncomment this config in PIPELINE_CONFIGS inside run_post_processing.sh
# ==============================================================================

# ==============================================================================
# MASTER REFERENCE
# ==============================================================================
#
# M1 uses genome-level references.
#
# Reference → Source files:
#   GPE001970_genome →  inputs/fasta/reference_genome/GPE001970_genome.fa
#                        inputs/gtf/reference/GPE001970_genome.gtf

MASTER_REFERENCES=(
    "GPE001970_genome"
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
)
