#!/bin/bash

# ==============================================================================
# CONCORDANCE DEFAULTS CLASS
# ==============================================================================
# Baseline single-run defaults for run_concordance.sh.
# Any values here can be overridden by env vars or other classes loaded later.

MASTER_REFERENCE="GPE001970_genome"

METHODS="M1_HISAT2_RefGuided M2_HISAT2_DeNovo M3_STAR_Align M4_Salmon_Saf M5_RSEM_Bowtie2"

GENE_GROUPS="SmelDMPs_v5_with_18s_and_HAP2,Selected_SmelGRF-GIF_with_two_GIF"

# Keep combination loops disabled by default
RUN_ALL_MASTER_REFERENCES="FALSE"
RUN_ALL_METHOD_COMBINATIONS="FALSE"
RUN_ALL_GENE_GROUP_COMBINATIONS="FALSE"
