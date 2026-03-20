#!/bin/bash

# ==============================================================================
# CLASS: GENES VS GENES
# ==============================================================================
# Run concordance using multiple gene-group sets on one reference/method panel.

MASTER_REFERENCE="GPE001970_genome"
METHODS="M1_HISAT2_RefGuided M2_HISAT2_DeNovo M3_STAR_Align M4_Salmon_Saf M5_RSEM_Bowtie2"

GENE_GROUP_COMBINATIONS=(
    "SmelDMPs_v5_with_18s_and_HAP2"
    "Selected_SmelGRF-GIF_with_two_GIF"
    "SmelDMPs_v5_with_18s_and_HAP2,Selected_SmelGRF-GIF_with_two_GIF"
)

RUN_ALL_GENE_GROUP_COMBINATIONS="TRUE"
