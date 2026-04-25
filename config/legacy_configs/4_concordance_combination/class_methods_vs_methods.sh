#!/bin/bash

# ==============================================================================
# CLASS: METHODS VS METHODS
# ==============================================================================
# Run concordance for multiple method subsets on one reference.

MASTER_REFERENCE="GPE001970_genome"

METHOD_COMBINATIONS=(
    "M1_HISAT2_RefGuided,M3_STAR_Align"
    "M1_HISAT2_RefGuided,M4_Salmon_Saf"
    "M3_STAR_Align,M5_RSEM_Bowtie2"
    "M2_HISAT2_DeNovo,M4_Salmon_Saf,M5_RSEM_Bowtie2"
)

RUN_ALL_METHOD_COMBINATIONS="TRUE"

GENE_GROUPS="SmelDMPs_v5_with_18s_and_HAP2,Selected_SmelGRF-GIF_with_two_GIF"
