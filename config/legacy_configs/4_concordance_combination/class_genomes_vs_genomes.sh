#!/bin/bash

# ==============================================================================
# CLASS: GENOMES VS GENOMES
# ==============================================================================
# Run concordance once per reference in MASTER_REFERENCES.
# This compares method concordance across different references by producing
# per-reference reports and figures.
# Includes both genome and transcript reference variants for the two primary
# master references.

MASTER_REFERENCES=(
    "GPE001970_genome"
    "GPE001970_transcripts"
    "Eggplant_V4.1_genome"
    "Eggplant_V4.1_transcripts.function"
)

RUN_ALL_MASTER_REFERENCES="TRUE"

# Use full method panel unless you want to restrict it
METHODS="M1_HISAT2_RefGuided M2_HISAT2_DeNovo M3_STAR_Align M4_Salmon_Saf M5_RSEM_Bowtie2"

# Keep default target gene groups unless overridden
GENE_GROUPS="SmelDMPs_v5_with_18s_and_HAP2,Selected_SmelGRF-GIF_with_two_GIF"
