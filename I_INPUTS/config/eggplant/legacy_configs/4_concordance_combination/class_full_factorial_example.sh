#!/bin/bash

# ==============================================================================
# CLASS: FULL FACTORIAL EXAMPLE
# ==============================================================================
# Demonstrates nested combination runs:
#   references x method_combinations x gene_group_combinations
#
# Keep list sizes small at first; total runs = R x M x G.

MASTER_REFERENCES=(
    "GPE001970_genome"
    "Eggplant_V4.1_genome"
)

METHOD_COMBINATIONS=(
    "M1_HISAT2_RefGuided,M3_STAR_Align"
    "M4_Salmon_Saf,M5_RSEM_Bowtie2"
)

GENE_GROUP_COMBINATIONS=(
    "SmelDMPs_v5_with_18s_and_HAP2"
    "Selected_SmelGRF-GIF_with_two_GIF"
)

RUN_ALL_MASTER_REFERENCES="TRUE"
RUN_ALL_METHOD_COMBINATIONS="TRUE"
RUN_ALL_GENE_GROUP_COMBINATIONS="TRUE"
