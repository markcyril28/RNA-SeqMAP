#!/bin/bash
# ===============================================
# BASH CONFIGURATIONS FOR METHOD 3 POST-PROCESSING
# ===============================================
# STAR Alignment + Salmon Quantification + tximport

# Method 3 specific paths
# The pipeline stores outputs at: M3_STAR_Align/
# Salmon quant: M3_STAR_Align/star_index/<fasta>_salmon_quant/
# Matrices: M3_STAR_Align/count_matrices_from_STAR/ (tximport output)

# Configure threads and paths
export THREADS=${THREADS:-8}
export METHOD_NAME="M3_STAR_Salmon"
