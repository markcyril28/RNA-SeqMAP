#!/bin/bash

# ========================================================================================================================================
# SWITCHES
# ========================================================================================================================================
WORKSTATION_GRF_GIF_TO_HPC=OFF
HPC_TO_WORKSTATION_GRF_GIF=ON
HPC_TO_DRIVE=OFF
timestamp=$(date +"%Y%m%d_%H%M%S")
# ========================================================================================================================================
# MAIN FOLDERS
# ========================================================================================================================================
HPC_MAIN_WORKING_DIR="admontecillo@10.0.9.31:/home/admontecillo/MCRM_pipeline_HPC/GRF-GIF_private_backup"
WORKSTATION_GRF_GIF_WORKING_DIR="/mnt/c/Users/Admon/Documents/z_Pipeline/GRF-GIF_private_backup"

# ========================================================================================================================================
# WORKSTATION TO HPC TRANSFER
# ========================================================================================================================================
#mega-cc_12.0.14-1_amd64_beta.deb
FILE_TO_TRANSFER_from_WORKSTATION_file_n="4_Phylogenetic_Anaylsis/1_CONFIG_FILES/mega-cc_12.0.14-1_amd64_beta.deb"
FILE_TO_TRANSFER_from_WORKSTATION_file_1="mega-cc_12.0.14-1_amd64_beta.deb"
HPC_TARGET_FOLDER="4_Phylogenetic_Anaylsis/1_CONFIG_FILES"

if [ "$WORKSTATION_GRF_GIF_TO_HPC" = "ON" ]; then
    echo "Transferring files from WORKSTATION to HPC..."
    scp -r "$WORKSTATION_GRF_GIF_WORKING_DIR/$FILE_TO_TRANSFER_from_WORKSTATION_file_1" "$HPC_MAIN_WORKING_DIR/$HPC_TARGET_FOLDER"
else
    echo "WORKSTATION to HPC transfer is OFF"
fi

# ========================================================================================================================================
# HPC TO WORKSTATION TRANSFER
# ========================================================================================================================================

# Preprocessing
TARGET_FOLDER="4_Phylogenetic_Anaylsis/2_PHYLOGENETIC_TREE_RESULTS"
tar -czf "$HPC_MAIN_WORKING_DIR/${TARGET_FOLDER}_${timestamp}.tar.gz" "$TARGET_FOLDER"

# Main Variables
#FILE_TO_TRANSFER_from_HPC_file_1="$TARGET_FOLDER_$timestamp.tar.gz"
FILE_TO_TRANSFER_from_HPC_file_1="2_PHYLOGENETIC_TREE_RESULTS_HPC_20251008_030148.tar.gz"
WORKSTATION_GRF_GIF_TARGET_FOLDER="4_Phylogenetic_Anaylsis/3_FROM_HPC_Phylo_Results"

if [ "$HPC_TO_WORKSTATION_GRF_GIF" = "ON" ]; then
    echo "Transferring files from HPC to WORKSTATION..."
    scp -r "$HPC_MAIN_WORKING_DIR/4_Phylogenetic_Anaylsis/$FILE_TO_TRANSFER_from_HPC_file_1" "$WORKSTATION_GRF_GIF_WORKING_DIR/$WORKSTATION_GRF_GIF_TARGET_FOLDER"
else
    echo "HPC to WORKSTATION transfer is OFF"
fi

scp -r admontecillo@10.0.9.31:/home/admontecillo/MCRM_pipeline_HPC/HeatSeq/HISAT2_DE_NOVO_ROOT_HPC_20251009_093227.tar.gz $PWD

scp admontecillo@10.0.9.31:/home/admontecillo/MCRM_pipeline_HPC/HeatSeq/4b_Method_2_HISAT2_De_Novo_20251009_095825.tar.gz $PWD
# ========================================================================================================================================
# EXTRACT TAR.GZ FILES
# ========================================================================================================================================
echo "Extracting tar.gz files..."
if [ -f "${FILE_TO_TRANSFER_from_HPC_file_1}" ]; then
    tar -xzf "${WORKSTATION_GRF_GIF_WORKING_DIR}/${WORKSTATION_GRF_GIF_TARGET_FOLDER}/${FILE_TO_TRANSFER_from_HPC_file_1}" -C "${WORKSTATION_GRF_GIF_WORKING_DIR}/${WORKSTATION_GRF_GIF_TARGET_FOLDER}"
    echo "Extracted ${FILE_TO_TRANSFER_from_HPC_file_1}"
fi

if [ -f "HISAT2_DE_NOVO_ROOT_HPC_20251009_093227.tar.gz" ]; then
    tar -xzf "HISAT2_DE_NOVO_ROOT_HPC_20251009_093227.tar.gz"
    echo "Extracted HISAT2_DE_NOVO_ROOT_HPC_20251009_093227.tar.gz"
fi


scp admontecillo@10.0.9.31:/home/admontecillo/MCRM_pipeline_HPC/HeatSeq/CMSC244_{}.tar.gz $PWD

scp admontecillo@10.0.9.31:/home/admontecillo/MCRM_pipeline_HPC/HeatSeq/CMSC244_last.zip $PWD

7z a -t7z -mx=9 last_transfer.7z 1_SRRs/C_FastQC/ 3_POST_PROC/ logs/ logs_archive/

scp admontecillo@10.0.9.31:/home/admontecillo/MCRM_pipeline_HPC/HeatSeq/last_transfer.7z $PWD 

scp Dataset_2_OPTIMIZATION.tar.xz mcmercado@10.0.9.31:/home/mcmercado/Training_MicroSpore/TRAINING_WD 

scp mcmercado@10.0.9.31:~/Training_Microspore/class_balancing_1280_extended_training.zip $PWD