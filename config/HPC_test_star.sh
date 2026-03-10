#!/bin/bash

# ==============================================================================
# IMPORTANT PARAMETERS
# ==============================================================================

# Runtime Configuration
THREADS=48                              # Threads for parallel operations
JOBS=3									# Parallel jobs for GNU Parallel
USE_GNU_PARALLEL="TRUE"                 # TRUE/FALSE for GNU Parallel
keep_bam_global="n"                     # y=keep BAM files, n=delete after

# Pipeline Stages (comment/uncomment to enable/disable)
PIPELINE_STAGES=(
	#"MAMBA_INSTALLATION"
	
	# Option A: Separate download and trim (keeps raw files)
	#"DOWNLOAD_SRR"
	#"TRIM_SRR"
	
	# Option B: Combined download+trim+cleanup (auto-deletes raw after trim)
	#"DOWNLOAD_TRIM_and_DELETE_RAW_SRR"
	
	#"GZIP_TRIMMED_FILES"
	#"QUALITY_CONTROL"

	#"DELETE_RAW_SRR"				# Manually delete raw SRR files
	#"DELETE_TRIMMED_FASTQ_FILES"	# Manually delete trimmed files

	#"METHOD_1_HISAT2_REF_GUIDED"
	#"METHOD_2_HISAT2_DE_NOVO"
	"METHOD_3_STAR_ALIGNMENT"
	#"METHOD_4_SALMON_SAF"
	#"METHOD_5_BOWTIE2_RSEM"

	#"HEATMAP_WRAPPER"
	#"ZIP_RESULTS"
)

# ==============================================================================
# CONDA ENVIRONMENT
# ==============================================================================

eval "$(conda shell.bash hook)"
conda activate gea

# ==============================================================================
# SOURCE MODULES
# ==============================================================================
# Structure: logging/, a_preprocessing/, b_main_methods/, 0_input_information/
# ==============================================================================

source "modules/modules_loader.sh"
#bash init_setup.sh

# Calculate threads per job
if [[ "$USE_GNU_PARALLEL" == "TRUE" ]]; then
	THREADS_PER_JOB=$((THREADS / JOBS))
	[[ $THREADS_PER_JOB -lt 1 ]] && THREADS_PER_JOB=1
else
	THREADS_PER_JOB=$THREADS
fi
export THREADS JOBS USE_GNU_PARALLEL THREADS_PER_JOB keep_bam_global

# ==============================================================================
# INPUT FILES AND DATA SOURCES
# ==============================================================================

# Legacy GTF and FASTA references (commented out)
#All_SmelGIF_GTF_FILE="inputs/All_SmelDMP_Head_Gene_Name_v4.gtf"
#Eggplant_V4_1_transcripts_function_FASTA_FILE="inputs/fasta/reference_genomes/Eggplant_V4.1_transcripts.function.fa"
#ALL_Smel_Genes_Full_Name_reformatted_GTF_FILE="inputs/gtf/reference/Eggplant_V4.1_function_IPR_reformatted_v4.gtf"
#ALL_Smel_Genes_Full_Name_reformatted_GTF_FILE="inputs/gtf/reference/GPE001970_transcripts.gtf"
ALL_Smel_Genes_Full_Name_reformatted_GTF_FILE="inputs/gtf/reference/GPE001970.gtf"

#decoy="inputs/fasta/experimental/TEST.fasta"
#gtf_file="${Eggplant_V4_1_transcripts_function_FASTA_FILE}"
gtf_file="${ALL_Smel_Genes_Full_Name_reformatted_GTF_FILE}"

# FASTA Files for Analysis
ALL_FASTA_FILES=(
	# List of FASTA files to process
	#"inputs/fasta/reference_genomes/Eggplant_V4.1_transcripts.function.fa"
	#"inputs/fasta/reference_genomes/GPE001970_transcripts.fa"
	"inputs/fasta/reference_genomes/GPE001970.fa"
)

# ==============================================================================
# RNA-SEQ DATA SOURCES (SRR LISTS)
# ==============================================================================

# SRA Run Selector Tool link: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA328564&o=acc_s%3Aa

SRR_LIST_PRJNA328564=(
	# Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA328564&o=acc_s%3Aa
	# Developmental stages arranged from early to late
	#SRR3884685	# Radicles (earliest - germination)
	#SRR3884677	# Cotyledons (seed leaves)
	#SRR3884675	# Roots (root development)
	#SRR3884690	# Stems (vegetative growth)
	#SRR3884689	# Leaves (vegetative growth)
	#SRR3884684	# Senescent_leaves (leaf aging)
	SRR3884686	# Buds_0.7cm (flower bud initiation) [MAIN INTEREST]
	SRR3884687	# Opened_Buds (flower development) 	 [MAIN INTEREST]
	SRR3884597	# Flowers (anthesis)/				 [MAIN INTEREST]	
	#SRR3884679	# Pistils (female reproductive parts)
	#SRR3884608	# Fruits_1cm (early fruit development)
	#SRR3884620	# Fruits_Stage_1 (early fruit stage)
	#SRR3884631	# Fruits_6cm (fruit enlargement)
	#SRR3884642	# Fruits_Skin_Stage_2 (mid fruit development)
	#SRR3884653	# Fruits_Flesh_Stage_2 (mid fruit development)
	#SRR3884664	# Fruits_Calyx_Stage_2 (mid fruit development)
	#SRR3884680	# Fruits_Skin_Stage_3 (late fruit development)
	#SRR3884681	# Fruits_Flesh_Stage_3 (late fruit development)
	#SRR3884678	# Fruits_peduncle (fruit attachment)
)

SRR_LIST_SAMN28540077=(
	# Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN28540077&o=acc_s%3Aa
	#SRR20722232	# Mature_fruits (10 GB file); corrected.
	#SRR20722226 # Young_fruits
	#SRR20722228	# sepals (too large; not included)
	SRR20722234	# Flowers 
	SRR4243802 # Buds, Adopted Dataset from ID: PRJNA341784 
	#SRR20722233	# leaf_buds 
	#SRR20722230	# mature_leaves (14 GB file; not included)
	#SRR20722227	# stems
	#SRR20722229	# roots
)

SRR_LIST_SAMN28540068=(
	#Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN28540068&o=acc_s%3Aa
	#SRR20722387 # mature_fruits
	SRR3884597 	# Flower
	SRR20722297 # flower_buds 
	#SRR20722385 # sepals (not included)
	#SRR20722296 # leaf_buds 
	#SRR20722386 # mature_leaves (not included)
	#SRR20722383 # young_leaves (not included)
	#SRR20722384 # stems
	#SRR31755282 # Roots (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP552204&o=acc_s%3Aa)
)

SRR_LIST_PRJNA865018=(
# Set_1: A Good Dataset for SmelDMP GEA: 
# 	https://www.ncbi.nlm.nih.gov/Traces/study/?acc=%20%20PRJNA865018&o=acc_s%3Aa PRJNA865018
	SRR21010466	# buds_1
	SRR21010456	# buds_2
	SRR21010454	# buds_3
	SRR21010462	# flowers_1
	SRR21010460	# flowers_2
	SRR21010458	# flowers_3
	SRR21010452	# fruits_1
	SRR21010450	# fruits_2
	SRR21010464	# fruits_3
)

SRR_LIST_PRJNA941250=(
# Set_2: A Good Dataset for SmelDMP GEA:
#	https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA941250&o=acc_s%3Aa PRJNA941250 # Buds, Opened Buds
	SRR23909869 # 8DBF_1
	SRR23909870 # 8DBF_2
	SRR23909871 # 8DBF_3
	SRR23909866 # 5DBF_1
	SRR23909867 # 5DBF_2
	SRR23909868 # 5DBF_3
	SRR23909863 # Fully Develop (FD) 1
	SRR23909864 # Fully Develop (FD) 2
	SRR23909865 # Fully Develop (FD) 3
)

OTHER_SRR_LIST=(
	# Possible Source: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP390977&o=acc_s%3Aa
	SRR34564302	# Fruits (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR34564302&display=metadata)
	SRR34848077 # Leaves (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR34848077&display=metadata)
	PRJNA613773 # Leaf (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA613773&o=acc_s%3Aa)
		# Cotyledons (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR3884677&display=metadata)
	SRR3479277 # Pistil (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR3479277&display=metadata)
	SRR3884597 # Flowers (https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR3884597&display=metadata)

	#PRJNA341784 # Flower buds lang (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA341784&o=acc_s%3Aa)
	#PRJNA477924 # Leaf and Root (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA477924&o=acc_s%3Aa)
	 # Leaves
	 # Stems
	 # Radicles
	
	# Add other SRR IDs here if needed
)

SRR_COMBINED_LIST=(
	"${SRR_LIST_PRJNA328564[@]}"	# Main Dataset for GEA. 
	#"${SRR_LIST_SAMN28540077[@]}"	# Chinese Dataset for replicability. 
	#"${SRR_LIST_SAMN28540068[@]}"	# Chinese Dataset for replicability. 
	#"${SRR_LIST_PRJNA865018[@]}"	# SET_1: Good Dataset for SmelDMP GEA.
	#"${SRR_LIST_PRJNA941250[@]}"	# SET_2: Good Dataset for SmelDMP GEA.
)

# ==============================================================================
# DIRECTORY STRUCTURE AND OUTPUT PATHS
# ==============================================================================

POST_PROC_ROOT="3_POST_PROC"
export POST_PROC_ROOT

# Create required directories
mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT" \
	"$HISAT2_REF_GUIDED_ROOT" "$HISAT2_REF_GUIDED_INDEX_DIR" "$STRINGTIE_HISAT2_REF_GUIDED_ROOT" \
	"$HISAT2_DE_NOVO_ROOT" "$HISAT2_DE_NOVO_INDEX_DIR" "$STRINGTIE_HISAT2_DE_NOVO_ROOT" \
	"$STAR_ALIGN_ROOT" "$STAR_INDEX_ROOT" \
	"$SALMON_SAF_ROOT" "$SALMON_INDEX_ROOT" "$SALMON_QUANT_ROOT" "$SALMON_SAF_MATRIX_ROOT" \
	"$BOWTIE2_RSEM_ROOT" "$RSEM_INDEX_ROOT" "$RSEM_QUANT_ROOT" "$RSEM_MATRIX_ROOT" \
	"$HISAT2_REF_GUIDED_MATRIX_ROOT" "$HISAT2_DE_NOVO_MATRIX_ROOT" "$STAR_MATRIX_ROOT"

# ==============================================================================
# CLEANUP OPTIONS AND TESTING ESSENTIALS
# ==============================================================================
# Uncomment lines below to remove previous results before re-running.
# WARNING: These are destructive operations — verify before uncommenting.
# ==============================================================================

ACTIVATE_RM=FALSE

# --- Preprocessing ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$RAW_DIR_ROOT"                          # Raw SRR downloads
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$TRIM_DIR_ROOT"                         # Trimmed FASTQ files
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$FASTQC_ROOT"                           # FastQC reports

# --- Method 1: HISAT2 Reference-Guided ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_REF_GUIDED_ROOT"                # HISAT2 ref-guided alignments
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_REF_GUIDED_INDEX_DIR"           # HISAT2 ref-guided index
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STRINGTIE_HISAT2_REF_GUIDED_ROOT"      # StringTie (ref-guided) assemblies
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_REF_GUIDED_MATRIX_ROOT"         # Count matrices (ref-guided)

# --- Method 2: HISAT2 De Novo ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_DE_NOVO_ROOT"                   # HISAT2 de novo alignments
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_DE_NOVO_INDEX_DIR"              # HISAT2 de novo index
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STRINGTIE_HISAT2_DE_NOVO_ROOT"         # StringTie (de novo) assemblies
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$HISAT2_DE_NOVO_MATRIX_ROOT"            # Count matrices (de novo)

# --- Method 3: STAR Alignment ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STAR_ALIGN_ROOT"                       # STAR alignment results (removes index + BAMs)
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STAR_INDEX_ROOT"                       # STAR genome index only
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$STAR_MATRIX_ROOT"                      # Count matrices (STAR)

# --- Method 4: Salmon SAF ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$SALMON_SAF_ROOT"                       # Salmon SAF root
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$SALMON_INDEX_ROOT"                     # Salmon index
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$SALMON_QUANT_ROOT"                     # Salmon quantification
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$SALMON_SAF_MATRIX_ROOT"                # Count matrices (Salmon)

# --- Method 5: Bowtie2 + RSEM ---
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$BOWTIE2_RSEM_ROOT"                     # Bowtie2/RSEM root
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$RSEM_INDEX_ROOT"                       # RSEM index
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$RSEM_QUANT_ROOT"                       # RSEM quantification
[[ "$ACTIVATE_RM" == "TRUE" ]] && rm -rf "$RSEM_MATRIX_ROOT"                      # Count matrices (RSEM)
