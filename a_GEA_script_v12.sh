#!/bin/bash
# ==============================================================================
# GENE EXPRESSION ANALYSIS (GEA) PIPELINE
# ==============================================================================
# RNA-seq analysis pipeline using multiple alignment/quantification methods
# Author: Mark Cyril R. Mercado | Version: v12 | Date: December 2025
# ==============================================================================

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT" || exit 1

# ==============================================================================
# CONFIGURATION FILE SELECTION
# ==============================================================================
# Comment/uncomment to select which configuration to load
# Only ONE config file should be active at a time
# ==============================================================================

CONFIG_FILES=(
	"config/HPC_full_run_config.sh"		# Full HPC run configuration
	#"config/local_test_config.sh"		# Local testing configuration
)

# Source the selected configuration file(s)
for config_file in "${CONFIG_FILES[@]}"; do
	if [[ -f "$config_file" ]]; then
		echo "Loading configuration: $config_file"
		source "$config_file"
	else
		echo "ERROR: Configuration file not found: $config_file"
		exit 1
	fi
done

# ==============================================================================
# DIRECTORY STRUCTURE AND OUTPUT PATHS
# ==============================================================================

# Create required directories (including log directories)
mkdir -p "$RAW_DIR_ROOT" "$TRIM_DIR_ROOT" "$FASTQC_ROOT" \
	"$HISAT2_REF_GUIDED_ROOT" "$HISAT2_REF_GUIDED_INDEX_DIR" "$STRINGTIE_HISAT2_REF_GUIDED_ROOT" \
	"$HISAT2_DE_NOVO_ROOT" "$HISAT2_DE_NOVO_INDEX_DIR" "$STRINGTIE_HISAT2_DE_NOVO_ROOT" \
	"$STAR_ALIGN_ROOT" "$STAR_INDEX_ROOT" "$STAR_GENOME_DIR" \
	"$SALMON_SAF_ROOT" "$SALMON_INDEX_ROOT" "$SALMON_QUANT_ROOT" "$SALMON_SAF_MATRIX_ROOT" \
	"$BOWTIE2_RSEM_ROOT" "$RSEM_INDEX_ROOT" "$RSEM_QUANT_ROOT" "$RSEM_MATRIX_ROOT" \
	"logs/log_files"

# ==============================================================================
# CLEANUP OPTIONS AND TESTING ESSENTIALS
# ==============================================================================

#rm -rf "$RAW_DIR_ROOT"                   # Remove previous raw SRR files
#rm -rf "$FASTQC_ROOT"                   # Remove previous FastQC results
#rm -rf "$HISAT2_DE_NOVO_ROOT"           # Remove previous HISAT2 results
#rm -rf "$HISAT2_DE_NOVO_INDEX_DIR"      # Remove previous HISAT2 index
#rm -rf "$STRINGTIE_HISAT2_DE_NOVO_ROOT" # Remove previous StringTie results

# ==============================================================================
# MAIN EXECUTION FUNCTIONS
# ==============================================================================

# Convert array to boolean flags for backward compatibility
RUN_MAMBA_INSTALLATION=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^MAMBA_INSTALLATION$" && echo "TRUE" || echo "FALSE")
RUN_DOWNLOAD_SRR=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^DOWNLOAD_SRR$" && echo "TRUE" || echo "FALSE")
RUN_TRIM_SRR=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^TRIM_SRR$" && echo "TRUE" || echo "FALSE")
RUN_DOWNLOAD_TRIM_and_DELETE_RAW_SRR=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^DOWNLOAD_TRIM_and_DELETE_RAW_SRR$" && echo "TRUE" || echo "FALSE")
RUN_GZIP_TRIMMED_FILES=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^GZIP_TRIMMED_FILES$" && echo "TRUE" || echo "FALSE")
RUN_DELETE_RAW_SRR=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^DELETE_RAW_SRR$" && echo "TRUE" || echo "FALSE")
RUN_QUALITY_CONTROL=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^QUALITY_CONTROL$" && echo "TRUE" || echo "FALSE")
RUN_METHOD_1_HISAT2_REF_GUIDED=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^METHOD_1_HISAT2_REF_GUIDED$" && echo "TRUE" || echo "FALSE")
RUN_METHOD_2_HISAT2_DE_NOVO=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^METHOD_2_HISAT2_DE_NOVO$" && echo "TRUE" || echo "FALSE")
RUN_METHOD_3_STAR_ALIGNMENT=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^METHOD_3_STAR_ALIGNMENT$" && echo "TRUE" || echo "FALSE")
RUN_METHOD_4_SALMON_SAF=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^METHOD_4_SALMON_SAF$" && echo "TRUE" || echo "FALSE")
RUN_METHOD_5_BOWTIE2_RSEM=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^METHOD_5_BOWTIE2_RSEM$" && echo "TRUE" || echo "FALSE")
RUN_HEATMAP_WRAPPER=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^HEATMAP_WRAPPER$" && echo "TRUE" || echo "FALSE")
RUN_ZIP_RESULTS=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^ZIP_RESULTS$" && echo "TRUE" || echo "FALSE")
RUN_DELETE_TRIMMED_FASTQ_FILES=$(printf '%s\n' "${PIPELINE_STAGES[@]}" | grep -q "^DELETE_TRIMMED_FASTQ_FILES$" && echo "TRUE" || echo "FALSE")

run_all() {
	# Main pipeline entrypoint: runs all steps for each FASTA and RNA-seq list
	# Steps: Logging, Download and Trim, HISAT2 alignment to Stringtie, and Cleanup. 
	local fasta=""
	local rnaseq_list=()
	# Parse arguments
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA)
				fasta="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1")
					shift
				done
				;;
			*)
				shift;;
		esac
	done

	local start_time end_time elapsed formatted_elapsed
	start_time=$(date +%s)
	setup_logging
	log_configuration
	log_step "Script started at: $(date -d @$start_time)"
	
	# Show pipeline configuration
	#show_pipeline_configuration

	log_info "SRR samples to process:"
	for SRR in "${rnaseq_list[@]}"; do
		log_info "$SRR"
	done

	if [[ $RUN_DOWNLOAD_SRR == "TRUE" ]]; then
		if [[ $RUN_DOWNLOAD_TRIM_and_DELETE_RAW_SRR == "TRUE" ]]; then
			log_warn "DOWNLOAD_SRR skipped: DOWNLOAD_TRIM_and_DELETE_RAW_SRR is enabled (use one or the other)"
		else
			log_step "STEP 01a: Download RNA-seq data"
			#download_srrs "${rnaseq_list[@]}"
			download_srrs_parallel "${rnaseq_list[@]}"
			#download_srrs_kingfisher "${rnaseq_list[@]}"
		fi
	fi

	if [[ $RUN_TRIM_SRR == "TRUE" ]]; then
		if [[ $RUN_DOWNLOAD_TRIM_and_DELETE_RAW_SRR == "TRUE" ]]; then
			log_warn "TRIM_SRR skipped: DOWNLOAD_TRIM_and_DELETE_RAW_SRR is enabled (use one or the other)"
		else
			log_step "STEP 01b: Trim RNA-seq data"
			#trim_srrs_trimmomatic "${rnaseq_list[@]}"
			trim_srrs_trimmomatic_parallel "${rnaseq_list[@]}"
		fi
	fi

	if [[ $RUN_DOWNLOAD_TRIM_and_DELETE_RAW_SRR == "TRUE" ]]; then
		log_step "STEP 01ab: Download, Trim, and Delete Raw SRR data"
		# Enable automatic deletion of raw files after successful trimming
		export DELETE_RAW_SRR_AFTER_DOWNLOAD_and_TRIMMING="TRUE"
		download_and_trim_srrs_parallel "${rnaseq_list[@]}"
	fi

	if [[ $RUN_DELETE_RAW_SRR == "TRUE" ]]; then
		log_step "STEP 01d: Delete Raw SRR files"
		delete_raw_srr_by_srr_list "${rnaseq_list[@]}"
	fi

	if [[ $RUN_QUALITY_CONTROL == "TRUE" ]]; then
		log_step "STEP 01c: Quality Control analysis"
		run_quality_control_all "${rnaseq_list[@]}"
	fi

	# Method 1: HISAT2 Reference-Guided Pipeline
	if [[ $RUN_METHOD_1_HISAT2_REF_GUIDED == "TRUE" ]]; then
		log_step "STEP 02a: HISAT2 Reference-Guided Pipeline"
		if [[ -z "$gtf_file" || ! -f "$gtf_file" ]]; then
			log_error "GTF file required for reference-guided alignment: $gtf_file"
			log_error "Skipping Method 1 - configure gtf_file variable"
		elif hisat2_ref_guided_pipeline --FASTA "$fasta" --GTF "$gtf_file" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 1 completed successfully"
		else
			log_error "Method 1 failed (exit code: $?) - continuing with remaining methods"
		fi
	fi

	# Method 2: HISAT2 De Novo Pipeline (Main method)
	if [[ $RUN_METHOD_2_HISAT2_DE_NOVO == "TRUE" ]]; then
		log_step "STEP 02b: HISAT2 De Novo Pipeline"
		if hisat2_de_novo_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 2 completed successfully"
		else
			log_error "Method 2 failed (exit code: $?) - continuing with remaining methods"
		fi
	fi

	# Method 3: STAR Alignment Pipeline
	if [[ $RUN_METHOD_3_STAR_ALIGNMENT == "TRUE" ]]; then
		log_step "STEP 03: STAR Splice-Aware Alignment"
		if star_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 3 completed successfully"
		else
			log_error "Method 3 failed (exit code: $?) - continuing with remaining methods"
		fi
	fi

	# Method 4: Salmon SAF Quantification
	if [[ $RUN_METHOD_4_SALMON_SAF == "TRUE" ]]; then
		log_step "STEP 04: Salmon SAF Quantification"
		local genome_file="$decoy"
		if [[ -f "$genome_file" ]]; then
			if salmon_saf_pipeline --FASTA "$fasta" --GENOME "$genome_file" --RNASEQ_LIST "${rnaseq_list[@]}"; then
				log_info "Method 4 completed successfully"
			else
				log_error "Method 4 failed (exit code: $?) - continuing with remaining methods"
			fi
		else
			log_warn "Genome file '$genome_file' not found. Skipping Salmon SAF pipeline."
			log_warn "Please provide genome file for decoy-aware Salmon quantification."
		fi
	fi

	# Method 5: Bowtie2 + RSEM Quantification
	if [[ $RUN_METHOD_5_BOWTIE2_RSEM == "TRUE" ]]; then
		log_step "STEP 05: Bowtie2 + RSEM Quantification"
		if bowtie2_rsem_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"; then
			log_info "Method 5 completed successfully"
		else
			log_error "Method 5 failed (exit code: $?) - continuing with remaining methods"
		fi
	fi

	# Generate cross-method validation summary
	local fasta_base="$(basename "$fasta")"
	local fasta_tag="${fasta_base%.*}"
	compare_methods_summary "$fasta_tag"

	end_time=$(date +%s)
	log_step "Final timing"
	log_info "Script ended at: $(date -d @$end_time)"
	elapsed=$((end_time - start_time))
	formatted_elapsed=$(date -u -d @${elapsed} +%H:%M:%S)
	log_info "Elapsed time: $formatted_elapsed"
}

# ==============================================================================
# SCRIPT EXECUTION
# ==============================================================================

# Call the function if installation is enabled
if [[ $RUN_MAMBA_INSTALLATION == "TRUE" ]]; then
	mamba_install
fi

if [[ $RUN_GZIP_TRIMMED_FILES == "TRUE" ]]; then
	log_step "Gzipping trimmed FASTQ files to save space"
	gzip_trimmed_fastq_files
fi

# Execute the pipeline for each FASTA input file
for fasta_input in "${ALL_FASTA_FILES[@]}"; do
	# Run the complete pipeline for each FASTA file with all SRR samples
	run_all --FASTA "$fasta_input" --RNASEQ_LIST "${SRR_COMBINED_LIST[@]}"
done

# ==============================================================================
# POST-PROCESSING: HEATMAP WRAPPER EXECUTION for HISAT2 DE NOVO
# ==============================================================================

if [[ $RUN_HEATMAP_WRAPPER == "TRUE" ]]; then
	log_step "Heatmap Wrapper post-processing enabled"
	
	# Navigate to post-processing directory
	if [[ ! -d "$POST_PROC_ROOT" ]]; then
		log_error "Directory '$POST_PROC_ROOT' not found"
	else
		cd "$POST_PROC_ROOT" || {
			log_error "Failed to change to $POST_PROC_ROOT directory"
			exit 1
		}
		
		# Execute the Heatmap Wrapper script for post-processing
		if [[ -f "run_all_post_processing.sh" ]]; then
			log_step "Executing Heatmap Wrapper post-processing script"
			chmod +x ./*.sh
			chmod +x run_all_post_processing.sh
			
			if bash "run_all_post_processing.sh" 2>&1; then
				log_info "Heatmap Wrapper completed successfully"
			else
				exit_code=$?
				log_error "Heatmap Wrapper failed with exit code $exit_code"
			fi
		else
			log_warn "Heatmap Wrapper script 'run_all_post_processing.sh' not found - skipping"
		fi
		
		# Return to original directory
		cd - > /dev/null || log_warn "Failed to return to previous directory"
	fi
fi

# ==============================================================================
# POST-PROCESSING OPTIONS (COMMENTED OUT)
# ==============================================================================

if [[ $RUN_ZIP_RESULTS == "TRUE" ]]; then
	# Optional: Archive StringTie results for sharing or backup
	#tar -czvf "stringtie_results_$(date +%Y%m%d_%H%M%S).tar.gz" "$STRINGTIE_HISAT2_DE_NOVO_ROOT"
	#tar -czvf HISAT2_DE_NOVO_ROOT_HPC_$(date +%Y%m%d_%H%M%S).tar.gz $HISAT2_DE_NOVO_ROOT
	#tar -czvf 4b_Method_2_HISAT2_De_Novo_$(date +%Y%m%d_%H%M%S).tar.gz 4b_Method_2_HISAT2_De_Novo/
	log_step "Creating compressed archive for folders: $POST_PROC_ROOT and logs"
	TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
	tar -c \
		--exclude="${POST_PROC_ROOT}/M3_STAR_Align" \
		"$POST_PROC_ROOT" logs | pigz -p "$THREADS" > "CMSC244_${TIMESTAMP}.tar.gz"
	log_info "Archive created: CMSC244_${TIMESTAMP}.tar.gz"
fi

# ==============================================================================
# CLEANUP: DELETE TRIMMED FASTQ FILES
# ==============================================================================

if [[ $RUN_DELETE_TRIMMED_FASTQ_FILES == "TRUE" ]]; then
	log_step "Deleting trimmed FASTQ files for SRR_COMBINED_LIST"
	delete_trimmed_fastq_by_srr_list "${SRR_COMBINED_LIST[@]}"
fi

# ==============================================================================
# SOFTWARE CATALOG
# ==============================================================================
catalog_all_software

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
echo "END OF SCRIPT"
