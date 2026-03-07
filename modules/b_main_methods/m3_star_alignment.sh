#!/bin/bash
# ==============================================================================
# METHOD 3: STAR SPLICE-AWARE ALIGNMENT PIPELINE
# ==============================================================================
# STAR splice-aware alignment with Salmon quantification and tximport for DESeq2
# Uses 2-pass mode for better junction detection
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${M3_STAR_SOURCED:-}" == "true" ]] && return 0
export M3_STAR_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_method.sh"

# ==============================================================================
# STAR CONFIGURATION - IMPORTANT PARAMETERS AT TOP
# ==============================================================================

# STAR-specific settings (GPU-aware defaults)
STAR_GENOME_LOAD="${STAR_GENOME_LOAD:-$(get_star_genome_load 2>/dev/null || echo NoSharedMemory)}"
STAR_READ_LENGTH="${STAR_READ_LENGTH:-100}"
STAR_STRAND_SPECIFIC="${STAR_STRAND_SPECIFIC:-intronMotif}"

# Delete transient big files after alignment (saves disk space)
# Set to "true" to delete intermediate files (2-pass genome, unsorted BAM, etc.)
# Set to "false" to keep all files for debugging
STAR_DELETE_TRANSIENT="${STAR_DELETE_TRANSIENT:-true}"

# GTF annotation file for splice junction detection (required for --sjdbOverhang)
STAR_GTF_FILE="${STAR_GTF_FILE:-0_INPUT_FASTAs/gtf/reference/Eggplant_V4.1_function_IPR_final_formatted_v3.gtf}"

# STAR temp directory configuration
# Set to "system" to use /tmp, "local" to use output dir, "cwd" for current directory, "none" to let STAR manage, or a specific path
STAR_TMP_MODE="${STAR_TMP_MODE:-cwd}"

# ==============================================================================
# STAR TISSUE-SPECIFIC PIPELINE
# ==============================================================================

star_tissue_specific_pipeline() {
	local fasta="" rnaseq_list=() metadata_file=""
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--METADATA) metadata_file="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1"); shift
				done;;
			*) log_error "Unknown option: $1"; return 1;;
		esac
	done
	
	[[ -z "$fasta" || ! -f "$fasta" ]] && { log_error "Valid FASTA required"; return 1; }
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")
	
	declare -A sample_metadata
	if [[ -n "$metadata_file" && -f "$metadata_file" ]]; then
		load_sample_metadata "$metadata_file" sample_metadata || {
			log_warn "Metadata load failed - running pooled alignment"
			star_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
			return $?
		}
	else
		log_warn "No metadata - running pooled alignment"
		star_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		return $?
	fi
	
	# Group samples by tissue
	declare -A tissue_samples
	for SRR in "${rnaseq_list[@]}"; do
		local tissue="${sample_metadata[${SRR}_condition]:-unknown}"
		tissue_samples["$tissue"]+="$SRR "
	done
	
	local tissue_count=${#tissue_samples[@]}
	log_info "[STAR TISSUE] Found $tissue_count tissue types"
	
	if [[ $tissue_count -lt 2 ]]; then
		log_warn "Single tissue detected - using pooled alignment"
		star_alignment_pipeline --FASTA "$fasta" --RNASEQ_LIST "${rnaseq_list[@]}"
		return $?
	fi
	
	log_info "[STAR TISSUE] Running tissue-specific alignments for better isoform detection"
	
	# Run STAR per tissue
	for tissue in "${!tissue_samples[@]}"; do
		local tissue_srrs=(${tissue_samples[$tissue]})
		log_step "STAR alignment for tissue: $tissue (${#tissue_srrs[@]} samples)"
		
		star_alignment_pipeline \
			--FASTA "$fasta" \
			--RNASEQ_LIST "${tissue_srrs[@]}" \
			--TISSUE_TAG "$tissue"
	done
}

# ==============================================================================
# MAIN STAR ALIGNMENT PIPELINE
# ==============================================================================

star_alignment_pipeline() {
	local fasta="" rnaseq_list=() tissue_tag=""
	
	# Check for GNU parallel
	if ! command -v parallel >/dev/null 2>&1; then
		log_warn "[PERFORMANCE] GNU parallel not found - Salmon quantification will run sequentially"
	fi
	
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--TISSUE_TAG) tissue_tag="$2"; shift 2;;
			--RNASEQ_LIST)
				shift
				while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
					rnaseq_list+=("$1"); shift
				done;;
			*) log_error "Unknown option: $1"; return 1;;
		esac
	done
	
	[[ -z "$fasta" ]] && { log_error "No FASTA file specified. Use --FASTA <fasta_file>."; return 1; }
	[[ ! -f "$fasta" ]] && { log_error "FASTA file '$fasta' not found."; return 1; }
	[[ ${#rnaseq_list[@]} -eq 0 ]] && rnaseq_list=("${SRR_COMBINED_LIST[@]}")
	[[ ${#rnaseq_list[@]} -eq 0 ]] && { log_error "No RNA-seq samples provided."; return 1; }

	local fasta_base fasta_tag star_index_dir star_genome_dir
	fasta_base="$(basename "$fasta")"
	fasta_tag="${fasta_base%.*}"
	
	# Get absolute paths - using realpath for robustness, fallback to manual resolution
	local abs_star_index_root abs_star_align_root
	
	# First ensure the directories exist
	mkdir -p "$STAR_INDEX_ROOT" "$STAR_ALIGN_ROOT" 2>/dev/null || true
	
	# Get absolute paths using realpath if available, otherwise use cd/pwd
	if command -v realpath >/dev/null 2>&1; then
		abs_star_index_root="$(realpath -m "$STAR_INDEX_ROOT" 2>/dev/null)" || abs_star_index_root=""
		abs_star_align_root="$(realpath -m "$STAR_ALIGN_ROOT" 2>/dev/null)" || abs_star_align_root=""
	fi
	
	# Fallback: use cd/pwd method
	if [[ -z "$abs_star_index_root" ]]; then
		abs_star_index_root="$(cd "$STAR_INDEX_ROOT" 2>/dev/null && pwd)" || abs_star_index_root="$STAR_INDEX_ROOT"
	fi
	if [[ -z "$abs_star_align_root" ]]; then
		abs_star_align_root="$(cd "$STAR_ALIGN_ROOT" 2>/dev/null && pwd)" || abs_star_align_root="$STAR_ALIGN_ROOT"
	fi
	
	# If still relative, prepend PROJECT_ROOT
	[[ "$abs_star_index_root" != /* ]] && abs_star_index_root="${PROJECT_ROOT}/${abs_star_index_root}"
	[[ "$abs_star_align_root" != /* ]] && abs_star_align_root="${PROJECT_ROOT}/${abs_star_align_root}"
	
	# Clean up any double slashes
	abs_star_index_root="${abs_star_index_root//\/\//\/}"
	abs_star_align_root="${abs_star_align_root//\/\//\/}"
	
	# Set directories based on tissue tag (using absolute paths)
	# Structure: operation_type/fasta_tag[_tissue_tag]/
	# This keeps all alignments together, all salmon_quant together, etc.
	local ref_tag="${fasta_tag}"
	[[ -n "$tissue_tag" ]] && ref_tag="${fasta_tag}_${tissue_tag}"
	
	if [[ -n "$tissue_tag" ]]; then
		star_index_dir="${abs_star_index_root}/5_star/index/${ref_tag}"
		star_genome_dir="${abs_star_align_root}/5_star/alignments/${ref_tag}"
		log_info "[STAR] Tissue-specific alignment for: $tissue_tag (${#rnaseq_list[@]} samples)"
	else
		star_index_dir="${abs_star_index_root}/5_star/index/${fasta_tag}"
		star_genome_dir="${abs_star_align_root}/5_star/alignments/${fasta_tag}"
		log_info "[STAR] Pooled alignment: ${#rnaseq_list[@]} samples"
	fi
	
	# Clean up any double slashes in final paths
	star_index_dir="${star_index_dir//\/\//\/}"
	star_genome_dir="${star_genome_dir//\/\//\/}"
	
	# Log resolved paths for debugging
	log_info "[STAR] Index directory (absolute): $star_index_dir"
	log_info "[STAR] Output directory (absolute): $star_genome_dir"
	
	local star_overhang=$((STAR_READ_LENGTH - 1))
	log_info "[STAR] CPU allocation: Total=$THREADS threads"

	# STEP 1: BUILD STAR GENOME INDEX
	log_step "STAR genome index generation for $fasta_tag"
	
	if [[ -f "$star_index_dir/SAindex" ]]; then
		log_info "[STAR INDEX] Index for $fasta_tag already exists. Skipping."
	else
		log_info "[STAR INDEX] Building STAR genome index from transcriptome..."
		mkdir -p "$star_index_dir"
		
		log_file_size "$fasta" "Input transcriptome FASTA for STAR indexing"
		
		# Validate GTF file exists
		if [[ ! -f "$STAR_GTF_FILE" ]]; then
			log_error "GTF annotation file not found: $STAR_GTF_FILE"
			log_error "Set STAR_GTF_FILE to a valid GTF file path"
			return 1
		fi
		log_info "[STAR INDEX] Using GTF annotation: $STAR_GTF_FILE"
		
		run_with_space_time_log --input "$fasta" --output "$star_index_dir" \
			STAR --runMode genomeGenerate \
				--genomeDir "$star_index_dir" \
				--genomeFastaFiles "$fasta" \
				--sjdbGTFfile "$STAR_GTF_FILE" \
				--sjdbOverhang "$star_overhang" \
				--genomeSAindexNbases 12 \
				--runThreadN "$THREADS"
		
		[[ ! -f "$star_index_dir/SAindex" ]] && { log_error "STAR index generation failed"; return 1; }
		log_info "[STAR INDEX] Index built successfully: $star_index_dir"
	fi
	
	# STEP 2: ALIGN READS WITH STAR (2-PASS MODE)
	mkdir -p "$star_genome_dir"
	
	# Verify directory exists and is writable
	if [[ ! -d "$star_genome_dir" ]]; then
		log_error "Failed to create STAR output directory: $star_genome_dir"
		return 1
	fi
	if [[ ! -w "$star_genome_dir" ]]; then
		log_error "STAR output directory is not writable: $star_genome_dir"
		return 1
	fi
	log_info "[STAR] Output directory verified: $star_genome_dir"
	
	log_step "STAR splice-aware alignment for $fasta_tag samples"
	
	for SRR in "${rnaseq_list[@]}"; do
		local bam_output="$star_genome_dir/${SRR}_Aligned.sortedByCoord.out.bam"
		
		# Check if BAM exists AND has content (not 0 bytes from failed run)
		if [[ -f "$bam_output" ]]; then
			local bam_size
			bam_size=$(stat -c%s "$bam_output" 2>/dev/null || stat -f%z "$bam_output" 2>/dev/null || echo "0")
			if [[ "$bam_size" -gt 1000 ]]; then
				log_info "[STAR] Alignment for $SRR already exists (${bam_size} bytes). Skipping."
				continue
			else
				log_warn "[STAR] Found empty/corrupt BAM for $SRR (${bam_size} bytes) - removing and re-running"
				rm -f "$bam_output"
			fi
		fi
		
		# Clean up any stale files from previous failed STAR runs for this sample
		log_info "[STAR] Cleaning up stale files for $SRR..."
		rm -f "${star_genome_dir}/${SRR}_Log.out" 2>/dev/null || true
		rm -f "${star_genome_dir}/${SRR}_Log.progress.out" 2>/dev/null || true
		rm -f "${star_genome_dir}/${SRR}_Log.final.out" 2>/dev/null || true
		rm -f "${star_genome_dir}/${SRR}_SJ.out.tab" 2>/dev/null || true
		rm -rf "${star_genome_dir}/${SRR}__STARgenome" 2>/dev/null || true
		rm -rf "${star_genome_dir}/${SRR}__STARpass1" 2>/dev/null || true
		rm -rf "${star_genome_dir}/${SRR}_STARtmp" 2>/dev/null || true
		
		find_trimmed_fastq "$SRR"
		[[ -z "$trimmed1" ]] && { log_warn "Trimmed FASTQ for $SRR not found; skipping."; continue; }
		
		log_info "[STAR] Aligning: $SRR"
		
		local star_reads=""
		if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
			star_reads="$trimmed1 $trimmed2"
			log_info "[STAR] Processing paired-end reads for $SRR"
		else
			star_reads="$trimmed1"
			log_info "[STAR] Processing single-end reads for $SRR"
		fi
		
		# Clean up any stale STAR temporary files from previous failed runs
		rm -rf "${star_genome_dir}/${SRR}_STARtmp" 2>/dev/null || true
		rm -rf "${star_genome_dir}/${SRR}_"*.tmp 2>/dev/null || true
		rm -rf "${star_genome_dir}/_STARtmp_${SRR}" 2>/dev/null || true
		# Also clean up any _STARtmp in project root (STAR default location)
		rm -rf "${PROJECT_ROOT}/_STARtmp" 2>/dev/null || true
		
		# CRITICAL: Ensure output directory exists and is writable RIGHT BEFORE STAR runs
		# This fixes "could not create output file" errors on HPC/server environments
		mkdir -p "$star_genome_dir"
		sync  # Force filesystem sync on HPC/NFS
		sleep 1  # Brief pause for filesystem consistency
		
		if [[ ! -d "$star_genome_dir" ]]; then
			log_error "[STAR] FATAL: Cannot create output directory: $star_genome_dir"
			return 1
		fi
		
		# Test write permissions by creating a test file with the exact output name pattern
		local test_file="${star_genome_dir}/${SRR}_test_write_$$"
		if ! touch "$test_file" 2>/dev/null; then
			log_error "[STAR] FATAL: Cannot create output file: $test_file"
			log_error "[STAR] Check disk space and permissions for: $star_genome_dir"
			return 1
		fi
		rm -f "$test_file"
		
		# Configure STAR temp directory - ALWAYS use output directory for temp
		# This ensures all STAR operations happen on the same filesystem
		local star_tmp_dir="${star_genome_dir}/_STARtmp_${SRR}"
		rm -rf "$star_tmp_dir" 2>/dev/null || true
		# Also clean up any stray temp dirs
		rm -rf "${PROJECT_ROOT}/_STARtmp" 2>/dev/null || true
		rm -rf "${PROJECT_ROOT}/_STARtmp_${SRR}" 2>/dev/null || true
		local star_tmp_arg="--outTmpDir $star_tmp_dir"
		
		log_error "[STAR DEBUG] Using temp directory: $star_tmp_dir"
		
		# Construct output prefix - ensure no double slashes and path is clean
		local out_prefix="${star_genome_dir}/${SRR}_"
		# Remove any double slashes
		out_prefix="${out_prefix//\/\//\/}"
		
		# ======================================================================
		# DEBUG INFO - Output to error log for troubleshooting
		# ======================================================================
		log_error "[STAR DEBUG] ============================================"
		log_error "[STAR DEBUG] Sample: $SRR"
		log_error "[STAR DEBUG] Output prefix: $out_prefix"
		log_error "[STAR DEBUG] Expected BAM: ${out_prefix}Aligned.sortedByCoord.out.bam"
		log_error "[STAR DEBUG] Output directory: $star_genome_dir"
		log_error "[STAR DEBUG] Index directory: $star_index_dir"
		log_error "[STAR DEBUG] Temp mode: ${STAR_TMP_MODE:-cwd}"
		[[ -n "$star_tmp_dir" ]] && log_error "[STAR DEBUG] Temp directory: $star_tmp_dir"
		log_error "[STAR DEBUG] Input reads: $star_reads"
		log_error "[STAR DEBUG] Current working dir: $(pwd)"
		log_error "[STAR DEBUG] PROJECT_ROOT: $PROJECT_ROOT"
		
		# Check and log disk space
		local available_space disk_info
		disk_info=$(df -Ph "$star_genome_dir" 2>&1) || disk_info="df command failed"
		log_error "[STAR DEBUG] Disk info for output dir:"
		log_error "$disk_info"
		
		available_space=$(df -P "$star_genome_dir" 2>/dev/null | awk 'NR==2 {print $4}')
		if [[ -n "$available_space" ]]; then
			local available_gb=$((available_space / 1024 / 1024))
			log_error "[STAR DEBUG] Available disk space: ${available_gb} GB"
			if [[ $available_gb -lt 10 ]]; then
				log_error "[STAR] FATAL: Less than 10GB available - STAR needs more disk space"
				return 1
			fi
		fi
		
		# Verify paths exist and are accessible
		log_error "[STAR DEBUG] Checking paths..."
		log_error "[STAR DEBUG] star_genome_dir exists: $([[ -d "$star_genome_dir" ]] && echo YES || echo NO)"
		log_error "[STAR DEBUG] star_genome_dir writable: $([[ -w "$star_genome_dir" ]] && echo YES || echo NO)"
		log_error "[STAR DEBUG] star_index_dir exists: $([[ -d "$star_index_dir" ]] && echo YES || echo NO)"
		log_error "[STAR DEBUG] Index SAindex exists: $([[ -f "$star_index_dir/SAindex" ]] && echo YES || echo NO)"
		
		# List output directory contents
		log_error "[STAR DEBUG] Output dir contents:"
		ls -la "$star_genome_dir" 2>&1 | head -20 | while read line; do log_error "  $line"; done
		
		# Test creating the exact BAM filename
		local test_bam="${out_prefix}Aligned.sortedByCoord.out.bam.test"
		if touch "$test_bam" 2>/dev/null; then
			log_error "[STAR DEBUG] Test BAM creation: SUCCESS"
			rm -f "$test_bam"
		else
			log_error "[STAR DEBUG] Test BAM creation: FAILED - cannot create $test_bam"
		fi
		
		log_error "[STAR DEBUG] ============================================"
		# ======================================================================
		
		# Run STAR alignment - output UNSORTED BAM first, then sort with samtools
		# This is more reliable than STAR's internal BAM sorting which can fail
		local unsorted_bam="${out_prefix}Aligned.out.bam"
		
		log_error "[STAR DEBUG] Running STAR with unsorted BAM output..."
		log_error "[STAR DEBUG] Unsorted BAM will be: $unsorted_bam"
		log_error "[STAR DEBUG] Will sort to: $bam_output"
		
		run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$star_genome_dir" \
			STAR --runMode alignReads \
				--genomeDir "$star_index_dir" \
				--readFilesIn $star_reads \
				--readFilesCommand zcat \
				--outFileNamePrefix "$out_prefix" \
				$star_tmp_arg \
				--outSAMtype BAM Unsorted \
				--outSAMstrandField "$STAR_STRAND_SPECIFIC" \
				--outSAMunmapped Within \
				--outFilterType BySJout \
				--outFilterMultimapNmax 20 \
				--alignSJoverhangMin 8 \
				--alignSJDBoverhangMin 1 \
				--outFilterMismatchNmax 999 \
				--outFilterMismatchNoverReadLmax 0.04 \
				--alignIntronMin 20 \
				--alignIntronMax 1000000 \
				--alignMatesGapMax 1000000 \
				--twopassMode Basic \
				--runThreadN "$THREADS"
		
		# Check if unsorted BAM was created
		if [[ ! -f "$unsorted_bam" ]]; then
			log_error "[STAR] FATAL: Unsorted BAM file not created for $SRR"
			log_error "[STAR] Check STAR log: ${out_prefix}Log.out"
			# Show last 30 lines of STAR log
			log_error "[STAR] Last 30 lines of STAR log:"
			tail -30 "${out_prefix}Log.out" 2>/dev/null | while read line; do log_error "  $line"; done
			return 1
		fi
		
		local unsorted_size
		unsorted_size=$(stat -c%s "$unsorted_bam" 2>/dev/null || stat -f%z "$unsorted_bam" 2>/dev/null || echo "0")
		log_error "[STAR DEBUG] Unsorted BAM size: $unsorted_size bytes"
		
		if [[ "$unsorted_size" -lt 1000 ]]; then
			log_error "[STAR] FATAL: Unsorted BAM is empty/corrupt for $SRR (${unsorted_size} bytes)"
			log_error "[STAR] Last 30 lines of STAR log:"
			tail -30 "${out_prefix}Log.out" 2>/dev/null | while read line; do log_error "  $line"; done
			return 1
		fi
		
		# Sort BAM with samtools (more reliable than STAR's internal sorting)
		log_info "[STAR] Sorting BAM with samtools..."
		log_error "[STAR DEBUG] Running samtools sort..."
		
		if ! samtools sort -@ "$THREADS" -m 2G -o "$bam_output" "$unsorted_bam" 2>&1; then
			log_error "[STAR] FATAL: samtools sort failed for $SRR"
			return 1
		fi
		
		# Verify sorted BAM
		local final_bam_size
		final_bam_size=$(stat -c%s "$bam_output" 2>/dev/null || stat -f%z "$bam_output" 2>/dev/null || echo "0")
		if [[ "$final_bam_size" -lt 1000 ]]; then
			log_error "[STAR] FATAL: Sorted BAM file is empty/corrupt for $SRR (${final_bam_size} bytes)"
			return 1
		fi
		
		log_info "[STAR] BAM sorted successfully: $final_bam_size bytes"
		
		# Remove unsorted BAM to save space
		rm -f "$unsorted_bam"
		log_info "[STAR] Removed unsorted BAM to save space"
		
		# Index the BAM
		log_info "[STAR] Indexing BAM..."
		samtools index -@ "$THREADS" "$bam_output" 2>&1 || log_warn "[STAR] BAM indexing failed (non-fatal)"
		
		# Clean up temp directories after successful alignment
		[[ -n "$star_tmp_dir" ]] && rm -rf "$star_tmp_dir" 2>/dev/null || true
		rm -rf "${PROJECT_ROOT}/_STARtmp" 2>/dev/null || true
		
		# Delete transient big files if enabled (saves significant disk space)
		if [[ "${STAR_DELETE_TRANSIENT:-true}" == "true" ]]; then
			log_info "[STAR] Cleaning up transient files for $SRR..."
			
			# Remove 2-pass intermediate directories (can be several GB each)
			rm -rf "${star_genome_dir}/${SRR}__STARgenome" 2>/dev/null || true
			rm -rf "${star_genome_dir}/${SRR}__STARpass1" 2>/dev/null || true
			
			# Remove any leftover temp directories
			rm -rf "${star_genome_dir}/${SRR}_STARtmp" 2>/dev/null || true
			rm -rf "${star_genome_dir}/_STARtmp_${SRR}" 2>/dev/null || true
			
			# Remove progress log (Log.out and Log.final.out are kept for diagnostics)
			rm -f "${star_genome_dir}/${SRR}_Log.progress.out" 2>/dev/null || true
			
			log_info "[STAR] Transient files cleaned up for $SRR"
		fi
		
		log_info "[STAR] Successfully aligned: $SRR"
	done
	
	log_info "[STAR] All samples aligned successfully"
	
	# STEP 3: SALMON QUANTIFICATION
	# Structure: 6_salmon/index/fasta_tag[_tissue_tag]/ and 6_salmon/quant/fasta_tag[_tissue_tag]/SRR/
	local ref_suffix="${fasta_tag}"
	[[ -n "${tissue_tag:-}" ]] && ref_suffix="${fasta_tag}_${tissue_tag}"
	local salmon_idx="${abs_star_align_root}/6_salmon/index/${ref_suffix}"
	local quant_root="${abs_star_align_root}/6_salmon/quant/${ref_suffix}"
	# Clean up double slashes
	salmon_idx="${salmon_idx//\/\//\/}"
	quant_root="${quant_root//\/\//\/}"
	mkdir -p "$quant_root"
	
	# Build Salmon index
	if [[ -f "$salmon_idx/versionInfo.json" ]]; then
		log_info "[SALMON INDEX] STAR Salmon index exists - skipping build"
	else
		log_step "Building Salmon index for transcriptome"
		run_with_space_time_log --input "$fasta" --output "$salmon_idx" \
			salmon index -t "$fasta" -i "$salmon_idx" -k 31 --threads "$THREADS"
	fi
	
	# Quantify samples
	log_info "[SALMON] Starting quantification for ${#rnaseq_list[@]} samples"
	
	for SRR in "${rnaseq_list[@]}"; do
		local quant_dir="$quant_root/$SRR"
		
		[[ -f "$quant_dir/quant.sf" ]] && { log_info "[SALMON] Quantification for $SRR exists - skipping"; continue; }
		
		find_trimmed_fastq "$SRR"
		[[ -z "$trimmed1" ]] && { log_warn "No trimmed reads for $SRR - skipping"; continue; }
		
		log_info "[SALMON] Quantifying: $SRR"
		mkdir -p "$quant_dir"
		
		if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
			run_with_space_time_log salmon quant -p "$THREADS" -i "$salmon_idx" -o "$quant_dir" \
				--validateMappings -l A -1 "$trimmed1" -2 "$trimmed2"
		else
			run_with_space_time_log salmon quant -p "$THREADS" -i "$salmon_idx" -o "$quant_dir" \
				--validateMappings -l A -r "$trimmed1"
		fi
		
		[[ -f "$quant_dir/quant.sf" ]] && log_info "[SALMON] Successfully quantified: $SRR"
	done

	# STEP 4: PREPARE TXIMPORT FILES FOR DESEQ2
	log_step "Preparing tximport input for DESeq2 (STAR + Salmon)"
	
	local method3_root="${STAR_ALIGN_ROOT%/*}"
	local matrix_dir="$method3_root/6_matrices_from_stringtie"
	mkdir -p "$matrix_dir"
	
	local sample_metadata="$matrix_dir/sample_info.tsv"
	local tx2gene_file="$matrix_dir/tx2gene_${fasta_tag}${tag_suffix}.tsv"
	
	# Create tx2gene mapping
	if [[ ! -f "$tx2gene_file" ]]; then
		log_info "[TXIMPORT] Creating transcript-to-gene mapping"
		awk '/^>/{tx=$1; gsub(/^>/, "", tx); gene=tx; gsub(/\.[0-9]+$/, "", gene); print tx "\t" gene}' "$fasta" > "$tx2gene_file"
		log_info "[TXIMPORT] Created tx2gene mapping: $(wc -l < "$tx2gene_file") transcripts"
	fi
	
	# Verify quantifications
	local quant_count=0
	for SRR in "${rnaseq_list[@]}"; do
		[[ -f "$quant_root/$SRR/quant.sf" ]] && ((quant_count++))
	done
	
	[[ $quant_count -lt 2 ]] && { log_error "Insufficient quantifications: $quant_count (need ≥2)"; return 1; }
	
	# Create sample metadata
	[[ ! -f "$sample_metadata" ]] && create_sample_metadata "$sample_metadata" rnaseq_list[@]
	
	# Generate tximport R script
	local tximport_script="$matrix_dir/run_tximport_star_salmon.R"
	if [[ ! -f "$tximport_script" ]]; then
		generate_tximport_star_script "$quant_root" "$sample_metadata" "$tx2gene_file" "$matrix_dir" "$tximport_script"
	fi
	
	# Run tximport if R is available
	if command -v Rscript >/dev/null 2>&1; then
		log_step "Running tximport to import Salmon quantifications"
		if Rscript "$tximport_script" "$quant_root" "$sample_metadata" "$tx2gene_file" "$matrix_dir" 2>&1 | tee -a "$LOG_FILE"; then
			log_info "[TXIMPORT] Successfully imported counts for DESeq2"
			local count_matrix="$matrix_dir/gene_counts_tximport.tsv"
			[[ -f "$count_matrix" ]] && validate_count_matrix "$count_matrix" "gene" 2
		else
			log_warn "[TXIMPORT] tximport failed - check R dependencies"
		fi
	else
		log_warn "[TXIMPORT] Rscript not found - run manually: Rscript $tximport_script"
	fi
	
	log_step "STAR alignment pipeline completed for $fasta_tag"
	log_info "STAR alignments: $star_genome_dir"
	log_info "Salmon quantifications: $quant_root"
	log_info "Tximport outputs: $matrix_dir"
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Run tximport for STAR+Salmon using external R helper
# Usage: run_tximport_star <quant_dir> <metadata_file> <tx2gene_file> [output_dir]
run_tximport_star() {
	local quant_dir="$1"
	local metadata_file="$2"
	local tx2gene_file="$3"
	local output_dir="${4:-$(dirname "$metadata_file")}"
	local helper_script="$SCRIPT_DIR/helpers/tximport_star_helper.R"
	
	if [[ ! -f "$helper_script" ]]; then
		log_error "tximport_star_helper.R not found: $helper_script"
		return 1
	fi
	
	log_info "[TXIMPORT] Running STAR+Salmon import..."
	Rscript "$helper_script" "$quant_dir" "$metadata_file" "$tx2gene_file" "$output_dir"
}

# Generate tximport script (legacy - copies helper)
generate_tximport_star_script() {
	local quant_root="$1"
	local sample_metadata="$2"
	local tx2gene_file="$3"
	local matrix_dir="$4"
	local output_script="$5"
	local helper_script="$SCRIPT_DIR/helpers/tximport_star_helper.R"
	
	if [[ -f "$helper_script" ]]; then
		cp "$helper_script" "$output_script"
		chmod +x "$output_script"
		log_info "[TXIMPORT] Copied helper to: $output_script"
	else
		log_error "tximport_star_helper.R not found: $helper_script"
		return 1
	fi
}
