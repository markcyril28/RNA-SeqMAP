#!/bin/bash
# ==============================================================================
# METHOD 3: STAR SPLICE-AWARE ALIGNMENT PIPELINE
# ==============================================================================
# STAR splice-aware alignment with Salmon quantification and tximport for DESeq2
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${M3_STAR_SOURCED:-}" == "true" ]] && return 0
export M3_STAR_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/shared_utils_method.sh"

# ==============================================================================
# STAR CONFIGURATION - IMPORTANT PARAMETERS (tweak here)
# ==============================================================================

# Genome loading mode: NoSharedMemory (default/safe), LoadAndKeep (multi-run HPC)
STAR_GENOME_LOAD="${STAR_GENOME_LOAD:-$(get_star_genome_load 2>/dev/null || echo NoSharedMemory)}"

# sjdbOverhang = read_length - 1; set to actual read length for best splice detection
STAR_READ_LENGTH="${STAR_READ_LENGTH:-100}"

# Delete intermediate STAR files after alignment to save disk space (true/false)
# Set to "false" to keep 2-pass genome, unsorted BAM, etc. for debugging
STAR_DELETE_TRANSIENT="${STAR_DELETE_TRANSIENT:-true}"

# GTF annotation file for splice junction detection (resolved at runtime in pipeline)
# NOTE: STAR_GTF_FILE is resolved inside star_alignment_pipeline() from $gtf_file

# Transcriptome FASTA for Salmon quantification step (auto-detected if not set)
# When using a genome FASTA for STAR, point this to the transcript-level FASTA.
# e.g. STAR_TRANSCRIPTOME_FASTA="0_INPUTs/fasta/reference_genomes/GPE001970_transcripts.fa"
# Auto-detect: looks for <genome_basename>_transcripts.fa; falls back to --FASTA if absent.

# STAR temp directory: "system"=/tmp, "local"=output dir, "cwd"=current dir, "none"=STAR default
STAR_TMP_MODE="${STAR_TMP_MODE:-cwd}"

# ==============================================================================
# MAIN STAR ALIGNMENT PIPELINE
# ==============================================================================

star_alignment_pipeline() {
	# Resolve GTF at runtime so the config's gtf_file is available
	STAR_GTF_FILE="${STAR_GTF_FILE:-$gtf_file}"

	local fasta="" transcriptome_fasta="" rnaseq_list=() tissue_tag=""

	# Check for GNU parallel
	if ! command -v parallel >/dev/null 2>&1; then
		log_warn "[PERFORMANCE] GNU parallel not found - Salmon quantification will run sequentially"
	fi

	while [[ $# -gt 0 ]]; do
		case "$1" in
			--FASTA) fasta="$2"; shift 2;;
			--TRANSCRIPTOME) transcriptome_fasta="$2"; shift 2;;
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

	# Resolve transcriptome FASTA for Salmon: explicit arg > config var > auto-detect > fallback to genome
	if [[ -z "$transcriptome_fasta" ]]; then
		transcriptome_fasta="${STAR_TRANSCRIPTOME_FASTA:-}"
	fi
	if [[ -z "$transcriptome_fasta" ]]; then
		# Auto-detect: look for <basename>_transcripts.fa alongside the genome FASTA
		local auto_tx="$(dirname "$fasta")/$(basename "${fasta%.*}")_transcripts.fa"
		if [[ -f "$auto_tx" ]]; then
			transcriptome_fasta="$auto_tx"
			log_info "[STAR] Auto-detected transcriptome FASTA: $transcriptome_fasta"
		else
			transcriptome_fasta="$fasta"
			log_warn "[STAR] No transcriptome FASTA found. Salmon will index the genome FASTA."
			log_warn "[STAR] Set STAR_TRANSCRIPTOME_FASTA or use --TRANSCRIPTOME to fix this."
		fi
	else
		if [[ ! -f "$transcriptome_fasta" ]]; then
			log_error "Transcriptome FASTA '$transcriptome_fasta' not found."
			return 1
		fi
		log_info "[STAR] Using transcriptome FASTA for Salmon: $transcriptome_fasta"
	fi
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
	# fasta_tag is already embedded in STAR_ALIGN_ROOT/STAR_INDEX_ROOT via set_fasta_output_dirs
	# Only append tissue_tag subdirectory when doing tissue-specific runs
	if [[ -n "$tissue_tag" ]]; then
		star_index_dir="${abs_star_index_root}/5_star/index/${tissue_tag}"
		star_genome_dir="${abs_star_align_root}/5_star/alignments/${tissue_tag}"
		log_info "[STAR] Tissue-specific alignment for: $tissue_tag (${#rnaseq_list[@]} samples)"
	else
		star_index_dir="${abs_star_index_root}/5_star/index"
		star_genome_dir="${abs_star_align_root}/5_star/alignments"
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

	if [[ -f "$star_index_dir/SAindex" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[STAR INDEX] Index for $fasta_tag already exists. Skipping."
	else
		log_info "[STAR INDEX] Building STAR genome index from genome FASTA..."
		mkdir -p "$star_index_dir"

		log_file_size "$fasta" "Input genome FASTA for STAR indexing"

		# Validate GTF file exists
		if [[ ! -f "$STAR_GTF_FILE" ]]; then
			log_error "GTF annotation file not found: $STAR_GTF_FILE"
			log_error "Set STAR_GTF_FILE to a valid GTF file path"
			return 1
		fi
		log_info "[STAR INDEX] Using GTF annotation: $STAR_GTF_FILE"

		# Compute recommended SAindex size: min(14, floor(log2(GenomeLength)/2 - 1))
		local genome_sa_index=14
		local genome_chr_bin=18
		local genome_bp num_seqs
		genome_bp=$(awk '!/^>/{total+=length($0)} END{printf "%d", total+0}' "$fasta" 2>/dev/null)
		num_seqs=$(grep -c '^>' "$fasta" 2>/dev/null || echo "1")
		if [[ -n "$genome_bp" && "$genome_bp" -gt 0 ]]; then
			genome_sa_index=$(awk "BEGIN{v=int(log($genome_bp)/log(2)/2-1); print (v<14)?v:14}")
			[[ "$genome_sa_index" -lt 1 ]] && genome_sa_index=1
			# genomeChrBinNbits: min(18, floor(log2(GenomeLength/NumberOfSequences) - 1))
			# Required for genomes with many scaffolds (e.g., draft assemblies)
			genome_chr_bin=$(awk "BEGIN{v=int(log($genome_bp/$num_seqs)/log(2)-1); print (v<18)?v:18}")
			[[ "$genome_chr_bin" -lt 1 ]] && genome_chr_bin=1
		fi
		log_info "[STAR INDEX] Genome size: ${genome_bp:-unknown} bp, $num_seqs sequences -> genomeSAindexNbases=$genome_sa_index, genomeChrBinNbits=$genome_chr_bin"

		run_with_space_time_log --input "$fasta" --output "$star_index_dir" \
			STAR --runMode genomeGenerate \
				--genomeDir "$star_index_dir" \
				--genomeFastaFiles "$fasta" \
				--sjdbGTFfile "$STAR_GTF_FILE" \
				--sjdbOverhang "$star_overhang" \
				--genomeSAindexNbases "$genome_sa_index" \
				--genomeChrBinNbits "$genome_chr_bin" \
				--runThreadN "$THREADS"

		[[ ! -f "$star_index_dir/SAindex" ]] && { log_error "STAR index generation failed"; return 1; }
		log_info "[STAR INDEX] Index built successfully: $star_index_dir"
	fi

	# STEP 2: ALIGN READS WITH STAR
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

	local parallel_jobs="${PARALLEL_JOBS:-${JOBS:-2}}"
	local threads_per_job=$((THREADS / parallel_jobs))
	[[ $threads_per_job -lt 1 ]] && threads_per_job=1

	if command -v parallel >/dev/null 2>&1 && [[ "$parallel_jobs" -gt 1 ]]; then
		log_step "[PARALLEL] STAR alignment: ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		log_warn "[PARALLEL] STAR is memory-intensive (~30GB/instance). Ensure sufficient RAM for $parallel_jobs concurrent jobs."
		_prepare_parallel_env

		export fasta_tag threads_per_job
		export star_index_dir star_genome_dir
		export STAR_DELETE_TRANSIENT PROJECT_ROOT

		_m3_star_parallel_worker() {
			local SRR="$1"
			_init_parallel_worker "$SRR"
			[[ -z "$trimmed1" ]] && { _parallel_log STAR "$SRR" WARN "Trimmed FASTQ not found - skipping"; return 0; }

			local bam_output="$star_genome_dir/${SRR}_Aligned.sortedByCoord.out.bam"

			# Check if BAM exists and has content
			if [[ -f "$bam_output" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
				local bam_size
				bam_size=$(stat -c%s "$bam_output" 2>/dev/null || stat -f%z "$bam_output" 2>/dev/null || echo "0")
				if [[ "$bam_size" -gt 1000 ]]; then
					_parallel_log STAR "$SRR" INFO "BAM exists (${bam_size} bytes) - skipping"
					return 0
				fi
				rm -f "$bam_output"
			fi

			# Clean up stale files from previous runs
			rm -f "${star_genome_dir}/${SRR}_Log.out" "${star_genome_dir}/${SRR}_Log.progress.out" \
				"${star_genome_dir}/${SRR}_Log.final.out" "${star_genome_dir}/${SRR}_SJ.out.tab" 2>/dev/null || true
			rm -rf "${star_genome_dir}/${SRR}__STARgenome" "${star_genome_dir}/${SRR}__STARpass1" \
				"${star_genome_dir}/${SRR}_STARtmp" "${star_genome_dir}/_STARtmp_${SRR}" 2>/dev/null || true

			local star_reads_args=("$trimmed1")
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				star_reads_args+=("$trimmed2")
			fi

			local star_tmp_dir="${star_genome_dir}/_STARtmp_${SRR}"
			rm -rf "$star_tmp_dir" 2>/dev/null || true

			local out_prefix="${star_genome_dir}/${SRR}_"
			out_prefix="${out_prefix//\/\//\/}"
			local unsorted_bam="${out_prefix}Aligned.out.bam"

			_parallel_log STAR "$SRR" INFO "Aligning with $threads_per_job threads"
			_parallel_log STAR "$SRR" INFO "--- BEGIN STAR OUTPUT ---"

			STAR --runMode alignReads \
				--genomeDir "$star_index_dir" \
				--readFilesIn "${star_reads_args[@]}" \
				--readFilesCommand "zcat -f" \
				--outFileNamePrefix "$out_prefix" \
				--outTmpDir "$star_tmp_dir" \
				--outSAMtype BAM Unsorted \
				--runThreadN "$threads_per_job" 2>&1 || \
				{ _parallel_log STAR "$SRR" ERROR "--- END STAR OUTPUT (FAILED) ---"; _parallel_log STAR "$SRR" ERROR "STAR alignment failed"; return 1; }

			_parallel_log STAR "$SRR" INFO "--- END STAR OUTPUT ---"

			if [[ ! -f "$unsorted_bam" ]]; then
				_parallel_log STAR "$SRR" ERROR "Unsorted BAM not created"
				return 1
			fi

			local unsorted_size
			unsorted_size=$(stat -c%s "$unsorted_bam" 2>/dev/null || stat -f%z "$unsorted_bam" 2>/dev/null || echo "0")
			if [[ "$unsorted_size" -lt 1000 ]]; then
				_parallel_log STAR "$SRR" ERROR "Unsorted BAM is empty/corrupt (${unsorted_size} bytes)"
				return 1
			fi

			_parallel_log STAR "$SRR" INFO "Sorting BAM with samtools"
			samtools sort -@ "$threads_per_job" -m 2G -o "$bam_output" "$unsorted_bam" 2>&1 | sed 's/^/\t/'
			local sort_exit=${PIPESTATUS[0]}
			[[ $sort_exit -ne 0 ]] && { _parallel_log STAR "$SRR" ERROR "samtools sort failed"; return 1; }

			local final_bam_size
			final_bam_size=$(stat -c%s "$bam_output" 2>/dev/null || stat -f%z "$bam_output" 2>/dev/null || echo "0")
			if [[ "$final_bam_size" -lt 1000 ]]; then
				_parallel_log STAR "$SRR" ERROR "Sorted BAM is empty/corrupt (${final_bam_size} bytes)"
				return 1
			fi

			rm -f "$unsorted_bam"
			samtools index -@ "$threads_per_job" "$bam_output" 2>&1 || true

			# Clean up transient files
			rm -rf "$star_tmp_dir" 2>/dev/null || true
			if [[ "${STAR_DELETE_TRANSIENT:-true}" == "true" ]]; then
				rm -rf "${star_genome_dir}/${SRR}__STARgenome" "${star_genome_dir}/${SRR}__STARpass1" \
					"${star_genome_dir}/${SRR}_STARtmp" "${star_genome_dir}/_STARtmp_${SRR}" 2>/dev/null || true
				rm -f "${star_genome_dir}/${SRR}_Log.progress.out" 2>/dev/null || true
			fi

			_parallel_log STAR "$SRR" INFO "Completed successfully"
			return 0
		}
		export -f _m3_star_parallel_worker

		printf '%s\n' "${rnaseq_list[@]}" | parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env abs_trim_dir_root --env abs_error_warn_file --env keep_bam_global \
			--env fasta_tag --env threads_per_job \
			--env star_index_dir --env star_genome_dir \
			--env STAR_DELETE_TRANSIENT --env PROJECT_ROOT \
			--env OVERWRITE_MODE \
			-j "$parallel_jobs" \
			--halt soon,fail=1 \
			--joblog "$star_genome_dir/parallel_star_align.log" \
			_m3_star_parallel_worker {}

		local par_exit=$?
		log_info "[PARALLEL] STAR alignment complete (exit=$par_exit)"
		[[ $par_exit -ne 0 ]] && log_warn "[PARALLEL] Some jobs failed - check $star_genome_dir/parallel_star_align.log"
	else
		# Sequential fallback
		for SRR in "${rnaseq_list[@]}"; do
			local bam_output="$star_genome_dir/${SRR}_Aligned.sortedByCoord.out.bam"

			# Check if BAM exists AND has content (not 0 bytes from failed run)
			if [[ -f "$bam_output" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
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

			# Clean up stale files from previous failed STAR runs
			log_info "[STAR] Cleaning up stale files for $SRR..."
			rm -f "${star_genome_dir}/${SRR}_Log.out" "${star_genome_dir}/${SRR}_Log.progress.out" \
				"${star_genome_dir}/${SRR}_Log.final.out" "${star_genome_dir}/${SRR}_SJ.out.tab" 2>/dev/null || true
			rm -rf "${star_genome_dir}/${SRR}__STARgenome" "${star_genome_dir}/${SRR}__STARpass1" \
				"${star_genome_dir}/${SRR}_STARtmp" "${star_genome_dir}/${SRR}_"*.tmp \
				"${star_genome_dir}/_STARtmp_${SRR}" "${PROJECT_ROOT}/_STARtmp" 2>/dev/null || true

			find_trimmed_fastq "$SRR"
			[[ -z "$trimmed1" ]] && { log_warn "Trimmed FASTQ for $SRR not found; skipping."; continue; }

			log_info "[STAR] Aligning: $SRR"

			local star_reads_args=("$trimmed1")
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				star_reads_args+=("$trimmed2")
				log_info "[STAR] Processing paired-end reads for $SRR"
			else
				log_info "[STAR] Processing single-end reads for $SRR"
			fi

			# CRITICAL: Ensure output directory exists and is writable RIGHT BEFORE STAR runs
			# This fixes "could not create output file" errors on HPC/server environments
			mkdir -p "$star_genome_dir"
			# On NFS/HPC systems, sync and a brief pause help avoid race conditions.
			# Set STAR_NFS_SYNC=true in your config if running on NFS storage.
			if [[ "${STAR_NFS_SYNC:-false}" == "true" ]]; then
				sync
				sleep 1
			fi

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
			log_info "[STAR] Using temp directory: $star_tmp_dir"

			# Construct output prefix - ensure no double slashes and path is clean
			local out_prefix="${star_genome_dir}/${SRR}_"
			# Remove any double slashes
			out_prefix="${out_prefix//\/\//\/}"

			# Check disk space before running STAR (needs ~30GB per run)
			local available_space
			available_space=$(df -P "$star_genome_dir" 2>/dev/null | awk 'NR==2 {print $4}')
			if [[ -n "$available_space" ]]; then
				local available_gb=$((available_space / 1024 / 1024))
				if [[ $available_gb -lt 30 ]]; then
					log_error "[STAR] FATAL: Less than 30GB available in $star_genome_dir - STAR needs ~30GB per run"
					return 1
				fi
			fi

			# Run STAR alignment - output UNSORTED BAM first, then sort with samtools
			local unsorted_bam="${out_prefix}Aligned.out.bam"

			run_with_space_time_log --input "$TRIM_DIR_ROOT/$SRR" --output "$star_genome_dir" \
				STAR --runMode alignReads \
					--genomeDir "$star_index_dir" \
					--readFilesIn "${star_reads_args[@]}" \
					--readFilesCommand "zcat -f" \
					--outFileNamePrefix "$out_prefix" \
					--outTmpDir "$star_tmp_dir" \
					--outSAMtype BAM Unsorted \
					--runThreadN "$THREADS"

			# Check if unsorted BAM was created
			if [[ ! -f "$unsorted_bam" ]]; then
				log_error "[STAR] FATAL: Unsorted BAM file not created for $SRR"
				log_error "[STAR] Check STAR log: ${out_prefix}Log.out"
				# Show last 30 lines of STAR log
				log_error "[STAR] Last 30 lines of STAR log:"
				tail -30 "${out_prefix}Log.out" 2>/dev/null | while IFS= read -r line; do log_error "  $line"; done
				return 1
			fi

			local unsorted_size
			unsorted_size=$(stat -c%s "$unsorted_bam" 2>/dev/null || stat -f%z "$unsorted_bam" 2>/dev/null || echo "0")
			if [[ "$unsorted_size" -lt 1000 ]]; then
				log_error "[STAR] FATAL: Unsorted BAM is empty/corrupt for $SRR (${unsorted_size} bytes)"
				log_error "[STAR] Last 30 lines of STAR log:"
				tail -30 "${out_prefix}Log.out" 2>/dev/null | while IFS= read -r line; do log_error "  $line"; done
				return 1
			fi

			# Sort BAM with samtools
			log_info "[STAR] Sorting BAM with samtools..."

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
			rm -rf "$star_tmp_dir" "${PROJECT_ROOT}/_STARtmp" 2>/dev/null || true

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
	fi

	log_info "[STAR] All samples aligned successfully"

	# STEP 3: SALMON QUANTIFICATION
	# fasta_tag is already embedded in abs_star_align_root via set_fasta_output_dirs
	# Only append tissue_tag subdirectory when doing tissue-specific runs
	local ref_suffix=""
	[[ -n "${tissue_tag:-}" ]] && ref_suffix="${tissue_tag}"
	local salmon_idx="${abs_star_align_root}/6_salmon/index${ref_suffix:+/${ref_suffix}}"
	local quant_root="${abs_star_align_root}/6_salmon/quant${ref_suffix:+/${ref_suffix}}"
	# Clean up double slashes
	salmon_idx="${salmon_idx//\/\//\/}"
	quant_root="${quant_root//\/\//\/}"
	mkdir -p "$quant_root"

	# Build Salmon index (using transcriptome FASTA, not genome)
	if [[ -f "$salmon_idx/versionInfo.json" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "[SALMON INDEX] STAR Salmon index exists - skipping build"
	else
		log_step "Building Salmon index for transcriptome"
		log_info "[SALMON INDEX] Indexing: $transcriptome_fasta"
		run_with_space_time_log --input "$transcriptome_fasta" --output "$salmon_idx" \
			salmon index -t "$transcriptome_fasta" -i "$salmon_idx" -k 31 --threads "$THREADS"
	fi

	# Quantify samples
	log_info "[SALMON] Starting quantification for ${#rnaseq_list[@]} samples"

	if command -v parallel >/dev/null 2>&1 && [[ "$parallel_jobs" -gt 1 ]]; then
		log_step "[PARALLEL] Salmon quant (STAR): ${#rnaseq_list[@]} samples, $parallel_jobs jobs x $threads_per_job threads"
		_prepare_parallel_env

		export fasta_tag threads_per_job salmon_idx quant_root

		_m3_salmon_parallel_worker() {
			local SRR="$1"
			_init_parallel_worker "$SRR"
			[[ -z "$trimmed1" ]] && { _parallel_log SALMON_STAR "$SRR" WARN "Trimmed FASTQ not found - skipping"; return 0; }

			local quant_dir="$quant_root/$SRR"
			[[ -f "$quant_dir/quant.sf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]] && { _parallel_log SALMON_STAR "$SRR" INFO "quant.sf exists - skipping"; return 0; }

			mkdir -p "$quant_dir"
			_parallel_log SALMON_STAR "$SRR" INFO "Quantifying with $threads_per_job threads"

			local quant_exit=0
			# Strip ANSI escape codes and carriage returns (Salmon uses colored progress bars)
			if [[ -n "$trimmed2" && -f "$trimmed2" ]]; then
				salmon quant -p "$threads_per_job" -i "$salmon_idx" -o "$quant_dir" \
					--validateMappings -l A -1 "$trimmed1" -2 "$trimmed2" 2>&1 | \
					sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g'
				quant_exit=${PIPESTATUS[0]}
			else
				salmon quant -p "$threads_per_job" -i "$salmon_idx" -o "$quant_dir" \
					--validateMappings -l A -r "$trimmed1" 2>&1 | \
					sed 's/\x1B\[[0-9;]*[a-zA-Z]//g; s/\r//g'
				quant_exit=${PIPESTATUS[0]}
			fi
			[[ $quant_exit -ne 0 ]] && { _parallel_log SALMON_STAR "$SRR" ERROR "Salmon quant failed (exit=$quant_exit)"; return $quant_exit; }

			_parallel_log SALMON_STAR "$SRR" INFO "Completed successfully"
			return 0
		}
		export -f _m3_salmon_parallel_worker

		printf '%s\n' "${rnaseq_list[@]}" | parallel \
			--env PATH --env CONDA_PREFIX --env CONDA_DEFAULT_ENV --env CONDA_EXE \
			--env abs_trim_dir_root --env abs_error_warn_file --env keep_bam_global \
			--env fasta_tag --env threads_per_job --env salmon_idx --env quant_root \
			--env OVERWRITE_MODE \
			-j "$parallel_jobs" \
			--halt soon,fail=1 \
			--joblog "$quant_root/parallel_salmon_star_quant.log" \
			_m3_salmon_parallel_worker {}

		local par_exit_salmon=$?
		log_info "[PARALLEL] Salmon quant complete (exit=$par_exit_salmon)"
		[[ $par_exit_salmon -ne 0 ]] && log_warn "[PARALLEL] Some Salmon jobs failed - check $quant_root/parallel_salmon_star_quant.log"
	else
		# Sequential fallback
		for SRR in "${rnaseq_list[@]}"; do
			local quant_dir="$quant_root/$SRR"

			[[ -f "$quant_dir/quant.sf" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]] && { log_info "[SALMON] Quantification for $SRR exists - skipping"; continue; }

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
	fi

	# STEP 4: PREPARE TXIMPORT FILES FOR DESEQ2
	log_step "Preparing tximport input for DESeq2 (STAR + Salmon)"

	local matrix_dir="$STAR_MATRIX_ROOT"
	mkdir -p "$matrix_dir"

	local sample_metadata="$matrix_dir/sample_info.tsv"
	local tx2gene_file="$matrix_dir/tx2gene_${fasta_tag}${ref_suffix:+_${ref_suffix}}.tsv"

	# Create tx2gene mapping from GTF (transcript_id -> gene_id attributes)
	# This correctly handles multi-transcript genes; version-stripping FASTA headers is not reliable.
	if [[ ! -f "$tx2gene_file" ]]; then
		log_info "[TXIMPORT] Creating transcript-to-gene mapping from GTF: $STAR_GTF_FILE"
		awk '$3=="transcript" {
			match($0, /transcript_id "([^"]+)"/, t)
			match($0, /gene_id "([^"]+)"/, g)
			if (t[1] && g[1]) print t[1] "\t" g[1]
		}' "$STAR_GTF_FILE" | sort -u > "$tx2gene_file"
		local tx2gene_count
		tx2gene_count=$(wc -l < "$tx2gene_file")
		if [[ "$tx2gene_count" -eq 0 ]]; then
			log_error "[TXIMPORT] tx2gene mapping is empty - check GTF has 'transcript' features with transcript_id/gene_id attributes"
			return 1
		fi
		log_info "[TXIMPORT] Created tx2gene mapping: $tx2gene_count transcripts"
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
		# Note: stdout already goes through tee via exec redirect; do NOT pipe to tee -a "$LOG_FILE" (causes double-logging)
		if Rscript "$tximport_script" "$quant_root" "$sample_metadata" "$tx2gene_file" "$matrix_dir" 2>&1; then
			log_info "[TXIMPORT] Successfully imported counts for DESeq2"
			local count_matrix="$matrix_dir/gene_counts_tximport.tsv"
			[[ -f "$count_matrix" ]] && validate_count_matrix "$count_matrix" "gene" 2
		else
			log_warn "[TXIMPORT] tximport failed - check R dependencies"
			log_warn "[TXIMPORT] Install missing packages: Rscript -e \"BiocManager::install('tximport')\""
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

# ==============================================================================
# ALTERNATIVE FLOW: TISSUE-SPECIFIC ALIGNMENT
# ==============================================================================
# Splits samples by condition/tissue and runs star_alignment_pipeline() per group.
# Falls back to pooled alignment when only one tissue type is detected.

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
