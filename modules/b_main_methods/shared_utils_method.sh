#!/bin/bash
# ==============================================================================
# MAIN METHODS SHARED UTILITIES
# ==============================================================================
# Common helper functions used across all GEA analysis methods
# Sourced by: All method scripts in b_main_methods/
# ==============================================================================

#set -euo pipefail

# Guard against double-sourcing
[[ "${METHOD_SHARED_SOURCED:-}" == "true" ]] && return 0
export METHOD_SHARED_SOURCED="true"

# Source dependencies
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/global_config_method.sh"
source "$SCRIPT_DIR/../logging/logging_utils.sh"
source "$SCRIPT_DIR/../a_preprocessing/shared_utils_preproc.sh"

# ==============================================================================
# SHARED COMPRESSION DETECTION
# ==============================================================================
# Detect pigz once at module load; methods use $_SHARED_GZIP_DC instead of
# repeatedly spawning `command -v pigz` per sample.

if command -v pigz &>/dev/null; then
	_SHARED_HAS_PIGZ="true"
	_SHARED_GZIP_DC="pigz -dc"
	_SHARED_GZIP_C="pigz"
else
	_SHARED_HAS_PIGZ="false"
	_SHARED_GZIP_DC="gzip -dc"
	_SHARED_GZIP_C="gzip"
fi
export _SHARED_HAS_PIGZ _SHARED_GZIP_DC _SHARED_GZIP_C

# ==============================================================================
# VALIDATION FUNCTIONS
# ==============================================================================

# Validate count matrix quality
# Single-pass AWK replaces 4 separate file reads (head, tail|wc, tail|cut|sort|wc)
validate_count_matrix() {
	local matrix="$1"
	local matrix_type="${2:-gene}"
	local min_samples="${3:-2}"

	[[ ! -f "$matrix" ]] && { log_error "Matrix not found: $matrix"; return 1; }

	# Handle empty files gracefully
	[[ ! -s "$matrix" ]] && { log_error "Matrix file is empty: $matrix"; return 1; }

	log_step "Validating $matrix_type matrix: $(basename "$matrix")"

	local delim=$'\t'
	[[ "$matrix" == *.csv ]] && delim=","

	# Single AWK pass: count samples from header, count total rows, count unique IDs
	local num_samples total_rows unique_ids
	read -r num_samples total_rows unique_ids < <(awk -F"$delim" '
		NR == 1 { samples = NF - 1 }
		NR > 1  { total++; ids[$1]++ }
		END     { print samples+0, total+0, length(ids) }
	' "$matrix")

	[[ $num_samples -lt $min_samples ]] && { log_error "Insufficient samples: $num_samples (need ≥$min_samples)"; return 1; }

	[[ $total_rows -ne $unique_ids ]] && { log_error "Duplicate ${matrix_type} IDs detected!"; return 1; }

	log_info "[VALIDATION] Passed - Samples: $num_samples, ${matrix_type^}s: $unique_ids"
	return 0
}

# ==============================================================================
# METADATA FUNCTIONS
# ==============================================================================

# Load sample metadata from file
load_sample_metadata() {
	local metadata_file="${1:-${SAMPLE_CONDITIONS_FILE:-sample_conditions.txt}}"
	local -n metadata_array=$2
	
	[[ ! -f "$metadata_file" ]] && { log_warn "Metadata file not found: $metadata_file"; return 1; }
	
	# Count valid sample lines (exclude comments, blanks, and header) — single pass
	local line_count
	line_count=$(awk '!/^#/ && !/^$/ && !/^SRR_ID/ {n++} END{print n+0}' "$metadata_file")
	[[ $line_count -lt 2 ]] && { log_error "Metadata file must contain at least 2 samples (found: $line_count)"; return 1; }
	
	while IFS=$'\t' read -r srr condition batch; do
		[[ -z "$srr" || "$srr" =~ ^# || "$srr" == "SRR_ID" ]] && continue
		metadata_array["${srr}_condition"]="$condition"
		metadata_array["${srr}_batch"]="$batch"
	done < "$metadata_file"
	
	log_info "Loaded metadata: $(( ${#metadata_array[@]} / 2 )) samples from $metadata_file"
	return 0
}

# Create sample metadata CSV for DESeq2
create_sample_metadata() {
	local metadata_file="$1"
	local -a sample_list=("${@:2}")
	local metadata_source="${SAMPLE_CONDITIONS_FILE:-sample_conditions.txt}"
	
	local delim=","
	[[ "$metadata_file" == *.tsv ]] && delim=$'\t'
	
	declare -A sample_metadata
	local has_external=false
	load_sample_metadata "$metadata_source" sample_metadata 2>/dev/null && has_external=true
	
	# Build output in memory, write once (avoids N+1 file opens for N samples)
	{
		echo -e "sample${delim}condition${delim}batch"
		for SRR in "${sample_list[@]}"; do
			if [[ "$has_external" == "true" ]]; then
				local condition="${sample_metadata[${SRR}_condition]:-unknown}"
				local batch="${sample_metadata[${SRR}_batch]:-1}"
			else
				local condition="treatment"
				local batch="1"
			fi
			echo -e "$SRR${delim}$condition${delim}$batch"
		done
	} > "$metadata_file"
	
	[[ "$has_external" == "false" ]] && \
		log_warn "Sample conditions need manual specification in: $metadata_file"
}

# ==============================================================================
# TXIMPORT SCRIPT GENERATION
# ==============================================================================

# Get the helper scripts directory (tximport helpers live in c_post_processing/preprocessing)
HELPERS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../c_post_processing/preprocessing" && pwd)"

# Run tximport using external R helper script
# Usage: run_tximport <method> <quant_dir> <metadata_file> [output_dir]
run_tximport() {
	local method="$1"
	local quant_dir="$2"
	local metadata_file="$3"
	local output_dir="${4:-$(dirname "$metadata_file")}"
	local helper_script="$HELPERS_DIR/tximport_helper.R"
	
	if [[ ! -f "$helper_script" ]]; then
		log_error "tximport_helper.R not found: $helper_script"
		return 1
	fi
	
	log_info "[TXIMPORT] Running $method import..."
	Rscript "$helper_script" "$method" "$quant_dir" "$metadata_file" "$output_dir"
}

# Generate tximport R script (legacy compatibility - copies helper)
generate_tximport_script() {
	local method="$1"
	local quant_dir="$2"
	local output_script="$3"
	local metadata_file="$4"
	local helper_script="$HELPERS_DIR/tximport_helper.R"
	
	if [[ -f "$helper_script" ]]; then
		cp "$helper_script" "$output_script"
		chmod +x "$output_script"
		log_info "[TXIMPORT] Copied helper to: $output_script"
	else
		log_error "tximport_helper.R not found: $helper_script"
		return 1
	fi
}

# ==============================================================================
# GENE-TRANSCRIPT MAPPING
# ==============================================================================

# Create gene-transcript mapping from FASTA
create_gene_trans_map() {
	local fasta="$1"
	local output_file="${2:-${fasta}.gene_trans_map}"
	
	if [[ -f "$output_file" && "${OVERWRITE_MODE:-skip}" != "overwrite" ]]; then
		log_info "Gene-transcript map already exists: $output_file"
		return 0
	fi
	
	log_info "Creating gene-transcript mapping from FASTA..."
	
	# Detect if this is a Trinity assembly
	local is_trinity=false
	if grep -q "^>TRINITY_" "$fasta" 2>/dev/null; then
		is_trinity=true
		log_info "Detected Trinity assembly format"
	fi
	
	if [[ "$is_trinity" == "true" ]]; then
		awk '/^>/ {
			trans = $1
			gsub(/^>/, "", trans)
			gene = trans
			# POSIX-portable: sub replaces _i<digits> suffix (no gawk capture groups)
			sub(/_i[0-9]+$/, "", gene)
			print gene "\t" trans
		}' "$fasta" > "$output_file"
	else
		# Single awk pass reads FASTA directly (replaces grep|sed|awk chain — 2 fewer processes)
		# Uses POSIX-portable match()+substr() instead of gawk-only capture groups
		awk '/^>/ {
			sub(/^>/, "")
			trans=$1
			gene=""
			if (match($0, /gene=[^ ]+/)) {
				gene=substr($0, RSTART+5, RLENGTH-5)
			} else if (match($0, /gene_id[=:][^ ]+/)) {
				gene=substr($0, RSTART+8, RLENGTH-8)
			} else if (match(trans, /^[^|]+\|/)) {
				gene=substr(trans, 1, RLENGTH-1)
			} else {
				gene=trans
			}
			print gene "\t" trans
		}' "$fasta" > "$output_file"
	fi
	
	local unique_genes total_transcripts
	read -r unique_genes total_transcripts < <(awk -F'\t' '{ genes[$1]++; total++ } END { print length(genes), total }' "$output_file")
	log_info "Created gene-transcript map: $unique_genes genes, $total_transcripts transcripts"
}

# ==============================================================================
# NORMALIZATION FUNCTIONS
# ==============================================================================

normalize_expression_data() {
	local matrix_dir="$1"
	local method="$2"
	
	if [[ -f "$matrix_dir/genes.counts.matrix" ]]; then
		log_step "Applying TMM normalization for $method"
		if command -v normalize_matrix.pl >/dev/null 2>&1; then
			run_with_space_time_log normalize_matrix.pl "$matrix_dir/genes.counts.matrix" \
				--est_method "$method" \
				--out_prefix "$matrix_dir/genes.TMM" || \
				log_warn "TMM normalization failed for $method"
		else
			log_warn "normalize_matrix.pl not found. Skipping TMM normalization."
		fi
	fi
}

# ==============================================================================
# MEMORY DETECTION (shared by _samtools_sort_mem and _star_sort_ram)
# ==============================================================================

# Get available system RAM in MB. Caches result in _CACHED_AVAIL_MB to avoid
# repeated /proc/meminfo reads in tight loops.
# Usage: _get_available_ram_mb
_get_available_ram_mb() {
	if [[ -n "${_CACHED_AVAIL_MB:-}" ]]; then
		echo "$_CACHED_AVAIL_MB"
		return
	fi
	local avail_mb
	avail_mb=$(awk '/MemAvailable/ {printf "%d", $2/1024}' /proc/meminfo 2>/dev/null) \
		|| avail_mb=$(sysctl -n hw.memsize 2>/dev/null | awk '{printf "%d", $1/1048576}') \
		|| avail_mb=8192  # fallback: 8GB
	export _CACHED_AVAIL_MB="$avail_mb"
	echo "$avail_mb"
}
export -f _get_available_ram_mb

# ==============================================================================
# SAMTOOLS SORT MEMORY CALCULATION
# ==============================================================================

# Calculate safe per-thread memory for samtools sort.
# samtools sort -m VALUE is per-thread (not total), so total RAM = VALUE x (threads + 1).
# Queries available system RAM and divides by active sorting threads,
# reserving headroom for the aligner and other processes.
# Usage: _samtools_sort_mem <num_threads> [parallel_jobs]
_samtools_sort_mem() {
	local sort_threads="${1:-4}"
	local parallel_jobs="${2:-1}"
	local total_slots=$(( (sort_threads + 1) * parallel_jobs ))
	[[ $total_slots -lt 1 ]] && total_slots=1

	local avail_mb
	avail_mb=$(_get_available_ram_mb)

	# Reserve 25% for aligner, OS, and other processes
	local usable_mb=$(( avail_mb * 75 / 100 ))
	local per_thread_mb=$(( usable_mb / total_slots ))

	# Clamp between 256MB and configurable max per thread (default 4GB; override with
	# SAMTOOLS_SORT_MEM_MAX_MB for high-memory systems, e.g. 8192 for 128GB+ servers)
	local max_per_thread="${SAMTOOLS_SORT_MEM_MAX_MB:-4096}"
	[[ $per_thread_mb -lt 256 ]] && per_thread_mb=256
	[[ $per_thread_mb -gt $max_per_thread ]] && per_thread_mb=$max_per_thread
	# Safety: re-check total allocation doesn't exceed usable memory (handles high thread counts)
	local total_alloc=$(( per_thread_mb * total_slots ))
	if [[ $total_alloc -gt $usable_mb ]]; then
		per_thread_mb=$(( usable_mb / total_slots ))
		[[ $per_thread_mb -lt 256 ]] && per_thread_mb=256
	fi

	echo "${per_thread_mb}M"
}
export -f _samtools_sort_mem

# Check if samtools supports --write-index (requires samtools >= 1.10)
# Caches result in _SAMTOOLS_HAS_WRITE_INDEX for repeated calls.
_samtools_has_write_index() {
	if [[ -z "${_SAMTOOLS_HAS_WRITE_INDEX:-}" ]]; then
		# Check samtools version >= 1.10 (when --write-index was added)
		local _st_ver
		_st_ver=$(samtools --version 2>/dev/null | head -1 | grep -oP '[0-9]+\.[0-9]+' | head -1) || _st_ver="0.0"
		local _st_major=${_st_ver%%.*} _st_minor=${_st_ver#*.}
		_st_minor=${_st_minor%%.*}
		if [[ "${_st_major:-0}" -gt 1 ]] || [[ "${_st_major:-0}" -eq 1 && "${_st_minor:-0}" -ge 10 ]]; then
			export _SAMTOOLS_HAS_WRITE_INDEX="yes"
		else
			export _SAMTOOLS_HAS_WRITE_INDEX="no"
		fi
	fi
	[[ "$_SAMTOOLS_HAS_WRITE_INDEX" == "yes" ]]
}
export -f _samtools_has_write_index

# ==============================================================================
# GNU PARALLEL HELPER FUNCTIONS
# ==============================================================================

# Initialize parallel worker environment (call at start of any parallel worker)
# Sets: trimmed1, trimmed2
# Requires exported: abs_trim_dir_root, CONDA_PREFIX, CONDA_EXE, CONDA_DEFAULT_ENV
_init_parallel_worker() {
	local SRR="$1"

	# Reactivate conda in subshell if needed
	if [[ -n "${CONDA_PREFIX:-}" && -n "${CONDA_EXE:-}" ]]; then
		source "$(dirname "$CONDA_EXE")/../etc/profile.d/conda.sh" 2>/dev/null || true
		conda activate "${CONDA_DEFAULT_ENV:-base}" 2>/dev/null || true
	fi

	# Inline find_trimmed_fastq (avoids function export issues)
	local trim_dir="${abs_trim_dir_root}/$SRR"
	trimmed1="" trimmed2=""
	if [[ -f "$trim_dir/${SRR}_1_val_1.fq.gz" && -f "$trim_dir/${SRR}_2_val_2.fq.gz" ]]; then
		trimmed1="$trim_dir/${SRR}_1_val_1.fq.gz"
		trimmed2="$trim_dir/${SRR}_2_val_2.fq.gz"
	elif [[ -f "$trim_dir/${SRR}_1_val_1.fq" && -f "$trim_dir/${SRR}_2_val_2.fq" ]]; then
		trimmed1="$trim_dir/${SRR}_1_val_1.fq"
		trimmed2="$trim_dir/${SRR}_2_val_2.fq"
	elif [[ -f "$trim_dir/${SRR}_trimmed.fq.gz" ]]; then
		trimmed1="$trim_dir/${SRR}_trimmed.fq.gz"
	elif [[ -f "$trim_dir/${SRR}_trimmed.fq" ]]; then
		trimmed1="$trim_dir/${SRR}_trimmed.fq"
	else
		for f in "$trim_dir"/${SRR}*val_1*.fq* "$trim_dir"/${SRR}*val_1*.gz; do
			[[ -f "$f" ]] && { trimmed1="$f"; break; }
		done
		for f in "$trim_dir"/${SRR}*val_2*.fq* "$trim_dir"/${SRR}*val_2*.gz; do
			[[ -f "$f" ]] && { trimmed2="$f"; break; }
		done
	fi
}
export -f _init_parallel_worker

# Log helper for parallel workers
# Usage: _parallel_log METHOD SRR LEVEL message
# Uses bash built-in printf %(%T)T to avoid forking a date subshell on every call.
_parallel_log() {
	local method="$1" SRR="$2" level="$3"; shift 3
	local ts
	printf -v ts '%(%Y-%m-%d %H:%M:%S)T' -1
	printf '[%s] [%s] [%s-%s] %s\n' "$ts" "$level" "$method" "$SRR" "$*"
	[[ "$level" != "INFO" && -n "${abs_error_warn_file:-}" ]] && \
		printf '[%s] [%s] [%s-%s] %s\n' "$ts" "$level" "$method" "$SRR" "$*" >> "$abs_error_warn_file"
}
export -f _parallel_log

# Export common environment variables for parallel workers
# Call before dispatching parallel jobs
_prepare_parallel_env() {
	export PATH CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_EXE
	export keep_bam_global

	abs_trim_dir_root="$TRIM_DIR_ROOT"
	[[ "$abs_trim_dir_root" != /* ]] && abs_trim_dir_root="$(pwd)/$abs_trim_dir_root"
	export abs_trim_dir_root

	abs_error_warn_file="${ERROR_WARN_FILE:-}"
	[[ -n "$abs_error_warn_file" && "$abs_error_warn_file" != /* ]] && abs_error_warn_file="$(pwd)/$abs_error_warn_file"
	export abs_error_warn_file
}

