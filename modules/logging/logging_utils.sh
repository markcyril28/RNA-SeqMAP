#!/bin/bash
# ==============================================================================
# LOGGING UTILITIES
# ==============================================================================
# Six-component logging system for comprehensive pipeline tracking:
# 1. Full logs (logs/log_files/*.log) - Complete execution output
# 2. Time logs (logs/time_logs/*.csv) - Time/CPU/memory metrics only
# 3. Space logs (logs/space_logs/*.csv) - File/directory size metrics only
# 4. Combined logs (logs/space_time_logs/*.csv) - Time + space metrics together
# 5. Error/Warning logs (logs/error_warn_logs/*.log) - Errors and warnings only
# 6. Software catalog (logs/software_catalogs/*.csv) - Software versions used
# ==============================================================================
# ERROR CAPTURE:
# - Monitors for: error, exception, fatal, failed, command not found,
#   no such file, cannot find, not installed, permission denied, traceback,
#   segmentation fault, killed, out of memory, no space left, broken pipe
# - Captures from: stderr, stdout, time output, and exit codes
# - Use run_with_error_capture() for simple commands
# - Use run_with_space_time_log() for resource-intensive commands
# ==============================================================================

# Guard against double-sourcing
[[ "${LOGGING_UTILS_SOURCED:-}" == "true" ]] && return 0
export LOGGING_UTILS_SOURCED="true"

# ==============================================================================
# LOGGING CONFIGURATION - IMPORTANT PARAMETERS AT TOP
# ==============================================================================

# Run ID for unique log file naming
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"

# Log directory structure
LOG_DIR="${LOG_DIR:-logs/log_files}"
TIME_DIR="${TIME_DIR:-logs/time_logs}"
SPACE_DIR="${SPACE_DIR:-logs/space_logs}"
SPACE_TIME_DIR="${SPACE_TIME_DIR:-logs/space_time_logs}"
ERROR_WARN_DIR="${ERROR_WARN_DIR:-logs/error_warn_logs}"
SOFTWARE_CATALOG_DIR="${SOFTWARE_CATALOG_DIR:-logs/software_catalogs}"
GPU_LOG_DIR="${GPU_LOG_DIR:-logs/gpu_log}"

# Log file paths (derived from directories and RUN_ID)
LOG_FILE="${LOG_FILE:-$LOG_DIR/pipeline_${RUN_ID}_full_log.log}"
TIME_FILE="${TIME_FILE:-$TIME_DIR/pipeline_${RUN_ID}_time_metrics.csv}"
TIME_TEMP="${TIME_TEMP:-$TIME_DIR/.time_temp_${RUN_ID}.txt}"
SPACE_FILE="${SPACE_FILE:-$SPACE_DIR/pipeline_${RUN_ID}_space_metrics.csv}"
SPACE_TIME_FILE="${SPACE_TIME_FILE:-$SPACE_TIME_DIR/pipeline_${RUN_ID}_combined_metrics.csv}"
ERROR_WARN_FILE="${ERROR_WARN_FILE:-$ERROR_WARN_DIR/pipeline_${RUN_ID}_errors_warnings.log}"
SOFTWARE_FILE="${SOFTWARE_FILE:-$SOFTWARE_CATALOG_DIR/software_catalog_${RUN_ID}.csv}"
GPU_LOG_FILE="${GPU_LOG_FILE:-$GPU_LOG_DIR/gpu_${RUN_ID}.log}"

# Logging behavior
log_choice="${log_choice:-1}"  # 1 = tee to console, 2 = file only

# ==============================================================================
# CORE LOGGING FUNCTIONS
# ==============================================================================

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
log() { local level="$1"; shift; printf '[%s] [%s] %s\n' "$(timestamp)" "$level" "$*"; }
log_info() { log INFO "$@"; }
log_warn() { log WARN "$@"; [[ -n "${ERROR_WARN_FILE:-}" ]] && printf '[%s] [WARN] %s\n' "$(timestamp)" "$*" >> "$ERROR_WARN_FILE"; }
log_error() { log ERROR "$@"; [[ -n "${ERROR_WARN_FILE:-}" ]] && printf '[%s] [ERROR] %s\n' "$(timestamp)" "$*" >> "$ERROR_WARN_FILE"; }
log_step() { log INFO "=============== $* ==============="; }

strip_ansi_stream() {
	# Strip ANSI escape codes (colors, cursor moves, erase sequences) and
	# carriage returns from a stream so log files remain human-readable.
	# CR (\r) is converted to newline so progress-bar overwrites become
	# separate lines instead of one giant unreadable blob.
	tr '\r' '\n' | sed -u 's/\x1B\[[0-9;?]*[a-zA-Z]//g; s/\x1B[()][A-Z0-9]//g'
}

# ==============================================================================
# LOGGING SETUP
# ==============================================================================

setup_logging() {
	# Set up logging and output redirection with dual-format support
	# Usage: setup_logging [clear_logs_flag]
	# clear_logs_flag: "true" to clear existing logs, anything else to keep them
	local clear_logs="${1:-false}"
	
	# Skip if already initialized
	if [[ "${LOGGING_INITIALIZED:-}" == "true" ]]; then
		log_info "Logging already initialized, skipping setup"
		return 0
	fi

	# Re-derive file paths from directories (in case directories were set after sourcing)
	# This ensures absolute paths are used when directories are set with absolute paths
	LOG_FILE="${LOG_DIR}/pipeline_${RUN_ID}_full_log.log"
	TIME_FILE="${TIME_DIR}/pipeline_${RUN_ID}_time_metrics.csv"
	TIME_TEMP="${TIME_DIR}/.time_temp_${RUN_ID}.txt"
	SPACE_FILE="${SPACE_DIR}/pipeline_${RUN_ID}_space_metrics.csv"
	SPACE_TIME_FILE="${SPACE_TIME_DIR}/pipeline_${RUN_ID}_combined_metrics.csv"
	ERROR_WARN_FILE="${ERROR_WARN_DIR}/pipeline_${RUN_ID}_errors_warnings.log"
	SOFTWARE_FILE="${SOFTWARE_CATALOG_DIR}/software_catalog_${RUN_ID}.csv"
	GPU_LOG_FILE="${GPU_LOG_DIR}/gpu_${RUN_ID}.log"

	# Create all log directories
	mkdir -p "$LOG_DIR" "$TIME_DIR" "$SPACE_DIR" "$SPACE_TIME_DIR" "$ERROR_WARN_DIR" "$SOFTWARE_CATALOG_DIR" "$GPU_LOG_DIR" || {
		echo "ERROR: Failed to create logging directories" >&2
		return 1
	}
	
	# Clear previous logs if requested (case-insensitive check)
	if [[ "${clear_logs^^}" == "TRUE" ]]; then
		rm -f "$LOG_DIR"/*.log 2>/dev/null || true
		rm -f "$TIME_DIR"/*.csv 2>/dev/null || true
		rm -f "$SPACE_DIR"/*.csv 2>/dev/null || true
		rm -f "$SPACE_TIME_DIR"/*.csv 2>/dev/null || true
		rm -f "$ERROR_WARN_DIR"/*.log 2>/dev/null || true
		rm -f "$SOFTWARE_CATALOG_DIR"/*.csv 2>/dev/null || true
		rm -f "$GPU_LOG_DIR"/*.log 2>/dev/null || true
		echo "Previous logs cleared"
	fi
	
	# Initialize CSV headers
	[[ ! -f "$TIME_FILE" ]] && echo "Timestamp,Command,Elapsed_Time_sec,CPU_Percent,Max_RSS_KB,User_Time_sec,System_Time_sec,Exit_Status" > "$TIME_FILE"
	[[ ! -f "$SPACE_FILE" ]] && echo "Timestamp,Type,Path,Size_KB,Size_MB,Size_GB,File_Count,Description" > "$SPACE_FILE"
	[[ ! -f "$SPACE_TIME_FILE" ]] && echo "Timestamp,Command,Elapsed_Time_sec,CPU_Percent,Max_RSS_KB,User_Time_sec,System_Time_sec,Input_Size_MB,Output_Size_MB,Exit_Status" > "$SPACE_TIME_FILE"
	[[ ! -f "$ERROR_WARN_FILE" ]] && touch "$ERROR_WARN_FILE"
	[[ ! -f "$SOFTWARE_FILE" ]] && echo "Software/Tool,Version" > "$SOFTWARE_FILE"
	[[ ! -f "$GPU_LOG_FILE" ]] && echo "=== GPU Log Started: $(timestamp) ===" > "$GPU_LOG_FILE"
	
	# Rotate old logs to prevent unbounded growth
	rotate_old_logs "$(dirname "$LOG_DIR")"

	# Set up output redirection (strip ANSI escape codes from log files)
	if [[ "$log_choice" == "2" ]]; then
		exec > >(strip_ansi_stream >> "$LOG_FILE") 2>&1
	else
		exec > >(tee >(strip_ansi_stream >> "$LOG_FILE")) 2>&1
	fi

	export LOGGING_INITIALIZED="true"
	log_info "Logging to: $LOG_FILE"
	log_info "Time metrics to: $TIME_FILE"
	log_info "Space metrics to: $SPACE_FILE"
	log_info "Combined metrics to: $SPACE_TIME_FILE"
	log_info "Errors & Warnings to: $ERROR_WARN_FILE"
	log_info "Software catalog to: $SOFTWARE_FILE"
	log_info "GPU logs to: $GPU_LOG_FILE"
}

# ==============================================================================
# STAGE-BASED LOG ROUTING
# ==============================================================================

switch_log_stage() {
	# Switch all log output to a stage-specific directory.
	# Usage: switch_log_stage <base_dir>
	# Example: switch_log_stage "1_SRRs"
	#          switch_log_stage "2_ALIGNMENT_RESULTs"
	#          switch_log_stage "3_POST_PROC"
	local stage_base="$1"

	# Convert to absolute path if relative
	if [[ "$stage_base" != /* ]]; then
		stage_base="${PROJECT_ROOT:-$(pwd)}/$stage_base"
	fi

	# Update directory paths
	LOG_DIR="${stage_base}/logs/log_files"
	TIME_DIR="${stage_base}/logs/time_logs"
	SPACE_DIR="${stage_base}/logs/space_logs"
	SPACE_TIME_DIR="${stage_base}/logs/space_time_logs"
	ERROR_WARN_DIR="${stage_base}/logs/error_warn_logs"
	SOFTWARE_CATALOG_DIR="${stage_base}/logs/software_catalogs"
	GPU_LOG_DIR="${stage_base}/logs/gpu_log"

	# Update file paths
	LOG_FILE="${LOG_DIR}/pipeline_${RUN_ID}_full_log.log"
	TIME_FILE="${TIME_DIR}/pipeline_${RUN_ID}_time_metrics.csv"
	TIME_TEMP="${TIME_DIR}/.time_temp_${RUN_ID}.txt"
	SPACE_FILE="${SPACE_DIR}/pipeline_${RUN_ID}_space_metrics.csv"
	SPACE_TIME_FILE="${SPACE_TIME_DIR}/pipeline_${RUN_ID}_combined_metrics.csv"
	ERROR_WARN_FILE="${ERROR_WARN_DIR}/pipeline_${RUN_ID}_errors_warnings.log"
	SOFTWARE_FILE="${SOFTWARE_CATALOG_DIR}/software_catalog_${RUN_ID}.csv"
	GPU_LOG_FILE="${GPU_LOG_DIR}/gpu_${RUN_ID}.log"

	# Create directories
	mkdir -p "$LOG_DIR" "$TIME_DIR" "$SPACE_DIR" "$SPACE_TIME_DIR" \
		"$ERROR_WARN_DIR" "$SOFTWARE_CATALOG_DIR" "$GPU_LOG_DIR" || {
		echo "ERROR: Failed to create log directories for stage: $stage_base" >&2
		return 1
	}

	# Initialize CSV headers if files don't exist
	[[ ! -f "$TIME_FILE" ]] && echo "Timestamp,Command,Elapsed_Time_sec,CPU_Percent,Max_RSS_KB,User_Time_sec,System_Time_sec,Exit_Status" > "$TIME_FILE"
	[[ ! -f "$SPACE_FILE" ]] && echo "Timestamp,Type,Path,Size_KB,Size_MB,Size_GB,File_Count,Description" > "$SPACE_FILE"
	[[ ! -f "$SPACE_TIME_FILE" ]] && echo "Timestamp,Command,Elapsed_Time_sec,CPU_Percent,Max_RSS_KB,User_Time_sec,System_Time_sec,Input_Size_MB,Output_Size_MB,Exit_Status" > "$SPACE_TIME_FILE"
	[[ ! -f "$ERROR_WARN_FILE" ]] && touch "$ERROR_WARN_FILE"
	[[ ! -f "$SOFTWARE_FILE" ]] && echo "Software/Tool,Version" > "$SOFTWARE_FILE"
	[[ ! -f "$GPU_LOG_FILE" ]] && echo "=== GPU Log Started: $(timestamp) ===" > "$GPU_LOG_FILE"

	# Re-setup output redirection to the new log file (strip ANSI codes)
	if [[ "$log_choice" == "2" ]]; then
		exec > >(strip_ansi_stream >> "$LOG_FILE") 2>&1
	else
		exec > >(tee >(strip_ansi_stream >> "$LOG_FILE")) 2>&1
	fi

	log_info "Switched logging to stage: $stage_base"
}

# ==============================================================================
# ERROR HANDLING
# ==============================================================================

# Error handling trap (can be enabled/disabled by caller)
enable_error_trap() {
	trap 'log_error "Command failed (rc=$?) at line $LINENO: ${BASH_COMMAND:-unknown}"; exit 1' ERR
}

# Cleanup trap
enable_exit_trap() {
	trap 'log_info "Script finished. See log: $LOG_FILE"; log_info "Time metrics: $TIME_FILE"; log_info "Errors & Warnings: $ERROR_WARN_FILE"' EXIT
}

# Log pipeline configuration settings
log_configuration() {
	log_step "PIPELINE CONFIGURATION"
	log_info "Run ID: ${RUN_ID:-N/A}"
	log_info "Threads: ${THREADS:-N/A}"
	log_info "Jobs: ${JOBS:-N/A}"
	log_info "Threads per job: ${THREADS_PER_JOB:-N/A}"
	log_info "GNU Parallel: ${USE_GNU_PARALLEL:-FALSE}"
	log_info "Keep BAM files: ${keep_bam_global:-n}"
	log_info "Project root: ${PROJECT_ROOT:-N/A}"
	
	# Log active pipeline stages if defined
	if [[ -n "${PIPELINE_STAGES[*]:-}" ]]; then
		log_info "Active pipeline stages:"
		for stage in "${PIPELINE_STAGES[@]}"; do
			log_info "  - $stage"
		done
	fi
	
	log_step "END CONFIGURATION"
}

# Unified error/warning regex patterns (single source of truth, exported for GNU Parallel)
_ERROR_PATTERN='error|exception|fatal|failed|command not found|no such file|cannot find|not installed|permission denied|traceback|segmentation fault|segfault|killed|out of memory|cannot allocate memory|no space left on device|disk full|broken pipe|filenotfound|access denied'
_WARN_PATTERN='warning|warn'
export _ERROR_PATTERN _WARN_PATTERN

capture_stderr_errors() {
	# Monitor stderr/stdout stream and capture errors to error log
	# Usage: command 2>&1 | capture_stderr_errors
	while IFS= read -r line; do
		echo "$line"
		if echo "$line" | grep -qiE "$_ERROR_PATTERN"; then
			printf '[%s] [ERROR] %s\n' "$(timestamp)" "$line" >> "$ERROR_WARN_FILE"
		fi
		if echo "$line" | grep -qiE "$_WARN_PATTERN"; then
			printf '[%s] [WARN] %s\n' "$(timestamp)" "$line" >> "$ERROR_WARN_FILE"
		fi
	done
}

run_with_error_capture() {
	# Simple wrapper to run commands with error capture (without time logging)
	# Usage: run_with_error_capture COMMAND...
	local cmd_string="$*"
	local exit_code=0

	"$@" 2>&1 | capture_stderr_errors
	exit_code=${PIPESTATUS[0]}

	if [[ $exit_code -ne 0 ]]; then
		log_error "Command failed (exit=$exit_code): $cmd_string"
	fi

	return $exit_code
}

# ==============================================================================
# TIME AND RESOURCE LOGGING
# ==============================================================================

run_with_space_time_log() {
	# Run a command and log resource usage (tracks time and memory)
	# Usage: run_with_space_time_log [--input PATH] [--output PATH] COMMAND...
	
	local input_path=""
	local output_path=""
	
	# Parse optional space tracking arguments
	while [[ $# -gt 0 ]]; do
		case "$1" in
			--input) input_path="$2"; shift 2 ;;
			--output) output_path="$2"; shift 2 ;;
			*) break ;;
		esac
	done
	
	local cmd_string="$*"
	local start_ts="$(timestamp)"
	
	# Measure input size before running command
	local input_size_mb="0"
	if [[ -n "$input_path" && -e "$input_path" ]]; then
		local input_kb=$(du -sk "$input_path" 2>/dev/null | awk '{print $1}')
		input_size_mb=$(awk "BEGIN{printf \"%.2f\", $input_kb / 1024}")
	fi
	
	mkdir -p "$TIME_DIR" || { log_error "Failed to create TIME_DIR: $TIME_DIR"; return 1; }
	
	local exit_code=0

	# Log abbreviated command before running (full command saved in CSV)
	local cmd_abbrev="${cmd_string:0:120}"
	[[ ${#cmd_string} -gt 120 ]] && cmd_abbrev="${cmd_abbrev}..."
	log_info "[CMD] $cmd_abbrev"

	# Write begin marker directly to log file (preserves ordering with tool stdout)
	printf '[%s] [INFO] --- BEGIN TOOL OUTPUT: %s ---\n' "$(timestamp)" "${1##*/}" >> "$LOG_FILE"
	# Strip ANSI escape codes and carriage returns before writing to log (e.g. Salmon progress bars)
	/usr/bin/time -v "$@" 2>"$TIME_TEMP" | strip_ansi_stream >> "$LOG_FILE"
	exit_code=${PIPESTATUS[0]}
	printf '[%s] [INFO] --- END TOOL OUTPUT: %s (exit=%d) ---\n' "$(timestamp)" "${1##*/}" "$exit_code" >> "$LOG_FILE"

	# Log key resource metrics as a single summary line (replaces 22-line verbose dump)
	local elapsed_raw cpu_raw rss_raw
	elapsed_raw=$(grep "Elapsed (wall clock)" "$TIME_TEMP" 2>/dev/null | awk '{print $NF}')
	cpu_raw=$(grep "Percent of CPU" "$TIME_TEMP" 2>/dev/null | awk '{print $NF}')
	rss_raw=$(grep "Maximum resident set size" "$TIME_TEMP" 2>/dev/null | awk '{print $NF}')
	log_info "[RESOURCES] Elapsed: ${elapsed_raw:-N/A} | CPU: ${cpu_raw:-N/A} | MaxRSS: ${rss_raw:-0} KB | Exit: $exit_code"

	# On failure: dump full time output for debugging
	if [[ $exit_code -ne 0 ]]; then
		printf '[%s] [DEBUG] --- TIME OUTPUT (failure details) ---\n' "$(timestamp)" >> "$LOG_FILE"
		cat "$TIME_TEMP" >> "$LOG_FILE" 2>&1
	fi

	# Capture errors/exceptions to error log (uses unified pattern)
	if [[ $exit_code -ne 0 ]] || grep -qiE "$_ERROR_PATTERN" "$TIME_TEMP" 2>/dev/null; then
		{
			printf '[%s] [ERROR] Command failed (exit=%d): %s\n' "$(timestamp)" "$exit_code" "$cmd_string"
			grep -iE "$_ERROR_PATTERN" "$TIME_TEMP" 2>/dev/null || true
		} >> "$ERROR_WARN_FILE"
	fi

	# Extract key metrics from time output (for CSV logging)
	local elapsed_time=$(grep "Elapsed (wall clock)" "$TIME_TEMP" | awk '{print $NF}' | awk -F: '{if (NF==3) print ($1*3600)+($2*60)+$3; else if (NF==2) print ($1*60)+$2; else print $1}')
	local cpu_percent=$(grep "Percent of CPU" "$TIME_TEMP" | awk '{print $NF}' | tr -d '%')
	local max_rss=$(grep "Maximum resident set size" "$TIME_TEMP" | awk '{print $NF}')
	local user_time=$(grep "User time" "$TIME_TEMP" | awk '{print $NF}')
	local system_time=$(grep "System time" "$TIME_TEMP" | awk '{print $NF}')
	
	# Measure output size after running command
	local output_size_mb="0"
	if [[ -n "$output_path" && -e "$output_path" ]]; then
		local output_kb=$(du -sk "$output_path" 2>/dev/null | awk '{print $1}')
		output_size_mb=$(awk "BEGIN{printf \"%.2f\", $output_kb / 1024}")
	fi
	
	# Append to CSV files (escape internal double quotes for valid CSV)
	local csv_cmd="${cmd_string//\"/\"\"}"
	echo "${start_ts},\"${csv_cmd}\",${elapsed_time:-0},${cpu_percent:-0},${max_rss:-0},${user_time:-0},${system_time:-0},${exit_code}" >> "$TIME_FILE"
	echo "${start_ts},\"${csv_cmd}\",${elapsed_time:-0},${cpu_percent:-0},${max_rss:-0},${user_time:-0},${system_time:-0},${input_size_mb},${output_size_mb},${exit_code}" >> "$SPACE_TIME_FILE"
	
	rm -f "$TIME_TEMP"
	return $exit_code
}

# ==============================================================================
# SPACE LOGGING FUNCTIONS
# ==============================================================================

log_file_size() {
	# Log size of a single file or directory
	local file_path="$1"
	local description="${2:-}"
	local type="FILE"
	
	[[ ! -e "$file_path" ]] && { log_warn "Path does not exist: $file_path"; return 1; }
	
	[[ -d "$file_path" ]] && type="DIR"
	
	local size_kb=$(du -sk "$file_path" 2>/dev/null | awk '{print $1}')
	local size_mb=$(awk "BEGIN{printf \"%.2f\", $size_kb / 1024}")
	local size_gb=$(awk "BEGIN{printf \"%.2f\", $size_kb / 1048576}")
	
	local file_count="-"
	[[ -d "$file_path" ]] && file_count=$(find "$file_path" -type f 2>/dev/null | wc -l)
	
	local ts="$(timestamp)"
	echo "${ts},${type},\"${file_path}\",${size_kb},${size_mb},${size_gb},${file_count},\"${description}\"" >> "$SPACE_FILE"
	log_info "Space logged: $file_path = ${size_mb}MB"
}

log_input_output_size() {
	# Log sizes of input and output files/directories
	local input_path="$1"
	local output_path="$2"
	local step_description="${3:-}"
	
	[[ -e "$input_path" ]] && log_file_size "$input_path" "Input: $step_description"
	[[ -e "$output_path" ]] && log_file_size "$output_path" "Output: $step_description"
}

# ==============================================================================
# SOFTWARE CATALOG FUNCTIONS
# ==============================================================================

log_software_version() {
	# Log software version to catalog
	local software="$1"
	local version="$2"
	
	echo "${software},${version}" >> "$SOFTWARE_FILE"
	log_info "Recorded software: $software v$version"
}

catalog_all_software() {
	# Catalog versions of all bioinformatics tools, R packages, and provenance metadata
	log_step "Cataloging software versions"

	# Record pipeline git commit SHA for provenance
	if command -v git >/dev/null 2>&1; then
		local git_sha
		git_sha=$(git -C "$(dirname "${BASH_SOURCE[0]}")/../.." rev-parse --short HEAD 2>/dev/null || echo "not_a_git_repo")
		local git_branch
		git_branch=$(git -C "$(dirname "${BASH_SOURCE[0]}")/../.." rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
		local git_dirty=""
		if ! git -C "$(dirname "${BASH_SOURCE[0]}")/../.." diff --quiet HEAD 2>/dev/null; then
			git_dirty="-dirty"
		fi
		echo "pipeline_git_commit,${git_sha}${git_dirty}" >> "$SOFTWARE_FILE"
		echo "pipeline_git_branch,${git_branch}" >> "$SOFTWARE_FILE"
		log_info "Pipeline git: ${git_branch}@${git_sha}${git_dirty}"
	fi

	# Record conda environment name
	echo "conda_env,${CONDA_DEFAULT_ENV:-unknown}" >> "$SOFTWARE_FILE"

	local tools=(
		"hisat2:hisat2 --version"
		"stringtie:stringtie --version"
		"samtools:samtools --version"
		"star:STAR --version"
		"salmon:salmon --version"
		"rsem:rsem-calculate-expression --version"
		"bowtie2:bowtie2 --version"
		"trim_galore:trim_galore --version"
		"trimmomatic:trimmomatic -version"
		"fastqc:fastqc --version"
		"multiqc:multiqc --version"
		"gffread:gffread --version"
		"cutadapt:cutadapt --version"
		"sra-tools:prefetch --version"
		"infer_experiment.py:infer_experiment.py --version"
		"prepDE.py:prepDE.py --version"
		"python:python3 --version"
		"parallel:parallel --version"
		"R:R --version"
	)

	for tool_cmd in "${tools[@]}"; do
		local tool="${tool_cmd%%:*}"
		local cmd="${tool_cmd#*:}"

		if command -v "${cmd%% *}" >/dev/null 2>&1; then
			local version
			version=$(eval "$cmd" 2>&1 | head -n1 | awk '{print $NF}' || echo "unknown")
			log_software_version "$tool" "$version"
		else
			echo "${tool},not_installed" >> "$SOFTWARE_FILE"
			log_info "Software not found: $tool"
		fi
	done

	# Catalog key R/Bioconductor packages used by analysis modules
	if command -v Rscript >/dev/null 2>&1; then
		log_info "Cataloging R package versions..."
		local r_pkgs=(DESeq2 tximport tximeta WGCNA clusterProfiler ComplexHeatmap
			ballgown AnnotationDbi enrichplot DOSE fgsea
			pheatmap ggplot2 corrplot dendextend gridExtra scales)
		for pkg in "${r_pkgs[@]}"; do
			local ver
			ver=$(Rscript -e "tryCatch(cat(as.character(packageVersion('$pkg'))), error=function(e) cat('not_installed'))" 2>/dev/null || echo "unknown")
			echo "R/${pkg},${ver}" >> "$SOFTWARE_FILE"
		done
		log_info "R package versions cataloged"

		# Save full R sessionInfo for complete reproducibility record
		local session_info_file
		session_info_file="$(dirname "$SOFTWARE_FILE")/R_sessionInfo_${RUN_ID}.txt"
		Rscript -e "writeLines(capture.output(sessionInfo()), '$session_info_file')" 2>/dev/null \
			&& log_info "R sessionInfo saved to: $session_info_file" \
			|| log_warn "Failed to capture R sessionInfo"
	else
		echo "R,not_installed" >> "$SOFTWARE_FILE"
		log_warn "Rscript not found — R package catalog skipped"
	fi
}

# ==============================================================================
# LOG ROTATION
# ==============================================================================

rotate_old_logs() {
	# Remove logs older than MAX_LOG_AGE_DAYS (default 30) to prevent unbounded growth.
	# Usage: rotate_old_logs [base_log_dir]
	# Called automatically by setup_logging; can also be called manually.
	local base_dir="${1:-$(dirname "$LOG_DIR")}"
	local max_age="${MAX_LOG_AGE_DAYS:-30}"

	[[ ! -d "$base_dir" ]] && return 0

	local count
	count=$(find "$base_dir" -type f \( -name '*.log' -o -name '*.csv' \) -mtime +"$max_age" 2>/dev/null | wc -l)
	if [[ "$count" -gt 0 ]]; then
		find "$base_dir" -type f \( -name '*.log' -o -name '*.csv' \) -mtime +"$max_age" -delete 2>/dev/null || true
		log_info "Log rotation: removed $count files older than ${max_age} days from $base_dir"
	fi
}

# ==============================================================================
# GPU LOGGING FUNCTIONS
# ==============================================================================

log_gpu() {
	# Log GPU-related message to GPU log file
	# Usage: log_gpu "message"
	local message="$*"
	printf '[%s] %s\n' "$(timestamp)" "$message" >> "$GPU_LOG_FILE"
	log_info "[GPU] $message"
}

log_gpu_info() {
	# Log GPU information and status using nvidia-smi
	# Usage: log_gpu_info [description]
	local description="${1:-GPU Status}"
	
	if ! command -v nvidia-smi >/dev/null 2>&1; then
		log_gpu "nvidia-smi not found - GPU monitoring unavailable"
		return 1
	fi
	
	{
		printf '\n=== %s: %s ===\n' "$description" "$(timestamp)"
		nvidia-smi
		printf '\n'
	} >> "$GPU_LOG_FILE"
	
	log_info "GPU info logged: $description"
}

log_gpu_memory() {
	# Log GPU memory usage
	# Usage: log_gpu_memory [description]
	local description="${1:-GPU Memory}"
	
	if ! command -v nvidia-smi >/dev/null 2>&1; then
		log_gpu "nvidia-smi not found - GPU monitoring unavailable"
		return 1
	fi
	
	local gpu_mem=$(nvidia-smi --query-gpu=memory.used,memory.total,memory.free --format=csv,noheader,nounits 2>/dev/null)
	
	{
		printf '[%s] %s: %s\n' "$(timestamp)" "$description" "$gpu_mem"
	} >> "$GPU_LOG_FILE"
	
	log_info "GPU memory logged: $gpu_mem"
}

log_gpu_utilization() {
	# Log GPU utilization percentage
	# Usage: log_gpu_utilization [description]
	local description="${1:-GPU Utilization}"
	
	if ! command -v nvidia-smi >/dev/null 2>&1; then
		log_gpu "nvidia-smi not found - GPU monitoring unavailable"
		return 1
	fi
	
	local gpu_util=$(nvidia-smi --query-gpu=utilization.gpu,utilization.memory,temperature.gpu --format=csv,noheader 2>/dev/null)
	
	{
		printf '[%s] %s: %s\n' "$(timestamp)" "$description" "$gpu_util"
	} >> "$GPU_LOG_FILE"
	
	log_info "GPU utilization logged: $gpu_util"
}

run_with_gpu_log() {
	# Run a command with continuous GPU VRAM monitoring (captures peak usage)
	# Usage: run_with_gpu_log COMMAND...
	local cmd_string="$*"
	local gpu_monitor_pid=""
	local peak_file=""

	log_gpu "Starting command: $cmd_string"
	log_gpu_memory "Before: $cmd_string"

	# Start background GPU monitor (polls every 2s, writes peak to temp file)
	if command -v nvidia-smi >/dev/null 2>&1; then
		peak_file=$(mktemp "${GPU_LOG_DIR}/.gpu_peak_XXXXXX")
		echo "0" > "$peak_file"
		(
			trap 'exit 0' TERM
			local peak_used=0
			while true; do
				local used
				used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | head -1)
				if [[ -n "$used" ]]; then
					printf '[%s] [GPU-MONITOR] VRAM used: %s MB\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$used" >> "$GPU_LOG_FILE"
					if (( used > peak_used )); then
						peak_used=$used
						echo "$peak_used" > "$peak_file"
					fi
				fi
				sleep 2
			done
		) &
		gpu_monitor_pid=$!
	fi

	local exit_code=0
	"$@" 2>&1 | tee >(strip_ansi_stream >> "$GPU_LOG_FILE")
	exit_code=${PIPESTATUS[0]}

	# Stop GPU monitor and log peak
	if [[ -n "$gpu_monitor_pid" ]]; then
		kill "$gpu_monitor_pid" 2>/dev/null; wait "$gpu_monitor_pid" 2>/dev/null || true
	fi
	if [[ -n "$peak_file" && -f "$peak_file" ]]; then
		local peak_vram
		peak_vram=$(cat "$peak_file" 2>/dev/null)
		[[ -n "$peak_vram" && "$peak_vram" -gt 0 ]] 2>/dev/null && \
			log_gpu "Peak VRAM usage: ${peak_vram} MB"
		rm -f "$peak_file"
	fi

	log_gpu_memory "After: $cmd_string"
	log_gpu "Finished command (exit=$exit_code): $cmd_string"

	return $exit_code
}
