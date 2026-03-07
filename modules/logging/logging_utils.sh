#!/bin/bash
# ==============================================================================
# LOGGING UTILITIES
# ==============================================================================
# Five-version logging system for comprehensive pipeline tracking:
# 1. Full logs (logs/log_files/*.log) - Complete execution output
# 2. Time logs (logs/time_logs/*.csv) - Time/CPU/memory metrics only
# 3. Space logs (logs/space_logs/*.csv) - File/directory size metrics only
# 4. Combined logs (logs/space_time_logs/*.csv) - Time + space metrics together
# 5. Error/Warning logs (logs/error_warn_logs/*.log) - Errors and warnings only
# 6. Software catalog (logs/software_catalogs/*.csv) - Software versions used
# ==============================================================================
# ERROR CAPTURE:
# - Monitors for: error, exception, fatal, failed, command not found, 
#   no such file, cannot find, not installed, permission denied, traceback
# - Captures from: stderr, stdout, time output, and exit codes
# - Use run_with_error_capture() for simple commands
# - Use run_with_space_time_log() for resource-intensive commands
# ==============================================================================

#set -euo pipefail

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
	
	# Set up output redirection
	if [[ "$log_choice" == "2" ]]; then
		exec >"$LOG_FILE" 2>&1
	else
		exec > >(tee -a "$LOG_FILE") 2>&1
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

capture_stderr_errors() {
	# Monitor stderr/stdout stream and capture errors to error log
	# Usage: command 2>&1 | capture_stderr_errors
	while IFS= read -r line; do
		echo "$line"
		if echo "$line" | grep -qiE 'error|exception|fatal|failed|command not found|no such file|cannot find|not installed|permission denied|traceback'; then
			printf '[%s] [ERROR] %s\n' "$(timestamp)" "$line" >> "$ERROR_WARN_FILE"
		fi
		if echo "$line" | grep -qiE 'warning|warn'; then
			printf '[%s] [WARN] %s\n' "$(timestamp)" "$line" >> "$ERROR_WARN_FILE"
		fi
	done
}

run_with_error_capture() {
	# Simple wrapper to run commands with error capture (without time logging)
	# Usage: run_with_error_capture COMMAND...
	local cmd_string="$*"
	local exit_code=0
	
	"$@" 2>&1 | capture_stderr_errors || exit_code=$?
	
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
		input_size_mb=$(echo "scale=2; $input_kb / 1024" | bc)
	fi
	
	mkdir -p "$TIME_DIR" || { log_error "Failed to create TIME_DIR: $TIME_DIR"; return 1; }
	
	local exit_code=0
	/usr/bin/time -v "$@" >> "$LOG_FILE" 2> "$TIME_TEMP" || exit_code=$?
	
	cat "$TIME_TEMP" >> "$LOG_FILE" 2>&1
	
	# Capture errors/exceptions to error log
	if [[ $exit_code -ne 0 ]] || grep -qiE 'exception|error|fatal|failed|traceback|command not found|no such file|cannot find|not installed' "$TIME_TEMP" 2>/dev/null; then
		{
			printf '[%s] [ERROR] Command failed (exit=%d): %s\n' "$(timestamp)" "$exit_code" "$cmd_string"
			grep -iE 'exception|error|fatal|failed|traceback|filenotfound|no such file|command not found|cannot find|not installed|permission denied|access denied' "$TIME_TEMP" 2>/dev/null || true
		} >> "$ERROR_WARN_FILE"
	fi
	
	# Extract key metrics from time output
	local elapsed_time=$(grep "Elapsed (wall clock)" "$TIME_TEMP" | awk '{print $NF}' | awk -F: '{if (NF==3) print ($1*3600)+($2*60)+$3; else if (NF==2) print ($1*60)+$2; else print $1}')
	local cpu_percent=$(grep "Percent of CPU" "$TIME_TEMP" | awk '{print $NF}' | tr -d '%')
	local max_rss=$(grep "Maximum resident set size" "$TIME_TEMP" | awk '{print $NF}')
	local user_time=$(grep "User time" "$TIME_TEMP" | awk '{print $NF}')
	local system_time=$(grep "System time" "$TIME_TEMP" | awk '{print $NF}')
	
	# Measure output size after running command
	local output_size_mb="0"
	if [[ -n "$output_path" && -e "$output_path" ]]; then
		local output_kb=$(du -sk "$output_path" 2>/dev/null | awk '{print $1}')
		output_size_mb=$(echo "scale=2; $output_kb / 1024" | bc)
	fi
	
	# Append to CSV files
	echo "${start_ts},\"${cmd_string}\",${elapsed_time:-0},${cpu_percent:-0},${max_rss:-0},${user_time:-0},${system_time:-0},${exit_code}" >> "$TIME_FILE"
	echo "${start_ts},\"${cmd_string}\",${elapsed_time:-0},${cpu_percent:-0},${max_rss:-0},${user_time:-0},${system_time:-0},${input_size_mb},${output_size_mb},${exit_code}" >> "$SPACE_TIME_FILE"
	
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
	local size_mb=$(echo "scale=2; $size_kb / 1024" | bc)
	local size_gb=$(echo "scale=2; $size_kb / 1048576" | bc)
	
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
	# Catalog versions of all bioinformatics tools
	log_step "Cataloging software versions"
	
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
	)
	
	for tool_cmd in "${tools[@]}"; do
		local tool="${tool_cmd%%:*}"
		local cmd="${tool_cmd#*:}"
		
		if command -v "${cmd%% *}" >/dev/null 2>&1; then
			local version=$(eval "$cmd" 2>&1 | head -n1 || echo "unknown")
			log_software_version "$tool" "$version"
		fi
	done
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
	# Run a command and log GPU usage before and after
	# Usage: run_with_gpu_log COMMAND...
	local cmd_string="$*"
	
	log_gpu "Starting command: $cmd_string"
	log_gpu_memory "Before: $cmd_string"
	
	local exit_code=0
	"$@" 2>&1 | tee -a "$GPU_LOG_FILE" || exit_code=$?
	
	log_gpu_memory "After: $cmd_string"
	log_gpu "Finished command (exit=$exit_code): $cmd_string"
	
	return $exit_code
}
