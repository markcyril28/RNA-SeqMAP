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
LOGGING_UTILS_SOURCED="true"

# ==============================================================================
# LOGGING CONFIGURATION - IMPORTANT PARAMETERS AT TOP
# ==============================================================================

# Run ID for unique log file naming (prefer printf builtin over date subprocess)
if [[ -z "${RUN_ID:-}" ]]; then
	printf -v RUN_ID '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || RUN_ID=$(date +%Y%m%d_%H%M%S)
fi

# Log directory structure — use PROJECT_ROOT prefix when available so log paths are
# absolute even when the working directory changes (e.g., Nextflow work directories).
# Callers like run_post_processing.sh override these with explicit absolute paths.
# Warn if orchestrator is active but PROJECT_ROOT is unset — logs would go to workDir.
if [[ -n "${WF_MANAGED_ENV:-}" && -z "${PROJECT_ROOT:-}" && -z "${LOG_DIR:-}" ]]; then
	echo "[WARN] logging_utils.sh: WF_MANAGED_ENV is set but PROJECT_ROOT is unset. Logs will be written to CWD ($(pwd))." >&2
fi
_LOG_ROOT="${PROJECT_ROOT:-.}"
LOG_DIR="${LOG_DIR:-${_LOG_ROOT}/logs/log_files}"
TIME_DIR="${TIME_DIR:-${_LOG_ROOT}/logs/time_logs}"
SPACE_DIR="${SPACE_DIR:-${_LOG_ROOT}/logs/space_logs}"
SPACE_TIME_DIR="${SPACE_TIME_DIR:-${_LOG_ROOT}/logs/space_time_logs}"
ERROR_WARN_DIR="${ERROR_WARN_DIR:-${_LOG_ROOT}/logs/error_warn_logs}"
SOFTWARE_CATALOG_DIR="${SOFTWARE_CATALOG_DIR:-${_LOG_ROOT}/logs/software_catalogs}"
GPU_LOG_DIR="${GPU_LOG_DIR:-${_LOG_ROOT}/logs/gpu_log}"
unset _LOG_ROOT

# Log file paths (derived from directories and RUN_ID)
LOG_FILE="${LOG_FILE:-$LOG_DIR/pipeline_${RUN_ID}_full_log.log}"
TIME_FILE="${TIME_FILE:-$TIME_DIR/pipeline_${RUN_ID}_time_metrics.csv}"
TIME_TEMP="${TIME_TEMP:-$TIME_DIR/.time_temp_${RUN_ID}.txt}"
SPACE_FILE="${SPACE_FILE:-$SPACE_DIR/pipeline_${RUN_ID}_space_metrics.csv}"
SPACE_TIME_FILE="${SPACE_TIME_FILE:-$SPACE_TIME_DIR/pipeline_${RUN_ID}_combined_metrics.csv}"
ERROR_WARN_FILE="${ERROR_WARN_FILE:-$ERROR_WARN_DIR/pipeline_${RUN_ID}_errors_warnings.log}"
SOFTWARE_FILE="${SOFTWARE_FILE:-$SOFTWARE_CATALOG_DIR/software_catalog_${RUN_ID}.csv}"
GPU_LOG_FILE="${GPU_LOG_FILE:-$GPU_LOG_DIR/gpu_${RUN_ID}.log}"

# Detect GNU time binary at module load (O(1) per session, avoids per-call PATH lookup)
# Linux: /usr/bin/time; macOS (Homebrew): gtime; fallback: run command without time wrapper
if [[ -z "${_GNU_TIME_CMD:-}" ]]; then
	if /usr/bin/time --version &>/dev/null; then
		_GNU_TIME_CMD="/usr/bin/time"
	elif command -v gtime &>/dev/null && gtime --version &>/dev/null; then
		_GNU_TIME_CMD="gtime"
	elif command -v time &>/dev/null && time --version &>/dev/null; then
		_GNU_TIME_CMD="time"
	else
		_GNU_TIME_CMD=""
	fi
	export _GNU_TIME_CMD
fi

# Logging behavior
log_choice="${log_choice:-1}"  # 1 = tee to console, 2 = file only

# ==============================================================================
# CORE LOGGING FUNCTIONS
# ==============================================================================

# timestamp: use bash built-in printf %(%T)T (no subprocess) with date fallback for bash < 4.2
# Output directly via printf (saves 1 echo subprocess per call vs prior && echo pattern)
timestamp() { printf '%(%Y-%m-%d %H:%M:%S)T\n' -1 2>/dev/null || date '+%Y-%m-%d %H:%M:%S'; }
# Unified log function: O(1) timestamp via bash printf (no subprocess), single source of truth.
# log_warn/log_error previously duplicated timestamp logic; now they call _log_impl.
_log_impl() {
	local level="$1"; shift
	local _ts; printf -v _ts '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || _ts=$(date '+%Y-%m-%d %H:%M:%S')
	local _msg; printf -v _msg '[%s] [%s] %s' "$_ts" "$level" "$*"
	if [[ "$level" == "WARN" || "$level" == "ERROR" ]]; then
		echo "$_msg" >&2
		if [[ -n "${ERROR_WARN_FILE:-}" ]]; then
			# Lazy mkdir: ensure parent dir exists once (avoids fork per call via sentinel)
			if [[ "${_LOG_DIR_VERIFIED:-}" != "true" ]]; then
				if ! mkdir -p "${ERROR_WARN_FILE%/*}" 2>/dev/null; then
					echo "[WARN] Cannot create error log directory: ${ERROR_WARN_FILE%/*}" >&2
				fi
				_LOG_DIR_VERIFIED="true"
			fi
			echo "$_msg" >> "$ERROR_WARN_FILE"
		fi
	else
		echo "$_msg"
	fi
}
log() { _log_impl "$@"; }
log_info() { _log_impl INFO "$@"; }
log_warn() { _log_impl WARN "$@"; }
log_error() { _log_impl ERROR "$@"; }
log_step() { log INFO "=============== $* ==============="; }

strip_ansi_stream() {
	# Strip ANSI escape codes (colors, cursor moves, erase sequences) and
	# carriage returns from a stream so log files remain human-readable.
	# CR (\r) is converted to newline so progress-bar overwrites become
	# separate lines instead of one giant unreadable blob.
	# Single sed pass (GNU sed): \r→\n replacement + ANSI stripping in one process.
	# O(L) where L = number of input lines — saves 1 fork vs prior tr|sed pipeline.
	sed -u $'s/\r/\\\n/g; s/\x1B\\[[0-9;?]*[a-zA-Z]//g; s/\x1B[()][A-Z0-9]//g'
}

# Initialize CSV headers for all log files if they don't exist yet.
# Called by setup_logging() and switch_log_stage() — single source of truth.
_init_csv_headers() {
	[[ ! -f "$TIME_FILE" ]] && echo "Timestamp,Command,Elapsed_Time_sec,CPU_Percent,Max_RSS_KB,User_Time_sec,System_Time_sec,Exit_Status" > "$TIME_FILE"
	[[ ! -f "$SPACE_FILE" ]] && echo "Timestamp,Type,Path,Size_KB,Size_MB,Size_GB,File_Count,Description" > "$SPACE_FILE"
	[[ ! -f "$SPACE_TIME_FILE" ]] && echo "Timestamp,Command,Elapsed_Time_sec,CPU_Percent,Max_RSS_KB,User_Time_sec,System_Time_sec,Input_Size_MB,Output_Size_MB,Exit_Status" > "$SPACE_TIME_FILE"
	[[ ! -f "$ERROR_WARN_FILE" ]] && touch "$ERROR_WARN_FILE"
	[[ ! -f "$SOFTWARE_FILE" ]] && echo "Software/Tool,Version" > "$SOFTWARE_FILE"
	# Inline printf -v avoids $(timestamp) subshell fork
	if [[ ! -f "$GPU_LOG_FILE" ]]; then
		local _gpu_ts; printf -v _gpu_ts '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || _gpu_ts=$(date '+%Y-%m-%d %H:%M:%S')
		printf '=== GPU Log Started: %s ===\n' "$_gpu_ts" > "$GPU_LOG_FILE"
	fi
}

# ==============================================================================
# LOGGING SETUP
# ==============================================================================

# Track background PIDs from process substitutions so switch_log_stage() can
# kill them before spawning replacements (prevents process/fd leak).
_LOGGING_BG_PIDS=()

_logging_cleanup_bg() {
	for _pid in "${_LOGGING_BG_PIDS[@]}"; do
		kill "$_pid" 2>/dev/null && wait "$_pid" 2>/dev/null || true
	done
	_LOGGING_BG_PIDS=()
}

_logging_setup_redirect() {
	# Redirect stdout to /dev/tty (or /dev/null) BEFORE killing the old tee process.
	# Without this, fd 1 still points to the dead pipe after the kill, and any
	# subsequent write (including the exec setup) triggers SIGPIPE → silent exit.
	if [[ ${#_LOGGING_BG_PIDS[@]} -gt 0 ]]; then
		exec > /dev/tty 2>/dev/null || exec > /dev/null 2>&1
		_logging_cleanup_bg
	fi

	if [[ "$log_choice" == "2" ]]; then
		exec > >(strip_ansi_stream >> "$LOG_FILE") 2>&1
	else
		exec > >(tee >(strip_ansi_stream >> "$LOG_FILE")) 2>&1
	fi
	# Capture the PID of the outermost process substitution.
	# $! is set by exec > >(...) in bash.
	[[ -n "$!" ]] && _LOGGING_BG_PIDS+=("$!")
}

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
	# Single find replaces 7 separate rm -f glob expansions
	if [[ "${clear_logs^^}" == "TRUE" ]]; then
		find "$LOG_DIR" "$TIME_DIR" "$SPACE_DIR" "$SPACE_TIME_DIR" "$ERROR_WARN_DIR" \
			"$SOFTWARE_CATALOG_DIR" "$GPU_LOG_DIR" \
			-maxdepth 1 -type f \( -name '*.log' -o -name '*.csv' \) -delete 2>/dev/null || true
		echo "Previous logs cleared"
	fi
	
	# Initialize CSV headers (single source of truth: _init_csv_headers)
	_init_csv_headers

	# Rotate old logs to prevent unbounded growth
	rotate_old_logs "${LOG_DIR%/*}"

	# Set up output redirection via process substitution (strip ANSI from log files).
	# Background PIDs are tracked in _LOGGING_BG_PIDS and cleaned up on stage switch.
	_logging_setup_redirect

	LOGGING_INITIALIZED="true"  # shell-local sentinel — NOT exported to avoid poisoning child process logging
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
	#          switch_log_stage "II_RESULTS/2_ALIGNMENT_RESULTs"
	#          switch_log_stage "II_RESULTS/3_POST_PROC"
	local stage_base="${WF_LOG_BASE:-$1}"

	# Convert to absolute path if relative
	if [[ "$stage_base" != /* ]]; then
		stage_base="${PROJECT_ROOT:-$PWD}/$stage_base"
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

	# Initialize CSV headers (single source of truth: _init_csv_headers)
	_init_csv_headers

	# Reset lazy-mkdir sentinel so _log_impl creates the new stage's error log dir
	_LOG_DIR_VERIFIED=""

	# Re-setup output redirection to the new log file (reuses _logging_setup_redirect
	# to kill previous background processes before spawning new ones)
	_logging_setup_redirect

	log_info "Switched logging to stage: $stage_base"
}

# ==============================================================================
# ERROR HANDLING
# ==============================================================================

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
	# Big O: O(L) where L = number of output lines; timestamp is O(1) per line.
	# Performance: uses gawk systime()/strftime() when available (zero subprocess
	# spawns) with automatic fallback to date(1) for mawk/nawk (caches timestamp
	# per epoch-second, reducing spawns from L to ~1 for burst output).
	awk -v err_file="$ERROR_WARN_FILE" \
		-v err_pat="$_ERROR_PATTERN" \
		-v warn_pat="$_WARN_PATTERN" '
	BEGIN {
		ts = ""; last_epoch = 0
		# Detect gawk: PROCINFO is gawk-only; on mawk/nawk this is an undefined array
		# whose elements return "" -- so the check safely falls through to has_systime=0.
		has_systime = 0
		if (PROCINFO["version"] != "") has_systime = 1
	}
	function get_ts() {
		if (has_systime) {
			# gawk path: zero subprocess spawns — O(1) built-in call
			epoch = systime()
			if (epoch != last_epoch) {
				last_epoch = epoch
				ts = strftime("%Y-%m-%d %H:%M:%S", epoch)
			}
		} else {
			# POSIX fallback: spawn date only when epoch second changes
			cmd = "date +\"%Y-%m-%d %H:%M:%S %s\""
			cmd | getline raw_ts
			close(cmd)
			n = split(raw_ts, parts, " ")
			epoch = parts[n] + 0
			if (epoch != last_epoch) {
				last_epoch = epoch
				ts = parts[1]
				for (i = 2; i < n; i++) ts = ts " " parts[i]
			}
		}
	}
	{
		print
		fflush()
		low = tolower($0)
		is_err = match(low, err_pat)
		is_warn = match(low, warn_pat)
		if (is_err || is_warn) {
			get_ts()
			if (is_err) {
				printf "[%s] [ERROR] %s\n", ts, $0 >> err_file
			}
			if (is_warn) {
				printf "[%s] [WARN] %s\n", ts, $0 >> err_file
			}
			fflush(err_file)
		}
	}'
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
	local start_ts; printf -v start_ts '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || start_ts=$(date '+%Y-%m-%d %H:%M:%S')
	
	# Measure input size before running command
	# O(1) stat for single files (avoids du|awk 2-process pipeline); du for directories
	# O(1) bash arithmetic — replaces awk subprocess for simple division
	local input_size_mb="0"
	if [[ -n "$input_path" ]]; then
		if [[ -f "$input_path" ]]; then
			local _sz_bytes
			_sz_bytes=$(stat -c%s "$input_path" 2>/dev/null || stat -f%z "$input_path" 2>/dev/null || echo 0)
			# Scaled integer: (bytes * 100 / 1048576) then insert decimal point
			local _scaled=$(( _sz_bytes * 100 / 1048576 ))
			local _frac; printf -v _frac '%02d' $(( _scaled % 100 ))
			input_size_mb="$(( _scaled / 100 )).$_frac"
		elif [[ -d "$input_path" ]]; then
			local _du_kb
			# Process substitution avoids herestring buffering for large du output
			read -r _du_kb _ < <(du -sk "$input_path" 2>/dev/null)
			_du_kb="${_du_kb:-0}"
			local _scaled=$(( _du_kb * 100 / 1024 ))
			local _frac; printf -v _frac '%02d' $(( _scaled % 100 ))
			input_size_mb="$(( _scaled / 100 )).$_frac"
		fi
	fi

	mkdir -p "$TIME_DIR" || { log_error "Failed to create TIME_DIR: $TIME_DIR"; return 1; }
	
	local exit_code=0

	# Log abbreviated command before running (full command saved in CSV)
	local cmd_abbrev="${cmd_string:0:120}"
	[[ ${#cmd_string} -gt 120 ]] && cmd_abbrev="${cmd_abbrev}..."
	log_info "[CMD] $cmd_abbrev"

	# Write begin marker directly to log file (preserves ordering with tool stdout)
	# Inline printf -v avoids $(timestamp) subshell forks — O(1) each, saves 2 forks per call
	local _ts_begin; printf -v _ts_begin '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || _ts_begin=$(date '+%Y-%m-%d %H:%M:%S')
	printf '[%s] [INFO] --- BEGIN TOOL OUTPUT: %s ---\n' "$_ts_begin" "${1##*/}" >> "$LOG_FILE"
	# Strip ANSI escape codes and carriage returns before writing to log (e.g. Salmon progress bars)
	# Use detected GNU time (_GNU_TIME_CMD); fall back to running command without time wrapper
	if [[ -n "${_GNU_TIME_CMD:-}" ]]; then
		"$_GNU_TIME_CMD" -v "$@" 2>"$TIME_TEMP" | strip_ansi_stream >> "$LOG_FILE"
	else
		"$@" 2>"$TIME_TEMP" | strip_ansi_stream >> "$LOG_FILE"
	fi
	exit_code=${PIPESTATUS[0]}
	local _ts_end; printf -v _ts_end '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || _ts_end=$(date '+%Y-%m-%d %H:%M:%S')
	printf '[%s] [INFO] --- END TOOL OUTPUT: %s (exit=%d) ---\n' "$_ts_end" "${1##*/}" "$exit_code" >> "$LOG_FILE"

	# Single-pass extraction of all metrics from GNU time output (if available).
	# Replaces 6 separate grep|awk pipelines (12 process spawns) with 1 awk process.
	local elapsed_raw cpu_raw rss_raw elapsed_time cpu_percent max_rss user_time system_time
	eval "$(awk '
	/Elapsed \(wall clock\)/ {
		raw = $NF
		# Convert h:mm:ss or m:ss to seconds
		n = split(raw, t, ":")
		if (n == 3) sec = t[1]*3600 + t[2]*60 + t[3]
		else if (n == 2) sec = t[1]*60 + t[2]
		else sec = raw
		printf "elapsed_raw=%s elapsed_time=%s ", raw, sec
	}
	/Percent of CPU/ {
		v = $NF; gsub(/%/, "", v)
		printf "cpu_raw=%s cpu_percent=%s ", $NF, v
	}
	/Maximum resident set size/ { printf "rss_raw=%s max_rss=%s ", $NF, $NF }
	/User time/    { printf "user_time=%s ", $NF }
	/System time/  { printf "system_time=%s ", $NF }
	' "$TIME_TEMP" 2>/dev/null)"

	log_info "[RESOURCES] Elapsed: ${elapsed_raw:-N/A} | CPU: ${cpu_raw:-N/A} | MaxRSS: ${rss_raw:-0} KB | Exit: $exit_code"

	# On failure: dump full time output for debugging
	if [[ $exit_code -ne 0 ]]; then
		local _ts_dbg; printf -v _ts_dbg '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || _ts_dbg=$(date '+%Y-%m-%d %H:%M:%S')
		printf '[%s] [DEBUG] --- TIME OUTPUT (failure details) ---\n' "$_ts_dbg" >> "$LOG_FILE"
		cat "$TIME_TEMP" >> "$LOG_FILE" 2>&1
	fi

	# Capture errors/exceptions to error log (single grep pass instead of two)
	local _err_lines=""
	_err_lines=$(grep -iE "$_ERROR_PATTERN" "$TIME_TEMP" 2>/dev/null) || true
	if [[ $exit_code -ne 0 ]] || [[ -n "$_err_lines" ]]; then
		{
			printf '[%s] [ERROR] Command failed (exit=%d): %s\n' "$_ts_end" "$exit_code" "$cmd_string"
			[[ -n "$_err_lines" ]] && printf '%s\n' "$_err_lines"
		} >> "$ERROR_WARN_FILE"
	fi
	
	# Measure output size after running command
	# O(1) stat for single files; du for directories
	# O(1) bash arithmetic — replaces awk subprocess for simple division
	local output_size_mb="0"
	if [[ -n "$output_path" ]]; then
		if [[ -f "$output_path" ]]; then
			local _sz_bytes
			_sz_bytes=$(stat -c%s "$output_path" 2>/dev/null || stat -f%z "$output_path" 2>/dev/null || echo 0)
			local _scaled=$(( _sz_bytes * 100 / 1048576 ))
			local _frac; printf -v _frac '%02d' $(( _scaled % 100 ))
			output_size_mb="$(( _scaled / 100 )).$_frac"
		elif [[ -d "$output_path" ]]; then
			local _du_kb
			# Process substitution avoids herestring buffering for large du output
			read -r _du_kb _ < <(du -sk "$output_path" 2>/dev/null)
			_du_kb="${_du_kb:-0}"
			local _scaled=$(( _du_kb * 100 / 1024 ))
			local _frac; printf -v _frac '%02d' $(( _scaled % 100 ))
			output_size_mb="$(( _scaled / 100 )).$_frac"
		fi
	fi

	# Append to both CSV files in a single write (reduces 2 file opens to 1 process)
	# O(1) — single printf with two redirect targets via tee replacement
	local csv_cmd="${cmd_string//\"/\"\"}"
	local _time_row="${start_ts},\"${csv_cmd}\",${elapsed_time:-0},${cpu_percent:-0},${max_rss:-0},${user_time:-0},${system_time:-0},${exit_code}"
	local _st_row="${start_ts},\"${csv_cmd}\",${elapsed_time:-0},${cpu_percent:-0},${max_rss:-0},${user_time:-0},${system_time:-0},${input_size_mb},${output_size_mb},${exit_code}"
	printf '%s\n' "$_time_row" >> "$TIME_FILE"
	printf '%s\n' "$_st_row" >> "$SPACE_TIME_FILE"
	
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
	
	# O(1) for files via stat (avoids du fork); O(F) for directories via du.
	# Bash arithmetic computes MB/GB without awk subprocess.
	local size_kb=0 size_mb="0.00" size_gb="0.00"
	if [[ -f "$file_path" ]]; then
		local _sz_bytes
		_sz_bytes=$(stat -c%s "$file_path" 2>/dev/null || stat -f%z "$file_path" 2>/dev/null || echo 0)
		size_kb=$(( (_sz_bytes + 1023) / 1024 ))
	elif [[ -d "$file_path" ]]; then
		read -r size_kb _ < <(du -sk "$file_path" 2>/dev/null)
		size_kb="${size_kb:-0}"
	fi
	# Bash integer math with 2-decimal precision (avoids awk subprocess)
	local _mb_scaled=$(( size_kb * 100 / 1024 ))
	local _mb_frac; printf -v _mb_frac '%02d' $(( _mb_scaled % 100 ))
	size_mb="$(( _mb_scaled / 100 )).$_mb_frac"
	local _gb_scaled=$(( size_kb * 100 / 1048576 ))
	local _gb_frac; printf -v _gb_frac '%02d' $(( _gb_scaled % 100 ))
	size_gb="$(( _gb_scaled / 100 )).$_gb_frac"

	local file_count="-"
	# POSIX-compatible file counting (find -printf is GNU-only, fails on macOS/BSD)
	[[ -d "$file_path" ]] && file_count=$(find "$file_path" -type f 2>/dev/null | wc -l)
	
	local ts; printf -v ts '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || ts=$(date '+%Y-%m-%d %H:%M:%S')
	echo "${ts},${type},\"${file_path}\",${size_kb},${size_mb},${size_gb},${file_count},\"${description}\"" >> "$SPACE_FILE"
	log_info "Space logged: $file_path = ${size_mb}MB"
}

# log_input_output_size() — removed (dead code; run_with_space_time_log calls log_file_size directly)

# ==============================================================================
# SOFTWARE CATALOG FUNCTIONS
# ==============================================================================

# log_software_version() — removed (dead code; catalog_all_software writes directly)

catalog_all_software() {
	# Catalog versions of all bioinformatics tools, R packages, and provenance metadata
	log_step "Cataloging software versions"

	# Record pipeline git commit SHA for provenance
	if command -v git >/dev/null 2>&1; then
		# Prefer PROJECT_ROOT (always correct under orchestrators); fall back to
		# BASH_SOURCE-relative navigation for standalone bash mode.
		local _repo_root="${PROJECT_ROOT:-}"
		if [[ -z "$_repo_root" ]]; then
			_repo_root="$(cd "${BASH_SOURCE[0]%/*}/../.." 2>/dev/null && pwd)" || _repo_root=""
		fi
		[[ -z "$_repo_root" ]] && return 0
		local git_sha
		git_sha=$(git -C "$_repo_root" rev-parse --short HEAD 2>/dev/null || echo "not_a_git_repo")
		local git_branch
		git_branch=$(git -C "$_repo_root" rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
		local git_dirty=""
		if ! git -C "$_repo_root" diff --quiet HEAD 2>/dev/null; then
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
		"R:R --version 2>&1 | awk '/^R version/{print \$3;exit}'"
	)

	# O(F) wall-clock (parallel) vs O(F) sequential — F tool version subshells run concurrently;
	# each subshell is O(1) process spawn; wait loop is O(F) to reap all background pids
	# Batch version checks: collect all installed tools and run version commands
	# concurrently via background subshells (reduces ~18 sequential spawns to parallel)
	local _ver_tmpdir
	_ver_tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/ver_check.XXXXXX")
	trap 'rm -rf "$_ver_tmpdir" 2>/dev/null' RETURN
	local _ver_pids=()
	for tool_cmd in "${tools[@]}"; do
		local tool="${tool_cmd%%:*}"
		local cmd="${tool_cmd#*:}"

		if command -v "${cmd%% *}" >/dev/null 2>&1; then
			# Single awk replaces head -n1 | awk (2 processes → 1 per tool)
			( eval "$cmd" 2>&1 | awk -v t="$tool" 'NR==1{print t","$NF; exit}' > "$_ver_tmpdir/$tool" ) &
			_ver_pids+=($!)
		else
			echo "${tool},not_installed" >> "$SOFTWARE_FILE"
			log_info "Software not found: $tool"
		fi
	done
	# Wait for all background version checks
	for _pid in "${_ver_pids[@]}"; do wait "$_pid" 2>/dev/null; done
	# Append results in deterministic order
	for tool_cmd in "${tools[@]}"; do
		local tool="${tool_cmd%%:*}"
		if [[ -f "$_ver_tmpdir/$tool" ]]; then
			local _ver_line
			_ver_line=$(<"$_ver_tmpdir/$tool")
			local _ver="${_ver_line#*,}"
			echo "$_ver_line" >> "$SOFTWARE_FILE"
			log_info "Software version: $tool = ${_ver:-unknown}"
		fi
	done
	rm -rf "$_ver_tmpdir"

	# Catalog key R/Bioconductor packages used by analysis modules
	# Single Rscript call replaces 16 separate invocations (~30s → ~2s)
	if command -v Rscript >/dev/null 2>&1; then
		log_info "Cataloging R package versions..."
		local session_info_file
		session_info_file="${SOFTWARE_FILE%/*}/R_sessionInfo_${RUN_ID}.txt"
		Rscript --vanilla -e "
pkgs <- c('DESeq2','tximport','tximeta','WGCNA','clusterProfiler','ComplexHeatmap',
          'ballgown','AnnotationDbi','enrichplot','DOSE','fgsea',
          'pheatmap','ggplot2','corrplot','dendextend','gridExtra','scales')
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) 'not_installed')
  cat(paste0('R/', p, ',', v), sep = '\n')
}
tryCatch(writeLines(capture.output(sessionInfo()), '$session_info_file'),
         error = function(e) message('sessionInfo capture failed'))
" >> "$SOFTWARE_FILE" 2>/dev/null \
			&& log_info "R package versions cataloged; sessionInfo saved to: $session_info_file" \
			|| log_warn "Failed to catalog R packages"
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
	local base_dir="${1:-${LOG_DIR%/*}}"
	local max_age="${MAX_LOG_AGE_DAYS:-30}"

	[[ ! -d "$base_dir" ]] && return 0

	# POSIX-compatible: count matching files, then delete in a single batched exec.
	# (Replaces GNU-only find -printf/-delete combo unavailable on macOS/BSD.)
	local count
	count=$(find "$base_dir" -type f \( -name '*.log' -o -name '*.csv' \) -mtime +"$max_age" 2>/dev/null | wc -l)
	if [[ "$count" -gt 0 ]]; then
		find "$base_dir" -type f \( -name '*.log' -o -name '*.csv' \) -mtime +"$max_age" -exec rm -f {} + 2>/dev/null
		log_info "Log rotation: removed $count files older than ${max_age} days from $base_dir"
	fi
}

# ==============================================================================
# GPU LOGGING FUNCTIONS
# ==============================================================================

# Cache nvidia-smi availability once (avoids repeated PATH lookups in GPU log functions).
# Guard: skip if already cached by gpu_utils.sh to avoid overwriting its result.
if [[ -z "${_HAS_NVIDIA_SMI+x}" ]]; then
    _HAS_NVIDIA_SMI=false
    command -v nvidia-smi &>/dev/null && _HAS_NVIDIA_SMI=true
fi

log_gpu() {
	# Log GPU-related message to GPU log file
	# Usage: log_gpu "message"
	local message="$*"
	# Inline timestamp — avoids $(timestamp) subshell fork. O(1).
	local _ts; printf -v _ts '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || _ts=$(date '+%Y-%m-%d %H:%M:%S')
	printf '[%s] %s\n' "$_ts" "$message" >> "$GPU_LOG_FILE"
	log_info "[GPU] $message"
}

