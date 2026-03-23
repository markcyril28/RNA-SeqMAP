#!/bin/bash
# ==============================================================================
# GPU UTILITIES
# ==============================================================================
# GPU detection for pipeline acceleration (standalone; not sourced by modules_loader.sh).
# R-side GPU detection is in 0_shared_config.R; this file is only needed for
# standalone GPU setup workflows (e.g., z_archive/prep_gpu.sh).
# ==============================================================================

# Guard against double-sourcing
[[ "${GPU_UTILS_SOURCED:-}" == "true" ]] && return 0
export GPU_UTILS_SOURCED="true"

# ==============================================================================
# GPU CONFIGURATION
# ==============================================================================

GPU_AVAILABLE="false"
GPU_COUNT=0
GPU_MEMORY_MB=0
CUDA_VERSION=""
CUDA_READY="false"
GPU_VENDOR=""
DISTRO=""
DISTRO_VERSION=""

# ==============================================================================
# LOGGING CONFIGURATION
# ==============================================================================

GPU_LOG_DIR="${GPU_LOG_DIR:-${SCRIPT_DIR:-$PWD}/logs}"
# Prefer printf builtin over date subprocess for log filename
if [[ -z "${GPU_LOG_FILE:-}" ]]; then
	printf -v _gpu_ts_id '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || _gpu_ts_id=$(date +%Y%m%d_%H%M%S)
	GPU_LOG_FILE="$GPU_LOG_DIR/gpu_prep_${_gpu_ts_id}.log"
	unset _gpu_ts_id
fi
GPU_LOG_ENABLED="${GPU_LOG_ENABLED:-true}"
GPU_LOG_LEVEL="${GPU_LOG_LEVEL:-INFO}"  # DEBUG, INFO, WARN, ERROR

# Ensure log directory exists
mkdir -p "$GPU_LOG_DIR" 2>/dev/null || true

# ==============================================================================
# LOGGING FUNCTIONS
# ==============================================================================

# Get timestamp for logging (prefer bash printf to avoid date subprocess)
# Reuse logging_utils.sh timestamp() if available; define standalone fallback otherwise.
if ! declare -f timestamp &>/dev/null; then
	# O(1) — printf directly to stdout (avoids echo subprocess for bash ≥ 4.2; date fallback for older)
	_get_timestamp() {
		printf '%(%Y-%m-%d %H:%M:%S)T\n' -1 2>/dev/null || date "+%Y-%m-%d %H:%M:%S"
	}
else
	_get_timestamp() { timestamp; }
fi

# Write to log file
# O(1) — uses printf -v to avoid subshell fork for timestamp on every call
_write_log() {
	local level="$1"
	local message="$2"
	if [[ "$GPU_LOG_ENABLED" == "true" && -n "$GPU_LOG_FILE" ]]; then
		local _ts; printf -v _ts '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || _ts=$(_get_timestamp)
		printf '[%s] [%s] %s\n' "$_ts" "$level" "$message" >> "$GPU_LOG_FILE" 2>/dev/null || true
	fi
}

# Log level check (returns 0 if should log)
_should_log() {
	local level="$1"
	case "$GPU_LOG_LEVEL" in
		DEBUG) return 0 ;;
		INFO)  [[ "$level" != "DEBUG" ]] && return 0 ;;
		WARN)  [[ "$level" == "WARN" || "$level" == "ERROR" ]] && return 0 ;;
		ERROR) [[ "$level" == "ERROR" ]] && return 0 ;;
	esac
	return 1
}

# Logging stubs (only when gpu_utils.sh is sourced standalone without logging_utils.sh).
# Single declare -f check replaces 5 separate checks; stubs are minimal wrappers
# that delegate to _write_log for GPU log file output.
if ! declare -f log_info &>/dev/null; then
	log_debug() { ! _should_log "DEBUG" || { echo "[DEBUG] $*"; _write_log "DEBUG" "$*"; }; }
	log_info()  { ! _should_log "INFO"  || { echo "[INFO] $*";  _write_log "INFO"  "$*"; }; }
	log_warn()  { ! _should_log "WARN"  || { echo "[WARN] $*" >&2; _write_log "WARN"  "$*"; }; }
	log_error() { ! _should_log "ERROR" || { echo "[ERROR] $*" >&2; _write_log "ERROR" "$*"; }; }
	log_step()  { echo "=== $* ==="; _write_log "STEP" "$*"; }
fi

# Log system info (called lazily via _ensure_gpu_detected)
_log_system_info() {
	_write_log "INFO" "========================================"
	_write_log "INFO" "GPU Utils initialized"
	_write_log "INFO" "Host: ${HOSTNAME:-unknown}"
	_write_log "INFO" "User: ${USER:-unknown}"
	_write_log "INFO" "Kernel: $(uname -r 2>/dev/null || echo unknown)"
	_write_log "INFO" "WSL: $(is_wsl && echo 'yes' || echo 'no')"
	_write_log "INFO" "========================================"
}

# ==============================================================================
# BINARY AVAILABILITY CACHE
# ==============================================================================
# Cache command -v results at module load time (O(1) per subsequent check).
# Reuse _HAS_NVIDIA_SMI from logging_utils.sh if already cached (avoids duplicate PATH lookup).
[[ -z "${_HAS_NVIDIA_SMI+x}" ]] && { _HAS_NVIDIA_SMI=false; command -v nvidia-smi &>/dev/null && _HAS_NVIDIA_SMI=true; }
_HAS_NVCC=false; command -v nvcc &>/dev/null && _HAS_NVCC=true

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Check if running in WSL — O(1) bash built-in, avoids grep subprocess
is_wsl() {
	if [[ -f /proc/version ]]; then
		local _ver
		_ver=$(</proc/version)
		_ver="${_ver,,}"  # lowercase (bash 4+)
		[[ "$_ver" == *microsoft* || "$_ver" == *wsl* ]]
	else
		return 1
	fi
}

# ==============================================================================
# GPU DETECTION
# ==============================================================================

detect_gpu() {
	GPU_AVAILABLE="false"
	GPU_COUNT=0
	CUDA_READY="false"

	# Single nvidia-smi invocation: query-gpu for metrics, parse CUDA from header
	if $_HAS_NVIDIA_SMI; then
		local _smi_full
		_smi_full=$(nvidia-smi 2>/dev/null)
		if [[ -n "$_smi_full" ]]; then
			# Single awk pass extracts CUDA version, GPU count, and memory from nvidia-smi output
			# Replaces 3 separate echo|grep|head/tail pipelines (9 processes → 1)
			# NOTE: match(s, re, a) capture groups require gawk; mawk silently yields empty
			# values, triggering the query-gpu fallback at the cost of one extra nvidia-smi call.
			local _cuda_ver _gpu_cnt _gpu_mem
			read -r _cuda_ver _gpu_cnt _gpu_mem < <(awk '
				/CUDA Version:/ && !cuda { match($0, /CUDA Version: ([0-9.]+)/, a); if (a[1]) cuda=a[1] }
				/MiB \/ [0-9]+MiB/ { cnt++; match($0, /([0-9]+)MiB \|/, a); if (a[1]) mem=a[1] }
				END { print (cuda ? cuda : ""), cnt+0, (mem ? mem : 0) }
			' <<< "$_smi_full")
			CUDA_VERSION="${_cuda_ver}"
			GPU_AVAILABLE="true"
			GPU_COUNT="${_gpu_cnt}"
			GPU_MEMORY_MB="${_gpu_mem}"

			# Fallback: if gawk capture groups failed (mawk/nawk), recover CUDA version
			# from the already-captured nvidia-smi output without an extra call.
			# Fallback: mawk/nawk lack capture groups — extract CUDA version with gsub
			if [[ -z "$CUDA_VERSION" ]]; then
				CUDA_VERSION=$(awk '/CUDA Version:/ {for(i=1;i<=NF;i++) if($i=="Version:") {gsub(/[[:space:]]/,"",$((i+1))); print $(i+1); exit}}' <<< "$_smi_full")
			fi

			# Fallback: if GPU count/memory parsing failed, use query-gpu (one extra call)
			if [[ "$GPU_COUNT" -eq 0 || "$GPU_MEMORY_MB" -eq 0 ]] 2>/dev/null; then
				local _query_info
				_query_info=$(nvidia-smi --query-gpu=count,memory.total,name,driver_version --format=csv,noheader,nounits 2>/dev/null)
				if [[ -n "$_query_info" ]]; then
					# Bash builtin line count avoids grep subprocess
				local -a _gpu_lines; mapfile -t _gpu_lines <<< "$_query_info"
				GPU_COUNT="${#_gpu_lines[@]}"
					# Single awk replaces head|cut|tr 3-process pipe
					GPU_MEMORY_MB=$(awk -F',' 'NR==1 {gsub(/[[:space:]]/,"",$2); print $2}' <<< "$_query_info")
				fi
			fi

			# Check if CUDA toolkit is ready (cached at module load)
			if $_HAS_NVCC; then
				CUDA_READY="true"
			fi
		fi
	fi

	export GPU_AVAILABLE GPU_COUNT GPU_MEMORY_MB CUDA_VERSION CUDA_READY
}

# Check if GPU is available (triggers lazy detection on first call)
has_gpu() {
	_ensure_gpu_detected
	[[ "$GPU_AVAILABLE" == "true" ]]
}

# Check if CUDA is ready for use (triggers lazy detection on first call)
is_cuda_ready() {
	_ensure_gpu_detected
	[[ "$CUDA_READY" == "true" ]]
}

# Log GPU status (triggers lazy detection on first call)
log_gpu_status() {
	_ensure_gpu_detected
	if has_gpu; then
		log_info "[GPU] Detected $GPU_COUNT GPU(s), ${GPU_MEMORY_MB}MB memory, CUDA: $CUDA_VERSION"
		if is_cuda_ready; then
			log_info "[GPU] CUDA toolkit ready (nvcc available)"
		else
			log_warn "[GPU] CUDA toolkit not installed - run setup_cuda for GPU acceleration"
		fi
	else
		log_info "[GPU] No GPU detected - using CPU only"
	fi
}

# Dead installation/setup functions removed (check_root, detect_distro, detect_gpu_type,
# install_base_packages, check_cuda_prerequisites, setup_cuda, install_cuda_apt/yum/conda,
# configure_cuda_env, configure_cuda_env_standalone, install_nvidia_driver,
# install_cuda_toolkit, install_cudnn, install_nvidia_monitoring, setup_nvidia/amd/intel,
# verify_installation, has_gpu_star, check_gpu_tools, get_star_genome_load,
# get_optimal_threads, gpu_prep_main) — none called by any active pipeline code.
# Previously used by z_archive/prep_gpu.sh which is archived.

# ==============================================================================
# LAZY GPU DETECTION
# ==============================================================================
# GPU detection is deferred until first use (has_gpu / is_cuda_ready / log_gpu_status).
# This avoids ~1-2s of nvidia-smi + system-info overhead on every source of
# modules_loader.sh — which matters when GPU is never used (most pipeline runs).

_GPU_DETECTED=false

_ensure_gpu_detected() {
	[[ "$_GPU_DETECTED" == "true" ]] && return 0
	_GPU_DETECTED=true
	_log_system_info
	detect_gpu
	_write_log "INFO" "GPU_AVAILABLE=$GPU_AVAILABLE, GPU_COUNT=$GPU_COUNT, GPU_MEMORY_MB=$GPU_MEMORY_MB"
	_write_log "INFO" "CUDA_VERSION=$CUDA_VERSION, CUDA_READY=$CUDA_READY"
}
