#!/bin/bash
# ==============================================================================
# SHARED RUNTIME DEFAULTS
# ==============================================================================
# Common thread calculation and conda activation used by all Stage 1-2 configs.
# Sourced AFTER per-config THREADS/JOBS/USE_GNU_PARALLEL are set.
#
# Usage in config files:
#   THREADS=64; JOBS=2; USE_GNU_PARALLEL="TRUE"
#   source "config/shared/runtime_defaults.sh"
# ==============================================================================

[[ "${_RUNTIME_DEFAULTS_SOURCED:-}" == "true" ]] && return 0
export _RUNTIME_DEFAULTS_SOURCED="true"

# ==============================================================================
# CONDA ENVIRONMENT
# ==============================================================================

# Skip conda hook if already in the gea environment (~0.3-0.5s saved per invocation)
if [[ "${CONDA_DEFAULT_ENV:-}" != "gea" ]]; then
	eval "$(conda shell.bash hook 2>/dev/null)" 2>/dev/null || true
	conda activate gea 2>/dev/null || true
fi

# ==============================================================================
# SOURCE MODULES
# ==============================================================================

# Resolve relative to this file's location (two levels up from config/shared/)
_RUNTIME_DEFAULTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_RUNTIME_DEFAULTS_DIR}/../../modules/modules_loader.sh"
unset _RUNTIME_DEFAULTS_DIR

# ==============================================================================
# AUTO-JOBS RESOLUTION — O(1)
# ==============================================================================
# When JOBS="auto", calculates optimal parallel job count from THREADS and the
# optimal thread count per program.  Called from this file (Stage 1-2 default)
# and from each method script with a program-specific optimal thread count.
#
# Usage:
#   _resolve_auto_jobs [optimal_threads_per_job]
#   Result is stored in JOBS (numeric).
#
# Known optimal thread counts per program:
#   TrimGalore       4     (I/O-bound adapter trimming)
#   Salmon quant     8     (EM algorithm, diminishing returns beyond 8)
#   HISAT2          16     (scales well up to ~16)
#   StringTie       16     (scales well up to ~16)
#   STAR            16     (scales well up to ~16)
#   RSEM/Bowtie2     8     (memory-constrained, ~2.5GB/thread)
#   Samtools sort    4     (I/O-bound, capped at 4 everywhere)
# ==============================================================================

# Default optimal threads when no program-specific value is given
OPTIMAL_THREADS_PER_JOB="${OPTIMAL_THREADS_PER_JOB:-8}"

_resolve_auto_jobs() {
	local optimal="${1:-${OPTIMAL_THREADS_PER_JOB:-8}}"
	local threads="${THREADS:-$(nproc 2>/dev/null || echo 12)}"

	if [[ "${JOBS:-}" == "auto" || "${JOBS:-}" == "AUTO" ]]; then
		JOBS=$(( threads / optimal ))
		(( JOBS < 1 )) && JOBS=1
		log_info "[AUTO-JOBS] THREADS=$threads / optimal=$optimal → JOBS=$JOBS"
	fi
	# Ensure JOBS is numeric after resolution
	[[ "$JOBS" =~ ^[0-9]+$ ]] || JOBS=2
}

# ==============================================================================
# THREAD CALCULATION — O(1)
# ==============================================================================
# Derives THREADS_PER_JOB from THREADS and JOBS based on parallelism mode.
# Each config sets its own THREADS/JOBS values before sourcing this file.
# ==============================================================================

# Resolve JOBS="auto" before computing THREADS_PER_JOB
_resolve_auto_jobs "${OPTIMAL_THREADS_PER_JOB:-8}"

if [[ "${USE_GNU_PARALLEL:-FALSE}" == "TRUE" ]]; then
	THREADS_PER_JOB=$((${THREADS:-4} / ${JOBS:-1}))
	[[ $THREADS_PER_JOB -lt 1 ]] && THREADS_PER_JOB=1
else
	THREADS_PER_JOB="${THREADS:-4}"
fi
export THREADS JOBS OPTIMAL_THREADS_PER_JOB USE_GNU_PARALLEL THREADS_PER_JOB keep_bam_global
