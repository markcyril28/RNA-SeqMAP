#!/bin/bash

#===============================================================================
# COMBINED CONDA/MAMBA ENVIRONMENT SETUP SCRIPT
#===============================================================================
# Creates/updates the conda environment named by setup_conda_geaManifest.yml, with all
# dependencies for:
#   - RNA-seq alignment & preprocessing       (GEA pipeline)
#   - Post-processing & statistical analysis  (R / Bioconductor)
#
# Spec     : setup_conda_geaManifest.yml       (hand-authored — edit packages HERE)
#            A real conda environment file, and the single source of truth for
#            everything this script needs to build the environment:
#              name:                    -> which environment to build
#              dependencies: - python=  -> which interpreter to solve for
#              dependencies: - <pkg>=   -> the direct, pinned package set
#            This script keeps no copy of any of it.
# Lockfile : setup_conda_gea.lock.yml  (generated — never edit)
#            The full resolved set the spec above solved to, transitive deps
#            included. Rewritten after every successful run.
# Logs     : logs/installation/setup_conda_gea_<RUN_ID>_*.log
#
# Quick start:
#   1. Ensure Conda/Miniconda is installed
#   2. Run     : bash setup_conda_gea.sh
#   3. Activate: conda activate "$(sed -n 's/^name:[[:space:]]*//p' setup_conda_geaManifest.yml | cut -d'#' -f1 | tr -d ' ')"
#
# To add, drop or repin a package: edit the 'dependencies:' list in
# setup_conda_geaManifest.yml and re-run this script. Nothing here needs changing.
#
# Run 'bash setup_conda_gea.sh --help' for the full flag list.
#
#-------------------------------------------------------------------------------
# Known benign log lines — upstream packaging defects, not faults of this script
# and not fixable from this repository. The issue log is deliberately
# over-inclusive, so it lists them under its KNOWN-BENIGN heading and keeps them
# out of the actionable count (see _BENIGN_PATTERN). Check them off rather than
# chasing them. All four were verified against a built environment on 2026-08-27.
#
#   "[kingfisher-*] The following files were already present in the environment:
#    - README.md"   (conda phrasing: ClobberError ... 'README.md')
#       kingfisher and argparse-manpage-birdtools both ship a top-level README.md
#       into $PREFIX. Only that file is affected; kingfisher 0.4.1 runs.
#
#   "[ossuuid-*] The following files were already present in the environment:
#    - lib/libuuid.a, lib/libuuid.so, lib/pkgconfig/uuid.pc"
#       ossuuid (pulled in by sra-tools) vs libuuid (util-linux). Only the
#       unversioned developer names collide, and which package ends up owning
#       them depends on link order — ossuuid won in the 'gea' prefix, libuuid in
#       'gea_test'. Harmless either way: both runtime SONAMEs are always present
#       and correct (libuuid.so.1 from util-linux, libuuid.so.16 from ossuuid),
#       and sra-tools does not link uuid dynamically at all — prefetch and
#       fasterq-dump 3.2.1 both run. Solving in one transaction (see SOLVE_SPECS)
#       downgrades this from conda's ClobberError to a mamba warning; it does not
#       remove it, because the two packages genuinely claim the same paths.
#
#   SafetyError / "Invalid package cache" on
#   r-base-*/lib/R/doc/html/packages.html
#       conda hardlinks that file from the shared package cache into the prefix,
#       and R regenerates it in place on every library install — so the write
#       lands on the cache's own inode and the recorded size stops matching.
#       Harmless (both solvers still use the cache), but it recurs whenever a
#       fresh environment is built. To silence it for a while:
#           rm -rf "$(conda info --base)"/pkgs/r-base-4.3.3-*/
#       Existing environments are unaffected — they hold hardlinks to the same
#       inodes, so only the cache copy goes away and r-base is re-downloaded the
#       next time an environment is created from scratch.
#
#   "Failed to parse field 'noarch' ... perl-net-http-6.19-*"
#       Malformed record in bioconda's sharded repodata. The package is not even
#       installed here; the solver ignores the field and moves on.
#
# One clobber phrasing stays ACTIONABLE on purpose: conda's "This transaction has
# incompatible packages due to a shared path." names the packages on a following
# line, which the issue-log filter does not capture, so it cannot be attributed
# to a known-benign pair. A clobber between arbitrary packages should be looked at.
#-------------------------------------------------------------------------------
#===============================================================================

set -euo pipefail

#===============================================================================
# CONFIGURATION
#===============================================================================

# Parameter expansion avoids nested $(dirname) subshell fork
SCRIPT_DIR="${BASH_SOURCE[0]%/*}"
[[ "$SCRIPT_DIR" == "${BASH_SOURCE[0]}" ]] && SCRIPT_DIR="."
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)" || { echo "[ERROR] setup_conda_gea.sh: Failed to resolve script directory" >&2; exit 1; }

# Single source of truth for consumed/generated paths.
# ENV_SPEC is read and never written; LOCKFILE is written and never read.
ENV_SPEC="$SCRIPT_DIR/setup_conda_geaManifest.yml"
LOCKFILE="$SCRIPT_DIR/setup_conda_gea.lock.yml"
R_INSTALLER="$SCRIPT_DIR/modules_gea/c_post_processing/utilities/install_R_packages.R"
LOG_DIR="$SCRIPT_DIR/logs/installation"

#-------------------------------------------------------------------------------
# ENVIRONMENT SPEC — declared once, in setup_conda_geaManifest.yml
#-------------------------------------------------------------------------------
# Everything needed to build the environment used to be split across two files
# that each drifted from the other: a hardcoded name/python at the top of this
# script, and the package set in bash arrays below it, while the .yml carried a
# rival copy of the first two. All three now have exactly one home — the
# authored spec, which is a real conda environment file:
#
#   name: gea_test           -> ENV_NAME      (which environment to build)
#   dependencies: - python=  -> PYTHON_SPEC   (which interpreter to solve for)
#   dependencies: - <pkg>=   -> ALL_PACKAGES  (the direct, pinned package set)
#
# Nothing here writes to that file, so hand-authored grouping, per-package
# rationale comments and the commented-out `#name: gea` alternative simply stay
# put. The generated companion, setup_conda_gea.lock.yml, is the one this script
# overwrites (see EXPORT LOCKFILE).

read_env_name() {
    [[ -f "$ENV_SPEC" ]] || return 1
    local line
    line="$(grep -m1 -E '^name:[[:space:]]*[^[:space:]#]' "$ENV_SPEC" 2>/dev/null)" || return 1
    line="${line#name:}"                    # drop the key
    line="${line%%#*}"                      # drop any trailing comment
    line="${line//[[:space:]]/}"            # drop spaces/tabs/CR
    line="${line#\"}"; line="${line%\"}"    # drop optional quoting
    line="${line#\'}"; line="${line%\'}"
    [[ -n "$line" ]] || return 1
    printf '%s' "$line"
}

# Parse the `dependencies:` list into ALL_PACKAGES — the direct dependencies,
# and the only package source this script has.
#
# Deliberately plain bash: this script bootstraps a bare machine, so it cannot
# depend on python, yq or an activated environment to read its own input. The
# grammar it accepts is the subset the spec file actually uses, and it is strict
# about the rest rather than silently dropping packages:
#
#   * only entries under a column-0 `dependencies:` key are read, so the
#     `channels:` list cannot leak in as a package spec;
#   * `#` opens a comment anywhere, so a commented-out `#   - foo=1.0` line
#     stays commented out and per-package rationale comments are ignored;
#   * a nested sub-list (`- pip:`) is a hard error, not a warning — pip specs
#     silently vanishing from an environment build is not an acceptable
#     failure mode;
#   * trailing CR is stripped per line. The repo has a documented history of
#     CRLF contamination from Windows editors, and a spec carrying `\r` would
#     otherwise reach the solver as the package name `pigz=2.8<CR>`.
read_direct_specs() {
    ALL_PACKAGES=()
    _SPEC_PARSE_ERR=""
    [[ -f "$ENV_SPEC" ]] || { _SPEC_PARSE_ERR="file not found"; return 1; }

    local line spec in_deps=false lineno=0
    while IFS= read -r line || [[ -n "$line" ]]; do
        (( ++lineno ))
        line="${line%$'\r'}"

        # A key at column 0 opens or closes the dependencies block. Comment
        # lines and indented content start with '#' or whitespace and fall
        # through, leaving the current block state alone.
        if [[ "$line" =~ ^[^[:space:]#-] ]]; then
            if [[ "$line" =~ ^dependencies:[[:space:]]*(#.*)?$ ]]; then
                in_deps=true
            else
                in_deps=false
            fi
            continue
        fi
        [[ "$in_deps" == true ]] || continue

        [[ "$line" =~ ^[[:space:]]*-[[:space:]]+(.*)$ ]] || continue
        spec="${BASH_REMATCH[1]}"
        spec="${spec%%#*}"                                # drop trailing comment
        spec="${spec#"${spec%%[![:space:]]*}"}"           # ltrim
        spec="${spec%"${spec##*[![:space:]]}"}"           # rtrim
        spec="${spec#\"}"; spec="${spec%\"}"              # drop optional quoting
        spec="${spec#\'}"; spec="${spec%\'}"
        [[ -n "$spec" ]] || continue

        if [[ "$spec" == *: ]]; then
            # Cleared, not left half-filled: the caller's only failure test is
            # an empty ALL_PACKAGES, and a partial list would sail past it and
            # build an environment missing whatever followed this line.
            ALL_PACKAGES=()
            _SPEC_PARSE_ERR="line ${lineno}: nested sub-lists are unsupported ('- ${spec}')"
            return 1
        fi
        ALL_PACKAGES+=("$spec")
    done < "$ENV_SPEC"

    (( ${#ALL_PACKAGES[@]} > 0 )) || { _SPEC_PARSE_ERR="no entries under 'dependencies:'"; return 1; }
    return 0
}

# Picked out of the parsed list rather than re-read from the file, so the spec
# handed to the solver and the one reported by --help cannot disagree.
# The anchored alternation is what keeps python-dateutil, python_abi,
# python-gil and cpython out of the result.
read_python_spec() {
    local spec re='^python([=<>!~].*)?$'
    for spec in "${ALL_PACKAGES[@]}"; do
        if [[ "$spec" =~ $re ]]; then printf '%s' "$spec"; return 0; fi
    done
    return 1
}

# Resolved before usage()/argument parsing so --help can name the environment.
# A failure here is not fatal yet: --help still works, and the hard check runs
# once arguments are parsed.
ALL_PACKAGES=()
ENV_NAME="$(read_env_name || true)"
read_direct_specs || true
PYTHON_SPEC="$(read_python_spec || true)"

# Array form so every call site expands identically without word-splitting.
#
# --override-channels discards whatever ~/.condarc lists and uses ONLY the two
# channels named here. Without it `-c` merely *prepends* to the user's channel
# list, so `defaults` (repo.anaconda.com) stayed in every solve.
#
# This is not cosmetic — it decides whether the environment resolves at all.
# Measured on 2026-08-27, same specs, same machine, same package cache:
#   mamba create --override-channels -c conda-forge -c bioconda ...
#       -> solved in 10.3 s, 886 packages
#   mamba create                     -c conda-forge -c bioconda ...
#       -> "critical libmamba Could not solve for environment specs"
#
# `defaults` supplies nothing this environment uses (every package in the
# resolved set comes from conda-forge or bioconda) but it does publish rival
# records under names the R stack needs — pkgs/r carries its own
# r-futile.logger, r-getoptlong and friends — and the solver then wanders into
# their prerequisites (r-base 3.3.x, r 3.2.2, icu 54, zlib 1.2.8) and dies
# there. It also emits Anaconda commercial Terms-of-Service warnings on every
# run, and forces four flat-repodata downloads because pkgs/main and pkgs/r
# publish no shard index.
CHANNELS=(--override-channels -c conda-forge -c bioconda)

UPDATE_MODE=false
ENV_RESTART_MODE=false
DRY_RUN=false           # If true, only show what would be done without installing
SKIP_UPDATE_CHECK=false # If true, skip slow update availability check
RECLASSIFY_LOG=""       # If set, re-derive that transcript's issue log and exit

#===============================================================================
# USAGE
#===============================================================================

usage() {
    cat << EOF
Usage: bash setup_conda_gea.sh [OPTIONS]

Creates or updates the '${ENV_NAME:-<unset>}' conda environment for the GEA pipeline.

Environment name, python version and the whole direct package set come from
setup_conda_geaManifest.yml ('name:' and the 'dependencies:' list). Edit them there;
this script keeps no copy of any of them.
  spec     : $ENV_SPEC
  name     : ${ENV_NAME:-<unreadable — see setup_conda_geaManifest.yml>}
  python   : ${PYTHON_SPEC:-<unreadable — see setup_conda_geaManifest.yml>}
  packages : ${#ALL_PACKAGES[@]} direct (transitive deps resolved by the solver)
  lockfile : $LOCKFILE  (generated)

Options:
  --update              Update an existing environment
  --restart             Remove and recreate the environment
  --dry-run             Show what would be done without installing
  --skip-update-check   Skip the slow update availability check
  --reclassify=LOG      Re-derive LOG's issue log with the current classifier
                        and exit; installs nothing and does not touch LOG
  -h, --help            Show this help and exit

With no flags, the environment is created if missing, and packages are
(re)installed into it if it already exists.

Every run (except --help and --reclassify) writes a full transcript to:
  logs/installation/setup_conda_gea_<RUN_ID>_full_log.log
Warnings and errors are distilled from that transcript into:
  logs/installation/setup_conda_gea_<RUN_ID>_errors_warnings.log
split into ACTIONABLE and KNOWN-BENIGN (see the note at the top of this file).
EOF
}

#===============================================================================
# ARGUMENT PARSING
#===============================================================================

for arg in "$@"; do
    case "$arg" in
        --update)             UPDATE_MODE=true ;;
        --restart)            ENV_RESTART_MODE=true ;;
        --dry-run)            DRY_RUN=true ;;
        --skip-update-check)  SKIP_UPDATE_CHECK=true ;;
        --reclassify=*)       RECLASSIFY_LOG="${arg#*=}" ;;
        -h|--help)            usage; exit 0 ;;
        *)  echo "ERROR: Unknown argument: $arg" >&2; usage >&2; exit 1 ;;
    esac
done

# Deferred from the CONFIGURATION block so --help and a bad flag still behave.
#
# Handing an empty or half-parsed list to the solver is the one failure this
# script must never have: `conda create -n gea` with no specs succeeds, leaving
# a plausible-looking but empty environment behind. So a spec that does not
# parse is fatal here, before anything is created.
if [[ -z "$ENV_NAME" || -z "$PYTHON_SPEC" || ${#ALL_PACKAGES[@]} -eq 0 ]]; then
    echo "ERROR: setup_conda_geaManifest.yml is the single source of truth for the" >&2
    echo "       environment spec, and this script could not read it:" >&2
    echo "         $ENV_SPEC" >&2
    [[ -z "$ENV_NAME"    ]] && echo "       missing/unreadable: 'name: <env>'" >&2
    [[ ${#ALL_PACKAGES[@]} -eq 0 ]] && \
        echo "       missing/unreadable: 'dependencies:' list${_SPEC_PARSE_ERR:+ — $_SPEC_PARSE_ERR}" >&2
    [[ -z "$PYTHON_SPEC" ]] && echo "       missing/unreadable: 'dependencies:' entry '- python=<version>'" >&2
    echo "       Restore it with: git checkout -- setup_conda_geaManifest.yml" >&2
    exit 1
fi

#===============================================================================
# PACKAGE LISTS
#===============================================================================
# There are none here any more. The direct, pinned package set lives in the
# authored spec and is parsed into ALL_PACKAGES by read_direct_specs() above:
#
#     setup_conda_geaManifest.yml  ->  dependencies:
#
# Five bash arrays (PREPROCESSING_TOOLS, ALIGNMENT_TOOLS, R_BASE, BIOCONDUCTOR,
# CRAN_PACKAGES) used to sit here, which meant the one file a reader opens to
# see what the environment contains — the .yml — was the one file that did not
# declare it. Add or repin a package in the spec; nothing in this script needs
# to change.

#===============================================================================
# VERIFICATION LISTS
#===============================================================================

# Executables expected on PATH once the environment is active
VERIFY_CMDS=(
    R Rscript python3
    samtools salmon hisat2 STAR bowtie2 stringtie rsem-calculate-expression
    fastqc trim_galore trimmomatic
    prefetch fasterq-dump pigz
    gffread infer_experiment.py gtfToGenePred genePredToBed dos2unix
)

# R packages spot-checked after install.
# The authoritative, complete R package lists live in install_R_packages.R.
VERIFY_R_PKGS=(
    DESeq2 ComplexHeatmap WGCNA tximport clusterProfiler
    ggplot2 pheatmap igraph corrplot dendextend gridExtra scales
)

#===============================================================================
# UTILITY FUNCTIONS
#===============================================================================

# Message format mirrors modules_gea/logging/logging_utils.sh (_log_impl) so
# installation logs read the same as pipeline logs. WARN/ERROR go to stderr,
# which is merged into the transcript; ERROR_WARN_FILE is distilled from that
# transcript at teardown (see _write_issue_log) so tool output is captured too.
_log_impl() {
    local level="$1"; shift
    local _ts; printf -v _ts '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || _ts=$(date '+%Y-%m-%d %H:%M:%S')
    if [[ "$level" == "WARN" || "$level" == "ERROR" ]]; then
        echo "[$_ts] [$level] $*" >&2
    else
        echo "[$_ts] [$level] $*"
    fi
    return 0
}

log_info()  { _log_impl INFO  "$@"; }
log_warn()  { _log_impl WARN  "$@"; }
log_error() { _log_impl ERROR "$@"; }

# Banner used for the major section headings printed at runtime
log_section() {
    log_info "========================================"
    log_info "$*"
    log_info "========================================"
}

check_command() { command -v "$1" &> /dev/null; }

# Run a command, or print it in dry-run mode.
# NOTE: redirections are applied by the caller, so this cannot wrap a command
#       whose output is redirected into a file (see the lockfile export below).
run_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        echo "[DRY RUN] Would execute: $*"
    else
        "$@"
    fi
}

#===============================================================================
# LOGGING
#===============================================================================
# Deliberately self-contained instead of sourcing modules_gea/logging/logging_utils.sh:
# this script bootstraps the environment on a bare machine, while that system is
# built for pipeline runs — it would create six sibling metric directories
# (time_logs, space_logs, space_time_logs, software_catalogs, gpu_log) that carry
# no meaning for an env setup, and its rotate_old_logs() operates one level above
# LOG_DIR (i.e. all of logs/). Format and ANSI/CR stripping match it regardless.
#
# Set up after argument parsing so `--help` and bad flags do not litter the
# directory with empty logs.

printf -v RUN_ID '%(%Y%m%d_%H%M%S)T' -1 2>/dev/null || RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="$LOG_DIR/setup_conda_gea_${RUN_ID}_full_log.log"
ERROR_WARN_FILE="$LOG_DIR/setup_conda_gea_${RUN_ID}_errors_warnings.log"

# CR -> newline so conda/mamba progress-bar redraws become readable lines rather
# than one overwritten blob; ANSI colour/cursor codes are dropped from the file
# but left intact on the console.
_strip_ansi_stream() {
    sed -u $'s/\r/\\\n/g; s/\x1B\\[[0-9;?]*[a-zA-Z]//g; s/\x1B[()][A-Z0-9]//g'
}

# Lines worth pulling out of the transcript, in one ERE:
#   1. this script's own [WARN]/[ERROR] messages
#   2. libmamba's `error`/`warning` column output
#   3. Python-style exception headers — conda reports SafetyError, ClobberError,
#      CondaError, UnsatisfiableError and friends this way
#   4. conda's bare WARNING lines
#   5. the ✗ marks emitted by the executable / R-package verification blocks
_ISSUE_PATTERN='^\[[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9:]+\] \[(WARN|ERROR)\]'
_ISSUE_PATTERN+='|^(error|warning)[[:space:]]+libmamba'
_ISSUE_PATTERN+='|^[A-Za-z]*(Error|Exception):'
_ISSUE_PATTERN+='|^WARNING[[:space:]]'
_ISSUE_PATTERN+='|^[[:space:]]*✗'

# Of the lines _ISSUE_PATTERN catches, these are upstream packaging defects with
# no fix available from this repository. They are still written to the issue log
# — under a KNOWN-BENIGN heading, so nothing is hidden — but they do not count
# towards the actionable total, which is what makes that total worth reading.
#
# Every entry is deliberately narrow. A *different* package clobbering files, or
# a cache-integrity failure on some other path, must still read as actionable.
# Verified harmless against the built environment on 2026-08-27:
#   * perl-net-http is not even installed — the warning is about a malformed
#     record inside bioconda's sharded repodata, which the solver then ignores.
#   * ossuuid and libuuid ship the same unversioned lib/libuuid.{a,so} and
#     lib/pkgconfig/uuid.pc. Whichever links last wins that name, so it varies
#     between runs — but both runtime SONAMEs (libuuid.so.1 from util-linux,
#     libuuid.so.16 from ossuuid) are always present and correct, and sra-tools
#     turns out not to link uuid dynamically at all. prefetch/fasterq-dump 3.2.1
#     both run.
#   * kingfisher and argparse-manpage-birdtools both ship $PREFIX/README.md.
#     kingfisher 0.4.1 runs.
#   * r-base's lib/R/doc/html/packages.html — see the header block.
# NOT listed here on purpose: conda's other ClobberError phrasing, "This
# transaction has incompatible packages due to a shared path." That line carries
# no package name (the names sit on a continuation line _ISSUE_PATTERN does not
# capture), so it always reads as actionable — which is the right default for a
# clobber between arbitrary packages.
_BENIGN_PATTERN='libmamba Failed to parse field .* in shard package record'
_BENIGN_PATTERN+='|libmamba \[(ossuuid|kingfisher)-[^]]*\] The following files were already present'
_BENIGN_PATTERN+='|ClobberError.*(ossuuid|kingfisher|argparse-manpage-birdtools)'
_BENIGN_PATTERN+='|SafetyError: The package for r-base located at'
_BENIGN_PATTERN+='|libmamba Invalid package cache.*packages\.html'
_BENIGN_PATTERN+='|WARNING conda\.conda_pypi'

# Distil ERROR_WARN_FILE from the finished transcript instead of only from this
# script's own log_warn/log_error calls.
#
# conda and mamba write their failures straight to stdout/stderr, so the previous
# scheme — appending inside _log_impl — could only ever see messages the script
# itself produced. Run 20260827_113346 closed with "warnings: 1 | errors: 0" and
# a one-line issue log, while its transcript held
# "error libmamba Could not solve for environment specs", a SafetyError on the
# r-base package cache and four ClobberErrors. Every one of them was invisible in
# the file the usage text tells you to read.
#
# Must run *after* the log pipeline is reaped, or the tail of the transcript is
# not on disk yet and the trailing errors are the ones that go missing. Output is `grep -n`, so each entry doubles as a jump target into the
# full log, and the captured set is then split into ACTIONABLE / KNOWN-BENIGN.
# Echoes "<actionable> <benign>" for the caller.
_write_issue_log() {
    [[ -n "${LOG_FILE:-}" && -f "$LOG_FILE" && -n "${ERROR_WARN_FILE:-}" ]] || { echo "0 0"; return 0; }

    local _tmp="${ERROR_WARN_FILE}.tmp"
    grep -nE "$_ISSUE_PATTERN" "$LOG_FILE" > "$_tmp" 2>/dev/null || true

    if [[ ! -s "$_tmp" ]]; then
        rm -f "$_tmp" "$ERROR_WARN_FILE" 2>/dev/null || true
        echo "0 0"
        return 0
    fi

    # Partition what was already captured — one pass over the transcript, not two.
    # Every _BENIGN_PATTERN entry is unanchored, so the "NNN:" prefixes grep -n
    # added do not interfere.
    local _act _ben _nact=0 _nben=0
    _act=$(grep -vE "$_BENIGN_PATTERN" "$_tmp" 2>/dev/null) || true
    _ben=$(grep -E  "$_BENIGN_PATTERN" "$_tmp" 2>/dev/null) || true
    if [[ -n "$_act" ]]; then _nact=$(printf '%s\n' "$_act" | wc -l); fi
    if [[ -n "$_ben" ]]; then _nben=$(printf '%s\n' "$_ben" | wc -l); fi

    {
        echo "# setup_conda_gea.sh — run ${RUN_ID}"
        echo "# Distilled from: $LOG_FILE"
        echo "# Entries are <line-number>:<text> from that file."
        if [[ -n "${_ISSUE_LOG_NOTE:-}" ]]; then echo "# ${_ISSUE_LOG_NOTE}"; fi
        echo
        echo "#=== ACTIONABLE (${_nact}) ==================================================="
        if (( _nact > 0 )); then printf '%s\n' "$_act"; else echo "(none)"; fi
        echo
        echo "#=== KNOWN-BENIGN (${_nben}) — upstream packaging, no fix available here ===="
        echo "# See the 'Known benign log lines' block at the top of setup_conda_gea.sh."
        if (( _nben > 0 )); then printf '%s\n' "$_ben"; else echo "(none)"; fi
    } > "$ERROR_WARN_FILE"

    rm -f "$_tmp" 2>/dev/null || true
    echo "${_nact} ${_nben}"
    return 0
}

_teardown_logging() {
    local _rc=$?

    # Restore the original stdout/stderr so the log process substitution sees EOF,
    # then reap it. That subshell waits for its own ANSI stripper before exiting
    # (see the redirect below), so this one wait means the whole chain has drained:
    # every line is on disk before _write_issue_log greps the transcript, and
    # before the summary block appends to it. Without it the tail is lost at exit.
    exec 1>&3 2>&4
    if [[ -n "${_LOG_TEE_PID:-}" ]]; then
        wait "$_LOG_TEE_PID" 2>/dev/null || true
    fi

    [[ -n "${LOG_FILE:-}" ]] || return "$_rc"

    local _counts; _counts=$(_write_issue_log)
    local _nact="${_counts%% *}" _nben="${_counts##* }"
    local _ts; printf -v _ts '%(%Y-%m-%d %H:%M:%S)T' -1 2>/dev/null || _ts=$(date '+%Y-%m-%d %H:%M:%S')

    # fds are restored, so this reaches the console; tee -a puts it in the log too.
    {
        echo "[$_ts] [INFO] Exit status : ${_rc}  |  actionable: ${_nact}  |  known-benign: ${_nben}"
        echo "[$_ts] [INFO] Full log    : $LOG_FILE"
        if (( _nact > 0 )); then
            echo "[$_ts] [INFO] Issue log   : $ERROR_WARN_FILE"
        elif (( _nben > 0 )); then
            echo "[$_ts] [INFO] Issue log   : $ERROR_WARN_FILE  (known-benign only)"
        fi
    } | tee -a "$LOG_FILE" 2>/dev/null || true

    return "$_rc"
}

#===============================================================================
# RECLASSIFY MODE
#===============================================================================
# Re-derive an existing run's issue log with the current classifier, without
# installing anything.
#
# The issue log is a *derived* artifact — the `_full_log.log` transcript is the
# record, and it is never written here. So when the classifier improves, an old
# issue log can be regenerated from its own transcript rather than left in a
# stale format that contradicts what the script now produces. The regenerated
# file carries a provenance note saying exactly that.
#
#   bash setup_conda_gea.sh --reclassify=logs/installation/<run>_full_log.log
#
# Deliberately placed before the tee/trap setup below: this mode reads, reports
# and exits, so it should not open a transcript of its own.
if [[ -n "$RECLASSIFY_LOG" ]]; then
    if [[ ! -f "$RECLASSIFY_LOG" ]]; then
        echo "ERROR: no such transcript: $RECLASSIFY_LOG" >&2
        exit 1
    fi
    if [[ "$RECLASSIFY_LOG" != *_full_log.log ]]; then
        echo "ERROR: --reclassify expects a *_full_log.log transcript, got: $RECLASSIFY_LOG" >&2
        exit 1
    fi

    LOG_FILE="$RECLASSIFY_LOG"
    ERROR_WARN_FILE="${LOG_FILE%_full_log.log}_errors_warnings.log"
    RUN_ID="${LOG_FILE##*/setup_conda_gea_}"
    RUN_ID="${RUN_ID%_full_log.log}"
    _ISSUE_LOG_NOTE="Re-derived by 'setup_conda_gea.sh --reclassify'; the transcript itself is unmodified."

    _reclass_counts="$(_write_issue_log)"
    log_info "Reclassified : $LOG_FILE"
    log_info "Issue log    : $ERROR_WARN_FILE"
    log_info "actionable: ${_reclass_counts%% *}  |  known-benign: ${_reclass_counts##* }"
    exit 0
fi

if mkdir -p "$LOG_DIR" 2>/dev/null; then
    exec 3>&1 4>&2
    # The trailing `wait` is load-bearing. _strip_ansi_stream sits in an *inner*
    # process substitution, so its PID never reaches this shell — $! yields the
    # outer subshell only. The teardown therefore reaped tee while sed still held
    # the tail of the run in its buffers: on a 20k-line burst ~400 lines were
    # still in flight when _write_issue_log ran, which is precisely where a failed
    # run puts the errors that explain it, and the summary block then appended to
    # a file sed was concurrently writing. `wait` makes the subshell outlive its
    # own child, so waiting on $! below transitively means "chain fully drained".
    #
    # tee keeps the *inherited* stdout rather than reopening /dev/fd/3: fd 3 is a
    # dup, sharing one file offset with fd 1, and a reopen leaves that offset at 0
    # — the teardown's restored fd 1 would then overwrite the head of the console
    # transcript whenever this script is run with stdout redirected to a file.
    exec > >(tee >(_strip_ansi_stream >> "$LOG_FILE"); wait) 2>&1
    _LOG_TEE_PID=$!
    trap _teardown_logging EXIT
else
    echo "[WARN] Cannot create log directory: $LOG_DIR — continuing without file logging" >&2
    LOG_FILE=""
    ERROR_WARN_FILE=""
fi

log_section "setup_conda_gea.sh — run ${RUN_ID}"
log_info "Command   : bash setup_conda_gea.sh${*:+ $*}"
log_info "Host      : ${HOSTNAME:-unknown}"
log_info "User      : ${USER:-unknown}"
log_info "Directory : $SCRIPT_DIR"
if [[ -n "$LOG_FILE" ]]; then
    log_info "Log file  : $LOG_FILE"
fi

#===============================================================================
# DETECT PACKAGE MANAGER  (prefer mamba, fall back to conda)
#===============================================================================

if check_command mamba; then
    PKG_MGR="mamba"
    log_info "Using mamba for faster installation"
else
    PKG_MGR="conda"
    log_warn "Mamba not found, using conda. Install mamba for faster setup:"
    log_warn "  conda install -c conda-forge mamba"
fi

#===============================================================================
# CONDA INITIALIZATION
#===============================================================================

log_info "Initializing conda..."

if [[ -z "${CONDA_EXE:-}" ]]; then
    if   [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    elif [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
        source "/opt/conda/etc/profile.d/conda.sh"
    else
        log_error "Cannot find conda installation"
        exit 1
    fi
fi

if [[ "${WF_MANAGED_ENV:-}" != "true" ]]; then
    eval "$(conda shell.bash hook)"
fi

#===============================================================================
# HELPERS: ENVIRONMENT / PACKAGE INSPECTION
#===============================================================================

# True if the conda environment already exists.
# Matches on the name column rather than a line anchor: conda prints
# "gea  /path" while mamba 2.x prints an indented "  gea  Active  /path" table,
# so "^${ENV_NAME} " silently never matched under mamba.
# Kept as a single awk rather than `awk | grep -q`: grep -q exits on first match,
# which can SIGPIPE awk (exit 141); under `set -o pipefail` that made an existing
# environment look missing, and the resulting `create` then failed on "prefix
# already exists".
env_exists() {
    conda env list 2>/dev/null | awk -v want="${ENV_NAME}" '$1 == want { found = 1 } END { exit !found }'
}

# Returns a space-separated list of package names missing from the environment
check_packages_installed() {
    # Build associative array for O(1) lookup (avoids echo|grep fork per package)
    local -A _installed_set=()
    local _line
    while IFS= read -r _line; do
        [[ -n "$_line" ]] && _installed_set["${_line%%=*}"]=1
    done < <("${PKG_MGR}" list -n "${ENV_NAME}" --export 2>/dev/null)

    local missing=()
    for pkg in "${ALL_PACKAGES[@]}"; do
        local pkg_name="${pkg%%[><=]*}"   # strip version constraint for comparison
        if [[ -z "${_installed_set[$pkg_name]+x}" ]]; then
            missing+=("$pkg_name")
        fi
    done
    echo "${missing[*]:-}"
}

#===============================================================================
# SOLVER INVOCATION
#===============================================================================

# A create that dies *after* the solve leaves a partial prefix on disk, and the
# next create then fails with "prefix already exists" instead of retrying
# anything useful. Clear it first. No-op for `install` (the prefix is meant to
# exist) and in dry-run (nothing was created).
_reset_partial_env() {
    [[ "$1" == "create" && "$DRY_RUN" != true ]] || return 0
    if env_exists; then
        log_warn "Removing partially created environment '${ENV_NAME}' before retrying..."
        conda env remove -n "${ENV_NAME}" -y &> /dev/null || true
    fi
    return 0
}

# run_solver <create|install> <spec>...
#
# Single entry point for every environment-building call, so the retry policy is
# applied uniformly:
#   1. try the preferred manager (mamba when present);
#   2. on failure, drop the repodata cache and retry once. mamba 2.x serves
#      conda-forge and bioconda from *sharded* repodata; this is the generic
#      "the index on disk is stale or truncated" recovery. Cheap, and the only
#      thing a retry can usefully change;
#   3. on a second failure, fall back to conda, which reads flat repodata and
#      has repeatedly solved specs mamba rejected (that is what rescued the
#      20260827_113346 run).
# When PKG_MGR is already conda there is nothing to fall back to, so the first
# failure is fatal rather than a repeat of the identical solve.
run_solver() {
    local action="$1"; shift

    if run_cmd "${PKG_MGR}" "$action" -n "${ENV_NAME}" "${CHANNELS[@]}" -y "$@"; then
        return 0
    fi

    if [[ "${PKG_MGR}" != "conda" ]]; then
        log_warn "${PKG_MGR} ${action} failed — refreshing the repodata/shard cache and retrying..."
        run_cmd "${PKG_MGR}" clean --index-cache -y || true
        _reset_partial_env "$action"
        if run_cmd "${PKG_MGR}" "$action" -n "${ENV_NAME}" "${CHANNELS[@]}" -y "$@"; then
            return 0
        fi

        log_warn "${PKG_MGR} ${action} failed again, falling back to conda..."
        _reset_partial_env "$action"
        if run_cmd conda "$action" -n "${ENV_NAME}" "${CHANNELS[@]}" -y "$@"; then
            return 0
        fi
    fi

    log_error "Failed to ${action} environment '${ENV_NAME}' — see $LOG_FILE"
    return 1
}

#===============================================================================
# ENVIRONMENT CREATION / UPDATE
#===============================================================================

log_section "Setting up environment: $ENV_NAME"
log_info "Total packages: ${#ALL_PACKAGES[@]}"

# Python is solved TOGETHER with the tool list, never before it.
#
# The previous two-step shape — create the env with `python=3.11`, then install
# the 66 packages into it — makes the bare-Python solve commit to build variants
# the real dependency set cannot keep. It picks libsqlite-3.53.4-h13e7031_1,
# whose `depends` carries `icu >=78.3,<79.0a0`, so icu 78.3 lands in the prefix;
# the R 4.3 / Bioconductor 3.18 stack needs icu 75.1, and step two has to undo
# step one. Measured on 2026-08-27 (identical channels either way):
#   two-step  -> Change: libsqlite h13e7031_1 -> h0737f62_1
#                Downgrade: icu 78.3 -> 75.1
#   one solve -> icu 75.1 and libsqlite h0737f62_1 chosen outright
#
# It also removes the second-transaction file collisions. `ossuuid` (pulled in by
# sra-tools) and `libuuid` both ship lib/libuuid.{a,so} and pkgconfig/uuid.pc;
# when libuuid arrives in an earlier transaction conda reports three
# ClobberErrors against paths it already linked. Within one transaction it
# orders the links itself. The two libraries do coexist — their runtime SONAMEs
# are libuuid.so.1 and libuuid.so.16, and only the unversioned developer
# symlinks overlap — so this was never fatal, just noise that looked fatal.
#
# NOTE: this is not what fixed the "Could not solve for environment specs"
# failure; --override-channels above is. The two-step shape solves fine once
# `defaults` is out of the channel list. This change is about not doing work
# twice, and about a clean transaction log.
# ALL_PACKAGES already carries the python spec — it is just another entry in the
# spec file's dependency list — so it is NOT prepended here. Doing both would
# hand the solver the same spec twice.
SOLVE_SPECS=("${ALL_PACKAGES[@]}")

NEEDS_CREATE=true

if env_exists; then
    log_info "Environment '${ENV_NAME}' exists."

    if [[ "$ENV_RESTART_MODE" == true ]]; then
        log_info "Removing existing environment '${ENV_NAME}'..."
        run_cmd "${PKG_MGR}" env remove -n "${ENV_NAME}" -y
    else
        NEEDS_CREATE=false

        if [[ "$UPDATE_MODE" == true ]]; then
            MISSING_PKGS=$(check_packages_installed)
            if [[ -n "$MISSING_PKGS" ]]; then
                log_info "Missing packages detected: $MISSING_PKGS"
                log_info "Installing/updating all packages..."
            else
                if [[ "$SKIP_UPDATE_CHECK" == true ]]; then
                    log_info "All packages installed. Skipping update check (--skip-update-check)."
                    exit 0
                fi
                log_info "All packages installed. Running update..."
                run_cmd "${PKG_MGR}" update -n "${ENV_NAME}" "${CHANNELS[@]}" --all -y
                log_info "Update complete."
                exit 0
            fi
        fi
    fi
fi

#===============================================================================
# SOLVE AND INSTALL
#===============================================================================

if [[ "$NEEDS_CREATE" == true ]]; then
    log_info "Creating environment '${ENV_NAME}' (${#SOLVE_SPECS[@]} specs, single solve)..."
    run_solver create "${SOLVE_SPECS[@]}" || exit 1
else
    log_info "Installing packages into '${ENV_NAME}' (${#SOLVE_SPECS[@]} specs, single solve)..."
    run_solver install "${SOLVE_SPECS[@]}" || exit 1
fi

#===============================================================================
# CONFIGURE SRA TOOLS
#===============================================================================

log_info "Configuring SRA tools..."
run_cmd conda run -n "${ENV_NAME}" vdb-config --prefetch-to-cwd || \
    log_warn "vdb-config failed — SRA prefetch-to-cwd not set. Pipeline will still work but prefetch may use default cache location."

#===============================================================================
# ACTIVATE AND VERIFY
#===============================================================================

if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY-RUN] Would activate environment '${ENV_NAME}' and verify tooling"
else
    log_info "Activating environment..."
    if ! conda activate "${ENV_NAME}" 2>/dev/null; then
        # Not fatal: the environment is already built. Failing hard here would
        # skip the R-package installer and the lockfile export below.
        log_warn "Could not activate '${ENV_NAME}' in this shell — skipping verification."
        log_warn "Verify manually with: conda activate ${ENV_NAME}"
    else
        log_info "Verifying key executables..."
        for cmd in "${VERIFY_CMDS[@]}"; do
            if check_command "$cmd"; then
                log_info "  ✓ $cmd"
            else
                log_warn "  ✗ $cmd not found"
            fi
        done

        log_info "Checking R packages..."
        if check_command Rscript; then
            Rscript -e '
for (pkg in commandArgs(trailingOnly = TRUE)) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(paste0("  ✓ ", pkg, "\n"))
    } else {
        cat(paste0("  ✗ ", pkg, " NOT FOUND\n"))
    }
}
' "${VERIFY_R_PKGS[@]}"
        else
            log_warn "Rscript not available — skipping R package check"
        fi
    fi
fi

#===============================================================================
# INSTALL ADDITIONAL R PACKAGES
#===============================================================================
# install_R_packages.R holds the authoritative CRAN + Bioconductor lists and
# already installs only what is missing, so it is the only fallback needed here.

if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY-RUN] Would install missing R packages via $R_INSTALLER"
elif [[ ! -f "$R_INSTALLER" ]]; then
    log_warn "Not found: $R_INSTALLER — skipping R package fallback install"
elif ! check_command Rscript; then
    log_warn "Rscript not on PATH — skipping R package fallback install"
    log_warn "Run manually after activating: Rscript --no-save \"$R_INSTALLER\""
else
    log_info "Installing any missing R packages via install_R_packages.R..."
    Rscript --no-save "$R_INSTALLER" || log_warn "Some R packages may have failed"
fi

#===============================================================================
# EXPORT LOCKFILE
#===============================================================================
# Cannot use run_cmd here: the redirection happens in this shell, so a dry run
# would still truncate the lockfile.

if [[ "$DRY_RUN" == true ]]; then
    log_info "[DRY-RUN] Would export environment lockfile to: $LOCKFILE"
else
    log_info "Exporting environment lockfile..."
    # Sibling temp file + atomic move: a failed export must not destroy the
    # previous, known-good lockfile by truncating it to the header.
    LOCKFILE_TMP="${LOCKFILE}.tmp"

    # No name-block preservation any more. That machinery existed only because
    # hand-authored content (the `name:` line and its commented-out `#name: gea`
    # alternative) lived inside a file this step overwrites. It now lives in
    # setup_conda_geaManifest.yml, which nothing here writes to, so the lockfile can be
    # generated end to end.
    {
        echo "#==============================================================================="
        echo "# CONDA ENVIRONMENT LOCKFILE — '${ENV_NAME}'"
        echo "#==============================================================================="
        echo "# GENERATED FILE — do not edit by hand; setup_conda_gea.sh overwrites it after"
        echo "# every successful run."
        echo "#"
        echo "# This is the full resolved set: every package the solver pulled in, transitive"
        echo "# dependencies included, produced by 'conda env export --no-builds'. Two lines"
        echo "# are stripped:"
        echo "#   'prefix:'     — machine-specific, keeps the file portable"
        echo "#   '- defaults'  — see the channel note in setup_conda_gea.sh"
        echo "#"
        echo "# The DIRECT, hand-pinned dependencies — the packages this project actually"
        echo "# asks for — live in the authored spec beside this file:"
        echo "#"
        echo "#   setup_conda_geaManifest.yml"
        echo "#"
        echo "# Edit packages there. This file only records what that spec resolved to."
        echo "#"
        echo "# Regenerate            : bash setup_conda_gea.sh"
        echo "# Rebuild same versions : conda env create -f setup_conda_gea.lock.yml"
        echo "#                         (--no-builds pins name=version, not the build string)"
        echo "# Re-solve from spec    : conda env create -f setup_conda_geaManifest.yml"
        echo "#==============================================================================="
    } > "$LOCKFILE_TMP"

    # `conda env export` takes no --override-channels: it reports the channels
    # recorded against the prefix, so `defaults` survives into the export even
    # when nothing in the environment came from it. Left in, the documented
    # "conda env create -f setup_conda_gea.lock.yml" recreate path would put
    # repo.anaconda.com straight back into the channel list and reproduce the
    # unsolvable environment this script exists to avoid. Strip it here, the
    # same way the machine-specific `prefix:` line is stripped.
    if conda env export -n "$ENV_NAME" --no-builds \
        | grep -vE '^prefix:|^[[:space:]]*-[[:space:]]*defaults[[:space:]]*$' \
          >> "$LOCKFILE_TMP"; then
        mv -f "$LOCKFILE_TMP" "$LOCKFILE"
        log_info "Lockfile saved to: $LOCKFILE"
    else
        rm -f "$LOCKFILE_TMP"
        log_warn "Failed to export environment lockfile — kept existing: $LOCKFILE"
    fi
fi

#===============================================================================
# COMPLETION
#===============================================================================

log_section "Setup complete!"
log_info ""
log_info "Environment : $ENV_NAME"
log_info "Spec        : $ENV_SPEC  (edit packages here)"
log_info "Lockfile    : $LOCKFILE  (generated)"
log_info "Log dir     : $LOG_DIR"
log_info "Activate    : conda activate $ENV_NAME"
log_info "Deactivate  : conda deactivate"
log_info ""
log_info "Run post-processing with:"
log_info "  cd \"$SCRIPT_DIR\" && bash run_post_processing.sh"
log_info ""
log_info "HTML results viewer is auto-generated by run_post_processing.sh."
log_info "For standalone viewer setup, see: modules_gea/other_tools/setup_html_viewer_env.sh"
log_info "  bash modules_gea/other_tools/setup_html_viewer_env.sh --generate   # generate viewer only"
log_info "  bash modules_gea/other_tools/setup_html_viewer_env.sh --serve      # serve on localhost:8080"
log_info ""
