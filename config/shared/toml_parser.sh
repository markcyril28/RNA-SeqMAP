#!/bin/bash
# ==============================================================================
# TOML PARSER FOR BASH
# ==============================================================================
# Lightweight pure-bash TOML parser for pipeline configuration files.
# Supports: scalars, arrays, inline comments, and [section] headers.
#
# Usage:
#   source "config/shared/toml_parser.sh"
#   load_toml "config/2_alignment/HPC_full_ref_guided.toml"
#
# After loading, variables are set as bash variables/arrays:
#   [runtime]
#   threads = 64          →  THREADS=64
#   jobs = 2              →  JOBS=2
#
#   [pipeline_stages]
#   enabled = ["METHOD_1_HISAT2_REF_GUIDED", "METHOD_3_STAR_ALIGNMENT"]
#                         →  PIPELINE_STAGES=("METHOD_1_HISAT2_REF_GUIDED" "METHOD_3_STAR_ALIGNMENT")
#
# Section prefixes are NOT added to variable names — variables are exported
# as flat uppercase names matching the existing pipeline conventions.
#
# Special handling:
#   - genome_ref_pairs: pipe-delimited multi-field entries stored as array
#   - Boolean true/false → "TRUE"/"FALSE" strings (pipeline convention)
#   - Inline comments (#) stripped
#   - Multiline arrays supported (one element per line)
# ==============================================================================

[[ "${_TOML_PARSER_SOURCED:-}" == "true" ]] && return 0
export _TOML_PARSER_SOURCED="true"

# load_toml <file>
# Parses a TOML file and sets bash variables in the caller's scope.
load_toml() {
    local toml_file="$1"

    if [[ ! -f "$toml_file" ]]; then
        echo "[ERROR] TOML config not found: $toml_file" >&2
        return 1
    fi

    local line key value section=""
    local in_array=false
    local array_var=""
    local -a array_values=()

    while IFS= read -r line || [[ -n "$line" ]]; do
        # Strip leading/trailing whitespace
        line="${line#"${line%%[![:space:]]*}"}"
        line="${line%"${line##*[![:space:]]}"}"

        # Skip empty lines and pure comment lines
        [[ -z "$line" || "$line" == \#* ]] && continue

        # Close multiline array if we hit the closing bracket
        if [[ "$in_array" == true ]]; then
            # Strip inline comment: everything after an unquoted #
            local _data_part="$line"
            # If the line contains a quoted string, only strip # after the last quote
            if [[ "$_data_part" == *'"'* ]]; then
                local _after_last_quote="${_data_part##*\"}"
                _after_last_quote="${_after_last_quote%%\#*}"
                _data_part="${_data_part%\"*}\"${_after_last_quote}"
            else
                _data_part="${_data_part%%\#*}"
            fi
            # Trim trailing whitespace from data portion
            _data_part="${_data_part%"${_data_part##*[![:space:]]}"}"

            if [[ "$_data_part" == *"]"* ]]; then
                # Extract any remaining values before the closing bracket
                local before_bracket="${_data_part%%]*}"
                before_bracket="${before_bracket#"${before_bracket%%[![:space:]]*}"}"
                if [[ -n "$before_bracket" && "$before_bracket" != "]" ]]; then
                    _toml_parse_array_elements "$before_bracket" array_values
                fi
                # Set the array variable
                eval "${array_var}=(\"\${array_values[@]}\")"
                in_array=false
                array_var=""
                array_values=()
                continue
            fi
            # Accumulate array elements (strip trailing comma and inline comments)
            _toml_parse_array_elements "$line" array_values
            continue
        fi

        # Section header: [section_name]
        if [[ "$line" =~ ^\[([a-zA-Z0-9_.:-]+)\]$ ]]; then
            section="${BASH_REMATCH[1]}"
            continue
        fi

        # Key-value pair: key = value
        if [[ "$line" =~ ^([a-zA-Z_][a-zA-Z0-9_]*)\ *=\ *(.*) ]]; then
            key="${BASH_REMATCH[1]}"
            value="${BASH_REMATCH[2]}"

            # Convert key to UPPERCASE (pipeline convention)
            key="${key^^}"

            # Check if value starts an array
            if [[ "$value" == "["* ]]; then
                # Single-line array: ["val1", "val2"]
                if [[ "$value" == *"]"* ]]; then
                    local inner="${value#\[}"
                    inner="${inner%\]}"
                    # Strip inline comment after closing bracket
                    inner="${inner%%\#*}"
                    local -a single_array=()
                    _toml_parse_array_elements "$inner" single_array
                    eval "${key}=(\"\${single_array[@]}\")"
                else
                    # Multiline array — start accumulating
                    in_array=true
                    array_var="$key"
                    array_values=()
                    # Parse any elements on the opening line after [
                    local after_bracket="${value#\[}"
                    after_bracket="${after_bracket#"${after_bracket%%[![:space:]]*}"}"
                    [[ -n "$after_bracket" ]] && _toml_parse_array_elements "$after_bracket" array_values
                fi
            else
                # Scalar value
                value="$(_toml_parse_scalar "$value")"
                eval "${key}=\"\${value}\""
            fi
        fi
    done < "$toml_file"
}

# _toml_parse_scalar <raw_value>
# Strips quotes, inline comments, and converts booleans.
_toml_parse_scalar() {
    local val="$1"

    # Strip inline comment (only if # is outside quotes)
    if [[ "$val" == '"'* ]]; then
        # Quoted string — extract content between quotes
        val="${val#\"}"
        val="${val%%\"*}"
    elif [[ "$val" == "'"* ]]; then
        # Single-quoted string (literal)
        val="${val#\'}"
        val="${val%%\'*}"
    else
        # Unquoted — strip inline comment
        val="${val%%\#*}"
        # Trim trailing whitespace
        val="${val%"${val##*[![:space:]]}"}"

        # Convert TOML booleans to pipeline convention
        case "$val" in
            true)  val="TRUE" ;;
            false) val="FALSE" ;;
        esac
    fi

    printf '%s' "$val"
}

# _toml_parse_array_elements <line> <array_name_ref>
# Parses comma-separated values from an array line and appends to the named array.
# Handles commas inside quoted strings correctly.
# O(n) via single awk pass — replaces O(n) bash character-by-character loop which
# was ~50-100x slower due to per-character substring extraction and string append.
_toml_parse_array_elements() {
    local line="$1"
    local -n _arr_ref=$2

    # Fast path: empty or whitespace-only line
    local _trimmed="${line#"${line%%[![:space:]]*}"}"
    [[ -z "$_trimmed" ]] && return 0

    # Use awk to split comma-separated values, respecting quoted strings and # comments.
    # Each output line is one element (already stripped of quotes and trimmed).
    # This replaces the O(n) bash character loop which spawned no subprocesses but
    # performed expensive ${line:$i:1} substring extraction per character.
    local _elem
    while IFS= read -r _elem; do
        [[ -n "$_elem" ]] && _arr_ref+=("$_elem")
    done < <(awk '
    BEGIN { FS="" }
    {
        in_q = 0; cur = ""
        for (i = 1; i <= NF; i++) {
            c = $i
            if (in_q) {
                if (c == "\"") in_q = 0
                else cur = cur c
            } else if (c == "\"") {
                in_q = 1
            } else if (c == "#") {
                break
            } else if (c == ",") {
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", cur)
                if (cur != "") print cur
                cur = ""
            } else {
                cur = cur c
            }
        }
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", cur)
        if (cur != "") print cur
    }' <<< "$line")
}

# load_toml_srr_datasets <file>
# Special loader for srr_datasets.toml which has [dataset.NAME] sections
# that need to be loaded as SRR_LIST_<NAME> arrays and a combined list.
load_toml_srr_datasets() {
    local toml_file="$1"

    if [[ ! -f "$toml_file" ]]; then
        echo "[ERROR] SRR datasets TOML not found: $toml_file" >&2
        return 1
    fi

    local line section="" current_list_var=""
    local in_array=false
    local -a array_values=()

    while IFS= read -r line || [[ -n "$line" ]]; do
        # Strip leading/trailing whitespace
        line="${line#"${line%%[![:space:]]*}"}"
        line="${line%"${line##*[![:space:]]}"}"

        # Skip empty lines and pure comment lines
        [[ -z "$line" || "$line" == \#* ]] && continue

        # Close multiline array
        if [[ "$in_array" == true ]]; then
            # Strip inline comment: everything after an unquoted #
            local _data_part="$line"
            if [[ "$_data_part" == *'"'* ]]; then
                local _after_last_quote="${_data_part##*\"}"
                _after_last_quote="${_after_last_quote%%\#*}"
                _data_part="${_data_part%\"*}\"${_after_last_quote}"
            else
                _data_part="${_data_part%%\#*}"
            fi
            _data_part="${_data_part%"${_data_part##*[![:space:]]}"}"

            if [[ "$_data_part" == *"]"* ]]; then
                local before_bracket="${_data_part%%]*}"
                before_bracket="${before_bracket#"${before_bracket%%[![:space:]]*}"}"
                [[ -n "$before_bracket" && "$before_bracket" != "]" ]] && \
                    _toml_parse_array_elements "$before_bracket" array_values
                eval "${current_list_var}=(\"\${array_values[@]}\")"
                in_array=false
                current_list_var=""
                array_values=()
                continue
            fi
            _toml_parse_array_elements "$line" array_values
            continue
        fi

        # Section header
        if [[ "$line" =~ ^\[([a-zA-Z0-9_.:-]+)\]$ ]]; then
            section="${BASH_REMATCH[1]}"
            continue
        fi

        # Key-value pair
        if [[ "$line" =~ ^([a-zA-Z_][a-zA-Z0-9_]*)\ *=\ *(.*) ]]; then
            local key="${BASH_REMATCH[1]}"
            local value="${BASH_REMATCH[2]}"

            # For dataset sections like [dataset.PRJNA328564], the key "samples"
            # maps to SRR_LIST_PRJNA328564
            if [[ "$section" == dataset.* ]]; then
                local dataset_name="${section#dataset.}"
                dataset_name="${dataset_name^^}"
                # Sanitize: replace chars invalid in bash variable names with _
                dataset_name="${dataset_name//[^A-Z0-9_]/_}"
                if [[ "$key" == "samples" ]]; then
                    current_list_var="SRR_LIST_${dataset_name}"
                    if [[ "$value" == "["* ]]; then
                        if [[ "$value" == *"]"* ]]; then
                            local inner="${value#\[}"
                            inner="${inner%\]}"
                            local -a single_array=()
                            _toml_parse_array_elements "$inner" single_array
                            eval "${current_list_var}=(\"\${single_array[@]}\")"
                        else
                            in_array=true
                            array_values=()
                            local after_bracket="${value#\[}"
                            after_bracket="${after_bracket#"${after_bracket%%[![:space:]]*}"}"
                            [[ -n "$after_bracket" ]] && _toml_parse_array_elements "$after_bracket" array_values
                        fi
                    fi
                fi
            elif [[ "$section" == "combined" && "$key" == "active_datasets" ]]; then
                # Parse active_datasets (may be single-line or multiline array)
                if [[ "$value" == "["* ]]; then
                    if [[ "$value" == *"]"* ]]; then
                        # Single-line array
                        local inner="${value#\[}"
                        inner="${inner%\]}"
                        inner="${inner%%\#*}"
                        local -a _active_datasets=()
                        _toml_parse_array_elements "$inner" _active_datasets
                    else
                        # Multiline array — use the existing in_array mechanism
                        # but store into a special variable
                        in_array=true
                        current_list_var="__ACTIVE_DATASETS_TMP"
                        array_values=()
                        local after_bracket="${value#\[}"
                        after_bracket="${after_bracket#"${after_bracket%%[![:space:]]*}"}"
                        [[ -n "$after_bracket" ]] && _toml_parse_array_elements "$after_bracket" array_values
                        continue
                    fi
                    # Build SRR_COMBINED_LIST from active dataset arrays
                    SRR_COMBINED_LIST=()
                    for _ds in "${_active_datasets[@]}"; do
                        local _ds_var="SRR_LIST_${_ds^^}"
                        if [[ -n "${!_ds_var+x}" ]]; then
                            local -n _ds_ref="$_ds_var"
                            SRR_COMBINED_LIST+=("${_ds_ref[@]}")
                        fi
                    done
                fi
            fi
        fi
    done < "$toml_file"

    # Post-process: if active_datasets was a multiline array, build SRR_COMBINED_LIST now
    if [[ -n "${__ACTIVE_DATASETS_TMP+x}" && ${#__ACTIVE_DATASETS_TMP[@]} -gt 0 ]]; then
        SRR_COMBINED_LIST=()
        local _ds _ds_var
        for _ds in "${__ACTIVE_DATASETS_TMP[@]}"; do
            _ds_var="SRR_LIST_${_ds^^}"
            if [[ -n "${!_ds_var+x}" ]]; then
                local -n _ds_ref="$_ds_var"
                SRR_COMBINED_LIST+=("${_ds_ref[@]}")
            fi
        done
        unset __ACTIVE_DATASETS_TMP
    fi
}
