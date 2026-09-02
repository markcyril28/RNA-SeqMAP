#!/bin/bash
# cd to script's own directory so **/*.csv glob works from any cwd
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1

# CSVs are grouped in single_project/ and combined_projects/; globstar makes the
# **/*.csv patterns below reach into them. nullglob keeps an empty group quiet.
shopt -s globstar nullglob

# Function to validate a single CSV file
validate_csv() {
    local file="$1"
    local filename="${file#./}"
    
    # Extract header and check columns
    local header=$(head -1 "$file")
    local has_srr_id=0
    local has_organ=0
    local has_batch=0
    
    local has_class=0
    [[ "$header" =~ "SRR_ID" ]] && has_srr_id=1
    [[ "$header" =~ "Organ" ]] && has_organ=1
    [[ "$header" =~ "Batch" ]] && has_batch=1
    [[ "$header" =~ "Tissue_Class" ]] && has_class=1
    
    echo "=== FILE: $filename ==="
    echo "Headers: $header"
    echo "  $([ "$has_srr_id" -eq 1 ] && echo '✓' || echo '✗') SRR_ID column: $([ "$has_srr_id" -eq 1 ] && echo 'YES' || echo 'NO')"
    echo "  $([ "$has_organ" -eq 1 ] && echo '✓' || echo '✗') Organ column: $([ "$has_organ" -eq 1 ] && echo 'YES' || echo 'NO')"
    echo "  $([ "$has_batch" -eq 1 ] && echo '✓' || echo '✗') Batch column: $([ "$has_batch" -eq 1 ] && echo 'YES' || echo 'NO')"
    echo "  $([ "$has_class" -eq 1 ] && echo '✓' || echo '✗') Tissue_Class column: $([ "$has_class" -eq 1 ] && echo 'YES' || echo 'NO')"
    if [[ "$has_class" -eq 1 ]]; then
        local bad_class
        bad_class=$(awk -F',' 'NR>1 && $1!="" && $3!="Meristematic" && $3!="Non_meristematic" {print $1": "$3}' "$file")
        if [[ -n "$bad_class" ]]; then
            echo "  ✗ Invalid Tissue_Class values (expect Meristematic / Non_meristematic):"
            echo "$bad_class" | sed 's/^/      /'
        else
            local n_mer n_non
            n_mer=$(awk -F',' 'NR>1 && $3=="Meristematic"' "$file" | wc -l)
            n_non=$(awk -F',' 'NR>1 && $3=="Non_meristematic"' "$file" | wc -l)
            echo "  ✓ Tissue_Class valid: $n_mer meristematic / $n_non non-meristematic"
        fi
    fi
    
    # Count data lines (excluding header and comment lines)
    local data_lines=$(grep -v "^#" "$file" | grep -v "^SRR_ID" | grep -v "^$" | wc -l)
    echo "  Sample count: $data_lines"
    
    # Check for commented entries
    local commented=$(grep -c "^#" "$file")
    echo "  Commented entries: $commented"
    
    # Check for empty values in SRR_ID or Organ columns
    local empty_srr=$(awk -F',' 'NR>1 && !/^#/ {if ($1 == "" || $1 ~ /^[[:space:]]*$/) print NR": "$0}' "$file")
    local empty_organ=$(awk -F',' 'NR>1 && !/^#/ {if ($2 == "" || $2 ~ /^[[:space:]]*$/) print NR": "$0}' "$file")
    
    if [ -n "$empty_srr" ]; then
        echo "  ✗ Empty SRR_ID values:"
        echo "$empty_srr" | sed 's/^/    /'
    else
        echo "  ✓ No empty SRR_ID values"
    fi
    
    if [ -n "$empty_organ" ]; then
        echo "  ✗ Empty Organ values:"
        echo "$empty_organ" | sed 's/^/    /'
    else
        echo "  ✓ No empty Organ values"
    fi
    
    # Validate SRR format
    local invalid_srr=$(awk -F',' 'NR>1 && !/^#/ && $1 != "" {if ($1 !~ /^SRR[0-9]+$/ && $1 !~ /^[[:space:]]*$/) print NR": "$1}' "$file")
    if [ -n "$invalid_srr" ]; then
        echo "  ✗ Invalid SRR accession format:"
        echo "$invalid_srr" | sed 's/^/    /'
    else
        echo "  ✓ All SRR accessions match SRR[0-9]+ format"
    fi
    
    echo ""
}

# Process all CSV files
for csv in **/*.csv; do
    [[ -f "$csv" ]] || continue
    validate_csv "$csv"
done
