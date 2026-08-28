#!/bin/bash
# cd to script's own directory so **/*.csv glob works from any cwd
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1

# CSVs are grouped in single_project/ and combined_projects/; globstar makes the
# **/*.csv patterns below reach into them. nullglob keeps an empty group quiet.
shopt -s globstar nullglob

echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║           SRR CSV FILES VALIDATION REPORT - COMPREHENSIVE              ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""

# Create summary table data
declare -A file_data
declare -a file_list

for csv in **/*.csv; do
    [[ -f "$csv" ]] || continue
    file_list+=("$csv")
    
    # Count samples (data lines)
    samples=$(grep -v "^#" "$csv" | grep -v "^SRR_ID" | grep -v "^$" | wc -l)
    
    # Count commented entries
    comments=$(grep -c "^#" "$csv")
    
    # Check for batch column
    has_batch=$(head -1 "$csv" | grep -c "Batch")
    batch_status="NO"
    [ "$has_batch" -eq 1 ] && batch_status="YES"
    
    # Check for issues
    has_issues="NO"
    issues=""
    
    # Check for empty values
    empty_srr=$(awk -F',' 'NR>1 && !/^#/ {if ($1 == "" || $1 ~ /^[[:space:]]*$/) print}' "$csv" | wc -l)
    empty_organ=$(awk -F',' 'NR>1 && !/^#/ {if ($2 == "" || $2 ~ /^[[:space:]]*$/) print}' "$csv" | wc -l)
    
    if [ "$empty_srr" -gt 0 ] || [ "$empty_organ" -gt 0 ]; then
        has_issues="YES"
        issues="${issues}Empty values ($empty_srr SRR_ID, $empty_organ Organ) "
    fi
    
    # Store in associative array
    file_data["$csv"]="$samples|$comments|$batch_status|$has_issues|$issues"
    # Fields: samples|comments|batch_status|has_issues|issues (5 fields)
done

# Print detailed summary table
printf "%-70s | %8s | %8s | %8s | %10s\n" "File Name" "Samples" "Comments" "Batch" "Issues"
printf "%s\n" "$(printf '%.0s─' {1..149})"

for csv in "${file_list[@]}"; do
    IFS='|' read -r samples comments batch has_issues issues <<< "${file_data[$csv]}"
    issue_display=$([ "$has_issues" == "NO" ] && echo "None" || echo "$issues")
    printf "%-70s | %8d | %8d | %8s | %10s\n" "$csv" "$samples" "$comments" "$batch" "$issue_display"
done

echo ""
echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║                        GLOBAL ANALYSIS RESULTS                         ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""

# Count total samples
total_samples=0
for csv in "${file_list[@]}"; do
    IFS='|' read -r samples _ _ _ _ <<< "${file_data[$csv]}"
    ((total_samples += samples))
done

# Get unique SRR IDs across all files
unique_srrs=$(grep -h "^SRR[0-9]" **/*.csv | awk -F',' '{print $1}' | sort -u | wc -l)

echo "Total sample entries: $total_samples"
echo "Total unique SRR_IDs: $unique_srrs"
echo "Total CSV files: ${#file_list[@]}"
echo ""

echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║                     DUPLICATE SRR_ID ANALYSIS                          ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""

# Find all duplicates with context
all_srrs=$(grep -h "^SRR[0-9]" **/*.csv | awk -F',' '{print $1}')
global_dups=$(echo "$all_srrs" | sort | uniq -d)

dup_count=$(echo "$global_dups" | grep -c "SRR")

if [ "$dup_count" -gt 0 ]; then
    echo "✗ FOUND $dup_count DUPLICATE SRR_IDs ACROSS FILES"
    echo ""
    echo "Duplicates by originating file:"
    echo "$global_dups" | while read srr; do
        locations=$(grep -l "^$srr," **/*.csv | tr '\n' ',' | sed 's/,$//')
        count=$(grep -l "^$srr," **/*.csv | wc -l)
        printf "  • %-20s (appears in %d files: %s)\n" "$srr" "$count" "$locations"
    done
else
    echo "✓ NO DUPLICATE SRR_IDs FOUND ACROSS FILES (excluding headers)"
fi

echo ""
echo "NOTES ON CROSS-FILE DUPLICATES:"
echo "  • single_project/ holds one CSV per accession (incl. _selected / _test subsets)"
echo "  • combined_projects/ holds CSVs spanning more than one accession"
echo "  • PRJNA328564_selected.csv & PRJNA328564_test.csv are subsets of PRJNA328564.csv"
echo "  • PRJNA865018_and_PRJNA941250.csv is the union of PRJNA865018.csv + PRJNA941250.csv"
echo "  • PRJNA328564_selected_plus_PRJNA1417697_selected.csv is the union of those two"
echo "  • These overlaps are INTENTIONAL for subset/combination workflows"
echo ""

echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║                    TRANSCRIPTOMICS COMPLIANCE RESULTS                  ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""

# Dynamically check compliance instead of hardcoding results
_col_pass=0 _col_fail=0 _fmt_pass=0 _fmt_fail=0 _empty_pass=0 _empty_fail=0
_dup_pass=0 _dup_fail=0 _batch_count=0 _comment_files=0

for csv in "${file_list[@]}"; do
    header=$(head -1 "$csv")
    # Required columns
    if echo "$header" | grep -q "SRR_ID" && echo "$header" | grep -q "Organ"; then
        ((_col_pass++))
    else
        ((_col_fail++))
    fi
    # SRR format
    bad_fmt=$(grep -v "^#" "$csv" | grep -v "^SRR_ID" | grep -v "^$" | awk -F',' '$1 !~ /^SRR[0-9]+$/ {print}' | wc -l)
    [ "$bad_fmt" -eq 0 ] && ((_fmt_pass++)) || ((_fmt_fail++))
    # Empty values
    IFS='|' read -r _ _ _ has_issues _ <<< "${file_data[$csv]}"
    [ "$has_issues" == "NO" ] && ((_empty_pass++)) || ((_empty_fail++))
    # Within-file duplicates
    within_dups=$(grep -v "^#" "$csv" | grep -v "^SRR_ID" | grep -v "^$" | awk -F',' '{print $1}' | sort | uniq -d | wc -l)
    [ "$within_dups" -eq 0 ] && ((_dup_pass++)) || ((_dup_fail++))
    # Batch column
    echo "$header" | grep -q "Batch" && ((_batch_count++))
    # Comment lines
    [ "$(grep -c "^#" "$csv")" -gt 0 ] && ((_comment_files++))
done

_total=${#file_list[@]}
[ "$_col_fail" -eq 0 ]   && echo "✓ Required Columns (SRR_ID, Organ): ALL FILES PASS"   || echo "✗ Required Columns (SRR_ID, Organ): $_col_fail/$_total FAILED"
[ "$_fmt_fail" -eq 0 ]   && echo "✓ SRR Accession Format (SRR[0-9]+): ALL FILES PASS"   || echo "✗ SRR Accession Format (SRR[0-9]+): $_fmt_fail/$_total FAILED"
[ "$_empty_fail" -eq 0 ] && echo "✓ No Empty/Missing Values: ALL FILES PASS"             || echo "✗ No Empty/Missing Values: $_empty_fail/$_total FAILED"
[ "$_dup_fail" -eq 0 ]   && echo "✓ No Duplicates Within Files: ALL FILES PASS"          || echo "✗ No Duplicates Within Files: $_dup_fail/$_total FAILED"
[ "$_batch_count" -gt 0 ] && echo "✓ Batch Metadata Column: FOUND in $_batch_count/$_total files" || echo "✗ Batch Metadata Column: NOT FOUND (optional)"
echo ""

_pass_count=$(( (_col_fail == 0) + (_fmt_fail == 0) + (_empty_fail == 0) + (_dup_fail == 0) ))
echo "Compliance Summary:"
echo "  • Format compliance: $_pass_count/4 checks passed"
echo "  • Data integrity: $(( _total - _empty_fail ))/$_total files clean"
[ "$_comment_files" -gt 0 ] && echo "  • Commented entries: Present in $_comment_files file(s)"
echo ""
