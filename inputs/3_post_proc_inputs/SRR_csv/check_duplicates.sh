#!/bin/bash
# cd to script's own directory so *.csv glob works from any cwd
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1

echo "=== GLOBAL DUPLICATE CHECK ==="
echo ""

# Extract all non-commented SRR_IDs from all files (^SRR[0-9] excludes header line)
all_srrs=$(grep -h "^SRR[0-9]" *.csv | awk -F',' '{print $1}' | sort)

# Find global duplicates
global_dups=$(echo "$all_srrs" | sort | uniq -d)

if [ -n "$global_dups" ]; then
    echo "✗ DUPLICATE SRR_IDs ACROSS ALL FILES:"
    echo "$global_dups" | while read srr; do
        echo "  $srr found in:"
        grep -l "^$srr," *.csv | sed 's/^/    /'
    done
else
    echo "✓ No duplicate SRR_IDs across all files"
fi

echo ""
echo "=== PER-FILE DUPLICATE CHECK ==="
echo ""

has_within_dups=false
for csv in *.csv; do
    [[ -f "$csv" ]] || continue
    file_dups=$(grep "^SRR[0-9]" "$csv" | awk -F',' '{print $1}' | sort | uniq -d)
    if [ -n "$file_dups" ]; then
        has_within_dups=true
        echo "✗ $csv: Duplicate SRR_IDs within file:"
        echo "$file_dups" | sed 's/^/    /'
    fi
done

if [ "$has_within_dups" = false ]; then
    echo "✓ No duplicates within individual files"
fi
echo ""

# Show statistics
echo "=== SAMPLE STATISTICS ==="
echo ""
total_unique=$(echo "$all_srrs" | sort -u | grep -c .)
total_entries=$(echo "$all_srrs" | grep -c .)
echo "Total unique SRR_IDs: $total_unique"
echo "Total SRR entries: $total_entries"
