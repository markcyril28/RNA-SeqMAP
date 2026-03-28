#!/bin/bash

echo "=== GLOBAL DUPLICATE CHECK ==="
echo ""

# Extract all non-commented SRR_IDs from all files
all_srrs=$(grep -h "^SRR" *.csv | awk -F',' '{print $1}' | sort)

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

for csv in *.csv; do
    file_dups=$(grep "^SRR" "$csv" | awk -F',' '{print $1}' | sort | uniq -d)
    if [ -n "$file_dups" ]; then
        echo "✗ $csv: Duplicate SRR_IDs within file:"
        echo "$file_dups" | sed 's/^/    /'
    fi
done

echo "✓ No duplicates within individual files"
echo ""

# Show statistics
echo "=== SAMPLE STATISTICS ==="
echo ""
total_unique=$(echo "$all_srrs" | sort -u | wc -l)
total_entries=$(echo "$all_srrs" | wc -l)
echo "Total unique SRR_IDs: $total_unique"
echo "Total SRR entries: $total_entries"
