#!/bin/bash

echo "=== COMMENTED ENTRIES AND NOTES ==="
echo ""

for csv in *.csv; do
    commented=$(grep "^#" "$csv")
    if [ -n "$commented" ]; then
        echo "File: $csv"
        echo "$commented" | sed 's/^/  /'
        echo ""
    fi
done
