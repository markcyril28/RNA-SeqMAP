#!/bin/bash
# cd to script's own directory so *.csv glob works from any cwd
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1

echo "=== COMMENTED ENTRIES AND NOTES ==="
echo ""

for csv in *.csv; do
    [[ -f "$csv" ]] || continue
    commented=$(grep "^#" "$csv")
    if [ -n "$commented" ]; then
        echo "File: $csv"
        echo "$commented" | sed 's/^/  /'
        echo ""
    fi
done
