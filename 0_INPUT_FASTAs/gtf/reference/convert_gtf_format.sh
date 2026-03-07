#!/bin/bash
# ==============================================================================
# GTF FORMAT CONVERTER
# ==============================================================================
# Converts Eggplant_V4.1_function_IPR_final_formatted_v3.gtf format to
# All_Smel_Genes_Full_Name_reformatted.gtf format
#
# Key transformations:
# 1. Fix attribute quoting: ""value"" -> "value"
# 2. Remove outer quotes from attribute column
# 3. Optionally map transcript IDs to chromosome names using a mapping file
#
# Usage: ./convert_gtf_format.sh
# ==============================================================================

set -euo pipefail

# ==============================================================================
# CONFIGURATION - Edit these variables as needed
# ==============================================================================

# Input GTF file (v3 format with escaped quotes)
INPUT_GTF="Eggplant_V4.1_function_IPR_final_formatted_v3.gtf"

# Output GTF file (reformatted)
OUTPUT_GTF="Eggplant_V4.1_function_IPR_reformatted.gtf"

# Optional: Mapping file for transcript->chromosome conversion (leave empty to skip)
# Format: transcript_id<TAB>chromosome<TAB>offset
MAPPING_FILE=""

# ==============================================================================
# SCRIPT LOGIC - No need to edit below this line
# ==============================================================================

if [[ ! -f "$INPUT_GTF" ]]; then
    echo "ERROR: Input file not found: $INPUT_GTF"
    exit 1
fi

echo "Converting GTF format..."
echo "  Input:  $INPUT_GTF"
echo "  Output: $OUTPUT_GTF"

# Process the GTF file
if [[ -n "$MAPPING_FILE" && -f "$MAPPING_FILE" ]]; then
    echo "  Mapping: $MAPPING_FILE"
    
    # With chromosome mapping
    awk -F'\t' -v OFS='\t' '
    BEGIN {
        # Load mapping file
        while ((getline line < "'"$MAPPING_FILE"'") > 0) {
            split(line, arr, "\t")
            chr_map[arr[1]] = arr[2]
            offset_map[arr[1]] = arr[3]
        }
    }
    {
        # Get transcript ID from column 1
        tid = $1
        
        # Map to chromosome if available
        if (tid in chr_map) {
            $1 = chr_map[tid]
            # Adjust coordinates if offset provided
            if (tid in offset_map && offset_map[tid] != "") {
                offset = offset_map[tid] - 1
                $4 = $4 + offset
                $5 = $5 + offset
            }
        }
        
        # Fix attribute column (column 9)
        # Remove leading/trailing quotes and fix escaped quotes
        attr = $9
        gsub(/^"/, "", attr)           # Remove leading quote
        gsub(/"$/, "", attr)           # Remove trailing quote
        gsub(/""/, "\"", attr)         # Fix doubled quotes
        $9 = attr
        
        print
    }' "$INPUT_GTF" > "$OUTPUT_GTF"
else
    echo "  No mapping file - keeping original sequence names"
    
    # Without chromosome mapping - just fix attribute formatting
    awk -F'\t' -v OFS='\t' '
    {
        # Fix attribute column (column 9)
        # Remove leading/trailing quotes and fix escaped quotes
        attr = $9
        gsub(/^"/, "", attr)           # Remove leading quote
        gsub(/"$/, "", attr)           # Remove trailing quote  
        gsub(/""/, "\"", attr)         # Fix doubled quotes
        gsub(/";$/, ";", attr)         # Clean up trailing
        $9 = attr
        
        print
    }' "$INPUT_GTF" > "$OUTPUT_GTF"
fi

# Validate output
INPUT_LINES=$(wc -l < "$INPUT_GTF")
OUTPUT_LINES=$(wc -l < "$OUTPUT_GTF")

echo ""
echo "Conversion complete!"
echo "  Input lines:  $INPUT_LINES"
echo "  Output lines: $OUTPUT_LINES"

# Show sample of output
echo ""
echo "Sample output (first 3 lines):"
head -3 "$OUTPUT_GTF"

# Validate GTF format
echo ""
echo "Validating output format..."
if head -1 "$OUTPUT_GTF" | grep -q 'transcript_id "[^"]*"; gene_id "[^"]*"'; then
    echo "  ✓ Attribute format looks correct"
else
    echo "  ⚠ Attribute format may need review"
fi

echo ""
echo "Done! Output saved to: $OUTPUT_GTF"
