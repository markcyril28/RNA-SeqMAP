#!/bin/bash

echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║                      ORGAN DIVERSITY BY FILE                           ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""

for csv in $(ls -1 *.csv | sort); do
    organ_count=$(grep "^SRR" "$csv" | awk -F',' '{print $2}' | sort -u | wc -l)
    unique_organs=$(grep "^SRR" "$csv" | awk -F',' '{print $2}' | sort -u | tr '\n' ', ' | sed 's/,$//')
    sample_count=$(grep "^SRR" "$csv" | wc -l)
    
    printf "\n%-45s: %2d samples, %2d organs\n" "$csv" "$sample_count" "$organ_count"
    echo "  Organs: $unique_organs"
done

echo ""
echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║                    FILE RELATIONSHIP MAPPING                           ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "SUBSET RELATIONSHIPS:"
echo "  └─ PRJNA328564.csv (full dataset, 19 samples)"
echo "     ├─ PRJNA328564_selected.csv (12 samples - main interest tissues)"
echo "     └─ PRJNA328564_test.csv (3 samples - test subset)"
echo ""

echo "COMBINED/UNION RELATIONSHIPS:"
echo "  └─ PRJNA865018_and_PRJNA941250.csv (18 samples - union of:"
echo "     ├─ PRJNA865018.csv (9 samples)"
echo "     └─ PRJNA941250.csv (9 samples)"
echo ""

echo "INDEPENDENT/SUPPLEMENTARY FILES:"
echo "  ├─ SAMN28540068.csv (8 samples)"
echo "  ├─ SAMN28540077.csv (9 samples)"
echo "  └─ OTHER_SRR_LIST.csv (3 active samples + 7 commented)"
