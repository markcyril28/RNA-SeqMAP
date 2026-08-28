#!/bin/bash
# cd to script's own directory so **/*.csv glob works from any cwd
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1

# CSVs are grouped in single_project/ and combined_projects/; globstar makes the
# **/*.csv patterns below reach into them. nullglob keeps an empty group quiet.
shopt -s globstar nullglob

echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║                      ORGAN DIVERSITY BY FILE                           ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""

for csv in **/*.csv; do
    [[ -f "$csv" ]] || continue
    organ_count=$(grep "^SRR[0-9]" "$csv" | awk -F',' '{print $2}' | sort -u | wc -l)
    unique_organs=$(grep "^SRR[0-9]" "$csv" | awk -F',' '{print $2}' | sort -u | tr '\n' ', ' | sed 's/,$//')
    sample_count=$(grep "^SRR[0-9]" "$csv" | wc -l)
    
    printf "\n%-70s: %2d samples, %2d organs\n" "$csv" "$sample_count" "$organ_count"
    echo "  Organs: $unique_organs"
done

echo ""
echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║                    FILE RELATIONSHIP MAPPING                           ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "SUBSET RELATIONSHIPS:"
echo "  └─ single_project/PRJNA328564.csv (full dataset, 19 samples)"
echo "     ├─ single_project/PRJNA328564_selected.csv (12 samples - main interest tissues)"
echo "     └─ single_project/PRJNA328564_test.csv (3 samples - test subset)"
echo ""

echo "COMBINED/UNION RELATIONSHIPS:"
echo "  ├─ combined_projects/PRJNA865018_and_PRJNA941250.csv (18 samples - union of:"
echo "  │  ├─ single_project/PRJNA865018.csv (9 samples)"
echo "  │  └─ single_project/PRJNA941250.csv (9 samples)"
echo "  └─ combined_projects/PRJNA328564_selected_plus_PRJNA1417697_selected.csv (17 samples - union of:"
echo "     ├─ single_project/PRJNA328564_selected.csv (12 samples)"
echo "     └─ single_project/PRJNA1417697.csv (5 anther samples)"
echo ""

echo "INDEPENDENT/SUPPLEMENTARY FILES:"
echo "  ├─ single_project/SAMN28540068.csv (8 samples)"
echo "  ├─ single_project/SAMN28540077.csv (9 samples)"
echo "  └─ OTHER_SRR_LIST.csv (3 active samples + 7 commented; not project-scoped)"
