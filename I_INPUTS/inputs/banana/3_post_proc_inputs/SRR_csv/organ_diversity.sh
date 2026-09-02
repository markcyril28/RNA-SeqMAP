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
    
    mer_count=$(grep "^SRR[0-9]" "$csv" | awk -F',' '$3=="Meristematic"' | wc -l)
    non_count=$(grep "^SRR[0-9]" "$csv" | awk -F',' '$3=="Non_meristematic"' | wc -l)

    printf "\n%-45s: %2d samples, %2d organs\n" "$csv" "$sample_count" "$organ_count"
    echo "  Organs: $unique_organs"
    printf "  Tissue_Class: %d meristematic / %d non-meristematic\n" "$mer_count" "$non_count"
done

echo ""
echo "╔════════════════════════════════════════════════════════════════════════╗"
echo "║                    FILE RELATIONSHIP MAPPING                           ║"
echo "╚════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "FOLDER LAYOUT:"
echo "  single_project/     one BioProject per CSV"
echo "  combined_projects/  CSVs spanning two or more BioProjects"
echo ""

echo "SUBSET RELATIONSHIPS:"
echo "  └─ single_project/PRJNA1237608.csv (full atlas, 34 samples / 12 tissues, cv. Grand Naine)"
echo "     ├─ single_project/PRJNA1237608_selected.csv (12 samples - one deepest rep per tissue)"
echo "     └─ single_project/PRJNA1237608_test.csv (3 samples - smallest runs, distinct tissues)"
echo ""

echo "COMBINED/UNION RELATIONSHIPS:"
echo "  └─ combined_projects/PRJNA1237608_selected_plus_PRJNA606623_PRJNA935075.csv (18 samples - union of:"
echo "     ├─ single_project/PRJNA1237608_selected.csv (12 samples)"
echo "     └─ combined_projects/PRJNA606623_and_PRJNA935075.csv (6 samples)"
echo ""

echo "INDEPENDENT/SUPPLEMENTARY FILES:"
echo "  └─ combined_projects/PRJNA606623_and_PRJNA935075.csv (6 samples - CROSS-SPECIES M. itinerans:"
echo "     ├─ PRJNA606623 male flower bud (3 wild-type reps)"
echo "     └─ PRJNA935075 pollen (3 reps; cold-stress study CONTROL arm, not untreated)"
echo "     NOTE: different species from the M. acuminata reference - expect lower mapping."
