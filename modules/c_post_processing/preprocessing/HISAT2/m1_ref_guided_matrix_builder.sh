#!/bin/bash

# ===============================================
# M1 HISAT2 Ref-Guided Matrix Builder (thin wrapper)
# ===============================================
# Delegates to the unified stringtie_matrix_builder.sh with STRINGTIE_METHOD=M1.
# All M1-specific defaults (paths, column indices, abundance suffix) are set
# automatically by the unified script when STRINGTIE_METHOD=M1.
#
# Called by: prepde_matrix_linker.sh (as part of M1 preprocessing)
# ===============================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export STRINGTIE_METHOD="M1"
exec bash "$SCRIPT_DIR/stringtie_matrix_builder.sh"
