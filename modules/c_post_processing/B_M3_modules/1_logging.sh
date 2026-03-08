#!/bin/bash
# ===============================================
# LOGGING FUNCTIONS FOR METHOD 3 POST-PROCESSING
# ===============================================

log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_warn() {
    echo "[WARN] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_step() {
    echo "========================================"
    echo "[STEP] $1"
    echo "========================================"
}
