#!/usr/bin/env Rscript

# ===============================================
# CONCORDANCE BATCH DISPATCHER — R Session Batching
# ===============================================
# Runs concordance Steps 1-3 in a single R session, eliminating 2x R startup
# overhead (~2-3s per additional step: interpreter init, package loading,
# 0_concordance_config.R sourcing, 0_shared_config.R sourcing, GPU detection).
#
# Usage: Rscript concordance_batch_dispatcher.R \
#          --step1=1_load_matrices.R \
#          --step2=2_quantification_concordance.R \
#          --step3=4_generate_report.R
#
# Arguments are optional: omit any --stepN to skip that step.
# Step 1 must complete before Step 2; Step 2 must complete before Step 3.
#
# Backward compatibility: standalone Rscript calls still work (each concordance
# script sources its own config and has no execution guard to modify).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript concordance_batch_dispatcher.R --step1=<script> [--step2=<script>] [--step3=<script>]")
}

# Parse arguments: --stepN=script_name.R
steps <- list()
for (arg in args) {
  if (grepl("^--step[123]=", arg)) {
    key <- sub("^--(step[123])=.*", "\\1", arg)
    val <- sub("^--step[123]=", "", arg)
    steps[[key]] <- val
  }
}

if (length(steps) == 0) {
  stop("No valid --stepN=<script> arguments provided")
}

CONCORDANCE_SCRIPT_DIR <- Sys.getenv("CONCORDANCE_SCRIPT_DIR", {
  if (nzchar(Sys.getenv("WF_MANAGED_ENV", "")))
    stop("[CONCORDANCE BATCH] CONCORDANCE_SCRIPT_DIR is required under workflow manager (WF_MANAGED_ENV is set).")
  "."
})

# Source shared config ONCE (all 3 steps re-source this; we do it first)
.config_path <- file.path(CONCORDANCE_SCRIPT_DIR, "0_concordance_config.R")
if (!file.exists(.config_path)) {
  stop("[CONCORDANCE BATCH] Configuration file not found: ", .config_path,
       "\n  Ensure CONCORDANCE_SCRIPT_DIR points to the correct directory.")
}
source(.config_path)
# Set guard so child scripts skip redundant re-sourcing
.CONC_BATCH_CONFIG_LOADED <- TRUE

# Load ComplexHeatmap/circlize/grid ONCE (Steps 2 and 3 both load these)
# Lazy: only load if Step 2 or Step 3 is requested
if (!is.null(steps$step2) || !is.null(steps$step3)) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
  })
}

# Source method loaders ONCE (Step 1 sources this; has its own .METHOD_LOADERS_SOURCED guard)
if (!is.null(steps$step1)) {
  source(file.path(CONCORDANCE_SCRIPT_DIR, "0_method_loaders.R"))
}

# Execute each step with error isolation
total_requested <- length(steps)
completed <- 0L
failures <- 0L
ordered_steps <- c("step1", "step2", "step3")

for (step_key in ordered_steps) {
  script_name <- steps[[step_key]]
  if (is.null(script_name)) next

  script_path <- file.path(CONCORDANCE_SCRIPT_DIR, script_name)
  if (!file.exists(script_path)) {
    message("[CONC-BATCH] Script not found: ", script_path, " -- skipping")
    failures <- failures + 1L
    next
  }

  step_label <- switch(step_key,
    step1 = "Load & Harmonize Matrices",
    step2 = "Quantification Concordance",
    step3 = "Generate Report"
  )

  message("\n[CONC-BATCH] ", strrep("=", 50))
  message("[CONC-BATCH] Starting: ", step_label, " (", script_name, ")")
  message("[CONC-BATCH] ", strrep("=", 50))

  result <- tryCatch(
    {
      # Source the script into a child environment that inherits globalenv
      # so it can see all variables set by 0_concordance_config.R, but
      # its local assignments do not pollute globalenv across steps.
      #
      # EXCEPTION: Step 1 must write to globalenv because Steps 2/3 read
      # HARMONIZED_RDS and other variables it sets. Use globalenv() for Step 1.
      # Steps 2/3 are terminal (their outputs are files, not R objects).
      if (step_key == "step1") {
        source(script_path, local = globalenv())
      } else {
        source(script_path, local = new.env(parent = globalenv()))
      }
      TRUE
    },
    error = function(e) {
      message("[CONC-BATCH] FAILED: ", step_label, " -- ", conditionMessage(e))
      FALSE
    }
  )

  if (isTRUE(result)) {
    completed <- completed + 1L
  }

  if (!isTRUE(result)) {
    failures <- failures + 1L
    # Step dependency: if Step 1 fails, skip Steps 2 and 3
    if (step_key == "step1") {
      message("[CONC-BATCH] Step 1 failed -- skipping remaining steps (dependency)")
      break
    }
    # If Step 2 fails, skip Step 3 (report needs concordance results)
    if (step_key == "step2" && !is.null(steps$step3)) {
      message("[CONC-BATCH] Step 2 failed -- skipping Step 3 (dependency)")
      steps$step3 <- NULL
    }
  }

  # Check skip sentinel after Step 1 (e.g., zero common genes in cross_genome mode).
  # Step 1 writes this file to signal that further analysis is not meaningful.
  if (step_key == "step1" && isTRUE(result)) {
    .skip_sentinel <- file.path(OUTPUT_DIR, ".skip_sentinel")
    if (file.exists(.skip_sentinel)) {
      .reason <- readLines(.skip_sentinel, n = 1L, warn = FALSE)
      message("[CONC-BATCH] Skip sentinel detected: ", .reason)
      message("[CONC-BATCH] Skipping remaining steps (no data to analyze)")
      break
    }
  }

  # Free memory between steps
  gc(verbose = FALSE)
}

# Summary
message("\n[CONC-BATCH] ", strrep("=", 50))
message("[CONC-BATCH] Completed ", completed, " of ", total_requested, " step(s), ", failures, " failure(s)")
message("[CONC-BATCH] ", strrep("=", 50))

if (failures > 0) quit(status = 1)
