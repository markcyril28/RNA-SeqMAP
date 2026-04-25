#!/usr/bin/env Rscript

# ===============================================
# BATCH DISPATCHER — R Session Batching
# ===============================================
# Runs multiple analyses in a single R session, eliminating ~2-3s of startup
# overhead per additional analysis (R interpreter init, package loading,
# shared config sourcing, GPU detection, sample label loading).
#
# Usage: Rscript batch_dispatcher.R Matrix_Creation Basic_Heatmap Heatmap_with_CV
#
# Each analysis runs with tryCatch error isolation — one failure does not kill
# the batch. Memory is freed between tasks via gc().
#
# Backward compatibility: standalone Rscript calls still work (each analysis
# script retains its if (!interactive() && identical(environment(), globalenv()))
# guard for direct execution).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript batch_dispatcher.R <analysis1> [analysis2] ...")
}

SCRIPT_DIR <- Sys.getenv("ANALYSIS_MODULES_DIR", {
  if (nzchar(Sys.getenv("WF_MANAGED_ENV", "")))
    stop("[BATCH DISPATCHER] ANALYSIS_MODULES_DIR is required under workflow manager (WF_MANAGED_ENV is set).")
  "."
})

# ── Source shared modules ONCE ──────────────────────────────────────────────
# These are the common dependencies for all analysis scripts. Each analysis
# script also tries to source them, but guards prevent redundant execution:
#   - 0_shared_config.R: .SHARED_CONFIG_INITIALIZED guard (Section 5)
#   - 1_utility_functions.R: exists("CURRENT_METHOD") guard
#   - 2_processing_engine.R: exists("CURRENT_METHOD") guard (hard stop)
# Flag: tells analysis scripts not to auto-execute their main block when sourced.
# We source into globalenv() (so function defs and top-level vars are accessible),
# which would normally trigger the if (!interactive() && identical(environment(), globalenv()))
# guards. This flag suppresses that auto-execution — each script checks
# !isTRUE(get0(".BATCH_DISPATCHER_ACTIVE")) before running its main block.
.BATCH_DISPATCHER_ACTIVE <- TRUE

source(file.path(SCRIPT_DIR, "0_shared_config.R"))
source(file.path(SCRIPT_DIR, "1_utility_functions.R"))
source(file.path(SCRIPT_DIR, "3_Matrix_Creation_utils.R"))

# ── Lazy-load analysis modules on first use ─────────────────────────────────
# Avoids loading ComplexHeatmap/tximport when only Matrix_Creation is requested,
# and avoids loading tximport when only heatmaps are requested.
.batch_loaded <- new.env(parent = emptyenv())

.ensure_module <- function(module_name) {
  if (!is.null(.batch_loaded[[module_name]])) return(invisible())
  # source() inside a function defaults to the function's local env, not globalenv().
  # Use local = globalenv() to ensure function definitions (run_basic_heatmap, etc.)
  # and top-level variables (HEATMAP_OUT_DIR, etc.) land in the global environment
  # where the dispatch table and analysis functions can find them.
  switch(module_name,
    "matrix_creation" = {
      source(file.path(SCRIPT_DIR, "3_Matrix_Creation.R"), local = globalenv())
      source(file.path(SCRIPT_DIR, "3_Matrix_Creation_STAR.R"), local = globalenv())
      source(file.path(SCRIPT_DIR, "3_Matrix_Creation_Salmon.R"), local = globalenv())
      source(file.path(SCRIPT_DIR, "3_Matrix_Creation_RSEM.R"), local = globalenv())
    },
    "heatmaps" = {
      source(file.path(SCRIPT_DIR, "2_processing_engine.R"), local = globalenv())
      source(file.path(SCRIPT_DIR, "4_Basic_Heatmap.R"), local = globalenv())
      source(file.path(SCRIPT_DIR, "5_Heatmap_with_CV.R"), local = globalenv())
    }
  )
  .batch_loaded[[module_name]] <- TRUE
}

# ── Dispatch table ──────────────────────────────────────────────────────────
# Maps analysis names (as used in TOML configs and bash ANALYSES arrays) to
# their module group and entry-point function.
.dispatch <- list(
  Matrix_Creation = list(
    module = "matrix_creation",
    fn = function() {
      method_type <- get_method_type(CURRENT_METHOD)
      switch(method_type,
        "star"   = run_star_matrix_creation(),
        "salmon" = run_salmon_matrix_creation(),
        "rsem"   = run_rsem_matrix_creation(),
        run_matrix_creation_main()  # M1/M2 fallback
      )
    }
  ),
  Basic_Heatmap = list(
    module = "heatmaps",
    fn = function() run_basic_heatmap()
  ),
  Heatmap_with_CV = list(
    module = "heatmaps",
    fn = function() run_cv_heatmap()
  )
)

# ── Execute each task with error isolation ──────────────────────────────────
failures <- 0L
for (task in args) {
  entry <- .dispatch[[task]]
  if (is.null(entry)) {
    message("[BATCH] Unknown analysis: ", task, " \u2014 skipping")
    failures <- failures + 1L
    next
  }

  message("\n[BATCH] ", strrep("=", 50))
  message("[BATCH] Starting: ", task)
  message("[BATCH] ", strrep("=", 50))

  .ensure_module(entry$module)

  result <- tryCatch(
    { entry$fn(); TRUE },
    error = function(e) {
      message("[BATCH] FAILED: ", task, " \u2014 ", conditionMessage(e))
      FALSE
    }
  )

  if (!isTRUE(result)) failures <- failures + 1L

  # Free memory between tasks (R does not always collect between function calls)
  gc(verbose = FALSE)
}

# ── Summary ─────────────────────────────────────────────────────────────────
message("\n[BATCH] ", strrep("=", 50))
message("[BATCH] Completed ", length(args), " task(s), ", failures, " failure(s)")
message("[BATCH] ", strrep("=", 50))

if (failures > 0) quit(status = 1)