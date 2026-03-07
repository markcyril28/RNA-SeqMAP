#!/usr/bin/env Rscript
# Install necessary CRAN and Bioconductor packages used by the pipeline
# Usage: Rscript install_R_packages.R
# Optional: Rscript install_R_packages.R --with-gpu  (installs GPU packages)

args <- commandArgs(trailingOnly = TRUE)
INSTALL_GPU <- "--with-gpu" %in% args

required_cran <- c(
  "tidyverse",
  "dplyr",
  "tibble",
  "readr",
  "ggplot2",
  "RColorBrewer",
  "circlize",
  "getopt",
  "grid",
  "pheatmap",
  "ggrepel",
  "scales",
  "dynamicTreeCut",
  "fastcluster",
  "Rtsne",
  "umap",
  "factoextra",
  "igraph"
)

required_bioc <- c(
  "DESeq2",
  "ComplexHeatmap",
  "ballgown",
  "tximport",
  "tximeta",
  "AnnotationDbi",
  "org.Mm.eg.db",
  "WGCNA",
  "fgsea",
  "clusterProfiler",
  "enrichplot",
  "DOSE"
)

install_if_missing_cran <- function(pkgs) {
  to_install <- pkgs[!(pkgs %in% installed.packages()[, "Package"]) ]
  if (length(to_install) > 0) {
    message("Installing CRAN packages: ", paste(to_install, collapse = ", "))
    install.packages(to_install, repos = "https://cloud.r-project.org")
  } else {
    message("All CRAN packages already installed.")
  }
}

install_if_missing_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  to_install <- pkgs[!(pkgs %in% installed.packages()[, "Package"]) ]
  if (length(to_install) > 0) {
    message("Installing Bioconductor packages: ", paste(to_install, collapse = ", "))
    BiocManager::install(to_install, ask = FALSE, update = FALSE)
  } else {
    message("All Bioconductor packages already installed.")
  }
}

message("Checking / installing CRAN packages...")
install_if_missing_cran(required_cran)

message("Checking / installing Bioconductor packages...")
install_if_missing_bioc(required_bioc)

# ===============================================
# GPU PACKAGES (optional)
# ===============================================
# torch: PyTorch-based GPU acceleration for matrix operations
# Requires CUDA toolkit - uses conda CUDA 12.1 for compatibility

install_gpu_packages <- function() {
  message("\n=== GPU Package Installation ===")
  
  # Check for conda CUDA first (preferred)
  conda_cuda_version <- Sys.getenv("TORCH_CUDA_VERSION", unset = "")
  conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
  
  if (nzchar(conda_cuda_version)) {
    message("Using conda CUDA version: ", conda_cuda_version)
    message("CONDA_PREFIX: ", conda_prefix)
  }
  
  # Check for CUDA availability via nvidia-smi
  cuda_available <- tryCatch({
    result <- system("nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null", intern = TRUE)
    length(result) > 0
  }, error = function(e) FALSE)
  
  if (!cuda_available) {
    message("WARNING: No NVIDIA GPU detected. GPU packages may not work.")
    message("Continuing with installation anyway...")
  } else {
    # Check system CUDA version
    sys_cuda <- tryCatch({
      output <- system("nvidia-smi 2>/dev/null", intern = TRUE)
      cuda_line <- grep("CUDA Version", output, value = TRUE)
      if (length(cuda_line) > 0) {
        gsub(".*CUDA Version: ([0-9]+\\.[0-9]+).*", "\\1", cuda_line[1])
      } else NULL
    }, error = function(e) NULL)
    
    if (!is.null(sys_cuda)) {
      message("System CUDA version: ", sys_cuda)
      if (as.numeric(sys_cuda) >= 13.0 && !nzchar(conda_cuda_version)) {
        message("WARNING: System CUDA ", sys_cuda, " is not supported by R torch (max 12.1)")
        message("Activate conda environment with CUDA 12.1 for GPU support")
        message("Run: source activate gea")
      }
    }
  }
  
  # Install torch (recommended for GPU acceleration)
  if (!requireNamespace("torch", quietly = TRUE)) {
    message("Installing torch...")
    install.packages("torch", repos = "https://cloud.r-project.org")
    
    # Install CUDA-enabled version - use conda CUDA version if available
    install_cuda_version <- if (nzchar(conda_cuda_version)) conda_cuda_version else NULL
    
    tryCatch({
      if (!is.null(install_cuda_version) && install_cuda_version %in% c("12.1", "11.8", "11.7")) {
        message("Installing torch with CUDA ", install_cuda_version, "...")
        torch::install_torch(type = "cuda", version = install_cuda_version)
      } else {
        torch::install_torch(type = "cuda")
      }
      message("torch with CUDA support installed successfully")
    }, error = function(e) {
      message("Note: torch CUDA installation failed: ", e$message)
      message("Installing CPU version as fallback...")
      tryCatch({
        torch::install_torch()
        message("torch CPU version installed")
      }, error = function(e2) {
        message("Warning: torch installation failed: ", e2$message)
      })
    })
  } else {
    message("torch already installed")
    # Check CUDA availability
    tryCatch({
      if (torch::cuda_is_available()) {
        message("torch CUDA support: AVAILABLE (", torch::cuda_device_count(), " device(s))")
      } else {
        message("torch CUDA support: NOT AVAILABLE (CPU mode)")
        if (!nzchar(conda_cuda_version)) {
          message("TIP: Activate conda environment with CUDA 12.1: source activate gea")
        }
      }
    }, error = function(e) {
      message("torch CUDA check failed: ", e$message)
    })
  }
  
  # gpuR is an alternative but requires OpenCL setup
  # Uncomment if needed:
  # if (!requireNamespace("gpuR", quietly = TRUE)) {
  #   message("Installing gpuR...")
  #   install.packages("gpuR", repos = "https://cloud.r-project.org")
  # }
}

if (INSTALL_GPU) {
  install_gpu_packages()
} else {
  message("\nNote: GPU packages not installed. Run with --with-gpu flag to install:")
  message("  Rscript install_R_packages.R --with-gpu")
}

message("\nDone. You can now load required packages in R (e.g. library(ballgown)).")
