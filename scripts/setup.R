# ============================================================================ #
# Fast LLQR Project - Setup Script
# Description: Load required packages and source project functions
# Author: [Your Name]
# Date: [Date]
# ============================================================================ #

# Load required packages for the project
required_packages <- c(
  "ggplot2", "gridExtra", "microbenchmark", "quantreg", 
  "Rcpp", "RcppArmadillo", "dplyr", "tidyr", "purrr"
)

# Install missing packages if needed
missing_packages <- required_packages[!required_packages %in% installed.packages()]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Source project-specific functions
source_files <- c(
  "../R/simulation_helpers.R",
  "../R/llqr_functions.R", 
  "../R/tvcqr_functions.R",
  "../R/add_llqr_functions.R",
  "../R/visualization.R",
  "../R/performance_metrics.R"
)

for (file in source_files) {
  if (file.exists(file)) {
    source(file)
  } else {
    warning("File not found: ", file)
  }
}

# Compile C++ code
Rcpp::sourceCpp("../src/cpp/tvcqr_helpers_new.cpp")

# Load compiled Fortran libraries
# dyn.unload("../src/fortran/tvcqr_seq.so")
dyn.load("../src/fortran/tvcqr_seq.so")
# dyn.unload("../src/fortran/tvcqr_seq_M_acc.so")
dyn.load("../src/fortran/tvcqr_seq_M_acc.so")

cat("Project setup completed successfully!\n")