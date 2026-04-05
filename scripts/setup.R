# ============================================================================ #
# Fast LLQR Project - Setup Script
# Description: Load required packages and source project functions
# Author: Erkang
# Date: 2025-10-09
# ============================================================================ #

cat("=== Fast LLQR Project Setup ===\n")
cat("Working directory:", getwd(), "\n")

# Load required packages for the project
required_packages <- c(
  "ggplot2", "gridExtra", "microbenchmark", "quantreg", "quantdr",
  "Rcpp", "RcppArmadillo", "dplyr", "tidyr", "purrr", "writexl", "readr"
)

# Install missing packages if needed
missing_packages <- required_packages[!required_packages %in% installed.packages()]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Source project-specific functions - 
source_files <- c(
  "R/llqr_simulation_helpers.R",  
  "R/tvcqr_simulation_helpers.R",
  "R/llqr_functions.R",          
  "R/tvcqr_functions.R",         
  "R/add_llqr_functions.R",      
  "R/visualization.R",           
  "R/performance_measurement.R"      
)

for (file in source_files) {
  if (file.exists(file)) {
    source(file)
    cat("✓ Sourced:", file, "\n")
  } else {
    warning("File not found: ", file)
  }
}

# Compile C++ code 
if (file.exists("src/cpp/tvcqr_helpers_new.cpp")) {
  Rcpp::sourceCpp("src/cpp/tvcqr_helpers_new.cpp")  
  cat("✓ Compiled C++ code\n")
} else {
  warning("C++ file not found: src/cpp/tvcqr_helpers_new.cpp")
}

# Load compiled Fortran libraries - 
lib_files <- c(
  "src/fortran/tvcqr_seq.so",        
  "src/fortran/tvcqr_seq_M_acc.so",
  "src/fortran/llqr_seq.so",
  "src/fortran/llqr_ppro.so"
)

for (lib in lib_files) {
  if (file.exists(lib)) {
    dyn.load(lib)
    cat("✓ Loaded:", lib, "\n")
  } else {
    cat("⚠ Library not found:", lib, "\n")
  }
}

cat("Project setup completed successfully!\n")
