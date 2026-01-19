# scripts/setup_hpc.R
# HPC-safe setup: check only, never install

required_pkgs <- c(
  "microbenchmark", "quantreg", "doParallel", "foreach",
  "doRNG", "dplyr", "ggplot2", "gridExtra"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages on HPC: ", paste(missing_pkgs, collapse = ", "),
       "\nPlease install them on a login/test node before running jobs.")
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

# native libs / compiled code
dyn.load("src/fortran/llqr_seq.so")
dyn.load("src/fortran/llqr_ppro.so")
dyn.load("src/fortran/tvcqr_seq.so")
dyn.load("src/fortran/tvcqr_seq_M_acc.so")

# source project functions
source("R/llqr_functions.R")
source("R/llqr_simulation_helpers.R")
source("R/tvcqr_functions.R")
source("R/tvcqr_simulation_helpers.R")
source("R/performance_measurement.R")
