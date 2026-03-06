# Source setup script to load all dependencies
source("scripts/setup.R")

# ============================================================================ #
# Reproduce TVCQR Simulation Results
# ============================================================================ #

# Simulation parameters of TVCQR
sim_config_tvcqr <- list(
  case = 1,              # Data generation case
  tau = 0.5,             # Quantile level
  n = 200,               # Sample size
  num_rep = 10,        # Number of replications
  h = NULL,              # Bandwidth (NULL for automatic selection)
  h.factor = 1,          # Bandwidth factor
  tol = 1e-14,           # Convergence tolerance
  maxit = 1e6,           # Maximum iterations
  bland = FALSE,         # Bland's rule for simplex
  Mm.factor = c(1e-3, 1e-4),  # Multiple Mm factors to test
  eps = 1e-06,           # Epsilon for preprocessing
  cpp_helper = FALSE,    # Use C++ helper functions
  seed_base = 2025      # Base seed for reproducibility
)

# Run TVCQR Simulation
cat("Starting TVCQR simulation...\n")
cat("Configuration:\n")
print(sim_config_tvcqr)
cat("\n")
results_tvcqr <- run_tvcqr_simulation(sim_config_tvcqr)

# Print TVCQR Summary
print_tvcqr_simulation_summary(results_tvcqr)

# Save TVCQR Results
save_tvcqr_simulation_results(results_tvcqr)

# ============================================================================ #
# Reproduce LLQR Simulation Results
# ============================================================================ #

# Simulation parameters of LLQR
sim_config_llqr <- list(
  case = 1,              # Data generation case
  tau = 0.5,             # Quantile level
  n = 500,              # Sample size
  num_rep = 10,         # Number of replications
  h = NULL,              # Bandwidth (NULL for automatic selection)
  z = NULL,              # Evaluation points (NULL for using x)
  tol = 1e-14,           # Convergence tolerance
  maxit = 2e6,           # Maximum iterations
  bland = FALSE,         # Bland's rule for simplex
  track_order = TRUE,   # Track order of evaluation points 
  Mm.factor = c(1,1e-1,1e-2,1e-3, 1e-4),  # Multiple Mm factors to test
  seed_base = 2025       # Base seed for reproducibility
)

# Run LLQR Simulation
cat("Starting LLQR simulation...\n")
cat("Configuration:\n")
print(sim_config_llqr)
cat("\n")
results_llqr <- run_llqr_simulation(sim_config_llqr)

# Print LLQR Summary
print_llqr_simulation_summary(results_llqr)

# Save LLQR Results
save_llqr_simulation_results(results_llqr)

