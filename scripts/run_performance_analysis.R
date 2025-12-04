# ============================================================================ #
# Fast LLQR Project - Run Performance Analysis
# Description: Script to run performance analysis for simulation results
# Author: Erkang
# Date: 2025-12-04
# ============================================================================ #

# Source setup script
source("scripts/setup.R")

# ============================================================================ #
# Example 1: TVCQR Performance Analysis
# ============================================================================ #

if (TRUE) {
  # Run TVCQR analysis
  tvcqr_results <- run_performance_analysis(
    case = 1,
    tau = c(0.5),
    n = c(200, 500, 1000, 2000),
    rep = 500,
    model = "tvcqr",
    metrics = c("computation_time", "average_relative_bias"),
    probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
    output_format = "xlsx",
    output_dir = "results/table",
    filename_prefix = "tvcqr_result"
  )
  
  # View results
  print(head(tvcqr_results$computation_time, 20))
  print(head(tvcqr_results$average_relative_bias, 20))
}

# ============================================================================ #
# Example 2: LLQR Performance Analysis
# ============================================================================ #

if (FALSE) {
  # Run LLQR analysis
  llqr_results <- run_performance_analysis(
    case = 1,
    tau = c(0.5),
    n = c(200, 500, 1000, 2000),
    rep = 500,
    model = "llqr",
    metrics = c("computation_time", "average_relative_bias", "iteration_numbers"),
    probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
    output_format = "xlsx",
    output_dir = "results/table",
    filename_prefix = "llqr_result"
  )
  
  # View results
  print(head(llqr_results$computation_time, 20))
  print(head(llqr_results$average_relative_bias, 20))
  print(head(llqr_results$iteration_numbers, 20))
}

# ============================================================================ #
# Example 3: Batch Analysis for Multiple Cases
# ============================================================================ #

if (FALSE) {
  # Configuration
  cases <- c(1, 2, 3)
  tau_values <- c(0.2, 0.5, 0.8)
  sample_sizes <- c(200, 500, 1000, 2000)
  num_replications <- 500
  
  # Run for TVCQR
  for (case_num in cases) {
    tryCatch({
      cat(sprintf("\n\nProcessing TVCQR Case %d...\n", case_num))
      
      run_performance_analysis(
        case = case_num,
        tau = tau_values,
        n = sample_sizes,
        rep = num_replications,
        model = "tvcqr",
        metrics = c("computation_time", "average_relative_bias"),
        output_format = "xlsx",
        output_dir = "results/table",
        filename_prefix = sprintf("tvcqr_case%d_result", case_num)
      )
      
    }, error = function(e) {
      cat(sprintf("Error processing TVCQR Case %d: %s\n", case_num, e$message))
    })
  }
  
  # Run for LLQR
  for (case_num in cases) {
    tryCatch({
      cat(sprintf("\n\nProcessing LLQR Case %d...\n", case_num))
      
      run_performance_analysis(
        case = case_num,
        tau = tau_values,
        n = sample_sizes,
        rep = num_replications,
        model = "llqr",
        metrics = c("computation_time", "average_relative_bias", "iteration_numbers"),
        output_format = "xlsx",
        output_dir = "results/table",
        filename_prefix = sprintf("llqr_case%d_result", case_num)
      )
      
    }, error = function(e) {
      cat(sprintf("Error processing LLQR Case %d: %s\n", case_num, e$message))
    })
  }
}

# ============================================================================ #
# Example 4: Debug - Check data structure
# ============================================================================ #

if (FALSE) {
  # Load a single result to inspect structure
  results_list <- load_simulation_results(
    case = 1,
    tau = c(0.5),
    n = c(200),
    rep = 500,
    model = "tvcqr"
  )
  
  # Get first successful result
  for (result in results_list) {
    if (result$status == "success") {
      cat("\n=== Data Structure ===\n")
      cat("Names in result$data:\n")
      print(names(result$data))
      
      if (!is.null(result$data$estimates_list)) {
        cat("\nLength of estimates_list:", length(result$data$estimates_list), "\n")
        cat("Names in first replication:\n")
        print(names(result$data$estimates_list[[1]]))
        
        cat("\nFirst method estimate (first 10 values):\n")
        first_method <- names(result$data$estimates_list[[1]])[1]
        print(head(result$data$estimates_list[[1]][[first_method]], 10))
        
        cat("\nDimension of first method estimate:\n")
        print(dim(result$data$estimates_list[[1]][[first_method]]))
      }
      
      break
    }
  }
}
