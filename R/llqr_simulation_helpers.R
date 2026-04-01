# ============================================================================ #
# Fast LLQR Project - LLQR Simulation Helper Functions
# Description: Simulation study for local linear quantile regression
# Author: Erkang
# Date: 2025-11-22
# ============================================================================ #

#' Create LLQR Methods List
#' 
#' @param Mm.factor_vec Vector of Mm.factor values to test
#' @return Named list of method functions
#' @export
compute_llqr_rule_bandwidth <- function(x, y, tau, h = NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  m <- nrow(x)
  nvar <- ncol(x)

  if (!is.null(h)) {
    return(as.numeric(h))
  }

  red_dim <- floor(0.2 * m)
  index_y <- order(y)[red_dim:(m - red_dim)]
  h_val <- KernSmooth::dpill(x[index_y, , drop = FALSE], y[index_y])
  h_val <- 1.25 * h_val * (tau * (1 - tau) / (dnorm(qnorm(tau)))^2)^0.2
  if (is.nan(h_val)) {
    h_val <- 1.25 * max(m^(-1 / (nvar + 4)), min(2, sd(y)) * m^(-1 / (nvar + 4)))
  }

  as.numeric(h_val)
}

run_llqr_baseline_with_retry <- function(x, y, config) {
  retry_factors <- config$llqr_h_retry_factors
  if (is.null(retry_factors) || length(retry_factors) == 0) {
    retry_factors <- c(1, 1.25, 1.5, 2, 4, 8)
  }
  retry_factors <- as.numeric(retry_factors)
  retry_factors <- retry_factors[is.finite(retry_factors) & retry_factors > 0]
  if (length(retry_factors) == 0) {
    stop("config$llqr_h_retry_factors must contain positive numeric values.")
  }

  base_h <- compute_llqr_rule_bandwidth(x = x, y = y, tau = config$tau, h = config$h)
  last_error <- NULL

  for (i in seq_along(retry_factors)) {
    h_used <- as.numeric(base_h * retry_factors[i])
    fit <- tryCatch(
      quantdr::llqr(x = x, y = y, tau = config$tau, h = h_used,
                    method = "rule", x0 = NULL),
      error = function(e) e
    )

    if (!inherits(fit, "error")) {
      fit$h <- h_used
      fit$h_used <- h_used
      fit$h_retry_factor <- retry_factors[i]
      fit$llqr_attempts <- i
      return(fit)
    }

    last_error <- fit
    if (!grepl("Singular design matrix", conditionMessage(fit), fixed = TRUE)) {
      stop(fit)
    }
  }

  stop(last_error)
}

create_llqr_methods <- function(Mm.factor_vec) {
  methods <- list()
  include_ppro <- tolower(Sys.getenv("FASTQR_INCLUDE_LLQR_PPRO", unset = "0")) %in%
    c("1", "true", "t", "yes", "y")
  include_ppro_fortran <- tolower(Sys.getenv("FASTQR_INCLUDE_LLQR_PPRO_FORTRAN", unset = "1")) %in%
    c("1", "true", "t", "yes", "y")
  
  # Methods without Mm.factor parameter
  methods$llqr <- function(x, y, config) {
    run_llqr_baseline_with_retry(x = x, y = y, config = config)
  }
  
  methods$llqr_seq <- function(x, y, config) {
    llqr_seq(x = x, y = y, tau = config$tau, z = config$z, 
             h = config$h, tol = config$tol, 
             maxit = config$maxit, bland = config$bland,
             track_order = config$track_order)
  }
  
  methods$llqr_seq_fortran <- function(x, y, config) {
    llqr_seq_fortran_wrapper(x = x, y = y, z = config$z, 
                             tau = config$tau, h = config$h, 
                             tol = config$tol, maxit = config$maxit, 
                             bland = config$bland)
  }
  
  # Methods with Mm.factor parameter - create multiple versions
  for (i in seq_along(Mm.factor_vec)) {
    Mm_val <- Mm.factor_vec[i]
    
    # llqr_ppro is kept opt-in only until the R ppro path is stable.
    if (include_ppro) {
      method_name_ppro <- sprintf("llqr_ppro_%d", i)
      methods[[method_name_ppro]] <- local({
        Mm_factor_local <- Mm_val
        function(x, y, config) {
          llqr_ppro(x = x, y = y, tau = config$tau, z = config$z, 
                    case = config$case,
                    h = config$h, Mm.factor = Mm_factor_local, 
                    track_order = config$track_order)
        }
      })
    }
    
    # llqr_seq_ppro with different Mm.factor values
    method_name <- sprintf("llqr_seq_ppro_%d", i)
    methods[[method_name]] <- local({
      Mm_factor_local <- Mm_val
      function(x, y, config) {
        llqr_seq_ppro(x = x, y = y, tau = config$tau, z = config$z, 
                      case = config$case,
                      h = config$h, tol = config$tol, 
                      maxit = config$maxit, bland = config$bland,
                      Mm.factor = Mm_factor_local, 
                      track_order = config$track_order)
      }
    })
    
    # llqr_seq_ppro_fortran is kept opt-in only until the Fortran path is stable.
    if (include_ppro_fortran) {
      method_name_fortran <- sprintf("llqr_seq_ppro_fortran_%d", i)
      methods[[method_name_fortran]] <- local({
        Mm_factor_local <- Mm_val
        function(x, y, config) {
          llqr_seq_ppro_fortran_wrapper(x = x, y = y, z = config$z, 
                                        tau = config$tau, h = config$h,
                                        Mm.factor = Mm_factor_local,
                                        case = config$case,
                                        tol = config$tol, maxit = config$maxit, 
                                        bland = config$bland)
        }
      })
    }
  }
  
  return(methods)
}

#' Run Single LLQR Replication
#' 
#' @param rep_id Replication ID
#' @param config Configuration list
#' @param methods Named list of method functions
#' @return List containing estimates, H_seq, and timing results
#' @export
run_single_llqr_replication <- function(rep_id, config, methods) {
  seed_used <- if (!is.null(config$seed_used)) {
    as.integer(config$seed_used)
  } else {
    as.integer(config$seed_base + rep_id)
  }

  # Generate data
  data <- generate_data(n = config$n, case = config$case, 
                        seed = seed_used)
  x <- data$x
  y <- data$y
  
  # Initialize storage for this replication
  results <- list(
    estimates = vector("list", length(methods)),
    H_seq = vector("list", length(methods)),
    timing = numeric(length(methods)),
    method_metadata = vector("list", length(methods))
  )
  names(results$estimates) <- names(methods)
  names(results$H_seq) <- names(methods)
  names(results$timing) <- names(methods)
  names(results$method_metadata) <- names(methods)
  
  # Run each method
  for (method_name in names(methods)) {
    # Benchmark the method
    timing_result <- microbenchmark::microbenchmark(
      {
        fit <- methods[[method_name]](x, y, config)
      },
      times = 1,
      unit = "s"
    )
    
    # Store results - all methods return ll_est in the same format
    results$estimates[[method_name]] <- fit$ll_est
    
    # Store H_seq if it exists in the fit object
    if (!is.null(fit$H_seq)) {
      results$H_seq[[method_name]] <- fit$H_seq
    } else {
      results$H_seq[[method_name]] <- NULL
    }

    results$method_metadata[[method_name]] <- list(
      h_used = if (!is.null(fit$h_used)) as.numeric(fit$h_used) else if (!is.null(fit$h)) as.numeric(fit$h) else NA_real_,
      h_retry_factor = if (!is.null(fit$h_retry_factor)) as.numeric(fit$h_retry_factor) else NA_real_,
      llqr_attempts = if (!is.null(fit$llqr_attempts)) as.integer(fit$llqr_attempts) else NA_integer_
    )
    
    # Extract timing in seconds
    results$timing[[method_name]] <- summary(timing_result)$mean
  }

  results$seed_used <- seed_used
  
  return(results)
}

#' Run LLQR Simulation
#' 
#' @param config Configuration list containing simulation parameters
#' @return List containing simulation results
#' @export
run_llqr_simulation <- function(config) {
  cat("========================================\n")
  cat("Starting LLQR Simulation\n")
  cat("========================================\n")
  cat(sprintf("Case: %d\n", config$case))
  cat(sprintf("Tau: %.2f\n", config$tau))
  cat(sprintf("Sample size: %d\n", config$n))
  cat(sprintf("Number of replications: %d\n", config$num_rep))
  cat(sprintf("Mm.factor values: %s\n", paste(config$Mm.factor, collapse = ", ")))
  cat("========================================\n\n")
  
  # Create methods list with multiple Mm.factor values
  llqr_methods <- create_llqr_methods(config$Mm.factor)
  
  # Initialize storage for all replications
  num_methods <- length(llqr_methods)
  method_names <- names(llqr_methods)
  
  cat(sprintf("Total number of methods to compare: %d\n", num_methods))
  cat("Methods:\n")
  for (i in seq_along(method_names)) {
    cat(sprintf("  %d. %s\n", i, method_names[i]))
  }
  cat("\n")
  
  # Timing matrix: num_rep x num_methods
  timing_matrix <- matrix(NA, nrow = config$num_rep, ncol = num_methods)
  colnames(timing_matrix) <- method_names
  
  # Estimates list: one list per method, each containing num_rep vectors
  estimates_list <- vector("list", num_methods)
  names(estimates_list) <- method_names
  for (i in seq_along(method_names)) {
    estimates_list[[i]] <- vector("list", config$num_rep)
  }
  
  # H_seq list: one list per method, each containing num_rep matrices (or NULL)
  H_seq_list <- vector("list", num_methods)
  names(H_seq_list) <- method_names
  for (i in seq_along(method_names)) {
    H_seq_list[[i]] <- vector("list", config$num_rep)
  }
  
  # Run replications with progress tracking
  cat("Running replications...\n")
  pb <- txtProgressBar(min = 0, max = config$num_rep, style = 3)
  
  for (rep in 1:config$num_rep) {
    # Run single replication
    rep_results <- run_single_llqr_replication(rep, config, llqr_methods)
    
    # Store results
    for (method_name in method_names) {
      timing_matrix[rep, method_name] <- rep_results$timing[[method_name]]
      estimates_list[[method_name]][[rep]] <- rep_results$estimates[[method_name]]
      H_seq_list[[method_name]][[rep]] <- rep_results$H_seq[[method_name]]
    }
    
    # Update progress bar
    setTxtProgressBar(pb, rep)
  }
  close(pb)
  
  cat("\n\nSimulation completed!\n")
  
  # Prepare final results
  simulation_results <- list(
    config = config,
    timing_matrix = timing_matrix,
    estimates_list = estimates_list,
    H_seq_list = H_seq_list,
    method_names = method_names,
    Mm.factor_mapping = create_llqr_Mm_factor_mapping(method_names, config$Mm.factor),
    timestamp = Sys.time(),
    simulation_type = "llqr"  # Add identifier
  )
  
  return(simulation_results)
}

#' Create Mm.factor Mapping for LLQR Methods
#' 
#' @param method_names Vector of method names
#' @param Mm.factor_vec Vector of Mm.factor values
#' @return Data frame mapping methods to Mm.factor values
#' @export
create_llqr_Mm_factor_mapping <- function(method_names, Mm.factor_vec) {
  mapping <- data.frame(
    Method = method_names,
    Mm.factor = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(method_names)) {
    method <- method_names[i]
    # Check if method name contains ppro and ends with a number
    if (grepl("ppro.*_(\\d+)$", method)) {
      # Extract the index from method name (the last number)
      idx <- as.numeric(sub(".*_(\\d+)$", "\\1", method))
      mapping$Mm.factor[i] <- as.character(Mm.factor_vec[idx])
    } else {
      mapping$Mm.factor[i] <- "N/A"
    }
  }
  
  return(mapping)
}

#' Save LLQR Simulation Results
#' 
#' @param results Simulation results object
#' @param output_dir Directory to save results
#' @return File path where results were saved
#' @export
save_llqr_simulation_results <- function(results, output_dir = "data/llqr_simu_results") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Generate filename
  config <- results$config
  tau_str <- sprintf("tau%02d", round(config$tau * 100))
  filename <- sprintf("case%d_%s_n%d_rep%d.RData", 
                      config$case, tau_str, config$n, config$num_rep)
  filepath <- file.path(output_dir, filename)
  
  # Save results
  save(results, file = filepath)
  cat(sprintf("\nResults saved to: %s\n", filepath))
  
  return(filepath)
}

#' Print LLQR Simulation Summary
#' 
#' @param results Simulation results object
#' @export
print_llqr_simulation_summary <- function(results) {
  cat("\n========================================\n")
  cat("LLQR Simulation Summary\n")
  cat("========================================\n\n")
  
  timing_matrix <- results$timing_matrix
  
  cat("Timing Statistics (in seconds):\n")
  cat("----------------------------------------\n")
  
  timing_summary <- data.frame(
    Method = colnames(timing_matrix),
    Mm.factor = results$Mm.factor_mapping$Mm.factor,
    Mean = colMeans(timing_matrix, na.rm = TRUE),
    Median = apply(timing_matrix, 2, median, na.rm = TRUE),
    SD = apply(timing_matrix, 2, sd, na.rm = TRUE),
    Min = apply(timing_matrix, 2, min, na.rm = TRUE),
    Max = apply(timing_matrix, 2, max, na.rm = TRUE)
  )
  
  print(timing_summary, row.names = FALSE, digits = 4)
  cat("\n")
  
  # Calculate maximum average relative bias
  cat("\nMaximum Average Relative Bias (relative to llqr):\n")
  cat("----------------------------------------\n")
  
  max_avg_rel_bias <- compute_llqr_max_average_relative_bias(results)
  
  bias_summary <- data.frame(
    Method = names(max_avg_rel_bias),
    Mm.factor = results$Mm.factor_mapping$Mm.factor[
      match(names(max_avg_rel_bias), results$Mm.factor_mapping$Method)
    ],
    Max_Avg_Rel_Bias = max_avg_rel_bias,
    stringsAsFactors = FALSE
  )
  
  print(bias_summary, row.names = FALSE, digits = 6)
  cat("\n")
}

#' Compute Maximum Average Relative Bias for LLQR
#' 
#' @param results Simulation results object
#' @return Named vector of maximum average relative bias values
#' @export
compute_llqr_max_average_relative_bias <- function(results) {
  estimates_list <- results$estimates_list
  method_names <- results$method_names
  num_rep <- results$config$num_rep
  
  # Get llqr estimates as baseline
  if (!"llqr" %in% method_names) {
    stop("llqr method not found in results!")
  }
  
  baseline_estimates <- estimates_list[["llqr"]]
  
  # Initialize vector to store maximum average relative bias
  max_avg_rel_bias <- numeric(length(method_names))
  names(max_avg_rel_bias) <- method_names
  
  # For llqr itself, the bias is 0
  max_avg_rel_bias["llqr"] <- 0
  
  # Calculate for other methods
  for (method_name in setdiff(method_names, "llqr")) {
    method_estimates <- estimates_list[[method_name]]
    
    # Vector to store average relative bias for each replication
    avg_rel_bias <- numeric(num_rep)
    
    for (rep in 1:num_rep) {
      baseline_vec <- baseline_estimates[[rep]]
      
      # Check if method_estimates has this replication
      if (rep > length(method_estimates) || is.null(method_estimates[[rep]])) {
        warning(sprintf("Missing estimates for %s at replication %d", 
                        method_name, rep))
        avg_rel_bias[rep] <- NA
        next
      }
      
      method_vec <- method_estimates[[rep]]
      
      # Convert to numeric vectors if necessary
      baseline_vec <- as.numeric(baseline_vec)
      method_vec <- as.numeric(method_vec)
      
      # Check if vectors have the same length
      if (length(baseline_vec) != length(method_vec)) {
        warning(sprintf("Length mismatch for %s at replication %d: baseline=%d, method=%d", 
                        method_name, rep, length(baseline_vec), length(method_vec)))
        avg_rel_bias[rep] <- NA
        next
      }
      
      # Calculate element-wise relative bias
      abs_diff <- abs(method_vec - baseline_vec)
      abs_baseline <- abs(baseline_vec)
      
      # Avoid division by zero
      threshold <- 1e-10
      rel_bias <- ifelse(abs_baseline < threshold, 
                         0, 
                         abs_diff / abs_baseline)
      
      # Calculate average relative bias for this replication
      avg_rel_bias[rep] <- mean(rel_bias, na.rm = TRUE)
    }
    
    # Take the maximum across all replications
    max_avg_rel_bias[method_name] <- max(avg_rel_bias, na.rm = TRUE)
  }
  
  return(max_avg_rel_bias)
}

#' Find LLQR Replication with Maximum Bias
#' 
#' @param results Simulation results object
#' @param method_name Name of the method to analyze
#' @return List containing replication ID and bias information
#' @export
find_llqr_max_bias_replication <- function(results, method_name) {
  estimates_list <- results$estimates_list
  num_rep <- results$config$num_rep
  
  # Get baseline estimates
  if (!"llqr" %in% names(estimates_list)) {
    stop("llqr method not found in results!")
  }
  
  baseline_estimates <- estimates_list[["llqr"]]
  method_estimates <- estimates_list[[method_name]]
  
  # Vector to store average relative bias for each replication
  avg_rel_bias <- numeric(num_rep)
  
  for (rep in 1:num_rep) {
    baseline_vec <- baseline_estimates[[rep]]
    
    if (rep > length(method_estimates) || is.null(method_estimates[[rep]])) {
      avg_rel_bias[rep] <- NA
      next
    }
    
    method_vec <- method_estimates[[rep]]
    
    # Convert to numeric vectors
    baseline_vec <- as.numeric(baseline_vec)
    method_vec <- as.numeric(method_vec)
    
    # Check lengths
    if (length(baseline_vec) != length(method_vec)) {
      avg_rel_bias[rep] <- NA
      next
    }
    
    # Calculate element-wise relative bias
    abs_diff <- abs(method_vec - baseline_vec)
    abs_baseline <- abs(baseline_vec)
    
    threshold <- 1e-10
    rel_bias <- ifelse(abs_baseline < threshold, 
                       0, 
                       abs_diff / abs_baseline)
    
    # Calculate average relative bias for this replication
    avg_rel_bias[rep] <- mean(rel_bias, na.rm = TRUE)
  }
  
  # Find replication with maximum bias
  max_rep <- which.max(avg_rel_bias)
  max_bias_value <- avg_rel_bias[max_rep]
  
  cat("\n========================================\n")
  cat(sprintf("Maximum Bias Analysis for %s\n", method_name))
  cat("========================================\n")
  cat(sprintf("Replication with maximum bias: %d\n", max_rep))
  cat(sprintf("Maximum average relative bias: %.6e\n", max_bias_value))
  cat("========================================\n\n")
  
  return(list(
    replication_id = max_rep,
    max_bias = max_bias_value,
    all_biases = avg_rel_bias
  ))
}

#' Extract LLQR Data for Specific Replication
#' 
#' @param results Simulation results object
#' @param rep_id Replication ID
#' @return List containing x and y data
#' @export
extract_llqr_replication_data <- function(results, rep_id) {
  config <- results$config
  
  # Regenerate the same data using the stored seed
  data <- generate_data(
    n = config$n, 
    case = config$case, 
    seed = config$seed_base + rep_id
  )
  
  return(data)
}

#' Compare LLQR Estimates for Specific Replication
#' 
#' @param results Simulation results object
#' @param rep_id Replication ID
#' @param methods_to_compare Vector of method names to compare (optional)
#' @export
compare_llqr_replication_estimates <- function(results, rep_id, methods_to_compare = NULL) {
  estimates_list <- results$estimates_list
  method_names <- results$method_names
  
  # If no specific methods provided, use all methods
  if (is.null(methods_to_compare)) {
    methods_to_compare <- method_names
  }
  
  cat("\n========================================\n")
  cat(sprintf("LLQR Estimate Comparison for Replication %d\n", rep_id))
  cat("========================================\n\n")
  
  # Extract estimates for each method
  estimates_comparison <- list()
  for (method in methods_to_compare) {
    if (method %in% names(estimates_list)) {
      estimates_comparison[[method]] <- estimates_list[[method]][[rep_id]]
    }
  }
  
  # Print summary statistics
  cat("Summary Statistics:\n")
  cat("----------------------------------------\n")
  for (method in names(estimates_comparison)) {
    est <- estimates_comparison[[method]]
    est_vec <- as.numeric(est)
    
    cat(sprintf("\n%s:\n", method))
    cat(sprintf("  Length: %d\n", length(est_vec)))
    cat(sprintf("  Mean: %.6f\n", mean(est_vec, na.rm = TRUE)))
    cat(sprintf("  Median: %.6f\n", median(est_vec, na.rm = TRUE)))
    cat(sprintf("  SD: %.6f\n", sd(est_vec, na.rm = TRUE)))
    cat(sprintf("  Min: %.6f\n", min(est_vec, na.rm = TRUE)))
    cat(sprintf("  Max: %.6f\n", max(est_vec, na.rm = TRUE)))
    cat(sprintf("  NA count: %d\n", sum(is.na(est_vec))))
  }
  
  # Calculate pairwise differences with llqr
  if ("llqr" %in% names(estimates_comparison)) {
    cat("\n\nDifferences from llqr:\n")
    cat("----------------------------------------\n")
    baseline <- as.numeric(estimates_comparison[["llqr"]])
    
    for (method in setdiff(names(estimates_comparison), "llqr")) {
      method_est <- as.numeric(estimates_comparison[[method]])
      
      if (length(baseline) == length(method_est)) {
        diff <- method_est - baseline
        abs_diff <- abs(diff)
        rel_diff <- abs_diff / (abs(baseline) + 1e-10)
        
        cat(sprintf("\n%s vs llqr:\n", method))
        cat(sprintf("  Mean absolute diff: %.6e\n", mean(abs_diff, na.rm = TRUE)))
        cat(sprintf("  Max absolute diff: %.6e\n", max(abs_diff, na.rm = TRUE)))
        cat(sprintf("  Mean relative diff: %.6e\n", mean(rel_diff, na.rm = TRUE)))
        cat(sprintf("  Max relative diff: %.6e\n", max(rel_diff, na.rm = TRUE)))
      } else {
        cat(sprintf("\n%s vs llqr: Length mismatch (llqr=%d, %s=%d)\n", 
                    method, length(baseline), method, length(method_est)))
      }
    }
  }
  
  return(estimates_comparison)
}

#' Save Problematic LLQR Data
#' 
#' @param results Simulation results object
#' @param method_name Name of the method to analyze
#' @param output_dir Directory to save data
#' @return List containing filepath and problematic case data
#' @export
save_llqr_problematic_data <- function(results, method_name, 
                                       output_dir = "data/llqr_problematic_cases") {
  # Find replication with max bias
  max_bias_info <- find_llqr_max_bias_replication(results, method_name)
  rep_id <- max_bias_info$replication_id
  
  # Extract data
  data <- extract_llqr_replication_data(results, rep_id)
  
  # Compare estimates
  estimates_comparison <- compare_llqr_replication_estimates(results, rep_id)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Prepare data to save
  problematic_case <- list(
    config = results$config,
    replication_id = rep_id,
    max_bias = max_bias_info$max_bias,
    all_biases = max_bias_info$all_biases,
    data = data,
    estimates = estimates_comparison,
    method_name = method_name,
    timestamp = Sys.time(),
    simulation_type = "llqr"
  )
  
  # Generate filename
  config <- results$config
  tau_str <- sprintf("tau%02d", round(config$tau * 100))
  filename <- sprintf("llqr_problematic_case%d_%s_n%d_rep%d_%s.RData", 
                      config$case, tau_str, config$n, 
                      rep_id, gsub("[^[:alnum:]]", "_", method_name))
  filepath <- file.path(output_dir, filename)
  
  # Save
  save(problematic_case, file = filepath)
  cat(sprintf("\nProblematic data saved to: %s\n", filepath))
  
  return(list(
    filepath = filepath,
    problematic_case = problematic_case
  ))
}

#' Analyze LLQR H_seq Convergence
#' 
#' @param results Simulation results object
#' @param method_name Name of the method to analyze
#' @param rep_id Replication ID (optional, if NULL analyzes all replications)
#' @export
analyze_llqr_H_seq_convergence <- function(results, method_name, rep_id = NULL) {
  H_seq_list <- results$H_seq_list
  
  if (!method_name %in% names(H_seq_list)) {
    stop(sprintf("Method %s not found in H_seq_list!", method_name))
  }
  
  method_H_seq <- H_seq_list[[method_name]]
  
  if (is.null(rep_id)) {
    # Analyze all replications
    cat("\n========================================\n")
    cat(sprintf("H_seq Convergence Analysis for %s\n", method_name))
    cat("(All Replications)\n")
    cat("========================================\n\n")
    
    num_iterations <- numeric(length(method_H_seq))
    
    for (i in seq_along(method_H_seq)) {
      if (!is.null(method_H_seq[[i]])) {
        num_iterations[i] <- nrow(method_H_seq[[i]])
      } else {
        num_iterations[i] <- NA
      }
    }
    
    cat(sprintf("Number of replications: %d\n", length(method_H_seq)))
    cat(sprintf("Mean iterations: %.2f\n", mean(num_iterations, na.rm = TRUE)))
    cat(sprintf("Median iterations: %.2f\n", median(num_iterations, na.rm = TRUE)))
    cat(sprintf("Min iterations: %d\n", min(num_iterations, na.rm = TRUE)))
    cat(sprintf("Max iterations: %d\n", max(num_iterations, na.rm = TRUE)))
    cat(sprintf("SD iterations: %.2f\n", sd(num_iterations, na.rm = TRUE)))
    cat(sprintf("NA count: %d\n", sum(is.na(num_iterations))))
    
    return(num_iterations)
    
  } else {
    # Analyze specific replication
    cat("\n========================================\n")
    cat(sprintf("H_seq Convergence Analysis for %s\n", method_name))
    cat(sprintf("(Replication %d)\n", rep_id))
    cat("========================================\n\n")
    
    H_seq <- method_H_seq[[rep_id]]
    
    if (is.null(H_seq)) {
      cat("No H_seq data available for this replication.\n")
      return(NULL)
    }
    
    cat(sprintf("Number of iterations: %d\n", nrow(H_seq)))
    cat(sprintf("Number of points: %d\n", ncol(H_seq)))
    cat("\nH values at final iteration:\n")
    print(H_seq[nrow(H_seq), ])
    
    return(H_seq)
  }
}
