# ============================================================================ #
# TVCQR Simulation Helper Functions
# ============================================================================ #

complete_tvcqr_sim_config <- function(config) {
  case_info <- resolve_tvcqr_case(config$case)
  config$case <- case_info$case_id
  config$case_key <- case_info$case_key
  config$case_label <- case_info$case_label
  
  config$seed_base <- if (is.null(config$seed_base)) 2025L else as.integer(config$seed_base)
  config$J <- if (is.null(config$J)) 100L else as.integer(config$J)
  config$burn_in <- if (is.null(config$burn_in)) 500L else as.integer(config$burn_in)
  
  if (is.na(config$seed_base)) {
    stop("config$seed_base must be an integer.")
  }
  if (is.na(config$J) || config$J < 0L) {
    stop("config$J must be a non-negative integer.")
  }
  if (is.na(config$burn_in) || config$burn_in < 0L) {
    stop("config$burn_in must be a non-negative integer.")
  }
  
  config
}

#' Create TVCQR Methods List
#' 
#' @param Mm.factor_vec Vector of Mm.factor values to test
#' @return Named list of method functions
#' @export
create_tvcqr_methods <- function(Mm.factor_vec) {
  methods <- list()
  include_ppro <- tolower(Sys.getenv("FASTQR_INCLUDE_TVCQR_PPRO", unset = "0")) %in%
    c("1", "true", "t", "yes", "y")
  
  # Methods without Mm.factor parameter
  methods$tvc_rq <- function(x, y, config) {
    tvc_rq(x = x, y = y, tau = config$tau, h = config$h)
  }
  
  methods$tvcqr_seq <- function(x, y, config) {
    tvcqr_seq(x = x, y = y, tau = config$tau, h = config$h, 
              h.factor = config$h.factor, tol = config$tol, 
              maxit = config$maxit, bland = config$bland)
  }
  
  methods$tvcqr_seq_fortran <- function(x, y, config) {
    tvcqr_seq_fortran_wrapper(x = x, y = y, tau = config$tau, 
                              h = config$h, tol = config$tol, 
                              maxit = config$maxit, bland = config$bland)
  }
  
  # Methods with Mm.factor parameter - create multiple versions
  for (i in seq_along(Mm.factor_vec)) {
    Mm_val <- Mm.factor_vec[i]
    
    # tvc_rq_ppro is kept opt-in only until the R ppro path is stable.
    if (include_ppro) {
      method_name_ppro <- sprintf("tvc_rq_ppro_%d", i)
      methods[[method_name_ppro]] <- local({
        Mm_factor_local <- Mm_val
        function(x, y, config) {
          tvc_rq_ppro(x = x, y = y, tau = config$tau, h = config$h, 
                      Mm.factor = Mm_factor_local, pmethod = NULL)
        }
      })
    }
    
    # tvcqr_seq_ppro with different Mm.factor values
    method_name <- sprintf("tvcqr_seq_ppro_%d", i)
    methods[[method_name]] <- local({
      Mm_factor_local <- Mm_val
      function(x, y, config) {
        tvcqr_seq_ppro(x = x, y = y, tau = config$tau, h = config$h, 
                       h.factor = config$h.factor, tol = config$tol, 
                       maxit = config$maxit, bland = config$bland,
                       Mm.factor = Mm_factor_local, eps = config$eps,
                       cpp_helper = config$cpp_helper)
      }
    })
    
    # tvcqr_seq_ppro_fortran with different Mm.factor values
    method_name_fortran <- sprintf("tvcqr_seq_ppro_fortran_%d", i)
    methods[[method_name_fortran]] <- local({
      Mm_factor_local <- Mm_val
      function(x, y, config) {
        tvcqr_seq_ppro_fortran_wrapper(x = x, y = y, tau = config$tau, 
                                       h = config$h, h.factor = config$h.factor,
                                       tol = config$tol, maxit = config$maxit, 
                                       bland = config$bland, 
                                       Mm.factor = Mm_factor_local, 
                                       eps = config$eps,
                                       store_residual = FALSE)
      }
    })
  }
  
  return(methods)
}

#' Run Single TVCQR Replication
#' 
#' @param rep_id Replication ID
#' @param config Configuration list
#' @param methods Named list of method functions
#' @return List containing estimates, H_seq, and timing results
#' @export
run_single_tvcqr_replication <- function(rep_id, config, methods) {
  config <- complete_tvcqr_sim_config(config)
  seed_used <- if (!is.null(config$seed_used)) {
    as.integer(config$seed_used)
  } else {
    as.integer(config$seed_base + rep_id)
  }
  
  # Generate data
  data <- generate_ts(
    n = config$n,
    case = config$case,
    seed = seed_used,
    J = config$J,
    burn_in = config$burn_in
  )
  x <- data$x
  y <- data$y
  
  # Initialize storage for this replication
  results <- list(
    estimates = vector("list", length(methods)),
    H_seq = vector("list", length(methods)),
    timing = numeric(length(methods))
  )
  names(results$estimates) <- names(methods)
  names(results$H_seq) <- names(methods)
  names(results$timing) <- names(methods)
  
  # Run each method
  for (method_name in names(methods)) {
    # Benchmark the method
    timing_result <- tryCatch(
      microbenchmark::microbenchmark(
        {
          fit <- methods[[method_name]](x, y, config)
        },
        times = 1,
        unit = "s"
      ),
      error = function(e) {
        prefix <- if (identical(method_name, "tvc_rq")) "baseline_error:tvc_rq" else sprintf("method_error:%s", method_name)
        stop(sprintf("%s:%s", prefix, conditionMessage(e)), call. = FALSE)
      }
    )
    
    # Store results
    results$estimates[[method_name]] <- fit$theta_ll_est
    
    # Store H_seq if it exists in the fit object
    if (!is.null(fit$H_seq)) {
      results$H_seq[[method_name]] <- fit$H_seq
    } else {
      results$H_seq[[method_name]] <- NULL
    }
    
    # Extract timing in seconds
    results$timing[[method_name]] <- summary(timing_result)$mean
  }

  results$seed_used <- seed_used
  
  return(results)
}

#' Run TVCQR Simulation
#' 
#' @param config Configuration list containing simulation parameters
#' @return List containing simulation results
#' @export
run_tvcqr_simulation <- function(config) {
  config <- complete_tvcqr_sim_config(config)
  
  cat("========================================\n")
  cat("Starting TVCQR Simulation\n")
  cat("========================================\n")
  cat(sprintf("Case: %d (%s)\n", config$case, config$case_label))
  cat(sprintf("Tau: %.2f\n", config$tau))
  cat(sprintf("Sample size: %d\n", config$n))
  cat(sprintf("Number of replications: %d\n", config$num_rep))
  cat(sprintf("Seed base: %d\n", config$seed_base))
  if (config$case == 2L) {
    cat(sprintf("Case 2 truncation J: %d\n", config$J))
    cat(sprintf("Case 2 burn-in: %d\n", config$burn_in))
  }
  cat(sprintf("Mm.factor values: %s\n", paste(config$Mm.factor, collapse = ", ")))
  cat("========================================\n\n")
  
  # Create methods list with multiple Mm.factor values
  tvcqr_methods <- create_tvcqr_methods(config$Mm.factor)
  
  # Initialize storage for all replications
  num_methods <- length(tvcqr_methods)
  method_names <- names(tvcqr_methods)
  
  cat(sprintf("Total number of methods to compare: %d\n", num_methods))
  cat("Methods:\n")
  for (i in seq_along(method_names)) {
    cat(sprintf("  %d. %s\n", i, method_names[i]))
  }
  cat("\n")
  
  # Timing matrix: num_rep x num_methods
  timing_matrix <- matrix(NA, nrow = config$num_rep, ncol = num_methods)
  colnames(timing_matrix) <- method_names
  
  # Estimates list: one list per method, each containing num_rep matrices
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
    rep_results <- run_single_tvcqr_replication(rep, config, tvcqr_methods)
    
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
    Mm.factor_mapping = create_tvcqr_Mm_factor_mapping(method_names, config$Mm.factor),
    timestamp = Sys.time(),
    case_label = config$case_label,
    case_key = config$case_key,
    simulation_type = "tvcqr"  # Add identifier
  )
  
  return(simulation_results)
}

#' Create Mm.factor Mapping for TVCQR Methods
#' 
#' @param method_names Vector of method names
#' @param Mm.factor_vec Vector of Mm.factor values
#' @return Data frame mapping methods to Mm.factor values
#' @export
create_tvcqr_Mm_factor_mapping <- function(method_names, Mm.factor_vec) {
  mapping <- data.frame(
    Method = method_names,
    Mm.factor = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(method_names)) {
    method <- method_names[i]
    # Match both R and Fortran ppro variants, e.g.
    # tvcqr_seq_ppro_1 and tvcqr_seq_ppro_fortran_1.
    if (grepl("ppro.*_(\\d+)$", method)) {
      idx <- as.numeric(sub(".*_(\\d+)$", "\\1", method))
      if (!is.na(idx) && idx >= 1L && idx <= length(Mm.factor_vec)) {
        mapping$Mm.factor[i] <- as.character(Mm.factor_vec[idx])
      } else {
        mapping$Mm.factor[i] <- "N/A"
      }
    } else {
      mapping$Mm.factor[i] <- "N/A"
    }
  }
  
  return(mapping)
}

#' Save TVCQR Simulation Results
#' 
#' @param results Simulation results object
#' @param output_dir Directory to save results
#' @return File path where results were saved
#' @export
save_tvcqr_simulation_results <- function(results, output_dir = "data/tvcqr_simu_results") {
  results$config <- complete_tvcqr_sim_config(results$config)
  results$case_label <- results$config$case_label
  results$case_key <- results$config$case_key
  
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

#' Print TVCQR Simulation Summary
#' 
#' @param results Simulation results object
#' @export
print_tvcqr_simulation_summary <- function(results) {
  config <- complete_tvcqr_sim_config(results$config)
  
  cat("\n========================================\n")
  cat("TVCQR Simulation Summary\n")
  cat("========================================\n\n")
  cat(sprintf("Case: %d (%s)\n", config$case, config$case_label))
  cat(sprintf("Tau: %.2f\n", config$tau))
  cat(sprintf("Sample size: %d\n", config$n))
  cat(sprintf("Replications: %d\n\n", config$num_rep))
  
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
  cat("\nMaximum Average Relative Bias (relative to tvc_rq):\n")
  cat("----------------------------------------\n")
  
  max_avg_rel_bias <- compute_tvcqr_max_average_relative_bias(results)
  
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

#' Compute Maximum Average Relative Bias for TVCQR
#' 
#' @param results Simulation results object
#' @return Named vector of maximum average relative bias values
#' @export
compute_tvcqr_max_average_relative_bias <- function(results) {
  estimates_list <- results$estimates_list
  method_names <- results$method_names
  num_rep <- results$config$num_rep
  
  # Get tvc_rq estimates as baseline
  if (!"tvc_rq" %in% method_names) {
    stop("tvc_rq method not found in results!")
  }
  
  baseline_estimates <- estimates_list[["tvc_rq"]]
  
  # Initialize vector to store maximum average relative bias
  max_avg_rel_bias <- numeric(length(method_names))
  names(max_avg_rel_bias) <- method_names
  
  # For tvc_rq itself, the bias is 0
  max_avg_rel_bias["tvc_rq"] <- 0
  
  # Calculate for other methods
  for (method_name in setdiff(method_names, "tvc_rq")) {
    method_estimates <- estimates_list[[method_name]]
    
    # Vector to store average relative bias for each replication
    avg_rel_bias <- numeric(num_rep)
    
    for (rep in 1:num_rep) {
      baseline_mat <- baseline_estimates[[rep]]
      method_mat <- method_estimates[[rep]]
      
      # Check if matrices have the same dimensions
      if (!all(dim(baseline_mat) == dim(method_mat))) {
        warning(sprintf("Dimension mismatch for %s at replication %d", 
                        method_name, rep))
        avg_rel_bias[rep] <- NA
        next
      }
      
      # Calculate element-wise relative bias
      abs_diff <- abs(method_mat - baseline_mat)
      abs_baseline <- abs(baseline_mat)
      
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

#' Find TVCQR Replication with Maximum Bias
#' 
#' @param results Simulation results object
#' @param method_name Name of the method to analyze
#' @return List containing replication ID and bias information
#' @export
find_tvcqr_max_bias_replication <- function(results, method_name) {
  estimates_list <- results$estimates_list
  num_rep <- results$config$num_rep
  
  # Get baseline estimates
  if (!"tvc_rq" %in% names(estimates_list)) {
    stop("tvc_rq method not found in results!")
  }
  
  baseline_estimates <- estimates_list[["tvc_rq"]]
  method_estimates <- estimates_list[[method_name]]
  
  # Vector to store average relative bias for each replication
  avg_rel_bias <- numeric(num_rep)
  
  for (rep in 1:num_rep) {
    baseline_mat <- baseline_estimates[[rep]]
    
    if (rep > length(method_estimates) || is.null(method_estimates[[rep]])) {
      avg_rel_bias[rep] <- NA
      next
    }
    
    method_mat <- method_estimates[[rep]]
    
    # Check dimensions
    if (!all(dim(baseline_mat) == dim(method_mat))) {
      avg_rel_bias[rep] <- NA
      next
    }
    
    # Calculate element-wise relative bias
    abs_diff <- abs(method_mat - baseline_mat)
    abs_baseline <- abs(baseline_mat)
    
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

#' Extract TVCQR Data for Specific Replication
#' 
#' @param results Simulation results object
#' @param rep_id Replication ID
#' @return List containing x and y data
#' @export
extract_tvcqr_replication_data <- function(results, rep_id) {
  config <- complete_tvcqr_sim_config(results$config)
  
  # Regenerate the same data using the stored seed
  data <- generate_ts(
    n = config$n, 
    case = config$case, 
    seed = config$seed_base + rep_id,
    J = config$J,
    burn_in = config$burn_in
  )
  
  return(data)
}

#' Compare TVCQR Estimates for Specific Replication
#' 
#' @param results Simulation results object
#' @param rep_id Replication ID
#' @param methods_to_compare Vector of method names to compare (optional)
#' @export
compare_tvcqr_replication_estimates <- function(results, rep_id, methods_to_compare = NULL) {
  estimates_list <- results$estimates_list
  method_names <- results$method_names
  
  # If no specific methods provided, use all methods
  if (is.null(methods_to_compare)) {
    methods_to_compare <- method_names
  }
  
  cat("\n========================================\n")
  cat(sprintf("TVCQR Estimate Comparison for Replication %d\n", rep_id))
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
    cat(sprintf("\n%s:\n", method))
    cat(sprintf("  Dimensions: %d x %d\n", nrow(est), ncol(est)))
    cat(sprintf("  Mean: %.6f\n", mean(est, na.rm = TRUE)))
    cat(sprintf("  Median: %.6f\n", median(est, na.rm = TRUE)))
    cat(sprintf("  SD: %.6f\n", sd(est, na.rm = TRUE)))
    cat(sprintf("  Min: %.6f\n", min(est, na.rm = TRUE)))
    cat(sprintf("  Max: %.6f\n", max(est, na.rm = TRUE)))
    cat(sprintf("  NA count: %d\n", sum(is.na(est))))
  }
  
  # Calculate pairwise differences with tvc_rq
  if ("tvc_rq" %in% names(estimates_comparison)) {
    cat("\n\nDifferences from tvc_rq:\n")
    cat("----------------------------------------\n")
    baseline <- estimates_comparison[["tvc_rq"]]
    
    for (method in setdiff(names(estimates_comparison), "tvc_rq")) {
      method_est <- estimates_comparison[[method]]
      
      if (all(dim(baseline) == dim(method_est))) {
        diff <- method_est - baseline
        abs_diff <- abs(diff)
        rel_diff <- abs_diff / (abs(baseline) + 1e-10)
        
        cat(sprintf("\n%s vs tvc_rq:\n", method))
        cat(sprintf("  Mean absolute diff: %.6e\n", mean(abs_diff, na.rm = TRUE)))
        cat(sprintf("  Max absolute diff: %.6e\n", max(abs_diff, na.rm = TRUE)))
        cat(sprintf("  Mean relative diff: %.6e\n", mean(rel_diff, na.rm = TRUE)))
        cat(sprintf("  Max relative diff: %.6e\n", max(rel_diff, na.rm = TRUE)))
      } else {
        cat(sprintf("\n%s vs tvc_rq: Dimension mismatch\n", method))
      }
    }
  }
  
  return(estimates_comparison)
}

#' Save Problematic TVCQR Data
#' 
#' @param results Simulation results object
#' @param method_name Name of the method to analyze
#' @param output_dir Directory to save data
#' @return List containing filepath and problematic case data
#' @export
save_tvcqr_problematic_data <- function(results, method_name, 
                                        output_dir = "data/tvcqr_problematic_cases") {
  # Find replication with max bias
  max_bias_info <- find_tvcqr_max_bias_replication(results, method_name)
  rep_id <- max_bias_info$replication_id
  
  # Extract data
  data <- extract_tvcqr_replication_data(results, rep_id)
  
  # Compare estimates
  estimates_comparison <- compare_tvcqr_replication_estimates(results, rep_id)
  
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
    simulation_type = "tvcqr"
  )
  
  # Generate filename
  config <- results$config
  tau_str <- sprintf("tau%02d", round(config$tau * 100))
  filename <- sprintf("tvcqr_problematic_case%d_%s_n%d_rep%d_%s.RData", 
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
