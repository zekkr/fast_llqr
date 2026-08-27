# ============================================================================ #
# Fast LLQR Project - Performance Measurement Functions
# Description: Generic functions for summarizing simulation results
# Author: Erkang
# Date: 2025-12-04
# ============================================================================ #

#' Load simulation results from RData files
#'
#' @param case Integer, case number
#' @param tau Numeric vector, quantile levels
#' @param n Integer vector, sample sizes
#' @param rep Integer, number of replications
#' @param model Character, either "tvcqr" or "llqr"
#' @param data_dir Character, directory containing RData files (optional, will use default based on type)
#' @return A list containing loaded results or NULL for missing files
#' @export
load_simulation_results <- function(case, tau, n, rep, 
                                    model = c("tvcqr", "llqr"),
                                    data_dir = NULL) {
  # Match argument
  model <- match.arg(model)
  
  # Set default data directory if not provided
  if (is.null(data_dir)) {
    data_dir <- switch(model,
                       "tvcqr" = "data/tvcqr_simu_results",
                       "llqr" = "data/llqr_simu_results")
  }
  
  cat(sprintf("Loading %s simulation results from: %s\n", 
              toupper(model), data_dir))
  cat(strrep("-", 80), "\n")
  
  results_list <- list()
  case_label_default <- sprintf("Case %s", as.character(case))
  if (model == "tvcqr" && exists("resolve_tvcqr_case", mode = "function")) {
    case_label_default <- resolve_tvcqr_case(case)$case_label
  }
  
  for (t in tau) {
    for (sample_size in n) {
      # Convert tau to string format
      tau_str <- sprintf("tau%02d", as.integer(t * 100))
      
      # Construct filename based on model
      filename <- sprintf("case%d_%s_n%d_rep%d.RData", 
                          case, tau_str, sample_size, rep)
      filepath <- file.path(data_dir, filename)
      
      # Create unique key for this configuration
      key <- sprintf("case%d_tau%.1f_n%d", case, t, sample_size)
      
      # Try to load the file
      if (file.exists(filepath)) {
        tryCatch({
          # Load RData file
          env <- new.env()
          load(filepath, envir = env)
          
          # Extract the results object
          if (exists("results", envir = env)) {
            results_obj <- get("results", envir = env)
            
            # Verify model matches
            if (!is.null(results_obj$model)) {
              if (results_obj$model != model) {
                warning(sprintf("File %s has model '%s' but expecting '%s'",
                                filename, results_obj$model, model))
              }
            }
            
            results_list[[key]] <- list(
              case = case,
              case_label = if (!is.null(results_obj$config$case_label)) {
                results_obj$config$case_label
              } else {
                case_label_default
              },
              tau = t,
              n = sample_size,
              rep = rep,
              model = model,
              data = results_obj,
              status = "success"
            )
            cat(sprintf("✓ Loaded: %s\n", filename))
          } else {
            results_list[[key]] <- list(
              case = case,
              tau = t,
              n = sample_size,
              rep = rep,
              model = model,
              data = NULL,
              status = "error",
              error_msg = "No 'results' object found in file"
            )
            cat(sprintf("✗ No 'results' object in %s\n", filename))
          }
        }, error = function(e) {
          results_list[[key]] <- list(
            case = case,
            case_label = case_label_default,
            tau = t,
            n = sample_size,
            rep = rep,
            model = model,
            data = NULL,
            status = "error",
            error_msg = as.character(e)
          )
          cat(sprintf("✗ Error loading %s: %s\n", filename, e$message))
        })
      } else {
        results_list[[key]] <- list(
          case = case,
          case_label = case_label_default,
          tau = t,
          n = sample_size,
          rep = rep,
          model = model,
          data = NULL,
          status = "missing"
        )
        cat(sprintf("⚠ Missing: %s\n", filename))
      }
    }
  }
  
  cat(strrep("-", 80), "\n")
  cat(sprintf("Total configurations loaded: %d\n", length(results_list)))
  successful <- sum(sapply(results_list, function(x) x$status == "success"))
  cat(sprintf("Successful loads: %d\n", successful))
  cat(sprintf("Failed/Missing: %d\n", length(results_list) - successful))
  cat("\n")
  
  return(results_list)
}


#' Extract summary statistics from a numeric vector
#'
#' @param x Numeric vector
#' @param probs Numeric vector of probabilities for quantiles
#' @param na.rm Logical, whether to remove NA values
#' @return Named vector of summary statistics
#' @export
compute_summary_stats <- function(x, 
                                  probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                                  na.rm = TRUE) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    n_stats <- 5 + length(probs)
    quantile_names <- paste0("Q", sprintf("%02d", probs * 100))
    return(setNames(rep(NA_real_, n_stats), 
                    c("mean", "median", "sd", "min", "max", quantile_names)))
  }
  
  quantile_names <- paste0("Q", sprintf("%02d", probs * 100))
  quantile_values <- quantile(x, probs = probs, na.rm = na.rm)
  names(quantile_values) <- quantile_names
  
  stats <- c(
    mean = mean(x, na.rm = na.rm),
    median = median(x, na.rm = na.rm),
    sd = sd(x, na.rm = na.rm),
    min = min(x, na.rm = na.rm),
    max = max(x, na.rm = na.rm),
    quantile_values
  )
  
  return(stats)
}


#' Get all unique method names across all results
#'
#' @param results_list List of loaded simulation results
#' @return Character vector of unique method names
#' @export
get_all_method_names <- function(results_list) {
  all_methods <- character()
  
  for (result in results_list) {
    if (result$status == "success" && !is.null(result$data)) {
      # Get method names from timing_matrix or method_names
      if (!is.null(result$data$method_names)) {
        methods <- result$data$method_names
      } else if (!is.null(result$data$timing_matrix)) {
        methods <- colnames(result$data$timing_matrix)
      } else if (!is.null(result$data$estimates_list)) {
        # Try to get from estimates_list
        if (length(result$data$estimates_list) > 0) {
          first_est <- result$data$estimates_list[[1]]
          if (is.list(first_est)) {
            methods <- names(first_est)
          }
        }
      } else {
        methods <- NULL
      }
      
      if (!is.null(methods)) {
        all_methods <- union(all_methods, methods)
      }
    }
  }
  
  return(all_methods)
}


#' Extract computation time for a specific method
#'
#' @param results_obj The results object from loaded RData
#' @param method_name Character, name of the method
#' @return Numeric vector of computation times or NULL
#' @export
extract_computation_time <- function(results_obj, method_name) {
  if (is.null(results_obj) || is.null(method_name)) return(NULL)
  
  # Check if timing_matrix exists
  if (!is.null(results_obj$timing_matrix)) {
    timing <- results_obj$timing_matrix
    
    # Extract specific method's timing
    if (method_name %in% colnames(timing)) {
      return(as.numeric(timing[, method_name]))
    }
  }
  
  return(NULL)
}


#' Extract iteration numbers for a specific method
#'
#' @param results_obj The results object from loaded RData
#' @param method_name Character, name of the method
#' @return Numeric vector of iteration numbers or NULL
#' @export
extract_iteration_numbers <- function(results_obj, method_name) {
  if (is.null(results_obj) || is.null(method_name)) return(NULL)
  
  # Check if it_num_matrix exists (for sequential methods)
  if (!is.null(results_obj$it_num_matrix)) {
    it_num <- results_obj$it_num_matrix
    
    # Extract specific method's iteration numbers
    if (method_name %in% colnames(it_num)) {
      return(as.numeric(it_num[, method_name]))
    }
  }
  
  return(NULL)
}


#' Calculate average relative bias for TVCQR methods
#'
#' @param results_obj The results object from loaded RData
#' @param method_name Character, name of the method
#' @param reference_method Character, name of reference method (default: "tvc_rq")
#' @return Numeric vector of average relative bias values (one per replication)
#' @export
calculate_average_relative_bias_tvcqr <- function(results_obj, 
                                                  method_name,
                                                  reference_method = "tvc_rq") {
  if (is.null(results_obj) || is.null(method_name)) return(NULL)
  
  # Check if estimates_list exists
  if (is.null(results_obj$estimates_list)) {
    cat(sprintf("  Warning: estimates_list not found for method %s\n", method_name))
    return(NULL)
  }
  
  estimates_list <- results_obj$estimates_list
  
  # Check if both methods exist in estimates_list
  if (!reference_method %in% names(estimates_list)) {
    cat(sprintf("  Warning: Reference method '%s' not found in estimates_list\n", 
                reference_method))
    return(NULL)
  }
  
  if (!method_name %in% names(estimates_list)) {
    cat(sprintf("  Warning: Method '%s' not found in estimates_list\n", method_name))
    return(NULL)
  }
  
  # Get estimates for both methods
  reference_estimates <- estimates_list[[reference_method]]
  method_estimates <- estimates_list[[method_name]]
  
  # Check that both have the same number of replications
  n_rep <- length(reference_estimates)
  if (length(method_estimates) != n_rep) {
    cat(sprintf("  Warning: Number of replications mismatch for method %s\n", method_name))
    return(NULL)
  }
  
  # Initialize vector to store average relative bias for each replication
  avg_rel_bias <- numeric(n_rep)
  expected_size <- suppressWarnings(as.integer(results_obj$config$n))
  if (length(expected_size) != 1L || is.na(expected_size) || expected_size <= 0L) {
    expected_size <- NULL
  } else {
    expected_size <- 4L * expected_size
  }
  
  # Loop through each replication
  for (i in 1:n_rep) {
    # Extract theta_ll_est for reference and current method
    reference_theta <- reference_estimates[[i]]
    method_theta <- method_estimates[[i]]
    
    # Check if both exist
    if (is.null(reference_theta) || is.null(method_theta)) {
      avg_rel_bias[i] <- NA
      next
    }
    
    # Convert to matrix if needed
    if (!is.matrix(reference_theta)) {
      reference_theta <- as.matrix(reference_theta)
    }
    if (!is.matrix(method_theta)) {
      method_theta <- as.matrix(method_theta)
    }
    
    # Check dimensions match
    if (!all(dim(reference_theta) == dim(method_theta))) {
      cat(sprintf("  Warning: Dimension mismatch in replication %d for method %s\n", 
                  i, method_name))
      cat(sprintf("    Reference: %s, Method: %s\n", 
                  paste(dim(reference_theta), collapse = "x"),
                  paste(dim(method_theta), collapse = "x")))
      avg_rel_bias[i] <- NA
      next
    }

    if (!is.null(expected_size) && length(reference_theta) != expected_size) {
      cat(sprintf(
        "  Warning: Unexpected TVCQR estimate size in replication %d for method %s: expected %d, found %d\n",
        i, method_name, expected_size, length(reference_theta)
      ))
      avg_rel_bias[i] <- NA_real_
      next
    }

    if (!is.numeric(reference_theta) || !is.numeric(method_theta) ||
        any(!is.finite(reference_theta)) || any(!is.finite(method_theta))) {
      cat(sprintf("  Warning: Non-finite TVCQR estimates in replication %d for method %s\n",
                  i, method_name))
      avg_rel_bias[i] <- NA_real_
      next
    }
    
    # Calculate d_rk = |a_rk - b_rk| / max(|b_rk|, 1e-10).
    abs_diff <- abs(method_theta - reference_theta)
    abs_reference <- abs(reference_theta)

    rel_bias <- abs_diff / pmax(abs_reference, 1e-10)
    
    # Calculate average across all elements
    avg_rel_bias[i] <- mean(rel_bias)
  }
  
  return(avg_rel_bias)
}


#' Calculate average relative bias for LLQR methods
#'
#' @param results_obj The results object from loaded RData
#' @param method_name Character, name of the method
#' @param reference_method Character, name of reference method (default: "llqr")
#' @return Numeric vector of average relative bias values (one per replication)
#' @export
calculate_average_relative_bias_llqr <- function(results_obj, 
                                                 method_name,
                                                 reference_method = "llqr") {
  if (is.null(results_obj) || is.null(method_name)) return(NULL)
  
  # Check if estimates_list exists
  if (is.null(results_obj$estimates_list)) {
    cat(sprintf("  Warning: estimates_list not found for method %s\n", method_name))
    return(NULL)
  }
  
  estimates_list <- results_obj$estimates_list
  
  # Check if both methods exist in estimates_list
  if (!reference_method %in% names(estimates_list)) {
    cat(sprintf("  Warning: Reference method '%s' not found in estimates_list\n", 
                reference_method))
    return(NULL)
  }
  
  if (!method_name %in% names(estimates_list)) {
    cat(sprintf("  Warning: Method '%s' not found in estimates_list\n", method_name))
    return(NULL)
  }
  
  # Get estimates for both methods
  reference_estimates <- estimates_list[[reference_method]]
  method_estimates <- estimates_list[[method_name]]
  
  # Check that both have the same number of replications
  n_rep <- length(reference_estimates)
  if (length(method_estimates) != n_rep) {
    cat(sprintf("  Warning: Number of replications mismatch for method %s\n", method_name))
    return(NULL)
  }
  
  # Initialize vector to store average relative bias for each replication
  avg_rel_bias <- numeric(n_rep)
  expected_size <- suppressWarnings(as.integer(results_obj$config$n))
  if (length(expected_size) != 1L || is.na(expected_size) || expected_size <= 0L) {
    expected_size <- NULL
  }
  
  # Loop through each replication
  for (i in 1:n_rep) {
    # Extract ll_est for reference and current method
    reference_ll <- reference_estimates[[i]]
    method_ll <- method_estimates[[i]]
    
    # Check if both exist
    if (is.null(reference_ll) || is.null(method_ll)) {
      avg_rel_bias[i] <- NA
      next
    }
    
    # Convert to numeric vectors
    reference_ll <- as.numeric(reference_ll)
    method_ll <- as.numeric(method_ll)
    
    # Check lengths match
    if (length(reference_ll) != length(method_ll)) {
      cat(sprintf("  Warning: Length mismatch in replication %d for method %s\n", 
                  i, method_name))
      cat(sprintf("    Reference: %d, Method: %d\n", 
                  length(reference_ll), length(method_ll)))
      avg_rel_bias[i] <- NA
      next
    }

    if (!is.null(expected_size) && length(reference_ll) != expected_size) {
      cat(sprintf(
        "  Warning: Unexpected LLQR estimate size in replication %d for method %s: expected %d, found %d\n",
        i, method_name, expected_size, length(reference_ll)
      ))
      avg_rel_bias[i] <- NA_real_
      next
    }

    if (any(!is.finite(reference_ll)) || any(!is.finite(method_ll))) {
      cat(sprintf("  Warning: Non-finite LLQR estimates in replication %d for method %s\n",
                  i, method_name))
      avg_rel_bias[i] <- NA_real_
      next
    }
    
    # Calculate d_rk = |a_rk - b_rk| / max(|b_rk|, 1e-10).
    abs_diff <- abs(method_ll - reference_ll)
    abs_reference <- abs(reference_ll)

    rel_bias <- abs_diff / pmax(abs_reference, 1e-10)
    
    # Calculate average across all elements
    avg_rel_bias[i] <- mean(rel_bias)
  }
  
  return(avg_rel_bias)
}



#' Extract metric value based on model and metric name
#'
#' @param results_obj The results object from loaded RData
#' @param method_name Character, name of the method
#' @param metric_name Character, name of the metric
#' @param model Character, either "tvcqr" or "llqr"
#' @return Numeric vector of metric values
#' @export
extract_metric_value <- function(results_obj, method_name, metric_name, model) {
  
  if (metric_name == "computation_time") {
    return(extract_computation_time(results_obj, method_name))
    
  } else if (metric_name == "iteration_numbers") {
    return(extract_iteration_numbers(results_obj, method_name))
    
  } else if (metric_name == "average_relative_bias") {
    if (model == "tvcqr") {
      return(calculate_average_relative_bias_tvcqr(results_obj, method_name))
    } else if (model == "llqr") {
      return(calculate_average_relative_bias_llqr(results_obj, method_name))
    } else {
      warning(sprintf("Unknown model: %s", model))
      return(NULL)
    }
    
  } else {
    warning(sprintf("Unknown metric: %s", metric_name))
    return(NULL)
  }
}


#' Create summary table for a specific metric (all methods)
#'
#' @param results_list List of loaded simulation results
#' @param metric_name Character, name of the metric
#' @param probs Numeric vector of probabilities for quantiles
#' @return Data frame with summary statistics for all methods
#' @export
create_summary_table <- function(results_list, 
                                 metric_name,
                                 probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  
  # Get all unique method names
  all_methods <- get_all_method_names(results_list)
  
  if (length(all_methods) == 0) {
    cat(sprintf("Warning: No methods found in the data\n"))
    return(data.frame())
  }
  
  cat(sprintf("  Methods found: %s\n", paste(all_methods, collapse = ", ")))
  
  # Get model from first successful result
  model <- NULL
  for (result in results_list) {
    if (result$status == "success") {
      model <- result$model
      break
    }
  }
  
  if (is.null(model)) {
    stop("Could not determine model from results")
  }
  
  cat(sprintf("  Model: %s\n", model))
  
  # Initialize list to store all rows
  all_rows <- list()
  
  # Iterate through each result
  for (key in names(results_list)) {
    result <- results_list[[key]]
    
    # Process each method
    for (method in all_methods) {
      
      # Extract metric value
      metric_value <- NULL
      if (result$status == "success" && !is.null(result$data)) {
        tryCatch({
          metric_value <- extract_metric_value(
            results_obj = result$data,
            method_name = method,
            metric_name = metric_name,
            model = model
          )
        }, error = function(e) {
          cat(sprintf("  Warning: Error extracting %s for %s in %s: %s\n", 
                      metric_name, method, key, e$message))
        })
      }
      
      # Compute summary statistics
      if (!is.null(metric_value) && length(metric_value) > 0 && !all(is.na(metric_value))) {
        stats <- compute_summary_stats(metric_value, probs = probs)
      } else {
        # Create NA statistics
        n_stats <- 5 + length(probs)
        quantile_names <- paste0("Q", sprintf("%02d", probs * 100))
        stats <- setNames(rep(NA_real_, n_stats), 
                          c("mean", "median", "sd", "min", "max", quantile_names))
      }
      
      # Create row data frame
      row_df <- data.frame(
        case = result$case,
        case_label = if (!is.null(result$case_label)) result$case_label else sprintf("Case %s", result$case),
        tau = result$tau,
        n = result$n,
        rep = result$rep,
        model = result$model,
        method = method,
        metric = metric_name,
        stringsAsFactors = FALSE
      )
      
      # Add statistics columns
      for (stat_name in names(stats)) {
        row_df[[stat_name]] <- stats[stat_name]
      }
      
      all_rows[[length(all_rows) + 1]] <- row_df
    }
  }
  
  # Combine all rows
  if (length(all_rows) == 0) {
    return(data.frame())
  }
  
  summary_df <- do.call(rbind, all_rows)
  rownames(summary_df) <- NULL
  
  return(summary_df)
}


#' Create summary tables for multiple metrics
#'
#' @param results_list List of loaded simulation results
#' @param metrics Character vector of metric names
#' @param probs Numeric vector of probabilities for quantiles
#' @return List of data frames, one for each metric
#' @export
create_multi_metric_summary <- function(results_list,
                                        metrics = c("computation_time", 
                                                    "average_relative_bias"),
                                        probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  
  summary_tables <- lapply(metrics, function(metric_name) {
    cat(sprintf("\nProcessing metric: %s\n", metric_name))
    cat(strrep("-", 60), "\n")
    
    table <- create_summary_table(
      results_list = results_list,
      metric_name = metric_name,
      probs = probs
    )
    
    cat(sprintf("  Rows generated: %d\n", nrow(table)))
    
    return(table)
  })
  
  names(summary_tables) <- metrics
  return(summary_tables)
}


#' Save summary tables to files
#'
#' @param summary_tables List of summary tables
#' @param output_dir Character, output directory
#' @param filename_prefix Character, prefix for output files
#' @param format Character, either "csv" or "xlsx"
#' @export
save_summary_tables <- function(summary_tables,
                                output_dir = "results/table",
                                filename_prefix = "simulation_result",
                                format = "xlsx") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Created directory: %s\n", output_dir))
  }
  
  if (format == "xlsx") {
    # Save all tables to a single Excel file with multiple sheets
    output_file <- file.path(output_dir, paste0(filename_prefix, ".xlsx"))
    writexl::write_xlsx(summary_tables, path = output_file)
    cat(sprintf("\n✓ Saved Excel file: %s\n", output_file))
    cat(sprintf("  Sheets: %s\n", paste(names(summary_tables), collapse = ", ")))
    
  } else if (format == "csv") {
    # Save each table as a separate CSV file
    for (metric_name in names(summary_tables)) {
      output_file <- file.path(output_dir, 
                               sprintf("%s_%s.csv", filename_prefix, metric_name))
      readr::write_csv(summary_tables[[metric_name]], file = output_file)
      cat(sprintf("✓ Saved CSV file: %s\n", output_file))
    }
  } else {
    stop("Format must be either 'csv' or 'xlsx'")
  }
  
  invisible(TRUE)
}


#' Main function to run complete performance analysis
#'
#' @param case Integer, case number
#' @param tau Numeric vector, quantile levels
#' @param n Integer vector, sample sizes
#' @param rep Integer, number of replications
#' @param model Character, either "tvcqr" or "llqr"
#' @param metrics Character vector of metrics to analyze
#' @param probs Numeric vector of probabilities for quantiles
#' @param output_format Character, output file format ("csv" or "xlsx")
#' @param output_dir Character, output directory
#' @param filename_prefix Character, prefix for output files (optional, will use default based on type)
#' @return Invisibly returns summary tables
#' @export
run_performance_analysis <- function(case = 1,
                                     tau = c(0.2, 0.5, 0.8),
                                     n = c(200, 500, 1000, 2000),
                                     rep = 500,
                                     model = c("tvcqr", "llqr"),
                                     metrics = c("computation_time", 
                                                 "average_relative_bias"),
                                     probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                                     output_format = "xlsx",
                                     output_dir = "results/table",
                                     filename_prefix = NULL) {
  
  # Match argument
  model <- match.arg(model)
  
  # Set default filename prefix if not provided
  if (is.null(filename_prefix)) {
    filename_prefix <- paste0(model, "_result")
  }
  
  cat(strrep("=", 80), "\n")
  cat(sprintf("%s Performance Analysis\n", toupper(model)))
  cat(strrep("=", 80), "\n\n")
  
  cat("Configuration:\n")
  cat(sprintf("  Model: %s\n", model))
  case_label <- sprintf("Case %s", as.character(case))
  if (model == "tvcqr" && exists("resolve_tvcqr_case", mode = "function")) {
    case_label <- resolve_tvcqr_case(case)$case_label
  }
  cat(sprintf("  Case: %s\n", case_label))
  cat(sprintf("  Tau: %s\n", paste(tau, collapse = ", ")))
  cat(sprintf("  Sample sizes: %s\n", paste(n, collapse = ", ")))
  cat(sprintf("  Replications: %d\n", rep))
  cat(sprintf("  Metrics: %s\n", paste(metrics, collapse = ", ")))
  cat(sprintf("  Quantile probs: %s\n", paste(probs, collapse = ", ")))
  cat("\n")
  
  # Step 1: Load simulation results
  cat("Step 1: Loading simulation results...\n")
  cat(strrep("-", 80), "\n")
  results_list <- load_simulation_results(
    case = case, 
    tau = tau, 
    n = n, 
    rep = rep,
    model = model
  )
  
  # Show methods detected
  all_methods <- get_all_method_names(results_list)
  cat(sprintf("\nMethods detected: %s\n", paste(all_methods, collapse = ", ")))
  
  # Step 2: Create summary tables
  cat("\n\nStep 2: Creating summary tables...\n")
  cat(strrep("-", 80), "\n")
  summary_tables <- create_multi_metric_summary(results_list, metrics, probs)
  cat("\n")
  
  # Step 3: Save results
  cat("\nStep 3: Saving results...\n")
  cat(strrep("-", 80), "\n")
  save_summary_tables(summary_tables, output_dir, filename_prefix, output_format)
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("Analysis complete!\n")
  cat(strrep("=", 80), "\n")
  
  return(invisible(summary_tables))
}
