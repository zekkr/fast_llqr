# ============================================================================ #
# Helper Function: Find Replication with Maximum Bias
# ============================================================================ #
find_max_bias_replication <- function(results, method_name) {
  estimates_list <- results$estimates_list
  num_rep <- results$config$num_rep
  
  # Get llqr estimates as baseline
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
    
    # Check length match
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

# ============================================================================ #
# Helper Function: Extract Data for Specific Replication
# ============================================================================ #
extract_replication_data <- function(results, rep_id) {
  config <- results$config
  
  # Regenerate the same data using the stored seed
  data <- generate_data(
    n = config$n, 
    case = config$case, 
    seed = config$seed_base + rep_id
  )
  
  return(data)
}

# ============================================================================ #
# Helper Function: Compare Estimates for Specific Replication
# ============================================================================ #
compare_replication_estimates <- function(results, rep_id, methods_to_compare = NULL) {
  estimates_list <- results$estimates_list
  method_names <- results$method_names
  
  # If no specific methods provided, use all methods
  if (is.null(methods_to_compare)) {
    methods_to_compare <- method_names
  }
  
  cat("\n========================================\n")
  cat(sprintf("Estimate Comparison for Replication %d\n", rep_id))
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
    cat(sprintf("  Length: %d\n", length(est)))
    cat(sprintf("  Mean: %.6f\n", mean(est, na.rm = TRUE)))
    cat(sprintf("  Median: %.6f\n", median(est, na.rm = TRUE)))
    cat(sprintf("  SD: %.6f\n", sd(est, na.rm = TRUE)))
    cat(sprintf("  Min: %.6f\n", min(est, na.rm = TRUE)))
    cat(sprintf("  Max: %.6f\n", max(est, na.rm = TRUE)))
    cat(sprintf("  NA count: %d\n", sum(is.na(est))))
  }
  
  # Calculate pairwise differences with llqr
  if ("llqr" %in% names(estimates_comparison)) {
    cat("\n\nDifferences from llqr:\n")
    cat("----------------------------------------\n")
    baseline <- estimates_comparison[["llqr"]]
    
    for (method in setdiff(names(estimates_comparison), "llqr")) {
      method_est <- estimates_comparison[[method]]
      
      if (length(baseline) == length(method_est)) {
        diff <- method_est - baseline
        abs_diff <- abs(diff)
        rel_diff <- abs_diff / (abs(baseline) + 1e-10)
        
        cat(sprintf("\n%s vs llqr:\n", method))
        cat(sprintf("  Mean absolute diff: %.6e\n", mean(abs_diff, na.rm = TRUE)))
        cat(sprintf("  Max absolute diff: %.6e\n", max(abs_diff, na.rm = TRUE)))
        cat(sprintf("  Mean relative diff: %.6e\n", mean(rel_diff, na.rm = TRUE)))
        cat(sprintf("  Max relative diff: %.6e\n", max(rel_diff, na.rm = TRUE)))
        
        # Find indices with largest differences
        top_diff_idx <- order(abs_diff, decreasing = TRUE)[1:min(5, length(abs_diff))]
        cat("  Top 5 largest absolute differences (index: baseline, method, diff):\n")
        for (idx in top_diff_idx) {
          cat(sprintf("    %d: %.6f, %.6f, %.6e\n", 
                      idx, baseline[idx], method_est[idx], diff[idx]))
        }
      } else {
        cat(sprintf("\n%s vs llqr: Length mismatch (%d vs %d)\n", 
                    method, length(method_est), length(baseline)))
      }
    }
  }
  
  return(estimates_comparison)
}

# ============================================================================ #
# Helper Function: Save Problematic Data
# ============================================================================ #
save_problematic_data <- function(results, method_name, output_dir = "data/problematic_cases") {
  # Find replication with max bias
  max_bias_info <- find_max_bias_replication(results, method_name)
  rep_id <- max_bias_info$replication_id
  
  # Extract data
  data <- extract_replication_data(results, rep_id)
  
  # Compare estimates
  estimates_comparison <- compare_replication_estimates(results, rep_id)
  
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
    timestamp = Sys.time()
  )
  
  # Generate filename
  config <- results$config
  tau_str <- sprintf("tau%02d", round(config$tau * 100))
  filename <- sprintf("problematic_case%d_%s_n%d_rep%d_%s.RData", 
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

# ============================================================================ #
# Helper Function: Visualize Problematic Case
# ============================================================================ #
visualize_problematic_case <- function(problematic_case, save_plot = TRUE, 
                                       output_dir = "plots") {
  # Extract data
  data <- problematic_case$data
  estimates <- problematic_case$estimates
  
  # Create plots directory if needed
  if (save_plot && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set up plot layout
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  # Plot 1: Scatter plot of data
  plot(data$x, data$y, pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3),
       xlab = "x", ylab = "y", 
       main = sprintf("Data (Replication %d)", problematic_case$replication_id))
  
  # Plot 2: Compare estimates
  if ("llqr" %in% names(estimates)) {
    # Sort by x for better visualization
    sorted_idx <- order(data$x)
    x_sorted <- data$x[sorted_idx]
    
    plot(x_sorted, estimates[["llqr"]][sorted_idx], type = "l", 
         col = "black", lwd = 2,
         xlab = "x", ylab = "Estimate",
         main = "Estimate Comparison")
    
    colors <- rainbow(length(estimates) - 1)
    color_idx <- 1
    for (method in setdiff(names(estimates), "llqr")) {
      lines(x_sorted, estimates[[method]][sorted_idx], 
            col = colors[color_idx], lwd = 2, lty = 2)
      color_idx <- color_idx + 1
    }
    
    legend("topright", legend = names(estimates), 
           col = c("black", colors[1:(length(estimates)-1)]),
           lty = c(1, rep(2, length(estimates)-1)), lwd = 2, cex = 0.7)
  }
  
  # Plot 3: Absolute differences from llqr
  if ("llqr" %in% names(estimates)) {
    baseline <- estimates[["llqr"]]
    max_diff <- 0
    
    for (method in setdiff(names(estimates), "llqr")) {
      method_est <- estimates[[method]]
      if (length(baseline) == length(method_est)) {
        diff <- abs(method_est - baseline)
        max_diff <- max(max_diff, max(diff, na.rm = TRUE))
      }
    }
    
    plot(x_sorted, rep(0, length(x_sorted)), type = "n",
         ylim = c(0, max_diff * 1.1),
         xlab = "x", ylab = "Absolute Difference",
         main = "Absolute Differences from llqr")
    
    color_idx <- 1
    for (method in setdiff(names(estimates), "llqr")) {
      method_est <- estimates[[method]]
      if (length(baseline) == length(method_est)) {
        diff <- abs(method_est - baseline)
        lines(x_sorted, diff[sorted_idx], 
              col = colors[color_idx], lwd = 2)
        color_idx <- color_idx + 1
      }
    }
    
    legend("topright", legend = setdiff(names(estimates), "llqr"),
           col = colors[1:(length(estimates)-1)], lwd = 2, cex = 0.7)
  }
  
  # Plot 4: Relative differences from llqr
  if ("llqr" %in% names(estimates)) {
    baseline <- estimates[["llqr"]]
    max_rel_diff <- 0
    
    for (method in setdiff(names(estimates), "llqr")) {
      method_est <- estimates[[method]]
      if (length(baseline) == length(method_est)) {
        rel_diff <- abs(method_est - baseline) / (abs(baseline) + 1e-10)
        max_rel_diff <- max(max_rel_diff, max(rel_diff, na.rm = TRUE))
      }
    }
    
    plot(x_sorted, rep(0, length(x_sorted)), type = "n",
         ylim = c(0, max_rel_diff * 1.1),
         xlab = "x", ylab = "Relative Difference",
         main = "Relative Differences from llqr")
    
    color_idx <- 1
    for (method in setdiff(names(estimates), "llqr")) {
      method_est <- estimates[[method]]
      if (length(baseline) == length(method_est)) {
        rel_diff <- abs(method_est - baseline) / (abs(baseline) + 1e-10)
        lines(x_sorted, rel_diff[sorted_idx], 
              col = colors[color_idx], lwd = 2)
        color_idx <- color_idx + 1
      }
    }
    
    legend("topright", legend = setdiff(names(estimates), "llqr"),
           col = colors[1:(length(estimates)-1)], lwd = 2, cex = 0.7)
  }
  
  # Save plot if requested
  if (save_plot) {
    config <- problematic_case$config
    tau_str <- sprintf("tau%02d", round(config$tau * 100))
    filename <- sprintf("problematic_case%d_%s_n%d_rep%d_%s.pdf", 
                        config$case, tau_str, config$n, 
                        problematic_case$replication_id,
                        gsub("[^[:alnum:]]", "_", problematic_case$method_name))
    filepath <- file.path(output_dir, filename)
    
    dev.copy(pdf, filepath, width = 10, height = 10)
    dev.off()
    cat(sprintf("\nPlot saved to: %s\n", filepath))
  }
  
  par(mfrow = c(1, 1))
}

# ============================================================================ #
# Execute Simulation and Diagnostic (Add to end of main script)
# ============================================================================ #

# After running the simulation and printing summary:
# results <- run_llqr_simulation(sim_config)
# print_simulation_summary(results)
# save_simulation_results(results)

# Diagnose problematic method
cat("\n========================================\n")
cat("Diagnosing Problematic Methods\n")
cat("========================================\n\n")

# Find methods with high bias (you can adjust the threshold)
max_avg_rel_bias <- compute_max_average_relative_bias(results)
problematic_methods <- names(max_avg_rel_bias[max_avg_rel_bias > 0.01])  # 1% threshold

if (length(problematic_methods) > 0) {
  cat("Methods with bias > 1%:\n")
  print(max_avg_rel_bias[problematic_methods])
  cat("\n")
  
  # Analyze each problematic method
  for (method in problematic_methods) {
    cat(sprintf("\n--- Analyzing %s ---\n", method))
    
    # Save problematic data
    saved_data <- save_problematic_data(results, method)
    
    # Visualize
    visualize_problematic_case(saved_data$problematic_case)
  }
} else {
  cat("No methods with bias > 1% found.\n")
}

# 运行模拟后
results <- run_llqr_simulation(sim_config)

# 诊断特定方法（例如 llqr_seq_ppro_fortran_1）
saved_data <- save_problematic_data(results, "llqr_seq_ppro_fortran_1")

# 查看可视化
visualize_problematic_case(saved_data$problematic_case)

# 或者自动检测所有有问题的方法
max_avg_rel_bias <- compute_max_average_relative_bias(results)
problematic_methods <- names(max_avg_rel_bias[max_avg_rel_bias > 0.01])
for (method in problematic_methods) {
  saved_data <- save_problematic_data(results, method)
  visualize_problematic_case(saved_data$problematic_case)
}
