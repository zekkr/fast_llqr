# ============================================================================ #
# Diagnostic script: Check if memory continuously grows
# ============================================================================ #
diagnose_memory_leak <- function(config) {
  config$num_rep <- 100  # Test 100 times
  
  tvcqr_methods <- create_tvcqr_methods(config$Mm.factor)
  
  # Test only one method
  test_method <- tvcqr_methods[[1]]
  
  memory_trace <- numeric(config$num_rep)
  
  cat("Testing for memory leak...\n")
  pb <- txtProgressBar(min = 0, max = config$num_rep, style = 3)
  
  for (rep in 1:config$num_rep) {
    # Generate data
    data <- generate_ts(n = config$n, case = config$case, 
                        seed = config$seed_base + rep)
    x <- data$x
    y <- data$y
    
    # Run method
    fit <- test_method(x, y, config)
    
    # Record memory usage
    if (.Platform$OS.type == "windows") {
      memory_trace[rep] <- memory.size()
    } else {
      gc_info <- gc()
      memory_trace[rep] <- sum(gc_info[, 2])  # Total memory used
    }
    
    # Manual cleanup
    rm(data, x, y, fit)
    
    # Force garbage collection every 10 iterations
    if (rep %% 10 == 0) {
      gc(full = TRUE)
    }
    
    setTxtProgressBar(pb, rep)
  }
  close(pb)
  
  # Plot memory usage trend
  plot(1:config$num_rep, memory_trace, 
       type = "l", 
       xlab = "Replication", 
       ylab = "Memory (MB)",
       main = "Memory Usage Over Replications")
  abline(lm(memory_trace ~ I(1:config$num_rep)), col = "red", lty = 2)
  
  # Check if there's a continuous growth trend
  mem_lm <- lm(memory_trace ~ I(1:config$num_rep))
  slope <- coef(mem_lm)[2]
  
  cat(sprintf("\nMemory growth rate: %.4f MB per replication\n", slope))
  
  if (slope > 0.1) {
    cat("⚠️  WARNING: Significant memory leak detected!\n")
    cat(sprintf("Estimated memory at rep 500: %.2f MB\n", 
                memory_trace[1] + slope * 500))
  } else {
    cat("✓ No significant memory leak detected\n")
  }
  
  return(memory_trace)
}

# ============================================================================ #
# Run diagnosis
# ============================================================================ #
mem_trace <- diagnose_memory_leak(sim_config)

# ============================================================================ #
# Test each method one by one
# ============================================================================ #
test_each_method <- function(config) {
  config$num_rep <- 100
  
  tvcqr_methods <- create_tvcqr_methods(config$Mm.factor)
  
  results <- list()
  
  for (method_name in names(tvcqr_methods)) {
    cat(sprintf("\n=== Testing %s ===\n", method_name))
    
    tryCatch({
      mem_before <- gc()[2, 2]
      
      for (rep in 1:config$num_rep) {
        data <- generate_ts(n = config$n, case = config$case, 
                            seed = config$seed_base + rep)
        x <- data$x
        y <- data$y
        
        fit <- tvcqr_methods[[method_name]](x, y, config)
        
        rm(data, x, y, fit)
        
        if (rep %% 20 == 0) {
          cat(sprintf("Rep %d completed\n", rep))
          gc()
        }
      }
      
      mem_after <- gc()[2, 2]
      
      results[[method_name]] <- list(
        success = TRUE,
        mem_increase = mem_after - mem_before
      )
      
      cat(sprintf("✓ %s completed successfully\n", method_name))
      cat(sprintf("Memory increase: %.2f MB\n", mem_after - mem_before))
      
    }, error = function(e) {
      results[[method_name]] <<- list(
        success = FALSE,
        error = e$message
      )
      cat(sprintf("✗ %s failed: %s\n", method_name, e$message))
    })
    
    # Cleanup
    gc(full = TRUE)
    Sys.sleep(1)
  }
  
  return(results)
}

# ============================================================================ #
# Run test
# ============================================================================ #
method_test_results <- test_each_method(sim_config)

# ============================================================================ #
# Error Capture Function with Data Saving
# ============================================================================ #

test_method_with_data_capture <- function(config, method_name, method_func, 
                                          Mm.factor = NULL) {
  cat(sprintf("\n=== Testing %s ===\n", method_name))
  if (!is.null(Mm.factor)) {
    cat(sprintf("Mm.factor = %.4f\n", Mm.factor))
  }
  
  # Create directory for error data if it doesn't exist
  error_data_dir <- "data/error_cases"
  if (!dir.exists(error_data_dir)) {
    dir.create(error_data_dir, recursive = TRUE)
  }
  
  success_count <- 0
  
  for (rep in 1:config$num_rep) {
    # Generate data
    data <- generate_ts(n = config$n, case = config$case, 
                        seed = config$seed_base + rep)
    x <- data$x
    y <- data$y
    
    # Create filename for current replication data
    data_filename <- file.path(
      error_data_dir,
      sprintf("temp_data_n%d_tau%.2f_rep%04d.rds", 
              config$n, config$tau, rep)
    )
    
    # Save data and metadata before running
    temp_data <- list(
      x = x,
      y = y,
      rep = rep,
      config = list(
        n = config$n,
        tau = config$tau,
        case = config$case,
        h = config$h,
        h.factor = config$h.factor,
        seed = config$seed_base + rep
      ),
      method_name = method_name,
      Mm.factor = Mm.factor,
      timestamp = Sys.time()
    )
    
    saveRDS(temp_data, data_filename)
    
    # Try to run the method
    tryCatch({
      fit <- method_func(x, y, config)
      
      # If successful, delete the temporary data file
      file.remove(data_filename)
      
      success_count <- success_count + 1
      
      # Clean up
      rm(fit, data, x, y, temp_data)
      
      if (rep %% 10 == 0) {
        cat(sprintf("Rep %d completed successfully\n", rep))
        gc()
      }
      
    }, error = function(e) {
      # If error occurs, keep the data file and create error log
      error_filename <- file.path(
        error_data_dir,
        sprintf("ERROR_n%d_tau%.2f_rep%04d.txt", 
                config$n, config$tau, rep)
      )
      
      error_log <- sprintf(
        "=== ERROR REPORT ===\n\n
Method: %s\n
Mm.factor: %s\n
Sample size (n): %d\n
Tau: %.2f\n
Replication: %d\n
Seed: %d\n
Case: %d\n
Timestamp: %s\n\n
Error message:\n%s\n\n
Data file saved as: %s\n",
        method_name,
        ifelse(is.null(Mm.factor), "N/A", sprintf("%.4f", Mm.factor)),
        config$n,
        config$tau,
        rep,
        config$seed_base + rep,
        config$case,
        as.character(Sys.time()),
        e$message,
        basename(data_filename)
      )
      
      writeLines(error_log, error_filename)
      
      cat(sprintf("\n✗ ERROR at Rep %d\n", rep))
      cat(sprintf("Error message: %s\n", e$message))
      cat(sprintf("Data saved to: %s\n", data_filename))
      cat(sprintf("Error log saved to: %s\n", error_filename))
      
      stop(sprintf("Stopped at replication %d due to error", rep))
    })
  }
  
  cat(sprintf("\n✓ All %d replications completed successfully\n", success_count))
  return(list(success = TRUE, completed_reps = success_count))
}

# ============================================================================ #
# Test All Methods with Error Capture
# ============================================================================ #

test_all_methods_with_capture <- function(config) {
  tvcqr_methods <- create_tvcqr_methods(config$Mm.factor)
  
  results <- list()
  
  for (method_name in names(tvcqr_methods)) {
    # Extract Mm.factor if applicable
    Mm.factor <- NULL
    if (grepl("ppro", method_name)) {
      # Find corresponding Mm.factor from mapping
      idx <- which(names(tvcqr_methods) == method_name)
      if (grepl("_fortran_", method_name)) {
        # For fortran methods
        ppro_idx <- as.numeric(gsub(".*_(\\d+)$", "\\1", method_name))
        Mm.factor <- config$Mm.factor[ppro_idx]
      } else if (grepl("_\\d+$", method_name)) {
        # For regular ppro methods
        ppro_idx <- as.numeric(gsub(".*_(\\d+)$", "\\1", method_name))
        Mm.factor <- config$Mm.factor[ppro_idx]
      }
    }
    
    tryCatch({
      result <- test_method_with_data_capture(
        config = config,
        method_name = method_name,
        method_func = tvcqr_methods[[method_name]],
        Mm.factor = Mm.factor
      )
      
      results[[method_name]] <- result
      
    }, error = function(e) {
      results[[method_name]] <<- list(
        success = FALSE,
        error = e$message
      )
      
      cat(sprintf("\n⚠️  Stopping test for %s\n", method_name))
    })
    
    # Clean up between methods
    gc(full = TRUE)
    Sys.sleep(1)
  }
  
  return(results)
}

# ============================================================================ #
# Function to Reproduce Error from Saved Data
# ============================================================================ #

reproduce_error_from_file <- function(data_filename, method_name = NULL, 
                                      Mm.factor = NULL) {
  cat("=== Reproducing Error ===\n\n")
  
  # Load saved data
  temp_data <- readRDS(data_filename)
  
  cat("Loaded data information:\n")
  cat(sprintf("  Method: %s\n", temp_data$method_name))
  cat(sprintf("  Mm.factor: %s\n", 
              ifelse(is.null(temp_data$Mm.factor), "N/A", 
                     sprintf("%.4f", temp_data$Mm.factor))))
  cat(sprintf("  n: %d\n", temp_data$config$n))
  cat(sprintf("  tau: %.2f\n", temp_data$config$tau))
  cat(sprintf("  Replication: %d\n", temp_data$rep))
  cat(sprintf("  Seed: %d\n", temp_data$config$seed))
  cat("\n")
  
  # Use saved method name and Mm.factor if not provided
  if (is.null(method_name)) {
    method_name <- temp_data$method_name
  }
  if (is.null(Mm.factor)) {
    Mm.factor <- temp_data$Mm.factor
  }
  
  # Recreate config
  config <- list(
    n = temp_data$config$n,
    tau = temp_data$config$tau,
    case = temp_data$config$case,
    h = temp_data$config$h,
    h.factor = temp_data$config$h.factor,
    tol = 1e-6,
    maxit = 500,
    bland = TRUE,
    eps = 1e-8,
    cpp_helper = TRUE,
    Mm.factor = if (!is.null(Mm.factor)) Mm.factor else 0.01
  )
  
  # Try to run the method
  cat("Attempting to reproduce error...\n")
  
  # Determine which function to call based on method name
  if (grepl("fortran", method_name)) {
    if (grepl("ppro", method_name)) {
      cat("Running tvcqr_seq_ppro_fortran_wrapper...\n")
      fit <- tvcqr_seq_ppro_fortran_wrapper(
        x = temp_data$x, 
        y = temp_data$y, 
        tau = config$tau, 
        h = config$h,
        h.factor = config$h.factor, 
        tol = config$tol, 
        maxit = config$maxit, 
        bland = config$bland,
        Mm.factor = config$Mm.factor, 
        eps = config$eps
      )
    } else {
      cat("Running tvcqr_seq_fortran...\n")
      fit <- tvcqr_seq_fortran(
        x = temp_data$x, 
        y = temp_data$y, 
        tau = config$tau, 
        h = config$h,
        h.factor = config$h.factor, 
        tol = config$tol, 
        maxit = config$maxit, 
        bland = config$bland
      )
    }
  } else if (grepl("ppro", method_name)) {
    cat("Running tvcqr_seq_ppro...\n")
    fit <- tvcqr_seq_ppro(
      x = temp_data$x, 
      y = temp_data$y, 
      tau = config$tau, 
      h = config$h,
      h.factor = config$h.factor, 
      tol = config$tol, 
      maxit = config$maxit, 
      bland = config$bland,
      Mm.factor = config$Mm.factor, 
      eps = config$eps, 
      cpp_helper = config$cpp_helper
    )
  } else if (method_name == "tvcqr_seq") {
    cat("Running tvcqr_seq...\n")
    fit <- tvcqr_seq(
      x = temp_data$x, 
      y = temp_data$y, 
      tau = config$tau, 
      h = config$h,
      h.factor = config$h.factor, 
      tol = config$tol, 
      maxit = config$maxit, 
      bland = config$bland
    )
  } else if (method_name == "tvc_rq") {
    cat("Running tvc_rq...\n")
    fit <- tvc_rq(
      x = temp_data$x, 
      y = temp_data$y, 
      tau = config$tau, 
      h = config$h
    )
  }
  
  cat("\n✓ Method ran successfully (no error reproduced)\n")
  
  return(fit)
}

# ============================================================================ #
# Utility: List Error Data Files
# ============================================================================ #

list_error_files <- function() {
  error_data_dir <- "data/error_cases"
  
  if (!dir.exists(error_data_dir)) {
    cat("No error data directory found.\n")
    return(NULL)
  }
  
  data_files <- list.files(error_data_dir, pattern = "^temp_data.*\\.rds$", 
                           full.names = TRUE)
  error_logs <- list.files(error_data_dir, pattern = "^ERROR.*\\.txt$", 
                           full.names = TRUE)
  
  cat("=== Error Data Files ===\n\n")
  
  if (length(data_files) > 0) {
    cat("Data files with errors:\n")
    for (f in data_files) {
      cat(sprintf("  - %s\n", basename(f)))
    }
  } else {
    cat("No error data files found.\n")
  }
  
  cat("\n")
  
  if (length(error_logs) > 0) {
    cat("Error log files:\n")
    for (f in error_logs) {
      cat(sprintf("  - %s\n", basename(f)))
      # Print error log content
      cat("    Content:\n")
      log_content <- readLines(f)
      cat(paste("    ", log_content, collapse = "\n"))
      cat("\n\n")
    }
  } else {
    cat("No error log files found.\n")
  }
  
  return(list(data_files = data_files, error_logs = error_logs))
}

# ============================================================================ #
# Utility: Clean up error cases directory
# ============================================================================ #

clean_error_cases <- function() {
  error_data_dir <- "data/error_cases"
  
  if (!dir.exists(error_data_dir)) {
    cat("No error data directory found.\n")
    return(invisible(NULL))
  }
  
  files <- list.files(error_data_dir, full.names = TRUE)
  
  if (length(files) == 0) {
    cat("Error cases directory is already empty.\n")
    return(invisible(NULL))
  }
  
  cat(sprintf("Found %d files in error_cases directory.\n", length(files)))
  cat("Delete all files? (y/n): ")
  response <- readline()
  
  if (tolower(response) == "y") {
    file.remove(files)
    cat(sprintf("✓ Deleted %d files.\n", length(files)))
  } else {
    cat("Operation cancelled.\n")
  }
}

# 1. Run test with error capture (errors saved to /data/error_cases)
cat("Testing all methods with error capture...\n")
test_results <- test_all_methods_with_capture(sim_config)

# 2. If crash occurs, restart R and check error files
error_files <- list_error_files()

# 3. Reproduce the error using saved data
if (!is.null(error_files) && length(error_files$data_files) > 0) {
  # Use the first error file
  error_data_file <- error_files$data_files[1]
  
  cat("\nAttempting to reproduce error from:", basename(error_data_file), "\n")
  
  # This will either reproduce the error or run successfully
  fit <- reproduce_error_from_file(error_data_file)
}

# 4. Clean up error cases directory when done
clean_error_cases()

# ============================================================================ #
# Test tvcqr_seq_ppro_fortran_wrapper with Data Capture
# Simplified version for Fortran crashes
# ============================================================================ #

test_fortran_ppro_with_data_save <- function(config, Mm.factor_index) {
  
  Mm.factor_value <- config$Mm.factor[Mm.factor_index]
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Testing tvcqr_seq_ppro_fortran_wrapper\n"))
  cat(sprintf("Mm.factor[%d] = %.6f\n", Mm.factor_index, Mm.factor_value))
  cat(sprintf("n = %d, tau = %.2f\n", config$n, config$tau))
  cat(sprintf("========================================\n\n"))
  
  # Create error_cases directory if it doesn't exist
  error_dir <- "data/error_cases"
  if (!dir.exists(error_dir)) {
    dir.create(error_dir, recursive = TRUE)
  }
  
  for (rep in 1:config$num_rep) {
    
    cat(sprintf("Rep %d/%d... ", rep, config$num_rep))
    
    # Generate data
    data <- generate_ts(n = config$n, case = config$case, 
                        seed = config$seed_base + rep)
    x <- data$x
    y <- data$y
    
    # Create descriptive filename
    # Use tau without decimal point
    tau_int <- as.integer(config$tau * 100)
    
    rdata_filename <- file.path(
      error_dir,
      sprintf("CRASH_DATA_n%d_tau%d_Mm%d_rep%04d.RData",
              config$n, tau_int, Mm.factor_index, rep)
    )
    
    # Save everything BEFORE running the method
    crash_data <- list(
      x = x,
      y = y,
      rep = rep,
      n = config$n,
      tau = config$tau,
      case = config$case,
      h = config$h,
      h.factor = config$h.factor,
      seed = config$seed_base + rep,
      Mm.factor_index = Mm.factor_index,
      Mm.factor_value = Mm.factor_value,
      all_Mm.factors = config$Mm.factor,  # Save all Mm.factor values for reference
      tol = config$tol,
      maxit = config$maxit,
      bland = config$bland,
      eps = config$eps,
      method = "tvcqr_seq_ppro_fortran_wrapper",
      timestamp = Sys.time()
    )
    
    save(crash_data, file = rdata_filename)
    
    cat(sprintf("Data saved to: %s\n", rdata_filename))
    
    # Run the fortran function
    # If it crashes, the .RData file will remain
    fit <- tvcqr_seq_ppro_fortran_wrapper(
      x = x,
      y = y,
      tau = config$tau,
      h = config$h,
      h.factor = config$h.factor,
      tol = config$tol,
      maxit = config$maxit,
      bland = config$bland,
      Mm.factor = Mm.factor_value,
      eps = config$eps
    )
    
    # If we reach here, it succeeded - delete the .RData file
    file.remove(rdata_filename)
    
    cat("Success (file deleted)\n")
    
    # Clean up
    rm(data, x, y, fit, crash_data)
    
    # Garbage collection every 10 reps
    if (rep %% 10 == 0) {
      gc()
    }
  }
  
  cat(sprintf("\n✓ All %d replications completed successfully!\n", config$num_rep))
}

# ============================================================================ #
# Test all Mm.factor values
# ============================================================================ #

test_all_Mm_factors <- function(config) {
  
  cat("========================================\n")
  cat("Testing all Mm.factor values\n")
  cat("Mm.factor values:\n")
  for (i in seq_along(config$Mm.factor)) {
    cat(sprintf("  [%d] %.6f\n", i, config$Mm.factor[i]))
  }
  cat("========================================\n")
  
  for (i in seq_along(config$Mm.factor)) {
    cat(sprintf("\n\n*** Testing Mm.factor[%d] = %.6f ***\n", i, config$Mm.factor[i]))
    
    test_fortran_ppro_with_data_save(config, Mm.factor_index = i)
    
    # Clean up between Mm.factor tests
    gc(full = TRUE)
    Sys.sleep(1)
  }
  
  cat("\n\n========================================\n")
  cat("All tests completed!\n")
  cat("========================================\n")
}

# ============================================================================ #
# Load and inspect crash data
# ============================================================================ #

load_crash_data <- function(filename = NULL) {
  error_dir <- "data/error_cases"
  
  if (is.null(filename)) {
    # List all crash data files
    crash_files <- list.files(error_dir, pattern = "^CRASH_DATA.*\\.RData$", 
                              full.names = TRUE)
    
    if (length(crash_files) == 0) {
      cat("No crash data files found.\n")
      return(NULL)
    }
    
    cat("Found crash data files:\n")
    for (i in seq_along(crash_files)) {
      fname <- basename(crash_files[i])
      
      # Extract values from filename
      n_val <- as.numeric(gsub(".*_n(\\d+)_.*", "\\1", fname))
      tau_val <- as.numeric(gsub(".*_tau(\\d+)_.*", "\\1", fname)) / 100
      Mm_idx <- as.numeric(gsub(".*_Mm(\\d+)_.*", "\\1", fname))
      rep_val <- as.numeric(gsub(".*_rep(\\d+)\\.RData", "\\1", fname))
      
      cat(sprintf("%d. %s\n", i, fname))
      cat(sprintf("   (n=%d, tau=%.2f, Mm.factor[%d], rep=%d)\n", 
                  n_val, tau_val, Mm_idx, rep_val))
    }
    
    if (length(crash_files) == 1) {
      filename <- crash_files[1]
      cat(sprintf("\nLoading: %s\n", basename(filename)))
    } else {
      cat("\nPlease specify which file to load.\n")
      return(crash_files)
    }
  } else {
    # If filename doesn't include full path, assume it's in error_dir
    if (!file.exists(filename)) {
      filename <- file.path(error_dir, filename)
    }
  }
  
  # Load the crash data
  load(filename)
  
  cat("\n========================================\n")
  cat("Crash Data Information\n")
  cat("========================================\n")
  cat(sprintf("Method: %s\n", crash_data$method))
  cat(sprintf("Mm.factor[%d]: %.6f (%.0e)\n", 
              crash_data$Mm.factor_index, 
              crash_data$Mm.factor_value, 
              crash_data$Mm.factor_value))
  cat(sprintf("\nAll Mm.factor values:\n"))
  for (i in seq_along(crash_data$all_Mm.factors)) {
    marker <- if (i == crash_data$Mm.factor_index) " <-- CRASH" else ""
    cat(sprintf("  [%d] %.6f%s\n", i, crash_data$all_Mm.factors[i], marker))
  }
  cat(sprintf("\nSample size (n): %d\n", crash_data$n))
  cat(sprintf("Tau: %.2f\n", crash_data$tau))
  cat(sprintf("Replication: %d\n", crash_data$rep))
  cat(sprintf("Seed: %d\n", crash_data$seed))
  cat(sprintf("Case: %d\n", crash_data$case))
  cat(sprintf("h: %.4f\n", crash_data$h))
  cat(sprintf("h.factor: %.4f\n", crash_data$h.factor))
  cat(sprintf("Timestamp: %s\n", crash_data$timestamp))
  cat(sprintf("\nData dimensions:\n"))
  cat(sprintf("  x: %s\n", paste(dim(crash_data$x), collapse = " x ")))
  cat(sprintf("  y: length %d\n", length(crash_data$y)))
  cat("========================================\n\n")
  
  return(crash_data)
}

# ============================================================================ #
# Reproduce crash from saved data
# ============================================================================ #

reproduce_crash <- function(crash_data) {
  cat("Attempting to reproduce crash...\n\n")
  cat(sprintf("Running %s\n", crash_data$method))
  cat(sprintf("Mm.factor[%d] = %.6f\n", 
              crash_data$Mm.factor_index, crash_data$Mm.factor_value))
  cat("If this crashes, you've successfully reproduced the issue!\n\n")
  
  fit <- tvcqr_seq_ppro_fortran_wrapper(
    x = crash_data$x,
    y = crash_data$y,
    tau = crash_data$tau,
    h = crash_data$h,
    h.factor = crash_data$h.factor,
    tol = crash_data$tol,
    maxit = crash_data$maxit,
    bland = crash_data$bland,
    Mm.factor = crash_data$Mm.factor_value,
    eps = crash_data$eps
  )
  
  cat("\n✓ Method ran successfully (no crash)\n")
  return(fit)
}

# ============================================================================ #
# Clean up error_cases directory
# ============================================================================ #

clean_error_cases <- function(confirm = TRUE) {
  error_dir <- "data/error_cases"
  
  if (!dir.exists(error_dir)) {
    cat("Error cases directory does not exist.\n")
    return(invisible(NULL))
  }
  
  files <- list.files(error_dir, pattern = "^CRASH_DATA.*\\.RData$", 
                      full.names = TRUE)
  
  if (length(files) == 0) {
    cat("No crash data files to clean.\n")
    return(invisible(NULL))
  }
  
  cat(sprintf("Found %d crash data file(s):\n", length(files)))
  for (f in files) {
    cat(sprintf("  - %s\n", basename(f)))
  }
  
  if (confirm) {
    cat("\nDelete all these files? (y/n): ")
    response <- readline()
    if (tolower(response) != "y") {
      cat("Operation cancelled.\n")
      return(invisible(NULL))
    }
  }
  
  file.remove(files)
  cat(sprintf("✓ Deleted %d file(s).\n", length(files)))
}

# ============================================================================ #
# Step 1: Run the test
# ============================================================================ #

# Test all Mm.factor values
test_all_Mm_factors(sim_config)

# OR test just one specific Mm.factor index
# test_fortran_ppro_with_data_save(sim_config, Mm.factor_index = 1)

# ============================================================================ #
# Step 2: After R crashes and restarts, load the crash data
# ============================================================================ #

# List all crash files
crash_files <- load_crash_data()

# Load the crash data
crash_data <- load_crash_data(crash_files[1])

# ============================================================================ #
# Step 3: Try to reproduce the crash
# ============================================================================ #
fit <- reproduce_crash(crash_data)

# ============================================================================ #
# Step 4: Clean up when done
# ============================================================================ #

clean_error_cases()