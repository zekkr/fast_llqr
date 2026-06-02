# Source project setup and dependencies.
source("scripts/setup.R")

# ============================================================================ #
# Multivariate LLQR iteration-reduction simulation
# ============================================================================ #

simulate_llqr_multivar_data <- function(n, p, case = 1L, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (case == 1L) {
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    f <- 1 +
      1.2 * x[, 1L] +
      if (p >= 2L) -0.9 * x[, 2L] else 0 +
      if (p >= 3L) 0.5 * sin(x[, 1L] + x[, 3L]) else 0 +
      if (p >= 4L) 0.7 * x[, 3L] * x[, 4L] else 0
    eps <- 0.5 * rt(n, df = 4)
  } else if (case == 2L) {
    rho <- 0.5
    sigma <- outer(seq_len(p), seq_len(p), function(i, j) rho^abs(i - j))
    x <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol(sigma)
    hetero <- 0.4 + 0.25 * abs(x[, 1L]) + if (p >= 2L) 0.1 * abs(x[, 2L]) else 0
    f <- 0.5 +
      0.8 * x[, 1L] +
      if (p >= 2L) 0.6 * x[, 2L]^2 else 0 +
      if (p >= 3L) -0.4 * x[, 3L] else 0 +
      if (p >= 4L) 0.5 * x[, 1L] * x[, 4L] else 0
    eps <- hetero * rnorm(n)
  } else {
    stop("case must be 1 or 2.")
  }

  colnames(x) <- paste0("x", seq_len(p))
  y <- as.numeric(f + eps)
  list(x = x, y = y)
}

summarize_iteration_fit <- function(fit, method_label, case_id, n, p, rep_id) {
  data.frame(
    case = case_id,
    n = n,
    p = p,
    replication = rep_id,
    method = method_label,
    total_iterations = sum(fit$it_num),
    mean_iterations = mean(fit$it_num),
    median_iterations = median(fit$it_num),
    max_iterations = max(fit$it_num),
    elapsed = unname(fit$elapsed),
    stringsAsFactors = FALSE
  )
}

compute_iteration_reduction <- function(raw_metrics) {
  cold <- raw_metrics[raw_metrics$method == "cold", ]
  input <- raw_metrics[raw_metrics$method == "seq_input", ]
  mst <- raw_metrics[raw_metrics$method == "seq_mst", ]

  merge(
    merge(
      cold[, c("case", "n", "p", "replication", "total_iterations", "mean_iterations", "median_iterations", "max_iterations")],
      input[, c("case", "n", "p", "replication", "total_iterations", "mean_iterations", "median_iterations", "max_iterations")],
      by = c("case", "n", "p", "replication"),
      suffixes = c("_cold", "_input")
    ),
    mst[, c("case", "n", "p", "replication", "total_iterations", "mean_iterations", "median_iterations", "max_iterations")],
    by = c("case", "n", "p", "replication")
  )
}

read_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (identical(value, "")) {
    return(default)
  }

  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) {
    stop(sprintf("Environment variable %s must be a positive integer.", name))
  }
  parsed
}

config <- list(
  cases = c(1L, 2L),
  n_values = c(300L, 500L, 800L),
  p_values = c(4L),
  tau = 0.5,
  num_rep = read_env_int("FASTQR_LLQR_MULTI_NUM_REP", 50L),
  tol = 1e-14,
  maxit = 20000L,
  bland = FALSE,
  seed_base = 20260405L,
  check_accuracy = TRUE,
  accuracy_tol = 1e-8,
  output_dir = "results/llqr_multivar_iterations"
)

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

raw_rows <- list()
diagnostic_rows <- list()
row_id <- 1L
diag_id <- 1L

for (case_id in config$cases) {
  for (n in config$n_values) {
    for (p in config$p_values) {
      cat("========================================\n")
      cat(sprintf("Running case=%d, n=%d, p=%d\n", case_id, n, p))
      cat("========================================\n")

      for (rep_id in seq_len(config$num_rep)) {
        seed_used <- config$seed_base + 100000L * case_id + 1000L * n + 10L * p + rep_id
        dat <- simulate_llqr_multivar_data(n = n, p = p, case = case_id, seed = seed_used)
        x <- dat$x
        y <- dat$y
        z <- x

        h_common <- compute_llqr_multivar_bandwidth(x = x, y = y, tau = config$tau)

        cold_time <- system.time({
          fit_cold <- llqr_tau_multivar(
            x = x, y = y, tau = config$tau, z = z, h = h_common,
            tol = config$tol, maxit = config$maxit,
            bland = config$bland, track_order = TRUE
          )
        })
        fit_cold$elapsed <- cold_time["elapsed"]

        input_time <- system.time({
          fit_seq_input <- llqr_tau_seq_multivar(
            x = x, y = y, tau = config$tau, z = z, h = h_common,
            tol = config$tol, maxit = config$maxit,
            bland = config$bland, track_order = TRUE,
            order_method = "input"
          )
        })
        fit_seq_input$elapsed <- input_time["elapsed"]

        mst_time <- system.time({
          fit_seq_mst <- llqr_tau_seq_multivar(
            x = x, y = y, tau = config$tau, z = z, h = h_common,
            tol = config$tol, maxit = config$maxit,
            bland = config$bland, track_order = TRUE,
            order_method = "mst",
            root_method = "center",
            distance_scale = "bandwidth"
          )
        })
        fit_seq_mst$elapsed <- mst_time["elapsed"]

        if (isTRUE(config$check_accuracy)) {
          input_diff <- max(abs(as.numeric(fit_seq_input$ll_est) - as.numeric(fit_cold$ll_est)))
          mst_diff <- max(abs(as.numeric(fit_seq_mst$ll_est) - as.numeric(fit_cold$ll_est)))
          if (input_diff > config$accuracy_tol || mst_diff > config$accuracy_tol) {
            warning(
              sprintf(
                "Accuracy check failed at case=%d, n=%d, p=%d, rep=%d: input diff=%.3e, mst diff=%.3e",
                case_id, n, p, rep_id, input_diff, mst_diff
              )
            )
          }
        } else {
          input_diff <- NA_real_
          mst_diff <- NA_real_
        }

        raw_rows[[row_id]] <- summarize_iteration_fit(
          fit = fit_cold, method_label = "cold", case_id = case_id, n = n, p = p, rep_id = rep_id
        )
        row_id <- row_id + 1L
        raw_rows[[row_id]] <- summarize_iteration_fit(
          fit = fit_seq_input, method_label = "seq_input", case_id = case_id, n = n, p = p, rep_id = rep_id
        )
        row_id <- row_id + 1L
        raw_rows[[row_id]] <- summarize_iteration_fit(
          fit = fit_seq_mst, method_label = "seq_mst", case_id = case_id, n = n, p = p, rep_id = rep_id
        )
        row_id <- row_id + 1L

        diagnostic_rows[[diag_id]] <- data.frame(
          case = case_id,
          n = n,
          p = p,
          replication = rep_id,
          h_summary = paste(signif(as.numeric(h_common), 4), collapse = ","),
          input_max_abs_diff = input_diff,
          mst_max_abs_diff = mst_diff,
          mst_mean_edge_weight = mean(fit_seq_mst$edge_weight[-fit_seq_mst$root]),
          mst_max_edge_weight = max(fit_seq_mst$edge_weight[-fit_seq_mst$root]),
          input_mean_edge_weight = mean(fit_seq_input$edge_weight[-fit_seq_input$root]),
          input_max_edge_weight = max(fit_seq_input$edge_weight[-fit_seq_input$root]),
          stringsAsFactors = FALSE
        )
        diag_id <- diag_id + 1L
      }
    }
  }
}

raw_metrics <- do.call(rbind, raw_rows)
diagnostics <- do.call(rbind, diagnostic_rows)

summary_metrics <- raw_metrics |>
  dplyr::group_by(case, n, p, method) |>
  dplyr::summarise(
    mean_total_iterations = mean(total_iterations),
    median_total_iterations = median(total_iterations),
    mean_mean_iterations = mean(mean_iterations),
    mean_median_iterations = mean(median_iterations),
    mean_max_iterations = mean(max_iterations),
    mean_elapsed = mean(elapsed),
    .groups = "drop"
  )

contrast_base <- compute_iteration_reduction(raw_metrics)
reduction_metrics <- contrast_base |>
  dplyr::mutate(
    total_reduction_input_vs_cold = 1 - total_iterations_input / total_iterations_cold,
    total_reduction_mst_vs_cold = 1 - total_iterations / total_iterations_cold,
    total_reduction_mst_vs_input = 1 - total_iterations / total_iterations_input,
    mean_reduction_input_vs_cold = 1 - mean_iterations_input / mean_iterations_cold,
    mean_reduction_mst_vs_cold = 1 - mean_iterations / mean_iterations_cold,
    mean_reduction_mst_vs_input = 1 - mean_iterations / mean_iterations_input
  )

reduction_summary <- reduction_metrics |>
  dplyr::group_by(case, n, p) |>
  dplyr::summarise(
    avg_total_reduction_input_vs_cold = mean(total_reduction_input_vs_cold),
    avg_total_reduction_mst_vs_cold = mean(total_reduction_mst_vs_cold),
    avg_total_reduction_mst_vs_input = mean(total_reduction_mst_vs_input),
    avg_mean_reduction_input_vs_cold = mean(mean_reduction_input_vs_cold),
    avg_mean_reduction_mst_vs_cold = mean(mean_reduction_mst_vs_cold),
    avg_mean_reduction_mst_vs_input = mean(mean_reduction_mst_vs_input),
    .groups = "drop"
  )

timestamp_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
raw_path <- file.path(config$output_dir, paste0("llqr_multivar_iteration_raw_", timestamp_tag, ".csv"))
diag_path <- file.path(config$output_dir, paste0("llqr_multivar_iteration_diagnostics_", timestamp_tag, ".csv"))
summary_path <- file.path(config$output_dir, paste0("llqr_multivar_iteration_summary_", timestamp_tag, ".csv"))
reduction_path <- file.path(config$output_dir, paste0("llqr_multivar_iteration_reduction_", timestamp_tag, ".csv"))
rds_path <- file.path(config$output_dir, paste0("llqr_multivar_iteration_results_", timestamp_tag, ".Rds"))

readr::write_csv(raw_metrics, raw_path)
readr::write_csv(diagnostics, diag_path)
readr::write_csv(summary_metrics, summary_path)
readr::write_csv(reduction_summary, reduction_path)

saveRDS(
  list(
    config = config,
    raw_metrics = raw_metrics,
    diagnostics = diagnostics,
    summary_metrics = summary_metrics,
    reduction_metrics = reduction_metrics,
    reduction_summary = reduction_summary
  ),
  file = rds_path
)

cat("\n=== Iteration summary ===\n")
print(summary_metrics, n = Inf)

cat("\n=== Iteration reduction summary ===\n")
print(reduction_summary, n = Inf)

cat("\nSaved files:\n")
cat(raw_path, "\n")
cat(diag_path, "\n")
cat(summary_path, "\n")
cat(reduction_path, "\n")
cat(rds_path, "\n")
