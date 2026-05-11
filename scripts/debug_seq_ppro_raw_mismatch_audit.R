#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

source("R/llqr_functions.R")
source("R/tvcqr_functions.R")

parse_csv_num <- function(x) {
  if (is.null(x) || identical(x, "")) {
    return(numeric(0))
  }
  vals <- strsplit(x, ",", fixed = TRUE)[[1]]
  as.numeric(trimws(vals[nchar(trimws(vals)) > 0]))
}

parse_csv_int <- function(x) {
  if (is.null(x) || identical(x, "")) {
    return(integer(0))
  }
  vals <- strsplit(x, ",", fixed = TRUE)[[1]]
  as.integer(trimws(vals[nchar(trimws(vals)) > 0]))
}

parse_args <- function(args) {
  cfg <- list(
    rep = 500L,
    seed_base = 2026L,
    case_list = c(1L, 2L),
    tau_list = c(0.2, 0.5, 0.8),
    n_list = c(200L, 500L, 1000L),
    llqr_factors = c(1e-2, 1e-3, 1e-4),
    tvcqr_factors = c(1e-3, 1e-4),
    cores = max(1L, parallel::detectCores(logical = FALSE) - 1L),
    output_dir = "tmp/raw_backend_hseq_audit"
  )

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    value <- if (length(kv) > 1L) kv[2] else ""

    if (key == "rep") cfg$rep <- as.integer(value)
    if (key == "seed-base") cfg$seed_base <- as.integer(value)
    if (key == "case-list") cfg$case_list <- parse_csv_int(value)
    if (key == "tau-list") cfg$tau_list <- parse_csv_num(value)
    if (key == "n-list") cfg$n_list <- parse_csv_int(value)
    if (key == "llqr-factors") cfg$llqr_factors <- parse_csv_num(value)
    if (key == "tvcqr-factors") cfg$tvcqr_factors <- parse_csv_num(value)
    if (key == "cores") cfg$cores <- as.integer(value)
    if (key == "output-dir") cfg$output_dir <- value
  }

  cfg$seeds <- seq.int(cfg$seed_base + 1L, cfg$seed_base + cfg$rep)
  cfg
}

llqr_bandwidth <- function(x, y, tau, h = NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  m <- nrow(x)
  nvar <- ncol(x)

  if (!is.null(h) && !is.na(h)) {
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

tvcqr_bandwidth <- function(m, h = NULL, h.factor = 1) {
  if (!is.null(h) && !is.na(h)) {
    return(as.numeric(h))
  }
  as.numeric(m^(-0.2) * h.factor)
}

llqr_row_set_match <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  if (!all(dim(a) == dim(b))) {
    return(FALSE)
  }
  all(vapply(seq_len(nrow(a)), function(i) llqr_set_equal_int(a[i, ], b[i, ]), logical(1)))
}

tvcqr_row_set_match_vec <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  vapply(seq_len(nrow(a)), function(i) tvcqr_set_equal_int(a[i, ], b[i, ]), logical(1))
}

llqr_row_set_match_vec <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  vapply(seq_len(nrow(a)), function(i) llqr_set_equal_int(a[i, ], b[i, ]), logical(1))
}

tvcqr_row_set_match <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  if (!all(dim(a) == dim(b))) {
    return(FALSE)
  }
  all(vapply(seq_len(nrow(a)), function(i) tvcqr_set_equal_int(a[i, ], b[i, ]), logical(1)))
}

run_llqr_one <- function(case, tau, n, seed, mm_factor, method_label) {
  tryCatch({
    data <- generate_data(n = n, case = case, seed = seed)
    x <- as.matrix(data$x)
    y <- as.numeric(data$y)
    z <- x[order(x), , drop = FALSE]
    h <- llqr_bandwidth(x = x, y = y, tau = tau)

    seq_fit <- llqr_seq(
      x = x,
      y = y,
      tau = tau,
      z = z,
      h = h,
      track_order = FALSE
    )
    ppro_fit <- llqr_seq_ppro_fortran_wrapper(
      x = as.numeric(x),
      y = y,
      tau = tau,
      z = as.numeric(z),
      h = h,
      Mm.factor = mm_factor,
      case = case,
      bland = TRUE,
      track_order = FALSE,
      fallback = FALSE,
      return_raw_backend = TRUE
    )

    row_match <- llqr_row_set_match_vec(ppro_fit$H_seq, seq_fit$H_seq)
    h_match <- all(row_match)
    first_bad_round <- if (all(row_match)) NA_integer_ else which(!row_match)[1]
    data.frame(
      model = "llqr",
      method = method_label,
      Mm.factor = mm_factor,
      case = case,
      tau = tau,
      n = n,
      seed = seed,
      backend_ierr = if (!is.null(ppro_fit$backend_ierr)) ppro_fit$backend_ierr else 0L,
      returned_backend = if (!is.null(ppro_fit$returned_backend)) ppro_fit$returned_backend else NA_character_,
      cert_fail_detected = NA,
      first_bad_round = first_bad_round,
      H_set_match = h_match,
      H_exact_match = isTRUE(all.equal(ppro_fit$H_seq, seq_fit$H_seq, tolerance = 0)),
      status = if (isTRUE((if (!is.null(ppro_fit$backend_ierr)) ppro_fit$backend_ierr else 0L) == 0L)) "success" else "backend_ierr",
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      model = "llqr",
      method = method_label,
      Mm.factor = mm_factor,
      case = case,
      tau = tau,
      n = n,
      seed = seed,
      backend_ierr = NA_integer_,
      returned_backend = NA_character_,
      cert_fail_detected = NA,
      first_bad_round = NA_integer_,
      H_set_match = NA,
      H_exact_match = NA,
      status = "error",
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

run_tvcqr_one <- function(case, tau, n, seed, mm_factor, method_label) {
  tryCatch({
    data <- generate_ts(n = n, case = case, seed = seed)
    x <- as.matrix(data$x)
    y <- as.numeric(data$y)
    h <- tvcqr_bandwidth(nrow(x), h = NULL, h.factor = 1)

    seq_fit <- tvcqr_seq(
      x = x,
      y = y,
      tau = tau,
      h = h,
      h.factor = 1
    )
    ppro_fit <- tvcqr_seq_ppro_fortran_wrapper(
      x = x,
      y = y,
      tau = tau,
      h = h,
      h.factor = 1,
      Mm.factor = mm_factor,
      store_residual = FALSE,
      fallback = FALSE
    )

    row_match <- tvcqr_row_set_match_vec(ppro_fit$H_seq, seq_fit$H_seq)
    h_match <- all(row_match)
    first_bad_round <- if (all(row_match)) NA_integer_ else which(!row_match)[1]
    data.frame(
      model = "tvcqr",
      method = method_label,
      Mm.factor = mm_factor,
      case = case,
      tau = tau,
      n = n,
      seed = seed,
      backend_ierr = if (!is.null(ppro_fit$backend_ierr)) ppro_fit$backend_ierr else 0L,
      returned_backend = if (!is.null(ppro_fit$returned_backend)) ppro_fit$returned_backend else NA_character_,
      cert_fail_detected = NA,
      first_bad_round = first_bad_round,
      H_set_match = h_match,
      H_exact_match = isTRUE(all.equal(ppro_fit$H_seq, seq_fit$H_seq, tolerance = 0)),
      status = if (isTRUE((if (!is.null(ppro_fit$backend_ierr)) ppro_fit$backend_ierr else 0L) == 0L)) "success" else "backend_ierr",
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      model = "tvcqr",
      method = method_label,
      Mm.factor = mm_factor,
      case = case,
      tau = tau,
      n = n,
      seed = seed,
      backend_ierr = NA_integer_,
      returned_backend = NA_character_,
      cert_fail_detected = NA,
      first_bad_round = NA_integer_,
      H_set_match = NA,
      H_exact_match = NA,
      status = "error",
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

build_tasks <- function(cfg) {
  tasks <- list()
  idx <- 1L

  for (case in cfg$case_list) {
    for (tau in cfg$tau_list) {
      for (n in cfg$n_list) {
        if (length(cfg$llqr_factors) > 0L) {
          tasks[[idx]] <- list(
            model = "llqr",
            case = case,
            tau = tau,
            n = n,
            seeds = cfg$seeds,
            mm_factors = cfg$llqr_factors
          )
          idx <- idx + 1L
        }
        if (length(cfg$tvcqr_factors) > 0L) {
          tasks[[idx]] <- list(
            model = "tvcqr",
            case = case,
            tau = tau,
            n = n,
            seeds = cfg$seeds,
            mm_factors = cfg$tvcqr_factors
          )
          idx <- idx + 1L
        }
      }
    }
  }

  tasks
}

run_llqr_config <- function(task) {
  rows <- vector("list", length(task$seeds) * length(task$mm_factors))
  idx <- 1L

  for (seed in task$seeds) {
    data <- generate_data(n = task$n, case = task$case, seed = seed)
    x <- as.matrix(data$x)
    y <- as.numeric(data$y)
    z <- x[order(x), , drop = FALSE]
    h <- llqr_bandwidth(x = x, y = y, tau = task$tau)
    seq_fit <- llqr_seq(
      x = x,
      y = y,
      tau = task$tau,
      z = z,
      h = h,
      track_order = FALSE
    )

    for (i in seq_along(task$mm_factors)) {
      mm_factor <- task$mm_factors[i]
      method_label <- sprintf("llqr_seq_ppro_fortran_%d", i)
      rows[[idx]] <- tryCatch({
        ppro_fit <- llqr_seq_ppro_fortran_wrapper(
          x = as.numeric(x),
          y = y,
          tau = task$tau,
          z = as.numeric(z),
          h = h,
          Mm.factor = mm_factor,
          case = task$case,
          bland = TRUE,
          track_order = FALSE,
          fallback = FALSE,
          return_raw_backend = TRUE
        )

        row_match <- llqr_row_set_match_vec(ppro_fit$H_seq, seq_fit$H_seq)
        h_match <- all(row_match)
        first_bad_round <- if (all(row_match)) NA_integer_ else which(!row_match)[1]
        data.frame(
          model = "llqr",
          method = method_label,
          Mm.factor = mm_factor,
          case = task$case,
          tau = task$tau,
          n = task$n,
          seed = seed,
          backend_ierr = if (!is.null(ppro_fit$backend_ierr)) ppro_fit$backend_ierr else 0L,
          returned_backend = if (!is.null(ppro_fit$returned_backend)) ppro_fit$returned_backend else NA_character_,
          cert_fail_detected = NA,
          first_bad_round = first_bad_round,
          H_set_match = h_match,
          H_exact_match = isTRUE(all.equal(ppro_fit$H_seq, seq_fit$H_seq, tolerance = 0)),
          status = if (isTRUE((if (!is.null(ppro_fit$backend_ierr)) ppro_fit$backend_ierr else 0L) == 0L)) "success" else "backend_ierr",
          error_message = NA_character_,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        data.frame(
          model = "llqr",
          method = method_label,
          Mm.factor = mm_factor,
          case = task$case,
          tau = task$tau,
          n = task$n,
          seed = seed,
          backend_ierr = NA_integer_,
          returned_backend = NA_character_,
          cert_fail_detected = NA,
          first_bad_round = NA_integer_,
          H_set_match = NA,
          H_exact_match = NA,
          status = "error",
          error_message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      })
      idx <- idx + 1L
    }
  }

  do.call(rbind, rows)
}

run_tvcqr_config <- function(task) {
  rows <- vector("list", length(task$seeds) * length(task$mm_factors))
  idx <- 1L

  for (seed in task$seeds) {
    data <- generate_ts(n = task$n, case = task$case, seed = seed)
    x <- as.matrix(data$x)
    y <- as.numeric(data$y)
    h <- tvcqr_bandwidth(nrow(x), h = NULL, h.factor = 1)
    seq_fit <- tvcqr_seq(
      x = x,
      y = y,
      tau = task$tau,
      h = h,
      h.factor = 1
    )

    for (i in seq_along(task$mm_factors)) {
      mm_factor <- task$mm_factors[i]
      method_label <- sprintf("tvcqr_seq_ppro_fortran_%d", i)
      rows[[idx]] <- tryCatch({
        ppro_fit <- tvcqr_seq_ppro_fortran_wrapper(
          x = x,
          y = y,
          tau = task$tau,
          h = h,
          h.factor = 1,
          Mm.factor = mm_factor,
          store_residual = FALSE,
          fallback = FALSE
        )

        row_match <- tvcqr_row_set_match_vec(ppro_fit$H_seq, seq_fit$H_seq)
        h_match <- all(row_match)
        first_bad_round <- if (all(row_match)) NA_integer_ else which(!row_match)[1]
        data.frame(
          model = "tvcqr",
          method = method_label,
          Mm.factor = mm_factor,
          case = task$case,
          tau = task$tau,
          n = task$n,
          seed = seed,
          backend_ierr = if (!is.null(ppro_fit$backend_ierr)) ppro_fit$backend_ierr else 0L,
          returned_backend = if (!is.null(ppro_fit$returned_backend)) ppro_fit$returned_backend else NA_character_,
          cert_fail_detected = NA,
          first_bad_round = first_bad_round,
          H_set_match = h_match,
          H_exact_match = isTRUE(all.equal(ppro_fit$H_seq, seq_fit$H_seq, tolerance = 0)),
          status = if (isTRUE((if (!is.null(ppro_fit$backend_ierr)) ppro_fit$backend_ierr else 0L) == 0L)) "success" else "backend_ierr",
          error_message = NA_character_,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        data.frame(
          model = "tvcqr",
          method = method_label,
          Mm.factor = mm_factor,
          case = task$case,
          tau = task$tau,
          n = task$n,
          seed = seed,
          backend_ierr = NA_integer_,
          returned_backend = NA_character_,
          cert_fail_detected = NA,
          first_bad_round = NA_integer_,
          H_set_match = NA,
          H_exact_match = NA,
          status = "error",
          error_message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      })
      idx <- idx + 1L
    }
  }

  do.call(rbind, rows)
}

run_task <- function(task) {
  if (identical(task$model, "llqr")) {
    return(run_llqr_config(task))
  }
  run_tvcqr_config(task)
}

aggregate_summary <- function(results) {
  split_keys <- interaction(results$model, results$method, results$case, results$tau, results$n, drop = TRUE)
  chunks <- split(results, split_keys)

  rows <- lapply(chunks, function(df) {
    success_mask <- df$status == "success"
    mismatch_mask <- success_mask & !df$H_set_match
    data.frame(
      model = df$model[1],
      method = df$method[1],
      Mm.factor = df$Mm.factor[1],
      case = df$case[1],
      tau = df$tau[1],
      n = df$n[1],
      n_total_rep = nrow(df),
      n_success_rep = sum(success_mask, na.rm = TRUE),
      n_failed_rep = sum(df$status != "success", na.rm = TRUE),
      n_backend_ierr_rep = sum(df$status == "backend_ierr", na.rm = TRUE),
      n_error_rep = sum(df$status == "error", na.rm = TRUE),
      n_mismatch_rep = sum(mismatch_mask, na.rm = TRUE),
      mismatch_rep_ratio = sum(mismatch_mask, na.rm = TRUE) / nrow(df),
      n_cert_fail_rep = sum(success_mask & df$cert_fail_detected, na.rm = TRUE),
      median_first_bad_round = if (all(is.na(df$first_bad_round))) NA_real_ else stats::median(df$first_bad_round, na.rm = TRUE),
      min_first_bad_round = if (all(is.na(df$first_bad_round))) NA_integer_ else min(df$first_bad_round, na.rm = TRUE),
      returned_backend_values = paste(sort(unique(df$returned_backend[!is.na(df$returned_backend)])), collapse = ","),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

  tasks <- build_tasks(cfg)
  cat(sprintf("Running raw-backend audit: %d config tasks, %d cores\n", length(tasks), cfg$cores))

  results_list <- parallel::mclapply(tasks, run_task, mc.cores = cfg$cores)
  results <- do.call(rbind, results_list)
  summary_df <- aggregate_summary(results)
  summary_df <- summary_df[order(summary_df$model, summary_df$method, summary_df$case, summary_df$tau, summary_df$n), ]

  results_path <- file.path(cfg$output_dir, "raw_backend_rep_results.csv")
  summary_path <- file.path(cfg$output_dir, "raw_backend_summary.csv")
  rds_path <- file.path(cfg$output_dir, "raw_backend_results.rds")

  utils::write.csv(results, results_path, row.names = FALSE)
  utils::write.csv(summary_df, summary_path, row.names = FALSE)
  saveRDS(list(config = cfg, results = results, summary = summary_df), rds_path)

  print(summary_df)
  cat(sprintf("Saved rep-level results: %s\n", results_path))
  cat(sprintf("Saved summary results: %s\n", summary_path))
  cat(sprintf("Saved audit object: %s\n", rds_path))
}

main()
