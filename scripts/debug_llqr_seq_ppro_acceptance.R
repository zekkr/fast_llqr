#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

source("R/llqr_functions.R")

compute_bandwidth <- function(x, y, tau, h = NULL) {
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

parse_csv_int <- function(x) {
  if (is.null(x) || identical(x, "")) {
    return(integer(0))
  }
  vals <- strsplit(x, ",", fixed = TRUE)[[1]]
  as.integer(trimws(vals[nchar(trimws(vals)) > 0]))
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || identical(x, "")) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

parse_args <- function(args) {
  cfg <- list(
    mode = "round",
    backend = "r_ppro",
    n = 200L,
    case = 2L,
    seed = 2048L,
    rd = 2L,
    rd_max = NA_integer_,
    tau = 0.5,
    Mm.factor = 1,
    tol = 1e-14,
    maxit = 1e6,
    min_subsample_size = NA_integer_,
    debug_rounds = integer(0),
    fallback = TRUE,
    return_raw_backend = FALSE,
    seed_from = NA_integer_,
    seed_to = NA_integer_,
    seed_list = integer(0),
    n_list = integer(0),
    max_reject_print = 20L,
    output = NULL
  )

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    value <- if (length(kv) > 1) kv[2] else ""

    if (key == "mode") cfg$mode <- value
    if (key == "backend") cfg$backend <- value
    if (key == "n") cfg$n <- as.integer(value)
    if (key == "case") cfg$case <- as.integer(value)
    if (key == "seed") cfg$seed <- as.integer(value)
    if (key == "rd") cfg$rd <- as.integer(value)
    if (key == "rd-max") cfg$rd_max <- as.integer(value)
    if (key == "tau") cfg$tau <- as.numeric(value)
    if (key == "Mm.factor") cfg$Mm.factor <- as.numeric(value)
    if (key == "tol") cfg$tol <- as.numeric(value)
    if (key == "maxit") cfg$maxit <- as.numeric(value)
    if (key == "min_subsample_size") cfg$min_subsample_size <- as.integer(value)
    if (key == "debug-rounds") cfg$debug_rounds <- parse_csv_int(value)
    if (key == "fallback") cfg$fallback <- parse_bool(value, default = TRUE)
    if (key == "return_raw_backend") cfg$return_raw_backend <- parse_bool(value, default = FALSE)
    if (key == "seed-from") cfg$seed_from <- as.integer(value)
    if (key == "seed-to") cfg$seed_to <- as.integer(value)
    if (key == "seed-list") cfg$seed_list <- parse_csv_int(value)
    if (key == "n-list") cfg$n_list <- parse_csv_int(value)
    if (key == "max-reject-print") cfg$max_reject_print <- as.integer(value)
    if (key == "output") cfg$output <- value
  }

  cfg$mode <- match.arg(cfg$mode, c("round", "trace", "scan"))
  cfg$backend <- match.arg(cfg$backend, c("r_ppro", "fortran_ppro"))
  if (is.na(cfg$rd_max)) {
    cfg$rd_max <- cfg$rd
  }
  if (length(cfg$n_list) == 0L) {
    cfg$n_list <- cfg$n
  }
  if (length(cfg$seed_list) == 0L && !is.na(cfg$seed_from) && !is.na(cfg$seed_to)) {
    cfg$seed_list <- seq.int(cfg$seed_from, cfg$seed_to)
  }
  if (cfg$mode != "scan" && length(cfg$debug_rounds) == 0L) {
    if (cfg$mode == "round") {
      cfg$debug_rounds <- cfg$rd
    } else {
      cfg$debug_rounds <- seq_len(cfg$rd_max)
    }
  }
  if (is.null(cfg$output)) {
    if (cfg$mode == "scan") {
      cfg$output <- sprintf(
        "tmp/llqr_seq_ppro_scan_case%d_seed%s_n%s.rds",
        cfg$case,
        if (length(cfg$seed_list) > 0L) paste0(min(cfg$seed_list), "_", max(cfg$seed_list)) else as.character(cfg$seed),
        paste(cfg$n_list, collapse = "_")
      )
    } else {
      cfg$output <- sprintf(
        "tmp/llqr_seq_ppro_%s_case%d_n%d_seed%d_rd%d.rds",
        cfg$mode, cfg$case, cfg$n, cfg$seed, cfg$rd
      )
    }
  }

  cfg
}

estimate_from_H <- function(H, A, y) {
  as.numeric(solve(A[as.integer(H), , drop = FALSE], y[as.integer(H)]))
}

row_set_equal_matrix <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  if (!all(dim(a) == dim(b))) {
    return(FALSE)
  }
  all(vapply(seq_len(nrow(a)), function(i) llqr_set_equal_int(a[i, ], b[i, ]), logical(1)))
}

enrich_round_trace <- function(round_trace, rd, seq_fit, x, y, z, h, tau) {
  if (is.null(round_trace)) {
    return(NULL)
  }

  x <- as.matrix(x)
  y <- as.numeric(y)
  A <- cbind(1, x)
  w <- as.numeric(dnorm((z[rd] - x) / h))
  seq_H <- as.integer(seq_fit$H_seq[rd, ])
  seq_estimate <- estimate_from_H(seq_H, A, y)
  seq_objective <- llqr_full_objective(seq_estimate, A, y, w, tau)

  enriched_attempts <- lapply(round_trace$attempts, function(attempt) {
    if (is.null(attempt$estimate_candidate) || is.null(attempt$H_candidate)) {
      attempt$candidate_class <- NA_character_
      attempt$certification <- NULL
      return(attempt)
    }

    cert_record <- build_llqr_full_cert_record(
      round = rd,
      backend = "debug_trace",
      estimate_candidate = attempt$estimate_candidate,
      H_candidate = attempt$H_candidate,
      A = A,
      y = y,
      w = w,
      tau = tau,
      baseline_estimate = seq_estimate,
      baseline_H = seq_H
    )
    attempt$certification <- cert_record$certification
    attempt$candidate_class <- cert_record$candidate_class
    attempt$fallback_triggered <- cert_record$fallback_triggered
    attempt
  })

  round_trace$attempts <- enriched_attempts
  round_trace$seq_H <- seq_H
  round_trace$seq_estimate <- seq_estimate
  round_trace$seq_objective <- seq_objective
  round_trace$seq_ll_est <- seq_fit$ll_est[rd]
  round_trace
}

build_round_debug <- function(cfg) {
  data <- generate_data(n = cfg$n, case = cfg$case, seed = cfg$seed)
  x <- as.matrix(data$x)
  y <- as.numeric(data$y)
  z <- x[order(x), , drop = FALSE]
  h <- compute_bandwidth(x = x, y = y, tau = cfg$tau)

  min_subsample_size <- if (is.na(cfg$min_subsample_size)) NULL else cfg$min_subsample_size

  seq_fit <- llqr_seq(
    x = x,
    y = y,
    tau = cfg$tau,
    z = z,
    h = h,
    tol = cfg$tol,
    maxit = cfg$maxit,
    track_order = FALSE
  )
  if (cfg$backend == "r_ppro") {
    ppro_fit <- llqr_seq_ppro(
      x = x,
      y = y,
      tau = cfg$tau,
      z = z,
      h = h,
      tol = cfg$tol,
      maxit = cfg$maxit,
      Mm.factor = cfg$Mm.factor,
      case = cfg$case,
      track_order = FALSE,
      min_subsample_size = min_subsample_size,
      store_residual = FALSE,
      debug_trace = TRUE,
      debug_rounds = cfg$debug_rounds
    )
  } else {
    ppro_fit <- llqr_seq_ppro_fortran_wrapper(
      x = as.numeric(x),
      y = y,
      tau = cfg$tau,
      z = as.numeric(z),
      h = h,
      tol = cfg$tol,
      maxit = cfg$maxit,
      Mm.factor = cfg$Mm.factor,
      case = cfg$case,
      bland = TRUE,
      track_order = FALSE,
      fallback = cfg$fallback,
      return_raw_backend = cfg$return_raw_backend,
      debug_trace = TRUE,
      debug_rounds = cfg$debug_rounds
    )
  }

  traced_rounds <- sort(unique(as.integer(cfg$debug_rounds)))
  traced_rounds <- traced_rounds[traced_rounds >= 1L & traced_rounds <= nrow(z)]
  enriched_rounds <- vector("list", length(traced_rounds))
  names(enriched_rounds) <- as.character(traced_rounds)
  for (i in seq_along(traced_rounds)) {
    rd <- traced_rounds[i]
    enriched_rounds[[i]] <- enrich_round_trace(
      round_trace = ppro_fit$round_debug[[rd]],
      rd = rd,
      seq_fit = seq_fit,
      x = x,
      y = y,
      z = z,
      h = h,
      tau = cfg$tau
    )
  }

  list(
    config = cfg,
    x = x,
    y = y,
    z = z,
    h = h,
    seq_fit = seq_fit,
    ppro_fit = ppro_fit,
    traced_rounds = traced_rounds,
    round_debug = enriched_rounds
  )
}

summarize_round_debug <- function(debug_obj) {
  rows <- list()
  idx <- 1L

  for (nm in names(debug_obj$round_debug)) {
    rd <- as.integer(nm)
    entry <- debug_obj$round_debug[[nm]]
    last_attempt <- if (!is.null(entry) && length(entry$attempts) > 0L) {
      entry$attempts[[length(entry$attempts)]]
    } else {
      NULL
    }
    rows[[idx]] <- data.frame(
      round = rd,
      final_status = if (!is.null(entry)) entry$final_status else NA_character_,
      n_attempts = if (!is.null(entry)) length(entry$attempts) else NA_integer_,
      seq_H = if (!is.null(entry)) paste(entry$seq_H, collapse = ",") else NA_character_,
      ppro_H = paste(debug_obj$ppro_fit$H_seq[rd, ], collapse = ","),
      h_set_match_vs_seq = if (!is.null(entry)) llqr_set_equal_int(debug_obj$ppro_fit$H_seq[rd, ], entry$seq_H) else NA,
      ll_gap = debug_obj$ppro_fit$ll_est[rd] - debug_obj$seq_fit$ll_est[rd],
      last_attempt_class = if (!is.null(last_attempt)) last_attempt$candidate_class else NA_character_,
      last_attempt_obj_gap = if (!is.null(last_attempt) && !is.null(last_attempt$certification)) last_attempt$certification$obj_gap_vs_seq else NA_real_,
      last_attempt_H_valid = if (!is.null(last_attempt) && !is.null(last_attempt$certification)) last_attempt$certification$candidate_H_valid else NA,
      last_attempt_in_zero_set = if (!is.null(last_attempt) && !is.null(last_attempt$certification)) last_attempt$certification$candidate_in_zero_set else NA,
      last_attempt_fallback = if (!is.null(last_attempt) && !is.null(last_attempt$fallback_triggered)) last_attempt$fallback_triggered else NA,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  do.call(rbind, rows)
}

print_round_report <- function(debug_obj, rd) {
  entry <- debug_obj$round_debug[[as.character(rd)]]
  if (is.null(entry)) {
    cat(sprintf("Round %d was not traced.\n", rd))
    return(invisible(NULL))
  }

  cat(sprintf("Round %d\n", rd))
  cat(sprintf("backend: %s\n", debug_obj$config$backend))
  cat(sprintf("final_status: %s\n", entry$final_status))
  cat(sprintf("seq_H: %s\n", paste(entry$seq_H, collapse = ",")))
  cat(sprintf("ppro_H: %s\n", paste(debug_obj$ppro_fit$H_seq[rd, ], collapse = ",")))
  cat(sprintf("ll_gap_vs_seq: %.13f\n", debug_obj$ppro_fit$ll_est[rd] - debug_obj$seq_fit$ll_est[rd]))
  cat(sprintf("attempt_count: %d\n", length(entry$attempts)))

  for (attempt in entry$attempts) {
    cat(sprintf(
      "  attempt=%d action=%s local_accept=%s bad_signs=%s ms=%s M=%.12f\n",
      attempt$attempt_id,
      attempt$action,
      if (!is.null(attempt$local_accept)) attempt$local_accept else NA,
      if (!is.null(attempt$bad_signs)) attempt$bad_signs else NA,
      if (!is.null(attempt$ms)) attempt$ms else NA,
      if (!is.null(attempt$M)) attempt$M else NA_real_
    ))
    if (!is.null(attempt$H_candidate)) {
      cat(sprintf("    H_prev=%s | H_candidate=%s\n",
                  paste(attempt$H_prev, collapse = ","),
                  paste(attempt$H_candidate, collapse = ",")))
    }
    if (!is.null(attempt$certification)) {
      cert <- attempt$certification
      cat(sprintf(
        "    class=%s obj_gap=%.13f H_valid=%s in_zero_set=%s h_match=%s\n",
        attempt$candidate_class,
        cert$obj_gap_vs_seq,
        cert$candidate_H_valid,
        cert$candidate_in_zero_set,
        cert$h_set_match_vs_seq
      ))
    }
  }
}

scan_one_case <- function(seed, n, cfg) {
  data <- generate_data(n = n, case = cfg$case, seed = seed)
  x <- as.matrix(data$x)
  y <- as.numeric(data$y)
  z <- x[order(x), , drop = FALSE]
  h <- compute_bandwidth(x = x, y = y, tau = cfg$tau)

  min_subsample_size <- if (is.na(cfg$min_subsample_size)) NULL else cfg$min_subsample_size

  seq_fit <- llqr_seq(
    x = x,
    y = y,
    tau = cfg$tau,
    z = z,
    h = h,
    tol = cfg$tol,
    maxit = cfg$maxit,
    track_order = FALSE
  )
  if (cfg$backend == "r_ppro") {
    ppro_fit <- llqr_seq_ppro(
      x = x,
      y = y,
      tau = cfg$tau,
      z = z,
      h = h,
      tol = cfg$tol,
      maxit = cfg$maxit,
      Mm.factor = cfg$Mm.factor,
      case = cfg$case,
      track_order = FALSE,
      min_subsample_size = min_subsample_size,
      store_residual = FALSE
    )
  } else {
    ppro_fit <- llqr_seq_ppro_fortran_wrapper(
      x = as.numeric(x),
      y = y,
      tau = cfg$tau,
      z = as.numeric(z),
      h = h,
      tol = cfg$tol,
      maxit = cfg$maxit,
      Mm.factor = cfg$Mm.factor,
      case = cfg$case,
      bland = TRUE,
      track_order = FALSE,
      fallback = cfg$fallback,
      return_raw_backend = cfg$return_raw_backend
    )
  }

  reject_rounds <- which(vapply(ppro_fit$acceptance_diagnostics, Negate(is.null), logical(1)))
  first_reject_round <- if (length(reject_rounds) > 0L) reject_rounds[1] else NA_integer_
  candidate_class <- NA_character_
  obj_gap <- NA_real_
  candidate_H_valid <- NA
  candidate_in_zero_set <- NA
  bad_H <- NA_character_
  seq_H <- NA_character_
  fallback_triggered <- isTRUE(ppro_fit$fallback_triggered) || length(reject_rounds) > 0L

  if (!is.na(first_reject_round)) {
    diag_entry <- ppro_fit$acceptance_diagnostics[[first_reject_round]]
    if (!is.null(diag_entry$certification)) {
      candidate_class <- diag_entry$candidate_class
      obj_gap <- if (!is.null(diag_entry$obj_gap_vs_seq)) diag_entry$obj_gap_vs_seq else NA_real_
      candidate_H_valid <- diag_entry$certification$candidate_H_valid
      candidate_in_zero_set <- diag_entry$certification$candidate_in_zero_set
      bad_H <- paste(diag_entry$H_candidate, collapse = ",")
      seq_H <- if (!is.null(diag_entry$seq_H)) paste(diag_entry$seq_H, collapse = ",") else NA_character_
    }
  }

  data.frame(
    backend = cfg$backend,
    case = cfg$case,
    n = n,
    seed = seed,
    reject_round_count = length(reject_rounds),
    fallback_triggered = fallback_triggered,
    first_bad_round = first_reject_round,
    first_reject_round = first_reject_round,
    rd = first_reject_round,
    H_candidate = bad_H,
    seq_H = seq_H,
    candidate_class = candidate_class,
    obj_gap_vs_seq = obj_gap,
    candidate_H_valid = candidate_H_valid,
    candidate_in_zero_set = candidate_in_zero_set,
    backend_ierr = if (!is.null(ppro_fit$backend_ierr)) ppro_fit$backend_ierr else NA_integer_,
    returned_backend = if (!is.null(ppro_fit$returned_backend)) ppro_fit$returned_backend else NA_character_,
    cert_fail_detected = if (!is.null(ppro_fit$cert_fail_detected)) ppro_fit$cert_fail_detected else FALSE,
    ll_match = isTRUE(all.equal(ppro_fit$ll_est, seq_fit$ll_est, tolerance = 1e-10)),
    H_set_match = row_set_equal_matrix(ppro_fit$H_seq, seq_fit$H_seq),
    H_exact_match = isTRUE(all.equal(ppro_fit$H_seq, seq_fit$H_seq, tolerance = 0)),
    bad_H = bad_H,
    stringsAsFactors = FALSE
  )
}

run_scan <- function(cfg) {
  if (length(cfg$seed_list) == 0L) {
    stop("scan mode requires --seed-list or --seed-from/--seed-to.")
  }

  rows <- list()
  idx <- 1L
  for (n in cfg$n_list) {
    for (seed in cfg$seed_list) {
      rows[[idx]] <- scan_one_case(seed = seed, n = n, cfg = cfg)
      idx <- idx + 1L
    }
  }

  results <- do.call(rbind, rows)
  results
}

print_scan_report <- function(results, cfg) {
  cat(sprintf("Scan summary: backend=%s case=%d\n", cfg$backend, cfg$case))
  for (n in sort(unique(results$n))) {
    subset_n <- results[results$n == n, , drop = FALSE]
    cat(sprintf(
      "  n=%d | seeds=%d | reject_cases=%d | all_ll_match=%s | all_H_set_match=%s\n",
      n,
      nrow(subset_n),
      sum(subset_n$reject_round_count > 0L),
      all(subset_n$ll_match),
      all(subset_n$H_set_match)
    ))
  }

  rejected <- results[results$reject_round_count > 0L, , drop = FALSE]
  if (nrow(rejected) > 0L) {
    rejected <- rejected[order(rejected$n, rejected$seed), , drop = FALSE]
    cat("Rejected candidates intercepted:\n")
    print(utils::head(rejected, cfg$max_reject_print))
  } else {
    cat("No rejected candidates intercepted in this scan.\n")
  }
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(dirname(cfg$output), recursive = TRUE, showWarnings = FALSE)

  if (cfg$mode == "scan") {
    results <- run_scan(cfg)
    saveRDS(results, cfg$output)
    print_scan_report(results, cfg)
    cat(sprintf("Saved scan results: %s\n", cfg$output))
    return(invisible(NULL))
  }

  debug_obj <- build_round_debug(cfg)
  round_summary <- summarize_round_debug(debug_obj)
  output_obj <- list(
    config = cfg,
    round_summary = round_summary,
    round_debug = debug_obj$round_debug
  )
  saveRDS(output_obj, cfg$output)

  if (cfg$mode == "round") {
    print_round_report(debug_obj, cfg$rd)
  } else {
    print(round_summary)
  }
  cat(sprintf("Saved debug object: %s\n", cfg$output))
}

main()
