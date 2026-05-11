#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

source("R/tvcqr_functions.R")

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
    seed = 2027L,
    rd = 2L,
    rd_max = NA_integer_,
    tau = 0.5,
    h = NA_real_,
    h.factor = 1,
    Mm.factor = 1e-4,
    tol = 1e-14,
    maxit = 1e6,
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
    if (key == "h") cfg$h <- as.numeric(value)
    if (key == "h.factor") cfg$h.factor <- as.numeric(value)
    if (key == "Mm.factor") cfg$Mm.factor <- as.numeric(value)
    if (key == "tol") cfg$tol <- as.numeric(value)
    if (key == "maxit") cfg$maxit <- as.numeric(value)
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
  if (is.null(cfg$output)) {
    if (cfg$mode == "scan") {
      cfg$output <- sprintf(
        "tmp/tvcqr_seq_ppro_scan_case%d_seed%s_n%s.rds",
        cfg$case,
        if (length(cfg$seed_list) > 0L) paste0(min(cfg$seed_list), "_", max(cfg$seed_list)) else as.character(cfg$seed),
        paste(cfg$n_list, collapse = "_")
      )
    } else {
      cfg$output <- sprintf(
        "tmp/tvcqr_seq_ppro_%s_case%d_n%d_seed%d_rd%d.rds",
        cfg$mode, cfg$case, cfg$n, cfg$seed, cfg$rd
      )
    }
  }

  cfg
}

tvcqr_bandwidth <- function(m, h = NULL, h.factor = 1) {
  if (!is.null(h) && !is.na(h)) {
    return(as.numeric(h))
  }
  as.numeric(m^(-0.2) * h.factor)
}

tvcqr_design_matrix <- function(x) {
  x <- as.matrix(x)
  m <- nrow(x)
  base <- cbind(1, x)
  cbind(base, base * ((1:m) / m))
}

tvcqr_quantile_loss <- function(r, tau) {
  r <- as.numeric(r)
  r * (tau - (r < 0))
}

tvcqr_full_objective <- function(estimate, A, y, w, tau) {
  r <- as.numeric(y - A %*% as.numeric(estimate))
  sum(as.numeric(w) * tvcqr_quantile_loss(r, tau))
}

tvcqr_set_equal_int <- function(x, y) {
  x <- sort(unique(as.integer(x)))
  y <- sort(unique(as.integer(y)))
  identical(x, y)
}

tvcqr_estimate_from_H <- function(H, A, y) {
  H <- as.integer(H)
  fit <- try(solve(A[H, , drop = FALSE], as.numeric(y)[H]), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(rep(NA_real_, ncol(A)))
  }
  as.numeric(fit)
}

tvcqr_weights_at_round <- function(rd, m, h) {
  time_index <- (1:m) / m
  0.75 * (1 - ((rd / m - time_index) / h)^2) * (abs(rd / m - time_index) <= h)
}

row_set_equal_matrix <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  if (!all(dim(a) == dim(b))) {
    return(FALSE)
  }
  all(vapply(seq_len(nrow(a)), function(i) tvcqr_set_equal_int(a[i, ], b[i, ]), logical(1)))
}

build_tvcqr_full_cert_record <- function(round, backend = "debug_trace", estimate_candidate, H_candidate,
                                         A, y, w, tau,
                                         baseline_estimate = NULL,
                                         baseline_H = NULL,
                                         residual_tol = 1e-8,
                                         rank_tol = 1e-10,
                                         obj_tol = 1e-10) {
  estimate_candidate <- as.numeric(estimate_candidate)
  H_candidate <- as.integer(H_candidate)
  A <- as.matrix(A)
  y <- as.numeric(y)
  w <- as.numeric(w)
  p <- ncol(A)
  r <- as.numeric(y - A %*% estimate_candidate)
  zero_idx <- which(abs(r) <= residual_tol)
  candidate_H_in_range <- length(H_candidate) == p &&
    !anyNA(H_candidate) &&
    !any(H_candidate < 1L | H_candidate > nrow(A)) &&
    !anyDuplicated(H_candidate)
  candidate_rank_ok <- candidate_H_in_range &&
    (qr(A[H_candidate, , drop = FALSE], tol = rank_tol)$rank == p)
  candidate_zero_ok <- candidate_H_in_range &&
    all(abs(r[H_candidate]) <= residual_tol)
  candidate_H_valid <- candidate_rank_ok && candidate_zero_ok
  candidate_in_zero_set <- candidate_H_in_range && all(H_candidate %in% zero_idx)
  obj_full <- tvcqr_full_objective(estimate_candidate, A, y, w, tau)
  obj_full_seq <- if (!is.null(baseline_estimate)) {
    tvcqr_full_objective(baseline_estimate, A, y, w, tau)
  } else {
    NA_real_
  }
  obj_gap_vs_seq <- if (is.finite(obj_full_seq)) obj_full - obj_full_seq else NA_real_
  cert <- list(
    cert_ok = candidate_H_valid && candidate_in_zero_set,
    obj_full = obj_full,
    obj_full_seq = obj_full_seq,
    obj_gap_vs_seq = obj_gap_vs_seq,
    obj_ok_vs_seq = if (is.finite(obj_gap_vs_seq)) obj_gap_vs_seq <= obj_tol else NA,
    H_recovered_from_full = zero_idx,
    candidate_H_valid = candidate_H_valid,
    candidate_in_zero_set = candidate_in_zero_set,
    h_set_match_vs_seq = if (!is.null(baseline_H)) tvcqr_set_equal_int(H_candidate, baseline_H) else NA,
    residual = r
  )

  list(
    round = as.integer(round),
    estimate_candidate = as.numeric(estimate_candidate),
    H_candidate = as.integer(H_candidate),
    seq_H = if (!is.null(baseline_H)) as.integer(baseline_H) else NULL,
    certification = cert,
    candidate_class = if (!isTRUE(candidate_H_valid)) {
      "invalid_H"
    } else if (!isTRUE(candidate_in_zero_set)) {
      "H_not_in_zero_set"
    } else {
      "certified"
    },
    fallback_triggered = !isTRUE(cert$cert_ok),
    obj_full = obj_full,
    obj_full_seq = obj_full_seq,
    obj_gap_vs_seq = obj_gap_vs_seq,
    H_recovered_from_full = cert$H_recovered_from_full,
    h_set_match_vs_seq = cert$h_set_match_vs_seq
  )
}

enrich_round_trace <- function(round_trace, rd, seq_fit, x, y, h, tau) {
  if (is.null(round_trace)) {
    return(NULL)
  }

  x <- as.matrix(x)
  y <- as.numeric(y)
  A <- tvcqr_design_matrix(x)
  w <- as.numeric(tvcqr_weights_at_round(rd = rd, m = nrow(x), h = h))
  seq_H <- as.integer(seq_fit$H_seq[rd, ])
  seq_estimate <- tvcqr_estimate_from_H(seq_H, A, y)
  seq_objective <- tvcqr_full_objective(seq_estimate, A, y, w, tau)

  enriched_attempts <- lapply(round_trace$attempts, function(attempt) {
    if (is.null(attempt$estimate_candidate) || is.null(attempt$H_candidate)) {
      attempt$candidate_class <- NA_character_
      attempt$certification <- NULL
      return(attempt)
    }

    cert_record <- build_tvcqr_full_cert_record(
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
  round_trace$seq_theta <- seq_fit$theta_ll_est[rd, ]
  round_trace
}

build_round_debug <- function(cfg) {
  data <- generate_ts(n = cfg$n, case = cfg$case, seed = cfg$seed)
  x <- as.matrix(data$x)
  y <- as.numeric(data$y)
  h <- tvcqr_bandwidth(nrow(x), h = if (is.na(cfg$h)) NULL else cfg$h, h.factor = cfg$h.factor)

  seq_fit <- tvcqr_seq(
    x = x,
    y = y,
    tau = cfg$tau,
    h = h,
    h.factor = cfg$h.factor,
    tol = cfg$tol,
    maxit = cfg$maxit
  )
  if (cfg$backend == "r_ppro") {
    ppro_fit <- tvcqr_seq_ppro(
      x = x,
      y = y,
      tau = cfg$tau,
      h = h,
      h.factor = cfg$h.factor,
      tol = cfg$tol,
      maxit = cfg$maxit,
      Mm.factor = cfg$Mm.factor,
      store_residual = FALSE,
      debug_trace = TRUE,
      debug_rounds = if (cfg$mode == "round") cfg$rd else seq_len(cfg$rd_max)
    )
  } else {
    ppro_fit <- tvcqr_seq_ppro_fortran_wrapper(
      x = x,
      y = y,
      tau = cfg$tau,
      h = h,
      h.factor = cfg$h.factor,
      tol = cfg$tol,
      maxit = cfg$maxit,
      bland = FALSE,
      Mm.factor = cfg$Mm.factor,
      store_residual = FALSE,
      fallback = cfg$fallback
    )
  }

  traced_rounds <- if (cfg$mode == "round") cfg$rd else seq_len(cfg$rd_max)
  traced_rounds <- traced_rounds[traced_rounds >= 1L & traced_rounds <= nrow(x)]
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
      h = h,
      tau = cfg$tau
    )
  }

  list(
    config = cfg,
    x = x,
    y = y,
    h = h,
    seq_fit = seq_fit,
    ppro_fit = ppro_fit,
    round_debug = enriched_rounds
  )
}

summarize_round_debug <- function(debug_obj) {
  rows <- list()
  idx <- 1L
  for (nm in names(debug_obj$round_debug)) {
    rd <- as.integer(nm)
    entry <- debug_obj$round_debug[[nm]]
    last_attempt <- if (!is.null(entry) && length(entry$attempts) > 0L) entry$attempts[[length(entry$attempts)]] else NULL
    rows[[idx]] <- data.frame(
      round = rd,
      final_status = if (!is.null(entry)) entry$final_status else NA_character_,
      n_attempts = if (!is.null(entry)) length(entry$attempts) else NA_integer_,
      seq_H = if (!is.null(entry)) paste(entry$seq_H, collapse = ",") else NA_character_,
      ppro_H = paste(debug_obj$ppro_fit$H_seq[rd, ], collapse = ","),
      h_set_match_vs_seq = if (!is.null(entry)) tvcqr_set_equal_int(debug_obj$ppro_fit$H_seq[rd, ], entry$seq_H) else NA,
      theta_max_abs_gap = max(abs(debug_obj$ppro_fit$theta_ll_est[rd, ] - debug_obj$seq_fit$theta_ll_est[rd, ])),
      last_attempt_class = if (!is.null(last_attempt)) last_attempt$candidate_class else NA_character_,
      last_attempt_obj_gap = if (!is.null(last_attempt) && !is.null(last_attempt$certification)) last_attempt$certification$obj_gap_vs_seq else NA_real_,
      last_attempt_H_valid = if (!is.null(last_attempt) && !is.null(last_attempt$certification)) last_attempt$certification$candidate_H_valid else NA,
      last_attempt_in_zero_set = if (!is.null(last_attempt) && !is.null(last_attempt$certification)) last_attempt$certification$candidate_in_zero_set else NA,
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
  cat(sprintf("theta_max_abs_gap_vs_seq: %.13f\n", max(abs(debug_obj$ppro_fit$theta_ll_est[rd, ] - debug_obj$seq_fit$theta_ll_est[rd, ]))))
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
  data <- generate_ts(n = n, case = cfg$case, seed = seed)
  x <- as.matrix(data$x)
  y <- as.numeric(data$y)
  h <- tvcqr_bandwidth(nrow(x), h = if (is.na(cfg$h)) NULL else cfg$h, h.factor = cfg$h.factor)

  seq_fit <- tvcqr_seq(
    x = x,
    y = y,
    tau = cfg$tau,
    h = h,
    h.factor = cfg$h.factor,
    tol = cfg$tol,
    maxit = cfg$maxit
  )
  if (cfg$backend == "r_ppro") {
    ppro_fit <- tvcqr_seq_ppro(
      x = x,
      y = y,
      tau = cfg$tau,
      h = h,
      h.factor = cfg$h.factor,
      tol = cfg$tol,
      maxit = cfg$maxit,
      Mm.factor = cfg$Mm.factor,
      store_residual = FALSE
    )
  } else {
    ppro_fit <- tvcqr_seq_ppro_fortran_wrapper(
      x = x,
      y = y,
      tau = cfg$tau,
      h = h,
      h.factor = cfg$h.factor,
      tol = cfg$tol,
      maxit = cfg$maxit,
      bland = FALSE,
      Mm.factor = cfg$Mm.factor,
      store_residual = FALSE,
      fallback = cfg$fallback
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

  if (!is.na(first_reject_round)) {
    A <- tvcqr_design_matrix(x)
    w <- as.numeric(tvcqr_weights_at_round(first_reject_round, nrow(x), h))
    diag_entry <- ppro_fit$acceptance_diagnostics[[first_reject_round]]
    seq_H_row <- as.integer(seq_fit$H_seq[first_reject_round, ])
    seq_estimate <- tvcqr_estimate_from_H(seq_H_row, A, y)
    cert_record <- build_tvcqr_full_cert_record(
      round = first_reject_round,
      backend = cfg$backend,
      estimate_candidate = diag_entry$estimate_candidate,
      H_candidate = diag_entry$H_candidate,
      A = A,
      y = y,
      w = w,
      tau = cfg$tau,
      baseline_estimate = seq_estimate,
      baseline_H = seq_H_row
    )
    candidate_class <- cert_record$candidate_class
    obj_gap <- cert_record$obj_gap_vs_seq
    candidate_H_valid <- cert_record$certification$candidate_H_valid
    candidate_in_zero_set <- cert_record$certification$candidate_in_zero_set
    bad_H <- paste(diag_entry$H_candidate, collapse = ",")
    seq_H <- paste(seq_H_row, collapse = ",")
  }

  data.frame(
    backend = cfg$backend,
    case = cfg$case,
    n = n,
    seed = seed,
    reject_round_count = length(reject_rounds),
    fallback_triggered = isTRUE(ppro_fit$fallback_triggered) || length(reject_rounds) > 0L,
    first_bad_round = first_reject_round,
    candidate_class = candidate_class,
    obj_gap_vs_seq = obj_gap,
    candidate_H_valid = candidate_H_valid,
    candidate_in_zero_set = candidate_in_zero_set,
    H_candidate = bad_H,
    seq_H = seq_H,
    backend_ierr = if (!is.null(ppro_fit$backend_ierr)) ppro_fit$backend_ierr else NA_integer_,
    returned_backend = if (!is.null(ppro_fit$returned_backend)) ppro_fit$returned_backend else NA_character_,
    cert_fail_detected = if (!is.null(ppro_fit$cert_fail_detected)) ppro_fit$cert_fail_detected else FALSE,
    theta_max_abs_gap = max(abs(ppro_fit$theta_ll_est - seq_fit$theta_ll_est)),
    theta_match = isTRUE(all.equal(ppro_fit$theta_ll_est, seq_fit$theta_ll_est, tolerance = 1e-7)),
    H_set_match = row_set_equal_matrix(ppro_fit$H_seq, seq_fit$H_seq),
    H_exact_match = isTRUE(all.equal(ppro_fit$H_seq, seq_fit$H_seq, tolerance = 0)),
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
  do.call(rbind, rows)
}

print_scan_report <- function(results, cfg) {
  cat(sprintf("Scan summary: backend=%s case=%d\n", cfg$backend, cfg$case))
  for (n in sort(unique(results$n))) {
    subset_n <- results[results$n == n, , drop = FALSE]
    cat(sprintf(
      "  n=%d | seeds=%d | reject_cases=%d | fallback_cases=%d | all_theta_match=%s | all_H_set_match=%s\n",
      n,
      nrow(subset_n),
      sum(subset_n$reject_round_count > 0L),
      sum(subset_n$fallback_triggered),
      all(subset_n$theta_match),
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
