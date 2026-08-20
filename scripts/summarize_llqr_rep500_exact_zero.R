#!/usr/bin/env Rscript

parse_args <- function(args) {
  cfg <- list(
    base_dir = Sys.getenv("FASTQR_LLQR_BASE_DIR", unset = "data/llqr_simu_results"),
    output_dir = "results/table/llqr_rep500_exact_zero",
    rep = 500L,
    seed_base = 2025L,
    method = "llqr_seq_ppro_fortran_1",
    method_summary = NULL,
    paper_method_summary = "paper/data/rep500_seed2025_method_summary.csv",
    paper_screening = "paper/data/rep500_seed2025_screening_tau05.csv"
  )
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1L]
    value <- if (length(kv) > 1L) kv[2L] else ""
    if (key == "base-dir") cfg$base_dir <- value
    if (key == "output-dir") cfg$output_dir <- value
    if (key == "rep") cfg$rep <- as.integer(value)
    if (key == "seed-base") cfg$seed_base <- as.integer(value)
    if (key == "method") cfg$method <- value
    if (key == "method-summary") cfg$method_summary <- value
    if (key == "paper-method-summary") cfg$paper_method_summary <- value
    if (key == "paper-screening") cfg$paper_screening <- value
  }
  if (is.null(cfg$method_summary) || !nzchar(cfg$method_summary)) {
    cfg$method_summary <- file.path(cfg$output_dir, "method_summary.csv")
  }
  cfg
}

load_results <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists("results", envir = env)) stop("No results object in ", path)
  get("results", envir = env)
}

scalar_or <- function(x, default) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) default else x[[1L]]
}

bind_nonempty <- function(xs) {
  xs <- xs[vapply(xs, function(x) !is.null(x) && nrow(x) > 0L, logical(1L))]
  if (!length(xs)) data.frame() else do.call(rbind, xs)
}

aggregate_exact <- function(rows, case, tau, n, expected_rep, metadata) {
  audited <- nrow(rows)
  meta <- metadata[seq_len(min(length(metadata), expected_rep))]
  backend <- vapply(meta, function(x) as.character(scalar_or(x$returned_backend, NA_character_)), character(1L))
  fallback <- vapply(meta, function(x) isTRUE(x$fallback_triggered), logical(1L))
  ierr <- vapply(meta, function(x) as.integer(scalar_or(x$backend_ierr, NA_integer_)), integer(1L))
  independent <- vapply(meta, function(x) as.integer(scalar_or(x$independent_init_count, 0L)), integer(1L))
  recovery <- vapply(meta, function(x) as.integer(scalar_or(x$full_active_recovery_count, 0L)), integer(1L))

  data.frame(
    case = case, tau = tau, n = n,
    expected_replications = expected_rep,
    n_replications_audited = audited,
    expected_next_window_points = expected_rep * (n - 1L),
    n_next_window_points_audited = sum(rows$n_next_window_points),
    max_n_zero_exact_all = max(rows$max_n_zero_exact_all),
    max_n_zero_exact_current = max(rows$max_n_zero_exact_current),
    max_n_zero_exact_next_window = max(rows$max_n_zero_exact_next_window),
    n_literal_violation_points = sum(rows$n_literal_violation_next_window),
    n_replications_with_literal_violation = sum(rows$replication_literal_violation_next_window),
    max_n_zero_c10_next_window = max(rows$max_n_zero_c10_next_window),
    max_n_zero_c100_next_window = max(rows$max_n_zero_c100_next_window),
    max_n_zero_c1000_next_window = max(rows$max_n_zero_c1000_next_window),
    n_nearzero_exceeds_q_next_c10 = sum(rows$n_nearzero_exceeds_q_next_c10),
    n_nearzero_exceeds_q_next_c100 = sum(rows$n_nearzero_exceeds_q_next_c100),
    n_nearzero_exceeds_q_next_c1000 = sum(rows$n_nearzero_exceeds_q_next_c1000),
    audit_time_mean_sec = mean(rows$audit_elapsed_sec),
    audit_time_max_sec = max(rows$audit_elapsed_sec),
    backend_ppro_count = sum(backend == "ppro", na.rm = TRUE),
    fallback_count = sum(fallback, na.rm = TRUE),
    backend_ierr_nonzero_count = sum(ierr != 0L, na.rm = TRUE),
    independent_init_total = sum(independent, na.rm = TRUE),
    full_active_recovery_total = sum(recovery, na.rm = TRUE),
    audit_complete = audited == expected_rep && sum(rows$n_next_window_points) == expected_rep * (n - 1L),
    stringsAsFactors = FALSE
  )
}

screening_summary <- function(results, method, case, tau, n) {
  metadata <- results$method_metadata[[method]]
  if (is.null(metadata)) return(data.frame())
  first_pass <- fr <- repairs <- max_s <- rep(NA_real_, length(metadata))
  for (i in seq_along(metadata)) {
    item <- metadata[[i]]
    if (!is.list(item)) next
    first_pass[i] <- as.numeric(scalar_or(item$first_pass_rate, NA_real_))
    rc <- as.integer(item$repair_count)
    fs <- as.integer(item$first_n_sub)
    if (length(rc) >= 2L) {
      rc <- rc[-1L]
      fr[i] <- sum(rc > 0L, na.rm = TRUE)
      repairs[i] <- sum(rc, na.rm = TRUE)
    }
    if (length(fs) >= 2L) max_s[i] <- max(fs[-1L], na.rm = TRUE)
  }
  data.frame(
    case = case, tau = tau, n = n,
    ok_reps = sum(is.finite(first_pass)),
    gamma_mean = 0.1 * sqrt(log(n)) * n^(-2 / 5) * log(log(n)),
    first_pass_prop_mean = mean(first_pass, na.rm = TRUE),
    Fr_mean = mean(fr, na.rm = TRUE),
    total_repair_mean = mean(repairs, na.rm = TRUE),
    mean_max_Sj = mean(max_s, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

compare_methods <- function(new_summary, paper_summary, method) {
  new_screen <- new_summary[new_summary$model == "llqr" & new_summary$raw_method == method, ]
  new_base <- new_summary[new_summary$model == "llqr" & new_summary$raw_method == "llqr", ]
  old_screen <- paper_summary[paper_summary$model == "llqr" & paper_summary$raw_method == method, ]
  old_base <- paper_summary[paper_summary$model == "llqr" & paper_summary$raw_method == "llqr", ]
  key <- c("case", "tau", "n")
  out <- merge(new_screen, old_screen, by = key, suffixes = c("_new", "_paper"))
  nb <- new_base[, c(key, "time_mean_sec")]
  names(nb)[4L] <- "baseline_time_mean_sec_new"
  ob <- old_base[, c(key, "time_mean_sec")]
  names(ob)[4L] <- "baseline_time_mean_sec_paper"
  out <- merge(merge(out, nb, by = key), ob, by = key)
  out$runtime_reduction_new <- 100 * (1 - out$time_mean_sec_new / out$baseline_time_mean_sec_new)
  out$runtime_reduction_paper <- 100 * (1 - out$time_mean_sec_paper / out$baseline_time_mean_sec_paper)
  out$runtime_reduction_abs_diff <- abs(out$runtime_reduction_new - out$runtime_reduction_paper)
  out$discrepancy_warning_threshold <- pmax(1e-12, 100 * out$max_average_relative_bias_paper)
  out$discrepancy_hard_pass <- out$max_average_relative_bias_new <= 1e-10
  out$discrepancy_warning <- out$max_average_relative_bias_new > out$discrepancy_warning_threshold
  out$runtime_warning <- out$runtime_reduction_new <= 0 | out$runtime_reduction_abs_diff > 10
  out
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  cases <- 1:2
  taus <- c(0.2, 0.5, 0.8)
  ns <- c(1000L, 2000L, 5000L, 10000L)
  config_rows <- list()
  rep_rows <- list()
  witness_rows <- list()
  screening_rows <- list()
  pos <- 1L

  for (case in cases) for (tau in taus) for (n in ns) {
    tag <- sprintf("tau%02d", as.integer(round(100 * tau)))
    path <- file.path(cfg$base_dir, sprintf("case%d_%s_n%d_rep%d.RData", case, tag, n, cfg$rep))
    if (!file.exists(path)) next
    results <- load_results(path)
    audits <- results$exact_zero_audit[[cfg$method]]
    summaries <- bind_nonempty(lapply(audits, function(x) if (is.list(x)) x$summary else NULL))
    witnesses <- bind_nonempty(lapply(audits, function(x) if (is.list(x)) x$literal_violation_witnesses else NULL))
    if (!nrow(summaries)) next
    summaries$result_file <- path
    rep_rows[[pos]] <- summaries
    if (nrow(witnesses)) {
      witnesses$result_file <- path
      witness_rows[[pos]] <- witnesses
    }
    config_rows[[pos]] <- aggregate_exact(
      summaries, case, tau, n, cfg$rep, results$method_metadata[[cfg$method]]
    )
    if (case == 2L && tau == 0.5) {
      screening_rows[[pos]] <- screening_summary(results, cfg$method, case, tau, n)
    }
    pos <- pos + 1L
  }

  config_df <- bind_nonempty(config_rows)
  rep_df <- bind_nonempty(rep_rows)
  witness_df <- bind_nonempty(witness_rows)
  screening_df <- bind_nonempty(screening_rows)
  write.csv(config_df, file.path(cfg$output_dir, "exact_zero_config_summary.csv"), row.names = FALSE)
  write.csv(rep_df, file.path(cfg$output_dir, "exact_zero_replication_summary.csv"), row.names = FALSE)
  write.csv(witness_df, file.path(cfg$output_dir, "exact_zero_literal_violation_witnesses.csv"), row.names = FALSE)

  if (nrow(screening_df) && file.exists(cfg$paper_screening)) {
    old <- read.csv(cfg$paper_screening, check.names = FALSE)
    old <- old[old$model == "llqr" & old$case == 2 & old$tau == 0.5, ]
    cmp <- merge(screening_df, old[, c("case", "tau", "n", "gamma_mean", "first_pass_prop_mean",
                                      "Fr_mean", "total_repair_mean", "mean_max_Sj")],
                 by = c("case", "tau", "n"), suffixes = c("_new", "_paper"))
    cmp$gamma_config_error <- abs(cmp$gamma_mean_new - cmp$gamma_mean_paper) > 1e-8 * abs(cmp$gamma_mean_paper)
    cmp$first_pass_warning <- abs(cmp$first_pass_prop_mean_new - cmp$first_pass_prop_mean_paper) > 0.02
    cmp$retained_size_warning <- abs(cmp$mean_max_Sj_new / cmp$mean_max_Sj_paper - 1) > 0.10
    cmp$Fr_warning <- abs(cmp$Fr_mean_new - cmp$Fr_mean_paper) > pmax(1, 0.2 * cmp$Fr_mean_paper)
    cmp$total_repair_warning <- abs(cmp$total_repair_mean_new - cmp$total_repair_mean_paper) >
      pmax(1, 0.2 * cmp$total_repair_mean_paper)
    write.csv(cmp, file.path(cfg$output_dir, "screening_paper_comparison.csv"), row.names = FALSE)
  }

  if (file.exists(cfg$method_summary) && file.exists(cfg$paper_method_summary)) {
    comparison <- compare_methods(
      read.csv(cfg$method_summary, check.names = FALSE),
      read.csv(cfg$paper_method_summary, check.names = FALSE),
      cfg$method
    )
    write.csv(comparison, file.path(cfg$output_dir, "method_paper_comparison.csv"), row.names = FALSE)
  }

  cat("Exact-zero configurations:", nrow(config_df), "\n")
  cat("Audited replications:", sum(config_df$n_replications_audited), "\n")
  cat("Literal violation points:", sum(config_df$n_literal_violation_points), "\n")
  cat("Output directory:", cfg$output_dir, "\n")
}

main()
