#!/usr/bin/env Rscript

parse_args <- function(args) {
  cfg <- list(
    base_dir = Sys.getenv("FASTQR_TVCQR_BASE_DIR", unset = "data/tvcqr_simu_results"),
    output_dir = "results/table/tvcqr_exact_zero",
    rep = 100L,
    cases = 1:2,
    taus = 0.5,
    ns = 500L,
    method = "tvcqr_seq_ppro_fortran_1"
  )
  csv_int <- function(x) as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1L]]))
  csv_num <- function(x) as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1L]]))
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- kv[1L]
    value <- if (length(kv) > 1L) kv[2L] else ""
    if (key == "base-dir") cfg$base_dir <- value
    if (key == "output-dir") cfg$output_dir <- value
    if (key == "rep") cfg$rep <- as.integer(value)
    if (key == "cases") cfg$cases <- csv_int(value)
    if (key == "taus") cfg$taus <- csv_num(value)
    if (key == "ns") cfg$ns <- csv_int(value)
    if (key == "method") cfg$method <- value
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

aggregate_config <- function(rows, results, method, case, tau, n, expected_rep) {
  metadata <- results$method_metadata[[method]]
  backend <- vapply(metadata, function(x) {
    as.character(scalar_or(x$returned_backend, NA_character_))
  }, character(1L))
  fallback <- vapply(metadata, function(x) isTRUE(x$fallback_triggered), logical(1L))
  ierr <- vapply(metadata, function(x) {
    as.integer(scalar_or(x$backend_ierr, NA_integer_))
  }, integer(1L))
  independent <- vapply(metadata, function(x) {
    as.integer(scalar_or(x$independent_init_count, 0L))
  }, integer(1L))
  recovery <- vapply(metadata, function(x) {
    as.integer(scalar_or(x$full_active_recovery_count, 0L))
  }, integer(1L))

  data.frame(
    case = case, tau = tau, n = n,
    expected_replications = expected_rep,
    n_replications_audited = nrow(rows),
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
    audit_complete = nrow(rows) == expected_rep &&
      sum(rows$n_next_window_points) == expected_rep * (n - 1L),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  config_rows <- rep_rows <- witness_rows <- list()
  pos <- 1L
  for (case in cfg$cases) for (tau in cfg$taus) for (n in cfg$ns) {
    tag <- sprintf("tau%02d", as.integer(round(100 * tau)))
    path <- file.path(
      cfg$base_dir, sprintf("case%d_%s_n%d_rep%d.RData", case, tag, n, cfg$rep)
    )
    if (!file.exists(path)) next
    results <- load_results(path)
    audits <- results$exact_zero_audit[[cfg$method]]
    summaries <- bind_nonempty(lapply(audits, function(x) {
      if (is.list(x)) x$summary else NULL
    }))
    witnesses <- bind_nonempty(lapply(audits, function(x) {
      if (is.list(x)) x$literal_violation_witnesses else NULL
    }))
    if (!nrow(summaries)) next
    summaries$result_file <- path
    rep_rows[[pos]] <- summaries
    if (nrow(witnesses)) {
      witnesses$result_file <- path
      witness_rows[[pos]] <- witnesses
    }
    config_rows[[pos]] <- aggregate_config(
      summaries, results, cfg$method, case, tau, n, cfg$rep
    )
    pos <- pos + 1L
  }
  config_df <- bind_nonempty(config_rows)
  rep_df <- bind_nonempty(rep_rows)
  witness_df <- bind_nonempty(witness_rows)
  if (!nrow(config_df)) stop("No exact-zero audit results found for the requested grid.")
  write.csv(config_df, file.path(cfg$output_dir, "exact_zero_config_summary.csv"), row.names = FALSE)
  write.csv(rep_df, file.path(cfg$output_dir, "exact_zero_replication_summary.csv"), row.names = FALSE)
  write.csv(
    witness_df,
    file.path(cfg$output_dir, "exact_zero_literal_violation_witnesses.csv"),
    row.names = FALSE
  )
  cat("Exact-zero configurations:", nrow(config_df), "\n")
  cat("Audited replications:", sum(config_df$n_replications_audited), "\n")
  cat("Audited next-window transitions:", sum(config_df$n_next_window_points_audited), "\n")
  cat("Literal violation points:", sum(config_df$n_literal_violation_points), "\n")
  cat("Output directory:", cfg$output_dir, "\n")
}

main()
