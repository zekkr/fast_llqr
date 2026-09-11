#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

project <- normalizePath(Sys.getenv("SSQR_PROJECT_ROOT", getwd()), mustWork = TRUE)
setwd(project)
run_tag <- Sys.getenv("SSQR_RUN_TAG", "")
output_root <- Sys.getenv("SSQR_OUTPUT_ROOT", "")
num_rep <- as.integer(Sys.getenv("SSQR_NUM_REP", "500"))
seed_base <- as.integer(Sys.getenv("SSQR_SEED_BASE", "2025"))
stopifnot(nzchar(run_tag), nzchar(output_root), num_rep == 500L, seed_base == 2025L)

methods <- c("direct_baseline", "lean_seq", "unified_u11")
models <- c("llqr", "tvcqr")
cases <- 1:2
taus <- c(0.2, 0.5, 0.8)
ns <- c(1000L, 2000L, 5000L, 10000L)
run_dir <- file.path(output_root, run_tag)
tables <- file.path(run_dir, "tables")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)

metrics <- list()
attempts <- list()
integrity <- list()
missing_configs <- character()
for (model in models) for (case_id in cases) for (tau in taus) for (n in ns) {
  tag <- sprintf("case%d_tau%02d_n%d", case_id, round(100 * tau), n)
  directory <- file.path(run_dir, model, tag)
  paths <- file.path(directory, c("replication_metrics.rds", "attempt_metrics.rds", "integrity.csv"))
  if (any(!file.exists(paths))) {
    missing_configs <- c(missing_configs, paste(model, tag, sep = "/"))
    next
  }
  metrics[[length(metrics) + 1L]] <- readRDS(paths[[1L]])
  attempts[[length(attempts) + 1L]] <- readRDS(paths[[2L]])
  integrity[[length(integrity) + 1L]] <- read.csv(paths[[3L]], stringsAsFactors = FALSE)
}

m <- if (length(metrics)) do.call(rbind, metrics) else data.frame()
a <- if (length(attempts)) do.call(rbind, attempts) else data.frame()
ii <- if (length(integrity)) do.call(rbind, integrity) else data.frame()
sum_defined <- function(x) if (all(is.na(x))) NA_real_ else sum(x[!is.na(x)])
max_defined <- function(x) if (all(is.na(x))) NA_real_ else max(x[!is.na(x)])

rows <- list()
idx <- 0L
if (nrow(m)) for (model in models) for (case_id in cases) for (tau in taus) for (n in ns) for (method in methods) {
  d <- m[m$model == model & m$case == case_id & m$tau == tau & m$n == n & m$method == method, , drop = FALSE]
  idx <- idx + 1L
  solver <- nrow(d) == num_rep && all(d$solver_ok %in% TRUE)
  accepted <- nrow(d) == num_rep && all(d$accepted_ok %in% TRUE)
  time_valid <- solver && all(is.finite(d$elapsed_sec))
  discrepancy_valid <- nrow(d) == num_rep && all(d$discrepancy_status == "ok") && all(is.finite(d$discrepancy))
  positions <- table(factor(d$method_position, levels = seq_along(methods)))
  balanced <- nrow(d) == num_rep && sum(positions) == num_rep && diff(range(positions)) <= 1L
  is_u11 <- method == "unified_u11"
  decomposition_ok <- !is_u11 || (
    nrow(d) == num_rep && all(d$retained_decomposition_ok %in% TRUE)
  )
  rows[[idx]] <- data.frame(
    model, case = case_id, paper_case = if (model == "llqr") case_id else case_id + 2L,
    tau, n, num_rep, seed_base, method, n_present = nrow(d),
    n_solver_ok = sum(d$solver_ok %in% TRUE), n_accepted = sum(d$accepted_ok %in% TRUE),
    n_nonfinite = sum(!(d$finite_estimate %in% TRUE)),
    ierr_nonzero = sum(!is.na(d$ierr) & d$ierr != 0L),
    failed_eval_nonzero = sum(!is.na(d$failed_eval) & d$failed_eval != 0L),
    h_mismatch_count = sum(d$h_check_status == "row_set_mismatch"),
    fallback_count = sum(d$fallback_triggered %in% TRUE),
    retained_decomposition_failures = if (is_u11) sum(!(d$retained_decomposition_ok %in% TRUE)) else NA_integer_,
    repair_total = sum_defined(d$repair_total),
    repaired_points_total = sum_defined(d$repaired_points),
    internal_recovery_total = sum_defined(d$internal_recovery_total),
    certificate_hits = sum_defined(d$certificate_hits),
    residual_rows_checked = sum_defined(d$residual_rows_checked),
    threshold_expansion_total = sum_defined(d$threshold_expansion_total),
    threshold_expansion_replications = if (is_u11) sum(d$threshold_expansion_total > 0L) else NA_integer_,
    first_retained_size_mean_of_max = if (is_u11) mean(d$first_retained_size_max) else NA_real_,
    first_retained_size_max = if (is_u11) max(d$first_retained_size_max) else NA_integer_,
    time_mean_sec = if (time_valid) mean(d$elapsed_sec) else NA_real_,
    time_median_sec = if (time_valid) median(d$elapsed_sec) else NA_real_,
    time_min_sec = if (time_valid) min(d$elapsed_sec) else NA_real_,
    time_max_sec = if (time_valid) max(d$elapsed_sec) else NA_real_,
    max_average_relative_discrepancy = if (discrepancy_valid) max(d$discrepancy) else NA_real_,
    solver_stable = solver, path_accepted = accepted,
    timing_order_balanced = balanced, discrepancy_valid = discrepancy_valid,
    retained_decomposition_ok = decomposition_ok,
    stringsAsFactors = FALSE
  )
}

summary <- if (length(rows)) do.call(rbind, rows) else data.frame()
write.csv(summary, file.path(tables, "config_method_summary.csv"), row.names = FALSE, na = "")

key <- c("model", "case", "paper_case", "tau", "n")
wide_time <- reshape(
  summary[, c(key, "method", "time_mean_sec")], idvar = key,
  timevar = "method", direction = "wide"
)
names(wide_time) <- sub("^time_mean_sec[.]", "time_", names(wide_time))
wide_time$u11_over_direct <- wide_time$time_unified_u11 / wide_time$time_direct_baseline
wide_time$u11_over_lean <- wide_time$time_unified_u11 / wide_time$time_lean_seq
wide_time$lean_over_direct <- wide_time$time_lean_seq / wide_time$time_direct_baseline
write.csv(wide_time, file.path(tables, "time_ratios.csv"), row.names = FALSE, na = "")
write.csv(
  summary[, c(key, "method", "max_average_relative_discrepancy", "discrepancy_valid")],
  file.path(tables, "numerical_discrepancy.csv"), row.names = FALSE, na = ""
)

replacement_attempts <- if (nrow(a)) {
  a[a$attempt > 1L | a$baseline_retryable %in% TRUE, , drop = FALSE]
} else data.frame()
write.csv(replacement_attempts, file.path(tables, "baseline_regeneration_attempts.csv"), row.names = FALSE, na = "")
write.csv(ii, file.path(tables, "integrity.csv"), row.names = FALSE, na = "")

u11_rep <- if (nrow(m)) m[m$method == "unified_u11", c(
  "run_tag", "model", "case", "tau", "n", "rep_id", "seed",
  "first_retained_size_max", "first_tableau_rows_max", "initial_threshold_hits_max",
  "effective_threshold_hits_max", "basis_forced_rows_max", "first_aggregate_rows_max",
  "threshold_initial", "threshold_effective_max", "threshold_expansion_total",
  "threshold_expansion_points", "first_pass_count", "first_pass_total", "first_pass_rate",
  "repaired_points", "repair_total", "internal_recovery_total", "certificate_hits",
  "residual_rows_checked", "retained_decomposition_ok"
), drop = FALSE] else data.frame()
write.csv(u11_rep, file.path(tables, "retained_size_replication_metrics.csv"), row.names = FALSE, na = "")

screen_rows <- list()
sidx <- 0L
if (nrow(u11_rep)) {
  selected <- u11_rep[u11_rep$tau == 0.5 & (
    (u11_rep$model == "llqr" & u11_rep$case == 2L) |
      (u11_rep$model == "tvcqr" & u11_rep$case == 1L)
  ), , drop = FALSE]
  for (model in c("llqr", "tvcqr")) for (n in ns) {
    d <- selected[selected$model == model & selected$n == n, , drop = FALSE]
    if (nrow(d) != num_rep) next
    b_n <- n^(-0.2)
    ratio_nominal <- d$first_retained_size_max /
      (n * b_n * d$threshold_initial + log(n))
    ratio_effective_conservative <- d$first_retained_size_max /
      (n * b_n * d$threshold_effective_max + log(n))
    sidx <- sidx + 1L
    screen_rows[[sidx]] <- data.frame(
      model, case = unique(d$case), paper_case = if (model == "llqr") 2L else 3L,
      tau = 0.5, n, num_rep, seed_base, raw_method = "unified_u11",
      Mm.factor = if (model == "llqr") 0.1 else 1e-5,
      ok_reps = nrow(d), gamma_mean = mean(d$threshold_initial),
      gamma_min = min(d$threshold_initial), gamma_max = max(d$threshold_initial),
      effective_gamma_max_mean = mean(d$threshold_effective_max),
      mean_max_Sj = mean(d$first_retained_size_max),
      max_Sj_over_nb = mean(d$first_retained_size_max) / (n * b_n),
      max_Sj_over_nb_gamma_logn = mean(ratio_nominal),
      max_Sj_over_nb_effective_gamma_logn_conservative = mean(ratio_effective_conservative),
      first_pass_prop_mean = mean(d$first_pass_rate),
      first_pass_prop_min = min(d$first_pass_rate),
      first_pass_prop_max = max(d$first_pass_rate),
      Fr_mean = mean(d$repaired_points), Fr_median = median(d$repaired_points),
      Fr_q90 = as.numeric(quantile(d$repaired_points, 0.9, names = FALSE, type = 7)),
      Fr_max = max(d$repaired_points), total_repair_mean = mean(d$repair_total),
      total_repair_min = min(d$repair_total), total_repair_max = max(d$repair_total),
      uniform = mean(d$repaired_points == 0L),
      threshold_expansion_replications = sum(d$threshold_expansion_total > 0L),
      stringsAsFactors = FALSE
    )
  }
}
screening <- if (length(screen_rows)) do.call(rbind, screen_rows) else data.frame()
write.csv(screening, file.path(tables, "paper_screening_tau05_staging.csv"), row.names = FALSE, na = "")

paper_methods <- summary[summary$method %in% c("direct_baseline", "unified_u11"), , drop = FALSE]
lean_ablation <- summary[summary$method == "lean_seq", , drop = FALSE]
write.csv(paper_methods, file.path(tables, "paper_method_summary_staging.csv"), row.names = FALSE, na = "")
write.csv(lean_ablation, file.path(tables, "lean_seq_ablation_staging.csv"), row.names = FALSE, na = "")

grid_complete <- length(missing_configs) == 0L && nrow(ii) == 48L && all(ii$complete %in% TRUE)
all_methods_complete <- nrow(summary) == 144L && all(summary$n_present == num_rep)
all_solver_stable <- all_methods_complete && all(summary$solver_stable) && all(summary$path_accepted)
timing_ok <- all_methods_complete && all(summary$timing_order_balanced)
discrepancy_ok <- all_methods_complete && all(summary$discrepancy_valid)
u11 <- summary[summary$method == "unified_u11", , drop = FALSE]
u11_gate <- nrow(u11) == 48L && all(u11$ierr_nonzero == 0L) &&
  all(u11$failed_eval_nonzero == 0L) && all(u11$h_mismatch_count == 0L) &&
  all(u11$fallback_count == 0L) && all(u11$retained_decomposition_ok) &&
  all(is.finite(u11$max_average_relative_discrepancy)) &&
  all(u11$max_average_relative_discrepancy <= 1e-8)
paper_staging_ready <- grid_complete && all_solver_stable && timing_ok && discrepancy_ok && u11_gate

max_disc_row <- if (nrow(u11)) u11[which.max(u11$max_average_relative_discrepancy), , drop = FALSE] else NULL
lines <- c(
  sprintf("run_tag: %s", run_tag),
  sprintf("grid_complete: %s", grid_complete),
  sprintf("all_solver_stable_and_path_accepted: %s", all_solver_stable),
  sprintf("timing_order_balanced: %s", timing_ok),
  sprintf("all_discrepancies_valid: %s", discrepancy_ok),
  sprintf("u11_stability_gate: %s", u11_gate),
  sprintf("paper_staging_ready: %s", paper_staging_ready),
  sprintf("missing_configs: %s", if (length(missing_configs)) paste(missing_configs, collapse = ",") else "none"),
  if (is.null(max_disc_row)) "max_u11_discrepancy: unavailable" else sprintf(
    "max_u11_discrepancy: %.17g at %s case=%d tau=%.1f n=%d",
    max_disc_row$max_average_relative_discrepancy, max_disc_row$model,
    max_disc_row$case, max_disc_row$tau, max_disc_row$n
  )
)
writeLines(lines, file.path(tables, "run_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
if (!paper_staging_ready) quit(save = "no", status = 2L)
