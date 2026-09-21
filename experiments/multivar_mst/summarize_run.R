#!/usr/bin/env Rscript
run_tag <- Sys.getenv("MST_RUN_TAG", "")
output_root <- Sys.getenv("MST_OUTPUT_ROOT", "results/hpc/multivar_mst_runs")
expected_rep <- as.integer(Sys.getenv("MST_NUM_REP", "100"))
run_dir <- file.path(output_root, run_tag)
tables <- file.path(run_dir, "tables"); dir.create(tables, recursive = TRUE, showWarnings = FALSE)
configs <- expand.grid(case = 1:2, n = c(500L, 1000L, 2000L))
rows <- lapply(seq_len(nrow(configs)), function(i) {
  case <- configs$case[i]; n <- configs$n[i]
  dir <- file.path(run_dir, sprintf("case%d_n%d", case, n))
  if (!file.exists(file.path(dir, "CONFIG_PASS"))) stop("Missing CONFIG_PASS: ", dir)
  d <- read.csv(file.path(dir, "replication_metrics.csv"), stringsAsFactors = FALSE)
  data.frame(
    case = case, n = n, num_rep = nrow(d),
    direct_mean = mean(d$direct_time), direct_min = min(d$direct_time), direct_max = max(d$direct_time),
    screen_mean = mean(d$screen_time), screen_min = min(d$screen_time), screen_max = max(d$screen_time),
    reduction_percent = 100 * (1 - mean(d$screen_time) / mean(d$direct_time)),
    max_objective_gap = max(d$objective_gap), max_kkt_residual = max(d$kkt_residual),
    max_fit_discrepancy = max(d$fit_discrepancy),
    max_abs_fit_difference = max(d$max_abs_fit_difference),
    h_mismatch_count = sum(!d$h_set_match), audit_failure_count = sum(!d$audit_ok),
    fallback_count = 0L, recovery_total = sum(d$full_active_recovery_total),
    repair_total = sum(d$repair_total), repaired_points_total = sum(d$repaired_points),
    retained_size_mean = mean(d$retained_size_mean),
    retained_fraction_mean = mean(d$retained_fraction_mean),
    threshold_expansion_total = sum(d$threshold_expansion_total),
    underflow_total = sum(d$underflow_total), underflow_max = max(d$underflow_max),
    mst_time_mean = mean(d$mst_seconds), max_edge_mean = mean(d$max_edge),
    retry_count = sum(d$attempt - 1L), max_attempt = max(d$attempt))
})
summary <- do.call(rbind, rows)
write.csv(summary, file.path(tables, "runtime_summary.csv"), row.names = FALSE)
if (any(summary$num_rep != expected_rep | summary$h_mismatch_count != 0L |
        summary$audit_failure_count != 0L)) stop("Formal summary gate failed.")
writeLines("PASS", file.path(run_dir, "RESULTS_PASS"))
print(summary)
