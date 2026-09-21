#!/usr/bin/env Rscript
run_tag <- Sys.getenv("MST_RUN_TAG", "")
output_root <- Sys.getenv("MST_OUTPUT_ROOT", "results/hpc/multivar_mst_runs")
run_dir <- file.path(output_root, run_tag); tables <- file.path(run_dir, "tables")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
configs <- expand.grid(case = 1:2, n = c(500L, 1000L, 2000L))
all_rows <- do.call(rbind, lapply(seq_len(nrow(configs)), function(i) {
  dir <- file.path(run_dir, sprintf("case%d_n%d", configs$case[i], configs$n[i]))
  if (!file.exists(file.path(dir, "CONFIG_PASS"))) stop("Missing CONFIG_PASS: ", dir)
  read.csv(file.path(dir, "pilot_metrics.csv"), stringsAsFactors = FALSE)
}))
thresholds <- c(0.1, 0.2, 0.5, 1.0)
config_rows <- do.call(rbind, lapply(seq_len(nrow(configs)), function(i) {
  do.call(rbind, lapply(thresholds, function(factor) {
    d <- subset(all_rows, case == configs$case[i] & n == configs$n[i] &
                          abs(threshold_factor - factor) < 1e-12)
    data.frame(case = configs$case[i], n = configs$n[i], threshold_factor = factor,
      num_rep = nrow(d), direct_mean = mean(d$direct_time), screen_mean = mean(d$screen_time),
      time_ratio = mean(d$screen_time) / mean(d$direct_time),
      reduction_percent = 100 * (1 - mean(d$screen_time) / mean(d$direct_time)),
      repair_total = sum(d$repair_total), repaired_points_total = sum(d$repaired_points),
      recovery_total = sum(d$full_active_recovery_total),
      retained_size_mean = mean(d$retained_size_mean),
      retained_fraction_mean = mean(d$retained_fraction_mean),
      max_objective_gap = max(d$objective_gap), max_kkt_residual = max(d$kkt_residual),
      max_relative_fit_difference = max(d$max_relative_fit_difference),
      h_mismatch_count = sum(!d$h_set_match), audit_failure_count = sum(!d$audit_ok))
  }))
}))
aggregate <- do.call(rbind, lapply(thresholds, function(factor) {
  c <- subset(config_rows, abs(threshold_factor - factor) < 1e-12)
  d <- subset(all_rows, abs(threshold_factor - factor) < 1e-12)
  data.frame(threshold_factor = factor, eligible = all(d$audit_ok & d$h_set_match & d$direct_ok & d$screen_ok),
    equal_config_time_ratio = mean(c$time_ratio),
    equal_config_reduction_percent = 100 * (1 - mean(c$time_ratio)),
    repair_total = sum(d$repair_total), recovery_total = sum(d$full_active_recovery_total),
    retained_fraction_mean = mean(d$retained_fraction_mean),
    max_objective_gap = max(d$objective_gap), max_kkt_residual = max(d$kkt_residual),
    max_relative_fit_difference = max(d$max_relative_fit_difference),
    h_mismatch_count = sum(!d$h_set_match), audit_failure_count = sum(!d$audit_ok))
}))
eligible <- subset(aggregate, eligible)
if (!nrow(eligible)) stop("No threshold factor passed the pilot audit gate.")
best <- min(eligible$equal_config_time_ratio)
shortlist <- subset(eligible, equal_config_time_ratio <= best * 1.02)
shortlist <- shortlist[order(shortlist$recovery_total, shortlist$repair_total,
                             shortlist$retained_fraction_mean, shortlist$threshold_factor), ]
selected <- shortlist$threshold_factor[1L]
aggregate$within_two_percent <- aggregate$eligible & aggregate$equal_config_time_ratio <= best * 1.02
aggregate$selected <- abs(aggregate$threshold_factor - selected) < 1e-12
write.csv(all_rows, file.path(tables, "pilot_all_metrics.csv"), row.names = FALSE)
write.csv(config_rows, file.path(tables, "pilot_by_config.csv"), row.names = FALSE)
write.csv(aggregate, file.path(tables, "pilot_threshold_summary.csv"), row.names = FALSE)
writeLines(format(selected, scientific = FALSE, trim = TRUE), file.path(tables, "selected_threshold.txt"))
writeLines(c(
  "Eligibility: all 30 paired pilot data sets pass solver, H-set, objective, KKT and fit audits.",
  "Primary score: minimize the equal-weight mean of the six configuration mean screen/direct time ratios.",
  "Within 2% of the minimum: minimize total recovery, then repair, then retained fraction, then threshold factor."
), file.path(tables, "selection_rule.txt"))
print(config_rows); print(aggregate); cat("SELECTED_THRESHOLD=", selected, "\n", sep = "")
