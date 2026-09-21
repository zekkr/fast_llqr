#!/usr/bin/env Rscript
# Standalone tables only: this script never reads or writes manuscript files.
options(stringsAsFactors = FALSE)
run_dir <- file.path(Sys.getenv("SSQR_OUTPUT_ROOT"), Sys.getenv("SSQR_RUN_TAG"))
tables <- file.path(run_dir, "tables")
stopifnot(any(readLines(file.path(tables, "run_summary.txt")) == "paper_staging_ready: TRUE"))
read <- function(name) read.csv(file.path(tables, name), check.names = FALSE)
s <- read("config_method_summary.csv")
x <- read("screening_all_taus.csv")
r <- read("time_ratios.csv")
m <- read("all_replication_metrics.csv")
a <- read("baseline_regeneration_attempts.csv")
stopifnot(nrow(s) == 36L, nrow(x) == 12L, nrow(r) == 12L)
stage <- file.path(run_dir, "paper_staging")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
export <- function(d, name) write.csv(d, file.path(stage, paste0(name, ".csv")), row.names = FALSE, na = "")
runtime <- s[, c("tau", "n", "method", "time_mean_sec", "time_min_sec", "time_max_sec")]
runtime$time_range_sec <- runtime$time_max_sec - runtime$time_min_sec
export(runtime, "runtime_all")
export(x, "screening_all")
abl <- r[, c("tau", "n", "time_lean_seq", "time_unified_u11", "u11_over_lean", "u11_over_direct")]
abl$runtime_change_percent <- 100 * (abl$u11_over_lean - 1)
abl$direct_speedup <- 1 / abl$u11_over_direct
export(abl, "ablation_by_config")
pooled <- data.frame(seq_mean_sec = mean(abl$time_lean_seq),
                     screen_mean_sec = mean(abl$time_unified_u11))
pooled$runtime_change_percent <- 100 * (pooled$screen_mean_sec / pooled$seq_mean_sec - 1)
export(pooled, "ablation_case_mean")
disc <- s[s$method != "direct_baseline", c("tau", "n", "method", "max_average_relative_discrepancy")]
export(disc, "numerical_discrepancy")
u <- m[m$method == "unified_u11", ]
grid <- do.call(rbind, lapply(split(u, interaction(u$tau, u$n, drop = TRUE)), function(d) {
  data.frame(tau = d$tau[1L], n = d$n[1L], n_eval_mean = mean(d$n_eval),
             n_eval_min = min(d$n_eval), n_eval_max = max(d$n_eval),
             no_transition_reps = sum(d$n_transitions == 0L),
             placeholder_reps = sum(d$grid_placeholder))
}))
grid <- grid[order(grid$tau, grid$n), ]
export(grid, "evaluation_grid")
status <- s[, c("tau", "n", "method", "n_present", "n_solver_ok", "n_accepted",
                "n_nonfinite", "ierr_nonzero", "failed_eval_nonzero", "h_mismatch_count",
                "fallback_count", "internal_recovery_total", "threshold_expansion_total")]
export(status, "integrity_and_recovery")
regen <- if (nrow(a)) a[a$method == "direct_baseline" & a$baseline_retryable, ] else a
export(regen, "failed_direct_attempts")

fmt <- function(v) {
  if (is.numeric(v)) return(ifelse(is.na(v), "NA", formatC(v, digits = 5, format = "g")))
  ifelse(is.na(v), "NA", as.character(v))
}
md_table <- function(d) {
  vals <- lapply(d, fmt)
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |"),
    vapply(seq_len(nrow(d)), function(i) paste0("| ", paste(vapply(vals, `[`, character(1L), i), collapse = " | "), " |"), character(1L)))
}
lines <- c("# LLQR Case 2: standardized logistic errors, interior evaluation grid", "",
  paste("Run:", Sys.getenv("SSQR_RUN_TAG")), "",
  "X ~ U(0,1); Y = 1 + 2X^2 + logistic(0, sqrt(3)/pi). The fit uses all n observations; evaluation points are observed X in [0.1,0.9]. Epanechnikov kernel, h=n^(-1/5), threshold=0.1*sqrt(log(n)/(n*h))*log(log(n)).",
  "", "Each of 12 configurations has 500 final paired replications. Only direct-fit exceptions trigger regeneration (at most 20 attempts). Candidate failure never changes the seed. Tables below use final attempts only.", "",
  sprintf("Failed direct attempts: %d. Affected configuration-replication pairs: %d.", nrow(regen), if (nrow(regen)) nrow(unique(regen[,c("tau","n","rep_id")])) else 0L),
  "", "Runtime is seconds; the full fitting call is timed, with common grid construction/sorting and post-fit checks excluded. Methods rotate positions across replications.")
add <- function(title, d, note = NULL) {
  lines <<- c(lines, "", paste0("## ", title), "", note, "", md_table(d))
}
for (tau in c(.5, .2, .8)) {
 d <- x[x$tau == tau, ]; d <- d[order(d$n), ]
 first <- d[,c("n", "gamma_mean", "first_pass_prop_mean", "first_pass_prop_min", "first_pass_prop_max",
                "Fr_mean", "Fr_median", "Fr_q90", "Fr_max", "total_repair_mean", "total_repair_min", "total_repair_max", "uniform")]
 add(paste("First pass and repairs, tau =", tau), first,
     "First-pass rate is calculated within each replication over m_r-1 transitions, then averaged. The initial full fit is excluded. 'uniform' is the proportion with no repair anywhere on the grid; it is not a pooled pointwise rate.")
 kept <- d[,c("n", "gamma_mean", "mean_max_Sj", "max_Sj_over_nb", "max_Sj_over_nb_gamma_logn")]
 add(paste("Initial retained size, tau =", tau), kept,
     "Counts include interpolation-basis padding and exclude aggregate rows. nb is a reference scale, not the realized active-set size. Empty transition maxima are defined as zero.")
}
add("Runtime (Table S2 counterpart)", runtime, "Range = maximum minus minimum; direct_baseline = direct, lean_seq = seq, unified_u11 = screen-seq.")
add("Ablation by configuration", abl, "Positive runtime_change_percent means screen-seq is slower than seq; direct_speedup is direct/screen-seq.")
add("Case-level ablation (Table S3 counterpart)", pooled, "Means weight the 12 configuration-level mean runtimes equally; the percentage is the ratio of these means, not the mean of percentages.")
add("Numerical discrepancy (Tables S4-S6 counterparts)", disc,
    "For each replication, average |candidate-direct|/max(|direct|,1e-10) over its actual evaluation points. Report the maximum over 500 replications. Derivatives are excluded; this is numerical agreement, not statistical bias.")
add("Evaluation-grid sizes", grid)
add("Integrity and recovery", status,
    "H row sets are checked in memory against seq. Internal full-active recovery is distinct from sequential fallback. Full H paths are not saved.")
lines <- c(lines, "", "No manuscript, prior paper asset, or previous experiment result was modified.")
writeLines(lines, file.path(stage, "REPORT.md"))
cat("PASS standalone report and nine CSV tables written to", stage, "\n")
