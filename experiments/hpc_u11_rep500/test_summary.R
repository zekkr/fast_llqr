#!/usr/bin/env Rscript
root <- file.path(tempdir(), paste0("u11_summary_", Sys.getpid()))
run_tag <- "synthetic_rep500"
dir.create(root, recursive = TRUE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)
methods <- c("direct_baseline", "lean_seq", "unified_u11")

for (model in c("llqr", "tvcqr")) for (case_id in 1:2) for (tau in c(0.2, 0.5, 0.8)) for (n in c(1000L, 2000L, 5000L, 10000L)) {
  tag <- sprintf("case%d_tau%02d_n%d", case_id, round(100 * tau), n)
  out <- file.path(root, run_tag, model, tag)
  dir.create(out, recursive = TRUE)
  blocks <- lapply(methods, function(method) {
    is_u11 <- method == "unified_u11"
    pos <- ((match(method, methods) - 1L - (0:499) %% 3L) %% 3L) + 1L
    data.frame(
      run_tag, model, case = case_id, tau, n, rep_id = 1:500, seed = 2025L + 1:500,
      method, method_position = pos,
      elapsed_sec = c(direct_baseline = 2, lean_seq = 1, unified_u11 = 0.5)[[method]],
      solver_ok = TRUE, accepted_ok = TRUE, finite_estimate = TRUE,
      discrepancy = if (method == "direct_baseline") 0 else 1e-12,
      discrepancy_status = "ok", ierr = if (is_u11) 0L else NA_integer_,
      failed_eval = if (is_u11) 0L else NA_integer_,
      h_check_status = if (method == "direct_baseline") "not_applicable_direct" else "match",
      h_match_vs_lean = if (method == "direct_baseline") NA else TRUE,
      h_mismatch_eval_count = if (method == "direct_baseline") NA_integer_ else 0L,
      first_h_mismatch_eval = NA_integer_,
      first_tableau_rows_max = if (is_u11) 12L else NA_integer_,
      first_retained_size_max = if (is_u11) 10L else NA_integer_,
      initial_threshold_hits_max = if (is_u11) 8L else NA_integer_,
      effective_threshold_hits_max = if (is_u11) 8L else NA_integer_,
      threshold_expansion_total = if (is_u11) 0L else NA_integer_,
      threshold_expansion_points = if (is_u11) 0L else NA_integer_,
      threshold_expansion_steps_max = if (is_u11) 0L else NA_integer_,
      threshold_initial = if (is_u11) 0.02 else NA_real_,
      threshold_effective_max = if (is_u11) 0.02 else NA_real_,
      basis_forced_rows_max = if (is_u11) 2L else NA_integer_,
      first_aggregate_rows_max = if (is_u11) 2L else NA_integer_,
      retained_decomposition_ok = if (is_u11) TRUE else NA,
      first_pass_count = if (is_u11) n - 1L else NA_integer_,
      first_pass_total = if (is_u11) n - 1L else NA_integer_,
      first_pass_rate = if (is_u11) 1 else NA_real_,
      repaired_points = if (is_u11) 0L else NA_integer_,
      repair_total = if (is_u11) 0L else NA_integer_,
      internal_recovery_total = if (is_u11) 0L else NA_integer_,
      first_full_m_recovery = if (is_u11) 0L else NA_integer_,
      certificate_hits = if (is_u11) 0L else NA_integer_,
      residual_rows_checked = if (is_u11) n else NA_integer_,
      fallback_triggered = FALSE, error_message = NA_character_, threw_error = FALSE,
      error_stage = NA_character_, data_generation_sec = 0.01,
      node_name = "synthetic", cpu_model = "synthetic"
    )
  })
  d <- do.call(rbind, blocks)
  saveRDS(d, file.path(out, "replication_metrics.rds"))
  saveRDS(data.frame(), file.path(out, "attempt_metrics.rds"))
  write.csv(data.frame(model, case = case_id, tau, n, complete = TRUE), file.path(out, "integrity.csv"), row.names = FALSE)
}

experiment_dir <- normalizePath(Sys.getenv("SSQR_EXPERIMENT_DIR", "experiments/hpc_u11_rep500"))
env <- c(
  paste0("SSQR_PROJECT_ROOT=", normalizePath(".")),
  paste0("SSQR_OUTPUT_ROOT=", root), paste0("SSQR_RUN_TAG=", run_tag),
  "SSQR_NUM_REP=500", "SSQR_SEED_BASE=2025"
)
status <- system2("Rscript", file.path(experiment_dir, "summarize_run.R"), env = env)
stopifnot(status == 0L)
stage_status <- system2("Rscript", file.path(experiment_dir, "build_paper_staging.R"), env = env)
stopifnot(stage_status == 0L)
summary <- read.csv(file.path(root, run_tag, "tables", "config_method_summary.csv"))
ratios <- read.csv(file.path(root, run_tag, "tables", "time_ratios.csv"))
screen <- read.csv(file.path(root, run_tag, "tables", "paper_screening_tau05_staging.csv"))
lines <- readLines(file.path(root, run_tag, "tables", "run_summary.txt"))
stopifnot(nrow(summary) == 144L, nrow(ratios) == 48L, nrow(screen) == 8L)
stopifnot(all(ratios$u11_over_lean == 0.5), any(lines == "paper_staging_ready: TRUE"))
stage <- file.path(root, run_tag, "paper_staging")
stopifnot(
  file.exists(file.path(stage, "rep500_seed2025_method_summary_candidate.csv")),
  file.exists(file.path(stage, "runtime_tau05_preview.png")),
  file.exists(file.path(stage, "runtime_complete_candidate.tex")),
  nrow(read.csv(file.path(stage, "lean_seq_ablation_candidate.csv"))) == 48L
)
cat("PASS synthetic 48-config rep500 summary, stability gate, and paper staging\n")
