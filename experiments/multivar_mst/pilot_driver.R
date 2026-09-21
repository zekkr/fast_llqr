#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

env_int <- function(name, default = NA_integer_) {
  value <- Sys.getenv(name, "")
  if (nzchar(value)) as.integer(value) else as.integer(default)
}
project <- normalizePath(Sys.getenv("MST_PROJECT_ROOT", getwd()), mustWork = TRUE)
setwd(project)
case_id <- env_int("MST_CASE"); n <- env_int("MST_N")
num_rep <- env_int("MST_NUM_REP", 5L); seed_base <- env_int("MST_SEED_BASE", 9025L)
task_id <- env_int("SLURM_ARRAY_TASK_ID", 1L); run_tag <- Sys.getenv("MST_RUN_TAG", "")
output_root <- Sys.getenv("MST_OUTPUT_ROOT", "results/hpc/multivar_mst_runs")
stopifnot(case_id %in% 1:2, n %in% c(500L, 1000L, 2000L), num_rep == 5L,
          seed_base == 9025L, task_id %in% seq_len(num_rep), nzchar(run_tag))

source("experiments/multivar_mst/design.R")
source("experiments/multivar_mst/solver.R")
stopifnot(requireNamespace("quantreg", quietly = TRUE))
kernel <- load_multivar_mst_kernel("experiments/multivar_mst/build/multivar_mst_optimized.so")
thresholds <- c(0.1, 0.2, 0.5, 1.0)

cpu_model <- tryCatch({
  if (!file.exists("/proc/cpuinfo")) stop("cpuinfo unavailable")
  trimws(sub("^[^:]+:", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1L]))
}, error = function(e) unname(Sys.info()[["machine"]]))
thread_env <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"))
if (Sys.getenv("MST_REQUIRE_HPC", "0") == "1")
  stopifnot(getRversion() == "4.3.1", packageVersion("quantreg") == "5.94",
            grepl("6258R", cpu_model), all(thread_env == "1"))

run_dir <- file.path(output_root, run_tag)
config_dir <- file.path(run_dir, sprintf("case%d_n%d", case_id, n))
partial_dir <- file.path(config_dir, "partial")
dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
meta_dir <- file.path(run_dir, "_run_meta", "workers")
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(session = sessionInfo(), system = Sys.info(), cpu = cpu_model,
             thread_env = thread_env, git_sha = Sys.getenv("MST_PUSHED_SHA"),
             thresholds = thresholds, seed_base = seed_base),
        file.path(meta_dir, sprintf("pilot_case%d_n%d_task%d_%s.rds", case_id, n, task_id,
                                    unname(Sys.info()[["nodename"]]))))

row_set_match <- function(a, b) nrow(a) == nrow(b) &&
  all(vapply(seq_len(nrow(a)), function(j) setequal(a[j, ], b[j, ]), logical(1L)))
dataset_md5 <- function(dat) {
  path <- tempfile("mst-pilot-data-", fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(list(x = dat$x, y = dat$y), path, version = 3); unname(tools::md5sum(path))
}
timed <- function(expr) {
  invisible(gc()); started <- proc.time()[["elapsed"]]
  value <- tryCatch(force(expr), error = identity)
  list(value = value, elapsed = proc.time()[["elapsed"]] - started)
}

rep_id <- task_id; original_seed <- seed_base + rep_id
for (attempt in seq_len(20L)) {
  seed <- original_seed + (attempt - 1L) * 1000000L
  dat <- simulate_multivar_llqr(n, case_id, seed)
  labels <- c("direct", paste0("screen_", thresholds))
  shift <- (rep_id + case_id + match(n, c(500L, 1000L, 2000L)) - 3L) %% length(labels)
  order <- labels[c(seq_len(length(labels)) + shift - 1L) %% length(labels) + 1L]
  fit <- setNames(vector("list", length(labels)), labels)
  for (label in order) {
    if (label == "direct") fit[[label]] <- timed(llqr_direct_multivar(dat$x, dat$y, 0.5, dat$z))
    else {
      factor <- as.numeric(sub("screen_", "", label, fixed = TRUE))
      fit[[label]] <- timed(llqr_seq_screen_mst(dat$x, dat$y, 0.5, dat$z,
        threshold_factor = factor, kernel = kernel, diagnostics = TRUE))
    }
  }
  direct_error <- inherits(fit$direct$value, "error")
  if (direct_error && attempt < 20L) next
  if (direct_error) stop("Direct pilot fit failed after 20 attempts: ", conditionMessage(fit$direct$value))
  screen_errors <- vapply(labels[-1L], function(label) inherits(fit[[label]]$value, "error"), logical(1L))
  if (any(screen_errors)) {
    label <- labels[-1L][which(screen_errors)[1L]]
    stop(label, " pilot fit failed: ", conditionMessage(fit[[label]]$value))
  }
  direct <- fit$direct$value
  full <- llqr_seq_screen_mst(dat$x, dat$y, 0.5, dat$z, h = direct$h,
    threshold_factor = 1, kernel = kernel, diagnostics = TRUE, full_active = TRUE)
  rows <- lapply(thresholds, function(factor) {
    label <- paste0("screen_", factor); screen <- fit[[label]]$value
    os <- weighted_qr_objective(dat$x, dat$y, dat$z, screen$h, 0.5, screen$beta)
    of <- weighted_qr_objective(dat$x, dat$y, dat$z, screen$h, 0.5, full$beta)
    objective_gap <- max(abs(os - of) / pmax(1, abs(of)))
    kkt <- max(weighted_qr_kkt_residual(dat$x, dat$y, dat$z, screen$h, 0.5,
                                        screen$beta, screen$H_seq))
    h_match <- row_set_match(screen$H_seq, full$H_seq)
    relative <- abs(screen$ll_est - direct$ll_est) / pmax(abs(direct$ll_est), 1e-10)
    retained <- screen$diagnostics[, "first_retained_size"]
    retained <- retained[is.finite(retained)]
    audit_ok <- h_match && objective_gap <= 1e-8 && kkt <= 1e-7 && max(relative) <= 1e-6
    data.frame(case = case_id, n = n, tau = 0.5, rep_id = rep_id,
      original_seed = original_seed, seed = seed, attempt = attempt,
      dataset_md5 = dataset_md5(dat), method_order = paste(order, collapse = ","),
      threshold_factor = factor, direct_time = fit$direct$elapsed,
      screen_time = fit[[label]]$elapsed, direct_ok = TRUE, screen_ok = TRUE,
      audit_ok = audit_ok, h_set_match = h_match, objective_gap = objective_gap,
      kkt_residual = kkt, fit_discrepancy = mean(relative),
      max_relative_fit_difference = max(relative),
      max_abs_fit_difference = max(abs(screen$ll_est - direct$ll_est)),
      repair_total = sum(screen$diagnostics[, "repairs"]),
      repaired_points = sum(screen$diagnostics[, "repairs"] > 0L),
      full_active_recovery_total = sum(screen$diagnostics[, "full_recovery"]),
      independent_init_total = sum(screen$diagnostics[, "independent"]),
      threshold_expansion_total = sum(screen$diagnostics[, "threshold_expansion_steps"], na.rm = TRUE),
      retained_size_mean = mean(retained), retained_size_median = median(retained),
      retained_size_p95 = unname(quantile(retained, 0.95)), retained_size_max = max(retained),
      retained_fraction_mean = mean(retained) / n,
      underflow_total = sum(screen$solver_info$underflow_count),
      underflow_max = max(screen$solver_info$underflow_count),
      mst_seconds = screen$solver_info$mst_seconds, max_edge = max(screen$edge_weight))
  })
  rows <- do.call(rbind, rows)
  write.csv(rows, file.path(partial_dir, sprintf("task_%05d.csv", task_id)),
            row.names = FALSE, quote = TRUE)
  if (any(!rows$audit_ok | !rows$h_set_match)) quit(status = 2L)
  break
}
