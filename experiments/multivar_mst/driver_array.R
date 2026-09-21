#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

env_int <- function(name, default = NA_integer_) {
  value <- Sys.getenv(name, "")
  if (nzchar(value)) as.integer(value) else as.integer(default)
}
project <- normalizePath(Sys.getenv("MST_PROJECT_ROOT", getwd()), mustWork = TRUE)
setwd(project)
case_id <- env_int("MST_CASE")
n <- env_int("MST_N")
num_rep <- env_int("MST_NUM_REP", 100L)
seed_base <- env_int("MST_SEED_BASE", 2025L)
task_id <- env_int("SLURM_ARRAY_TASK_ID", 1L)
chunk <- env_int("MST_CHUNK_SIZE", 1L)
run_tag <- Sys.getenv("MST_RUN_TAG", "")
output_root <- Sys.getenv("MST_OUTPUT_ROOT", "results/hpc/multivar_mst_runs")
allow_smoke <- identical(Sys.getenv("MST_ALLOW_SMOKE", "0"), "1")
stopifnot(case_id %in% 1:2, n > 4L, nzchar(run_tag), chunk > 0L,
          allow_smoke || n %in% c(500L, 1000L, 2000L),
          allow_smoke || num_rep == 100L, allow_smoke || seed_base == 2025L)

first_id <- (task_id - 1L) * chunk + 1L
last_id <- min(task_id * chunk, num_rep)
if (first_id > last_id) quit(save = "no")
rep_ids <- first_id:last_id

source("experiments/multivar_mst/design.R")
source("experiments/multivar_mst/solver.R")
stopifnot(requireNamespace("quantreg", quietly = TRUE))
kernel <- load_multivar_mst_kernel("experiments/multivar_mst/build/multivar_mst_optimized.so")

run_dir <- file.path(output_root, run_tag)
config_dir <- file.path(run_dir, sprintf("case%d_n%d", case_id, n))
partial_dir <- file.path(config_dir, "partial")
dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
meta_dir <- file.path(run_dir, "_run_meta", "workers")
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

cpu_model <- tryCatch({
  if (!file.exists("/proc/cpuinfo")) stop("cpuinfo unavailable")
  line <- grep("^model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)[1L]
  trimws(sub("^[^:]+:", "", line))
}, error = function(e) unname(Sys.info()[["machine"]]))
thread_env <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"))
if (Sys.getenv("MST_REQUIRE_HPC", "0") == "1") {
  stopifnot(getRversion() == "4.3.1", packageVersion("quantreg") == "5.94",
            grepl("6258R", cpu_model), all(thread_env == "1"))
}
saveRDS(list(session = sessionInfo(), system = Sys.info(), cpu = cpu_model,
             thread_env = thread_env, git_sha = Sys.getenv("MST_PUSHED_SHA")),
        file.path(meta_dir, sprintf("case%d_n%d_task%d_%s.rds", case_id, n, task_id,
                                    unname(Sys.info()[["nodename"]]))))

row_set_match <- function(a, b) {
  nrow(a) == nrow(b) && all(vapply(seq_len(nrow(a)), function(j)
    setequal(a[j, ], b[j, ]), logical(1L)))
}
dataset_md5 <- function(dat) {
  path <- tempfile("mst-data-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(x = dat$x, y = dat$y), path, version = 3)
  unname(tools::md5sum(path))
}
timed_fit <- function(method, dat) {
  invisible(gc())
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(
    if (method == "direct_fit") llqr_direct_multivar(dat$x, dat$y, 0.5, dat$z) else
      llqr_seq_screen_mst(dat$x, dat$y, 0.5, dat$z, kernel = kernel, diagnostics = TRUE),
    error = identity)
  list(value = value, elapsed = proc.time()[["elapsed"]] - started)
}

run_replication <- function(rep_id) {
  original_seed <- seed_base + rep_id
  for (attempt in seq_len(20L)) {
    seed <- original_seed + (attempt - 1L) * 1000000L
    dat <- simulate_multivar_llqr(n, case_id, seed)
    order <- if ((rep_id + attempt) %% 2L == 0L)
      c("direct_fit", "seq_screen_mst") else c("seq_screen_mst", "direct_fit")
    fit <- setNames(vector("list", 2L), order)
    for (method in order) fit[[method]] <- timed_fit(method, dat)
    direct_error <- inherits(fit$direct_fit$value, "error")
    screen_error <- inherits(fit$seq_screen_mst$value, "error")
    if (direct_error && attempt < 20L) next
    if (direct_error || screen_error) {
      return(data.frame(
        case = case_id, n = n, tau = 0.5, rep_id = rep_id, original_seed = original_seed,
        seed = seed, attempt = attempt, dataset_md5 = dataset_md5(dat),
        method_order = paste(order, collapse = ","), direct_time = fit$direct_fit$elapsed,
        screen_time = fit$seq_screen_mst$elapsed, direct_ok = !direct_error,
        screen_ok = !screen_error, audit_ok = FALSE, h_set_match = FALSE,
        objective_gap = NA_real_, kkt_residual = NA_real_, fit_discrepancy = NA_real_,
        max_abs_fit_difference = NA_real_, repair_total = NA_integer_, repaired_points = NA_integer_,
        full_active_recovery_total = NA_integer_, independent_init_total = NA_integer_,
        threshold_expansion_total = NA_integer_, underflow_total = NA_integer_,
        underflow_max = NA_integer_, mst_seconds = NA_real_, max_edge = NA_real_,
        direct_error = if (direct_error) conditionMessage(fit$direct_fit$value) else "",
        screen_error = if (screen_error) conditionMessage(fit$seq_screen_mst$value) else ""))
    }
    direct <- fit$direct_fit$value; screen <- fit$seq_screen_mst$value
    full <- tryCatch(llqr_seq_screen_mst(dat$x, dat$y, 0.5, dat$z, h = screen$h,
                                         kernel = kernel, diagnostics = TRUE,
                                         full_active = TRUE), error = identity)
    if (inherits(full, "error")) {
      return(data.frame(
        case = case_id, n = n, tau = 0.5, rep_id = rep_id, original_seed = original_seed,
        seed = seed, attempt = attempt, dataset_md5 = dataset_md5(dat),
        method_order = paste(order, collapse = ","), direct_time = fit$direct_fit$elapsed,
        screen_time = fit$seq_screen_mst$elapsed, direct_ok = TRUE, screen_ok = TRUE,
        audit_ok = FALSE, h_set_match = FALSE, objective_gap = NA_real_,
        kkt_residual = NA_real_, fit_discrepancy = NA_real_, max_abs_fit_difference = NA_real_,
        repair_total = sum(screen$diagnostics[, "repairs"]),
        repaired_points = sum(screen$diagnostics[, "repairs"] > 0L),
        full_active_recovery_total = sum(screen$diagnostics[, "full_recovery"]),
        independent_init_total = sum(screen$diagnostics[, "independent"]),
        threshold_expansion_total = sum(screen$diagnostics[, "threshold_expansion_steps"], na.rm = TRUE),
        underflow_total = sum(screen$solver_info$underflow_count),
        underflow_max = max(screen$solver_info$underflow_count),
        mst_seconds = screen$solver_info$mst_seconds, max_edge = max(screen$edge_weight),
        direct_error = "", screen_error = paste("full audit:", conditionMessage(full))))
    }
    os <- weighted_qr_objective(dat$x, dat$y, dat$z, screen$h, 0.5, screen$beta)
    of <- weighted_qr_objective(dat$x, dat$y, dat$z, screen$h, 0.5, full$beta)
    objective_gap <- max(abs(os - of) / pmax(1, abs(of)))
    kkt <- max(weighted_qr_kkt_residual(dat$x, dat$y, dat$z, screen$h, 0.5,
                                        screen$beta, screen$H_seq))
    h_match <- row_set_match(screen$H_seq, full$H_seq)
    relative <- abs(screen$ll_est - direct$ll_est) / pmax(abs(direct$ll_est), 1e-10)
    fit_discrepancy <- mean(relative)
    max_abs <- max(abs(screen$ll_est - direct$ll_est))
    audit_ok <- h_match && objective_gap <= 1e-8 && kkt <= 1e-7 &&
      max(relative) <= 1e-6
    return(data.frame(
      case = case_id, n = n, tau = 0.5, rep_id = rep_id, original_seed = original_seed,
      seed = seed, attempt = attempt, dataset_md5 = dataset_md5(dat),
      method_order = paste(order, collapse = ","), direct_time = fit$direct_fit$elapsed,
      screen_time = fit$seq_screen_mst$elapsed, direct_ok = TRUE, screen_ok = TRUE,
      audit_ok = audit_ok, h_set_match = h_match, objective_gap = objective_gap,
      kkt_residual = kkt, fit_discrepancy = fit_discrepancy,
      max_abs_fit_difference = max_abs,
      repair_total = sum(screen$diagnostics[, "repairs"]),
      repaired_points = sum(screen$diagnostics[, "repairs"] > 0L),
      full_active_recovery_total = sum(screen$diagnostics[, "full_recovery"]),
      independent_init_total = sum(screen$diagnostics[, "independent"]),
      threshold_expansion_total = sum(screen$diagnostics[, "threshold_expansion_steps"], na.rm = TRUE),
      underflow_total = sum(screen$solver_info$underflow_count),
      underflow_max = max(screen$solver_info$underflow_count),
      mst_seconds = screen$solver_info$mst_seconds, max_edge = max(screen$edge_weight),
      direct_error = "", screen_error = ""))
  }
  stop("unreachable")
}

rows <- do.call(rbind, lapply(rep_ids, run_replication))
path <- file.path(partial_dir, sprintf("task_%05d.csv", task_id))
write.csv(rows, path, row.names = FALSE, quote = TRUE)
if (any(!rows$direct_ok | !rows$screen_ok | !rows$audit_ok)) quit(status = 2L)
