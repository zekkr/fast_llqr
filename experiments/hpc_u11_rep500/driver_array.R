#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

get_int <- function(key, default = NA_integer_) {
  value <- Sys.getenv(key, "")
  if (nzchar(value)) as.integer(value) else as.integer(default)
}
get_num <- function(key, default = NA_real_) {
  value <- Sys.getenv(key, "")
  if (nzchar(value)) as.numeric(value) else as.numeric(default)
}

project <- normalizePath(Sys.getenv("SSQR_PROJECT_ROOT", getwd()), mustWork = TRUE)
setwd(project)
experiment_dir <- Sys.getenv("SSQR_EXPERIMENT_DIR", "experiments/hpc_u11_rep500")
model <- tolower(Sys.getenv("SSQR_MODEL", ""))
case_id <- get_int("SSQR_CASE")
tau <- get_num("SSQR_TAU")
n <- get_int("SSQR_N")
num_rep <- get_int("SSQR_NUM_REP", 500L)
seed_base <- get_int("SSQR_SEED_BASE", 2025L)
task_id <- get_int("SLURM_ARRAY_TASK_ID", 1L)
chunk <- get_int("SSQR_CHUNK_SIZE", 1L)
cores <- get_int("SLURM_CPUS_PER_TASK", 1L)
run_tag <- Sys.getenv("SSQR_RUN_TAG", "")
output_root <- Sys.getenv("SSQR_OUTPUT_ROOT", "")
smoke <- identical(Sys.getenv("SSQR_ALLOW_SMOKE", "0"), "1")

stopifnot(
  model == "llqr", case_id == 2L,
  tau %in% c(0.2, 0.5, 0.8), n > 0L,
  smoke || n %in% c(1000L, 2000L, 5000L, 10000L),
  smoke || num_rep == 500L, smoke || seed_base == 2025L,
  nzchar(run_tag), nzchar(output_root), chunk > 0L, cores > 0L
)

all_ids <- seq_len(num_rep)
first_id <- (task_id - 1L) * chunk + 1L
last_id <- min(task_id * chunk, num_rep)
if (first_id > last_id) quit(save = "no")
rep_ids <- all_ids[first_id:last_id]

source(file.path(experiment_dir, "design.R"))
source(file.path(experiment_dir, "metrics.R"))
source(file.path(experiment_dir, "retry.R"))
source(file.path(experiment_dir, "model_methods.R"))
source(file.path(experiment_dir, "adapters.R"))
source("R/llqr_functions.R")
source("R/tvcqr_functions.R")

core <- ssqr_load(Sys.getenv("SSQR_BUILD_MODE", "optimized"))
lean_dll <- load_lean_kernels(experiment_dir)

node_name <- unname(Sys.info()[["nodename"]])
cpu_model <- tryCatch({
  if (.Platform$OS.type == "unix" && file.exists("/proc/cpuinfo")) {
    line <- grep("^model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)[1L]
    trimws(sub("^[^:]+:", "", line))
  } else {
    unname(Sys.info()[["machine"]])
  }
}, error = function(e) NA_character_)

safe_max_int <- function(x) {
  x <- as.integer(x)
  if (!length(x) || all(is.na(x))) NA_integer_ else max(x, na.rm = TRUE)
}
safe_sum_int <- function(x) {
  x <- as.integer(x)
  if (!length(x) || all(is.na(x))) NA_integer_ else as.integer(sum(x, na.rm = TRUE))
}

empty_compact_diagnostics <- function() {
  list(
    first_tableau_rows_max = NA_integer_, first_retained_size_max = NA_integer_,
    initial_threshold_hits_max = NA_integer_, effective_threshold_hits_max = NA_integer_,
    threshold_expansion_total = NA_integer_, threshold_expansion_points = NA_integer_,
    threshold_expansion_steps_max = NA_integer_, threshold_initial = NA_real_,
    threshold_effective_max = NA_real_,
    basis_forced_rows_max = NA_integer_, first_aggregate_rows_max = NA_integer_,
    retained_decomposition_ok = NA, first_pass_count = NA_integer_,
    first_pass_total = NA_integer_, first_pass_rate = NA_real_,
    repaired_points = NA_integer_, repair_total = NA_integer_,
    internal_recovery_total = NA_integer_, first_full_m_recovery = NA_integer_,
    certificate_hits = NA_integer_, residual_rows_checked = NA_integer_
  )
}

u11_diagnostics <- function(diag) {
  if (is.null(diag)) return(empty_compact_diagnostics())
  if (nrow(diag) == 1L) {
    out <- empty_compact_diagnostics()
    zero <- setdiff(names(out), c("first_pass_rate", "threshold_initial", "threshold_effective_max"))
    out[zero] <- 0L
    out$retained_decomposition_ok <- TRUE
    out$first_full_m_recovery <- as.integer(diag[1L, "first_full_m_recovery"])
    out$internal_recovery_total <- out$first_full_m_recovery
    return(out)
  }
  x <- diag[-1L, , drop = FALSE]
  repairs <- as.integer(x[, "repairs"])
  decomp <- x[, "first_tableau_rows"] ==
    x[, "first_retained_size"] + x[, "first_aggregate_rows"]
  basis <- x[, "basis_forced_rows"] ==
    x[, "first_retained_size"] - x[, "effective_threshold_hits"]
  list(
    first_tableau_rows_max = safe_max_int(x[, "first_tableau_rows"]),
    first_retained_size_max = safe_max_int(x[, "first_retained_size"]),
    initial_threshold_hits_max = safe_max_int(x[, "initial_threshold_hits"]),
    effective_threshold_hits_max = safe_max_int(x[, "effective_threshold_hits"]),
    threshold_expansion_total = safe_sum_int(x[, "threshold_expansion_steps"]),
    threshold_expansion_points = as.integer(sum(x[, "threshold_expansion_steps"] > 0L)),
    threshold_expansion_steps_max = safe_max_int(x[, "threshold_expansion_steps"]),
    basis_forced_rows_max = safe_max_int(x[, "basis_forced_rows"]),
    first_aggregate_rows_max = safe_max_int(x[, "first_aggregate_rows"]),
    retained_decomposition_ok = isTRUE(all(decomp & basis, na.rm = FALSE)),
    first_pass_count = as.integer(sum(repairs == 0L)),
    first_pass_total = as.integer(length(repairs)),
    first_pass_rate = mean(repairs == 0L),
    repaired_points = as.integer(sum(repairs > 0L)),
    repair_total = safe_sum_int(repairs),
    internal_recovery_total = safe_sum_int(c(
      x[, "full_recovery"], diag[1L, "first_full_m_recovery"]
    )),
    first_full_m_recovery = as.integer(diag[1L, "first_full_m_recovery"]),
    certificate_hits = safe_sum_int(x[, "certificate_hits"]),
    residual_rows_checked = safe_sum_int(x[, "residual_rows"])
  )
}

run_method <- function(method, position, dat) {
  ne <- length(dat$z)
  invisible(gc())
  stage <- "fit"
  started <- proc.time()[["elapsed"]]
  tryCatch({
    raw <- if (model == "llqr") {
      h <- llqr_default_bandwidth(dat$x, dat$y, tau, h = NULL, case = case_id, h.factor = 1)
      switch(method,
        direct_baseline = llqr_local_fit(dat$x, dat$y, tau, dat$z, h, case_id, 1),
        lean_seq = run_llqr_lean(lean_dll$llqr, dat$x, dat$y, dat$z, tau, h, case_id),
        unified_u11 = ssqr_llqr(core, dat$x, dat$y, tau, case_id, h, dat$z, Mm.factor = 0.1)
      )
    } else {
      switch(method,
        direct_baseline = tvc_rq(dat$x, dat$y, tau, h = NULL),
        lean_seq = run_tvcqr_lean(lean_dll$tvcqr, dat$x, dat$y, tau),
        unified_u11 = ssqr_tvcqr(core, dat$x, dat$y, tau, Mm.factor = 1e-5)
      )
    }
    elapsed <- proc.time()[["elapsed"]] - started
    stage <- "postprocess"

    if (model == "llqr") {
      estimate <- if (method == "unified_u11") raw$estimate else raw$ll_est
      derivative <- if (method == "unified_u11") raw$derivative else raw$d_ll_est
      h_path <- if (method == "direct_baseline") NULL else
        if (method == "unified_u11") raw$H else raw$H_mat
      state_finite <- length(derivative) == ne && all(is.finite(derivative))
    } else {
      estimate <- if (method == "unified_u11") raw$estimate else raw$theta_ll_est
      h_path <- if (method == "direct_baseline") NULL else
        if (method == "unified_u11") raw$H else raw$H_mat
      state_finite <- if (method == "unified_u11") {
        length(raw$beta) == 8L * n && all(is.finite(raw$beta))
      } else TRUE
    }

    expected_k <- if (model == "llqr") ne else 4L * n
    finite <- length(estimate) == expected_k && all(is.finite(estimate)) && state_finite
    ierr <- if (method == "unified_u11") as.integer(raw$ierr) else NA_integer_
    failed <- if (method == "unified_u11") as.integer(raw$failed_eval) else NA_integer_
    iteration_limit_hit <- if (method == "unified_u11") {
      any(raw$diagnostics[, "iterations"] >= 1000000L)
    } else if (method == "lean_seq") any(raw$it_num >= 1000000L) else FALSE
    solver_ok <- finite && !iteration_limit_hit && (is.na(ierr) || ierr == 0L) && (is.na(failed) || failed == 0L)
    compact <- if (method == "unified_u11") {
      u11_diagnostics(raw$diagnostics)
    } else {
      empty_compact_diagnostics()
    }
    if (method == "unified_u11") {
      compact$threshold_initial <- as.numeric(raw$threshold_initial)
      compact$threshold_effective_max <- compact$threshold_initial *
        1.5^compact$threshold_expansion_steps_max
    }

    c(list(
      method = method, position = as.integer(position), elapsed_sec = elapsed,
      estimate = estimate, H = h_path, solver_ok = solver_ok,
      finite_estimate = finite, ierr = ierr, failed_eval = failed, iteration_limit_hit = iteration_limit_hit,
      error_message = if (finite) NA_character_ else "nonfinite_or_dimension",
      threw_error = FALSE, error_stage = NA_character_
    ), compact)
  }, error = function(e) {
    elapsed <- proc.time()[["elapsed"]] - started
    empty_fit(method, position, elapsed, conditionMessage(e), stage)
  })
}

run_attempt <- function(rep_id, seed) {
  generated_at <- proc.time()[["elapsed"]]
  generated <- tryCatch({
    if (model == "llqr") {
      value <- generate_logistic_case2(n, seed)
      value
    } else {
      generate_ts(n, case_id, seed, J = 100L, burn_in = 500L)
    }
  }, error = function(e) e)
  generation_sec <- proc.time()[["elapsed"]] - generated_at
  ne <- if (inherits(generated, "error")) NA_integer_ else length(generated$z)

  shift <- (rep_id - 1L) %% length(experiment_methods)
  order_now <- experiment_methods[
    (seq_along(experiment_methods) + shift - 1L) %% length(experiment_methods) + 1L
  ]
  fits <- if (inherits(generated, "error")) {
    lapply(seq_along(order_now), function(pos) {
      empty_fit(order_now[[pos]], pos, message = conditionMessage(generated), stage = "data_generation")
    })
  } else {
    lapply(seq_along(order_now), function(pos) run_method(order_now[[pos]], pos, generated))
  }
  names(fits) <- vapply(fits, `[[`, character(1L), "method")
  baseline <- fits$direct_baseline
  lean <- fits$lean_seq
  expected_k <- if (model == "llqr") ne else 4L * n
  q <- if (model == "llqr") 2L else 8L

  do.call(rbind, lapply(experiment_methods, function(method) {
    fit <- fits[[method]]
    h_audit <- audit_h_path(method, fit$solver_ok, fit$H, lean$solver_ok, lean$H, n, q, ne)
    discrepancy <- if (!isTRUE(baseline$solver_ok)) {
      list(value = NA_real_, status = "baseline_failed")
    } else if (!isTRUE(fit$solver_ok)) {
      list(value = NA_real_, status = "method_failed")
    } else {
      paper_discrepancy(fit$estimate, baseline$estimate, expected_k)
    }
    data.frame(
      run_tag, model, case = case_id, tau, n, n_eval = ne,
      n_interior = if (inherits(generated, "error")) NA_integer_ else generated$n_interior,
      grid_placeholder = if (inherits(generated, "error")) NA else generated$placeholder,
      n_transitions = max(0L, ne - 1L), rep_id, seed, method,
      method_position = fit$position, elapsed_sec = fit$elapsed_sec,
      solver_ok = fit$solver_ok, accepted_ok = h_audit$accepted,
      finite_estimate = fit$finite_estimate, iteration_limit_hit = fit$iteration_limit_hit, discrepancy = discrepancy$value,
      discrepancy_status = discrepancy$status, ierr = fit$ierr,
      failed_eval = fit$failed_eval, h_check_status = h_audit$status,
      h_match_vs_lean = h_audit$match,
      h_mismatch_eval_count = h_audit$count,
      first_h_mismatch_eval = h_audit$first,
      first_tableau_rows_max = fit$first_tableau_rows_max,
      first_retained_size_max = fit$first_retained_size_max,
      initial_threshold_hits_max = fit$initial_threshold_hits_max,
      effective_threshold_hits_max = fit$effective_threshold_hits_max,
      threshold_expansion_total = fit$threshold_expansion_total,
      threshold_expansion_points = fit$threshold_expansion_points,
      threshold_expansion_steps_max = fit$threshold_expansion_steps_max,
      threshold_initial = fit$threshold_initial,
      threshold_effective_max = fit$threshold_effective_max,
      basis_forced_rows_max = fit$basis_forced_rows_max,
      first_aggregate_rows_max = fit$first_aggregate_rows_max,
      retained_decomposition_ok = fit$retained_decomposition_ok,
      first_pass_count = fit$first_pass_count,
      first_pass_total = fit$first_pass_total,
      first_pass_rate = fit$first_pass_rate,
      repaired_points = fit$repaired_points,
      repair_total = fit$repair_total,
      internal_recovery_total = fit$internal_recovery_total,
      first_full_m_recovery = fit$first_full_m_recovery,
      certificate_hits = fit$certificate_hits,
      residual_rows_checked = fit$residual_rows_checked,
      fallback_triggered = FALSE, error_message = fit$error_message,
      threw_error = fit$threw_error, error_stage = fit$error_stage,
      data_generation_sec = generation_sec, node_name, cpu_model,
      stringsAsFactors = FALSE
    )
  }))
}

run_rep <- function(rep_id) {
  value <- run_attempt_chain(model, rep_id, seed_base, run_attempt)
  config_tag <- sprintf("case%d_tau%02d_n%d", case_id, round(100 * tau), n)
  partial_dir <- file.path(output_root, run_tag, model, config_tag, "partials")
  dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(partial_dir, sprintf("rep%04d.rds", rep_id))
  temporary <- paste0(path, ".", Sys.getpid(), ".tmp")
  saveRDS(value, temporary, compress = FALSE)
  if (!file.rename(temporary, path)) stop("atomic save failed")
  data.frame(
    rep_id, ok = all(value$final$accepted_ok %in% TRUE),
    attempt = max(value$attempts$attempt), seed = value$final$seed[[1L]]
  )
}

cat(sprintf(
  "%s case=%d tau=%.1f n=%d reps=%d-%d workers=%d node=%s\n",
  model, case_id, tau, n, min(rep_ids), max(rep_ids),
  min(cores, length(rep_ids)), node_name
))
status <- if (.Platform$OS.type == "unix" && cores > 1L && length(rep_ids) > 1L) {
  parallel::mclapply(
    rep_ids, run_rep, mc.cores = min(cores, length(rep_ids)), mc.preschedule = FALSE
  )
} else {
  lapply(rep_ids, run_rep)
}
status <- do.call(rbind, status)
cat(sprintf("completed=%d accepted=%d\n", nrow(status), sum(status$ok)))
if (nrow(status) != length(rep_ids)) quit(save = "no", status = 2L)
