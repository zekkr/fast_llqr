#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

as_int <- function(name, default = NA_integer_) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  as.integer(value)
}

as_num <- function(name, default = NA_real_) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  as.numeric(value)
}

parse_rep_ids <- function(num_rep, task_id, chunk_size) {
  explicit <- trimws(Sys.getenv("FASTQR_REP_ID_LIST", unset = ""))
  if (nzchar(explicit)) {
    ids <- suppressWarnings(as.integer(strsplit(explicit, "[,[:space:]]+", perl = TRUE)[[1L]]))
    ids <- sort(unique(ids[is.finite(ids) & ids >= 1L & ids <= num_rep]))
    if (!length(ids)) stop("FASTQR_REP_ID_LIST contains no valid replication IDs.")
  } else {
    ids <- seq_len(num_rep)
  }
  first <- (task_id - 1L) * chunk_size + 1L
  last <- min(task_id * chunk_size, length(ids))
  if (first > last) integer() else ids[first:last]
}

paper_discrepancy <- function(estimate, baseline, expected_k) {
  a <- as.numeric(estimate)
  b <- as.numeric(baseline)
  if (length(a) != expected_k || length(b) != expected_k) {
    return(list(value = NA_real_, status = sprintf(
      "dimension_mismatch:candidate=%d:baseline=%d:expected=%d",
      length(a), length(b), expected_k
    )))
  }
  if (any(!is.finite(a)) || any(!is.finite(b))) {
    return(list(value = NA_real_, status = "nonfinite_estimate"))
  }
  value <- mean(abs(a - b) / pmax(abs(b), 1e-10))
  list(value = value, status = if (is.finite(value)) "ok" else "nonfinite_discrepancy")
}

empty_method_result <- function(method, position, elapsed = NA_real_, message = NA_character_) {
  list(
    method = method,
    position = as.integer(position),
    elapsed_sec = as.numeric(elapsed),
    estimate = NULL,
    method_ok = FALSE,
    finite_estimate = FALSE,
    ierr = NA_integer_,
    failed_eval = NA_integer_,
    first_point_rows = NA_integer_,
    repair_total = NA_integer_,
    independent_init_count = NA_integer_,
    full_active_recovery_count = NA_integer_,
    first_point_full_m = NA,
    returned_backend = NA_character_,
    error_message = as.character(message)
  )
}

summarize_fit <- function(method, position, elapsed, estimate, raw = NULL,
                          returned_backend) {
  ierr <- if (!is.null(raw$ierr)) as.integer(raw$ierr[[1L]]) else NA_integer_
  failed_eval <- if (!is.null(raw$failed_eval)) as.integer(raw$failed_eval[[1L]]) else NA_integer_
  first_n_sub <- if (!is.null(raw$first_n_sub)) as.integer(raw$first_n_sub) else integer()
  repair_count <- if (!is.null(raw$repair_count)) as.integer(raw$repair_count) else integer()
  init_mode <- if (!is.null(raw$init_mode)) as.integer(raw$init_mode) else integer()
  finite_estimate <- length(estimate) > 0L && all(is.finite(as.numeric(estimate)))
  kernel_ok <- (is.na(ierr) || ierr == 0L) && (is.na(failed_eval) || failed_eval == 0L)
  list(
    method = method,
    position = as.integer(position),
    elapsed_sec = as.numeric(elapsed),
    estimate = estimate,
    method_ok = isTRUE(kernel_ok && finite_estimate),
    finite_estimate = finite_estimate,
    ierr = ierr,
    failed_eval = failed_eval,
    first_point_rows = if (length(first_n_sub)) first_n_sub[[1L]] else NA_integer_,
    repair_total = if (length(repair_count)) sum(repair_count, na.rm = FALSE) else NA_integer_,
    independent_init_count = if (length(init_mode)) sum(init_mode == 2L, na.rm = TRUE) else NA_integer_,
    full_active_recovery_count = if (length(init_mode)) sum(init_mode == 3L, na.rm = TRUE) else NA_integer_,
    first_point_full_m = NA,
    returned_backend = returned_backend,
    error_message = NA_character_
  )
}

project_dir <- normalizePath(Sys.getenv("FASTQR_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
setwd(project_dir)

model <- tolower(Sys.getenv("FASTQR_MODEL", unset = ""))
if (!model %in% c("llqr", "tvcqr")) stop("FASTQR_MODEL must be llqr or tvcqr.")

case_id <- as_int("FASTQR_CASE")
tau <- as_num("FASTQR_TAU")
n <- as_int("FASTQR_N")
num_rep <- as_int("FASTQR_NUM_REP", 100L)
seed_base <- as_int("FASTQR_SEED_BASE", 2026L)
task_id <- as_int("SLURM_ARRAY_TASK_ID", 1L)
ncores <- as_int("SLURM_CPUS_PER_TASK", 1L)
chunk_size <- as_int("FASTQR_CHUNK_SIZE", ncores)
run_tag <- Sys.getenv("FASTQR_RUN_TAG", unset = "")
base_dir <- Sys.getenv(
  "FASTQR_V4142_BASE_DIR",
  unset = "data/v41_v42_rep100"
)

if (!case_id %in% c(1L, 2L)) stop("FASTQR_CASE must be 1 or 2.")
if (!is.finite(tau) || tau <= 0 || tau >= 1) stop("FASTQR_TAU must be in (0,1).")
if (!is.finite(n) || n <= 0L) stop("FASTQR_N must be positive.")
if (!is.finite(num_rep) || num_rep <= 0L) stop("FASTQR_NUM_REP must be positive.")
if (!is.finite(seed_base)) stop("FASTQR_SEED_BASE must be finite.")
if (!nzchar(run_tag)) stop("FASTQR_RUN_TAG is required.")

rep_ids <- parse_rep_ids(num_rep, task_id, chunk_size)
if (!length(rep_ids)) quit(save = "no", status = 0L)

experiment_dir <- file.path(project_dir, "experiments", "hpc_v41_v42_rep100")
build_dir <- file.path(experiment_dir, "build")

if (model == "llqr") {
  source(file.path(project_dir, "R", "llqr_functions.R"))
  dlls <- list(
    v41 = dyn.load(file.path(build_dir, "llqr_v41.so")),
    v42 = dyn.load(file.path(build_dir, "llqr_v42.so")),
    lean_seq = dyn.load(file.path(build_dir, "llqr_lean_seq.so"))
  )
} else {
  source(file.path(project_dir, "R", "tvcqr_functions.R"))
  dlls <- list(
    v41 = dyn.load(file.path(build_dir, "tvcqr_v41.so")),
    v42 = dyn.load(file.path(build_dir, "tvcqr_v42.so")),
    lean_seq = dyn.load(file.path(build_dir, "tvcqr_lean_seq.so"))
  )
}

run_llqr_ppro <- function(dll, x, y, z, tau, h, case_id) {
  m <- length(y)
  .Fortran(
    "llqr_ppro_fortran", PACKAGE = dll[["name"]],
    x = as.double(x), y = as.double(y), z = as.double(z),
    m = as.integer(m), nvar = 1L, rounds = as.integer(m),
    tau = as.double(tau), h = as.double(h), tol = 1e-14,
    maxit = 1000000L, Mm_factor = 0.1, case_int = as.integer(case_id),
    bland_int = 0L, min_subsample_size_in = 1L,
    ll_est = double(m), d_ll_est = double(m), it_num = integer(m),
    residual_est = double(m), H_mat = integer(2L * m),
    first_n_sub = integer(m), repair_count = integer(m),
    final_n_sub = integer(m), init_mode = integer(m),
    init_trigger = integer(m), ierr = integer(1L), failed_eval = integer(1L),
    always_same_h_refit_int = 1L, threshold_lower_bound_int = 0L,
    threshold_scale_mode_int = 1L
  )
}

run_llqr_seq <- function(dll, x, y, z, tau, h, case_id) {
  m <- length(y)
  .Fortran(
    "llqr_seq_fortran", PACKAGE = dll[["name"]],
    x = as.double(x), y = as.double(y), z = as.double(z),
    m = as.integer(m), nvar = 1L, rounds = as.integer(m),
    tau = as.double(tau), h = as.double(h), tol = 1e-14,
    maxit = 1000000L, case_int = as.integer(case_id), bland_int = 0L,
    ll_est = double(m), d_ll_est = double(m), it_num = integer(m),
    residual_est = double(m), H_mat = integer(2L * m)
  )
}

run_tvcqr_ppro <- function(dll, x, y, tau) {
  m <- nrow(x)
  nvar <- ncol(x)
  q <- 2L * (nvar + 1L)
  .Fortran(
    "tvcqr_seq_ppro_fortran", PACKAGE = dll[["name"]],
    x = as.double(x), y = as.double(y), m = as.integer(m),
    nvar = as.integer(nvar), tau = as.double(tau), h = 0.0,
    h_factor = 1.0, tol = 1e-14, maxit = 1000000L,
    bland_int = 0L, Mm_factor = 1e-5, eps = 1e-6,
    store_residual_int = 0L, debug_int = 0L,
    theta_ll_est = double(m * (nvar + 1L)),
    beta_full_est = double(m * q), it_num = integer(m),
    residual_est = double(1L), M_out = double(1L),
    first_n_sub = integer(m), repair_count = integer(m),
    final_n_sub = integer(m), init_mode = integer(m),
    init_trigger = integer(m), H_seq = integer(m * q),
    same_h_refit_attempted = integer(m),
    same_h_refit_recovered = integer(m),
    ierr = integer(1L), failed_eval = integer(1L),
    min_subsample_size_in = 1L, always_same_h_refit_int = 1L,
    threshold_lower_bound_int = 0L, threshold_scale_mode_int = 1L
  )
}

run_tvcqr_seq <- function(dll, x, y, tau) {
  m <- nrow(x)
  nvar <- ncol(x)
  q <- 2L * (nvar + 1L)
  .Fortran(
    "tvcqr_seq_fortran", PACKAGE = dll[["name"]],
    x = as.double(x), y = as.double(y), m = as.integer(m),
    nvar = as.integer(nvar), tau = as.double(tau), h = 0.0,
    tol = 1e-14, maxit = 1000000L, bland_int = 0L,
    theta_ll_est = double(m * (nvar + 1L)), it_num = integer(m),
    residual_est = double(1L), H_mat = integer(m * q)
  )
}

run_method <- function(method, position, data) {
  invisible(gc())
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    if (model == "llqr") {
      raw <- switch(
        method,
        direct_baseline = llqr_local_fit(
          x = data$x, y = data$y, tau = tau, z = data$z,
          h = data$h, case = case_id, h.factor = 1
        ),
        lean_seq = run_llqr_seq(dlls$lean_seq, data$x, data$y, data$z, tau, data$h, case_id),
        v41 = run_llqr_ppro(dlls$v41, data$x, data$y, data$z, tau, data$h, case_id),
        v42 = run_llqr_ppro(dlls$v42, data$x, data$y, data$z, tau, data$h, case_id)
      )
      estimate <- raw$ll_est
    } else {
      raw <- switch(
        method,
        direct_baseline = tvc_rq(x = data$x, y = data$y, tau = tau, h = NULL),
        lean_seq = run_tvcqr_seq(dlls$lean_seq, data$x, data$y, tau),
        v41 = run_tvcqr_ppro(dlls$v41, data$x, data$y, tau),
        v42 = run_tvcqr_ppro(dlls$v42, data$x, data$y, tau)
      )
      estimate <- if (method == "direct_baseline") {
        raw$theta_ll_est
      } else {
        matrix(raw$theta_ll_est, nrow = n, ncol = 4L)
      }
    }
    elapsed <- proc.time()[["elapsed"]] - started
    backend <- if (method %in% c("v41", "v42")) "ppro_kernel" else method
    summarized <- summarize_fit(method, position, elapsed, estimate, raw, backend)
    summarized$first_point_full_m <- if (method %in% c("v41", "v42")) {
      !is.na(summarized$first_point_rows) && summarized$first_point_rows == n
    } else {
      NA
    }
    summarized
  }, error = function(e) {
    elapsed <- proc.time()[["elapsed"]] - started
    empty_method_result(method, position, elapsed, conditionMessage(e))
  })
  result
}

run_replication <- function(rep_id) {
  seed <- as.integer(seed_base + rep_id)
  generated <- tryCatch({
    if (model == "llqr") {
      dat <- generate_data(n = n, case = case_id, seed = seed)
      dat$z <- sort(dat$x)
      dat$h <- llqr_default_bandwidth(
        x = dat$x, y = dat$y, tau = tau, h = NULL,
        case = case_id, h.factor = 1
      )
      dat
    } else {
      generate_ts(n = n, case = case_id, seed = seed, J = 100L, burn_in = 500L)
    }
  }, error = function(e) e)

  base_methods <- c("direct_baseline", "lean_seq", "v41", "v42")
  shift <- (rep_id - 1L) %% length(base_methods)
  method_order <- base_methods[((seq_along(base_methods) + shift - 1L) %% length(base_methods)) + 1L]

  if (inherits(generated, "error")) {
    method_results <- lapply(seq_along(method_order), function(position) {
      empty_method_result(
        method_order[[position]], position,
        message = paste0("data_generation:", conditionMessage(generated))
      )
    })
  } else {
    method_results <- lapply(seq_along(method_order), function(position) {
      run_method(method_order[[position]], position, generated)
    })
  }
  names(method_results) <- vapply(method_results, `[[`, character(1L), "method")

  baseline <- method_results[["direct_baseline"]]
  expected_k <- if (model == "llqr") n else 4L * n
  rows <- lapply(base_methods, function(method) {
    item <- method_results[[method]]
    discrepancy <- if (!isTRUE(baseline$method_ok)) {
      list(value = NA_real_, status = "baseline_failed")
    } else if (!isTRUE(item$method_ok)) {
      list(value = NA_real_, status = "method_failed")
    } else {
      paper_discrepancy(item$estimate, baseline$estimate, expected_k)
    }
    data.frame(
      run_tag = run_tag,
      model = model,
      case = case_id,
      tau = tau,
      n = n,
      rep_id = rep_id,
      seed = seed,
      method = method,
      method_position = item$position,
      elapsed_sec = item$elapsed_sec,
      method_ok = item$method_ok,
      finite_estimate = item$finite_estimate,
      discrepancy = discrepancy$value,
      discrepancy_status = discrepancy$status,
      ierr = item$ierr,
      failed_eval = item$failed_eval,
      first_point_rows = item$first_point_rows,
      first_point_full_m = item$first_point_full_m,
      repair_total = item$repair_total,
      independent_init_count = item$independent_init_count,
      full_active_recovery_count = item$full_active_recovery_count,
      returned_backend = item$returned_backend,
      fallback_triggered = FALSE,
      error_message = item$error_message,
      stringsAsFactors = FALSE
    )
  })
  rows <- do.call(rbind, rows)

  config_tag <- sprintf("case%d_tau%02d_n%d", case_id, as.integer(round(100 * tau)), n)
  partial_dir <- file.path(base_dir, run_tag, model, config_tag, "partials")
  dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
  output <- file.path(partial_dir, sprintf("rep%04d.rds", rep_id))
  temporary <- sprintf("%s.%d.tmp", output, Sys.getpid())
  saveRDS(rows, temporary, compress = FALSE)
  if (!file.rename(temporary, output)) stop("Failed to atomically install partial result: ", output)
  data.frame(rep_id = rep_id, ok = all(rows$method_ok), output = output)
}

cat(sprintf(
  "model=%s case=%d tau=%.1f n=%d reps=%s workers=%d run_tag=%s\n",
  model, case_id, tau, n, paste(range(rep_ids), collapse = "-"),
  min(ncores, length(rep_ids)), run_tag
))

if (.Platform$OS.type == "unix" && ncores > 1L && length(rep_ids) > 1L) {
  status <- parallel::mclapply(
    rep_ids, run_replication,
    mc.cores = min(ncores, length(rep_ids)), mc.preschedule = FALSE
  )
} else {
  status <- lapply(rep_ids, run_replication)
}

status <- do.call(rbind, status)
cat(sprintf("completed=%d all_methods_ok=%d\n", nrow(status), sum(status$ok)))
if (nrow(status) != length(rep_ids)) quit(save = "no", status = 2L)
