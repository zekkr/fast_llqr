#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

project <- normalizePath(Sys.getenv("SSQR_PROJECT_ROOT", getwd()), mustWork = TRUE)
setwd(project)
experiment_dir <- Sys.getenv("SSQR_EXPERIMENT_DIR", "experiments/hpc_u11_rep500")
preflight_dir <- Sys.getenv("SSQR_PREFLIGHT_DIR", "")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
stopifnot(nzchar(preflight_dir), task_id >= 1L, task_id <= 48L)

source(file.path(experiment_dir, "metrics.R"))
source(file.path(experiment_dir, "model_methods.R"))
source(file.path(experiment_dir, "adapters.R"))
source("R/llqr_functions.R")
source("R/tvcqr_functions.R")

grid <- expand.grid(
  n = c(1000L, 2000L, 5000L, 10000L),
  tau = c(0.2, 0.5, 0.8), case = 1:2,
  model = c("llqr", "tvcqr"), stringsAsFactors = FALSE
)
cfg <- grid[task_id, , drop = FALSE]
model <- cfg$model[[1L]]; case_id <- cfg$case[[1L]]
tau <- cfg$tau[[1L]]; n <- cfg$n[[1L]]

dll <- ssqr_load("optimized")
lean_dll <- load_lean_kernels(experiment_dir)
seed <- 2026L
dat <- if (model == "llqr") {
  x <- generate_data(n, case_id, seed)
  x$z <- sort(x$x)
  x
} else {
  generate_ts(n, case_id, seed, J = 100L, burn_in = 500L)
}

if (model == "llqr") {
  h <- llqr_default_bandwidth(dat$x, dat$y, tau, h = NULL, case = case_id, h.factor = 1)
  lean <- run_llqr_lean(lean_dll$llqr, dat$x, dat$y, dat$z, tau, h, case_id)
  run_u11 <- function(flags) {
    Sys.setenv(SSQR_CACHE_FLAGS = as.character(flags), SSQR_PROVIDER_FLAGS = "1")
    ssqr_llqr(dll, dat$x, dat$y, tau, case_id, h, dat$z, Mm.factor = 0.1)
  }
  expected_h <- lean$H_mat; q <- 2L
} else {
  lean <- run_tvcqr_lean(lean_dll$tvcqr, dat$x, dat$y, tau)
  run_u11 <- function(flags) {
    Sys.setenv(SSQR_CACHE_FLAGS = as.character(flags), SSQR_PROVIDER_FLAGS = "1")
    ssqr_tvcqr(dll, dat$x, dat$y, tau, Mm.factor = 1e-5)
  }
  expected_h <- lean$H_mat; q <- 8L
}

fits <- setNames(lapply(c(25L, 27L, 31L), run_u11), c("f25", "f27", "f31"))
ref <- fits$f25
max_abs <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))
h_audit <- lapply(fits, function(x) {
  audit_h_path("unified_u11", x$ierr == 0L && x$failed_eval == 0L,
               x$H, TRUE, expected_h, n, q)
})

rows <- do.call(rbind, lapply(names(fits), function(nm) {
  x <- fits[[nm]]; flags <- as.integer(sub("f", "", nm))
  data.frame(
    task_id, model, case = case_id, tau, n, seed, cache_flags = flags,
    ierr = x$ierr, failed_eval = x$failed_eval,
    beta_bitwise_equal_25 = identical(x$beta, ref$beta),
    beta_max_abs_diff_25 = max_abs(x$beta, ref$beta),
    h_exact_equal_25 = identical(x$H, ref$H),
    retained_diagnostics_equal_25 =
      identical(x$diagnostics[, 13:18, drop = FALSE], ref$diagnostics[, 13:18, drop = FALSE]),
    h_rowset_match_lean = isTRUE(h_audit[[nm]]$match),
    h_mismatch_eval_count = h_audit[[nm]]$count,
    certificate_hits = sum(x$diagnostics[-1L, "certificate_hits"]),
    residual_rows_checked = sum(x$diagnostics[-1L, "residual_rows"]),
    retained_decomposition_ok = all(
      x$diagnostics[-1L, "first_tableau_rows"] ==
        x$diagnostics[-1L, "first_retained_size"] +
        x$diagnostics[-1L, "first_aggregate_rows"] &
      x$diagnostics[-1L, "basis_forced_rows"] ==
        x$diagnostics[-1L, "first_retained_size"] -
        x$diagnostics[-1L, "effective_threshold_hits"]
    ), stringsAsFactors = FALSE
  )
}))

audit_dir <- file.path(preflight_dir, "cache_mode_audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
path <- file.path(audit_dir, sprintf("config%02d.rds", task_id))
tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
saveRDS(rows, tmp, compress = FALSE)
stopifnot(file.rename(tmp, path))
cat(sprintf("audit task=%d %s case=%d tau=%.1f n=%d PASS\n",
            task_id, model, case_id, tau, n))
