#!/usr/bin/env Rscript
source("R/llqr_functions.R")
source("experiments/hpc_u11_rep500/adapters.R")
dll <- ssqr_load("checked")

set.seed(1)
n <- 40L
x <- seq(0, 1, length.out = n)
y <- 1 + 2 * x + rnorm(n, sd = 0.01)
Sys.setenv(SSQR_CACHE_FLAGS = "31", SSQR_PROVIDER_FLAGS = "1")
out <- ssqr_call(
  dll, cbind(1, x), y, grid = c(0.1, 0.9), coordinate = x,
  kernel = 2L, h = 0.12, tau = 0.5,
  threshold = c(1e-12, 1e-12), min_keep = 1L
)
stopifnot(out$ierr == 0L, out$failed_eval == 0L)
stopifnot(all(is.na(out$diagnostics[1L, 13:18])))
d <- out$diagnostics[2L, ]
stopifnot(
  d[["threshold_expansion_steps"]] > 0L,
  d[["basis_forced_rows"]] > 0L,
  d[["first_aggregate_rows"]] >= 1L,
  d[["first_tableau_rows"]] == d[["first_retained_size"]] + d[["first_aggregate_rows"]],
  d[["basis_forced_rows"]] == d[["first_retained_size"]] - d[["effective_threshold_hits"]]
)
cat("PASS retained-size decomposition with expansion, basis padding, aggregates, and first-point exclusion\n")

dat <- generate_data(200L, case = 2L, seed = 2026L)
h <- llqr_default_bandwidth(dat$x, dat$y, 0.5, case = 2L)
fits <- lapply(c(25L, 27L, 31L), function(flags) {
  Sys.setenv(SSQR_CACHE_FLAGS = as.character(flags))
  ssqr_llqr(dll, dat$x, dat$y, 0.5, 2L, h, sort(dat$x), Mm.factor = 0.1)
})
reference <- fits[[1L]]
for (fit in fits[-1L]) {
  stopifnot(
    fit$ierr == 0L, fit$failed_eval == 0L,
    identical(fit$beta, reference$beta), identical(fit$H, reference$H),
    identical(fit$diagnostics[, 13:18], reference$diagnostics[, 13:18])
  )
}
cat("PASS cache flags 25/27/31 preserve estimates, H, and retained-size diagnostics\n")
