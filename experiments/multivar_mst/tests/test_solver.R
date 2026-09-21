#!/usr/bin/env Rscript
source("experiments/multivar_mst/design.R")
source("experiments/multivar_mst/solver.R")

mode <- Sys.getenv("MST_BUILD_MODE", "checked")
kernel <- load_multivar_mst_kernel(sprintf(
  "experiments/multivar_mst/build/multivar_mst_%s.so", mode))

row_set_match <- function(a, b) {
  nrow(a) == nrow(b) && all(vapply(seq_len(nrow(a)), function(j)
    setequal(a[j, ], b[j, ]), logical(1L)))
}

# The complete manuscript signals must be present, including the four terms
# that the historical unparenthesized generator accidentally skipped.
for (case in 1:2) {
  dat <- simulate_multivar_llqr(30L, case, 2026L)
  if (case == 1L) {
    expected <- 1 + 1.2 * dat$x[, 1L] - 0.9 * dat$x[, 2L] +
      0.5 * sin(dat$x[, 1L] + dat$x[, 3L]) + 0.7 * dat$x[, 3L] * dat$x[, 4L]
  } else {
    expected <- 0.5 + 0.8 * dat$x[, 1L] + 0.6 * dat$x[, 2L]^2 -
      0.4 * dat$x[, 3L] + 0.5 * dat$x[, 1L] * dat$x[, 4L]
  }
  stopifnot(max(abs(dat$signal - expected)) < 1e-14)
}

for (case in 1:2) {
  dat <- simulate_multivar_llqr(80L, case, 2030L + case)
  screened <- llqr_seq_screen_mst(dat$x, dat$y, z = dat$z, kernel = kernel,
                                  diagnostics = TRUE)
  full <- llqr_seq_screen_mst(dat$x, dat$y, z = dat$z, h = screened$h,
                              kernel = kernel, diagnostics = TRUE,
                              full_active = TRUE)
  direct <- llqr_direct_multivar(dat$x, dat$y, z = dat$z, h = screened$h)
  os <- weighted_qr_objective(dat$x, dat$y, dat$z, screened$h, 0.5, screened$beta)
  of <- weighted_qr_objective(dat$x, dat$y, dat$z, screened$h, 0.5, full$beta)
  od <- weighted_qr_objective(dat$x, dat$y, dat$z, screened$h, 0.5, direct$beta)
  kkt <- weighted_qr_kkt_residual(dat$x, dat$y, dat$z, screened$h, 0.5,
                                  screened$beta, screened$H_seq)
  stopifnot(
    row_set_match(screened$H_seq, full$H_seq),
    max(abs(os - of) / pmax(1, abs(of))) <= 1e-8,
    max(abs(os - od) / pmax(1, abs(od))) <= 1e-8,
    max(kkt) <= 1e-7,
    max(abs(screened$ll_est - direct$ll_est) /
          pmax(abs(direct$ll_est), 1e-10)) <= 1e-6,
    identical(screened$solver_info$fallback_triggered, FALSE),
    all(screened$parent[screened$eval_order[-1L]] %in%
          screened$eval_order[-length(screened$eval_order)])
  )
}

# Geometries that induce a chain and branching, including a parent that is not
# the node visited immediately before its child.
set.seed(99)
x <- matrix(rnorm(240L), 60L, 4L); y <- rnorm(60L)
chain_z <- cbind(seq(-1, 1, length.out = 9L), matrix(0, 9L, 3L))
chain <- llqr_seq_screen_mst(x, y, z = chain_z, h = rep(2, 4), kernel = kernel)
stopifnot(sum(chain$parent > 0L) == 8L)
branch_z <- rbind(c(0, 0, 0, 0), diag(4), -diag(4))
branch <- llqr_seq_screen_mst(x, y, z = branch_z, h = rep(2, 4), kernel = kernel)
previous <- setNames(c(NA_integer_, head(branch$eval_order, -1L)), branch$eval_order)
stopifnot(any(branch$parent[-branch$root] != previous[as.character(setdiff(seq_len(nrow(branch_z)), branch$root))]))

# Nondefault bandwidth, shuffled and duplicate evaluation points remain valid.
z <- rbind(x[10:1, , drop = FALSE], x[1:2, , drop = FALSE])
edge <- llqr_seq_screen_mst(x, y, z = z, h = rep(1.5, 4), kernel = kernel,
                            diagnostics = TRUE)
stopifnot(nrow(edge$beta) == nrow(z), all(is.finite(edge$beta)))

# Extreme locations can underflow to exact zero; ordinary fitting records these
# zeros, while a window with insufficient effective rank fails explicitly.
stopifnot(sum(gaussian_multivar_weights(x, rep(20, 4), rep(1, 4)) == 0) > 0L)
far_error <- try(llqr_seq_screen_mst(x, y, z = matrix(rep(1e3, 4), 1L),
                                     h = rep(1, 4), kernel = kernel), silent = TRUE)
stopifnot(inherits(far_error, "try-error"))

# A genuinely rank-deficient full design must fail explicitly.
rank_x <- matrix(rep(seq_len(30L), 4L), 30L, 4L)
rank_error <- try(llqr_seq_screen_mst(rank_x, rnorm(30L), h = rep(10, 4),
                                      kernel = kernel), silent = TRUE)
stopifnot(inherits(rank_error, "try-error"))

cat("PASS: DGP, tree ownership, full-path H, objective, KKT and edge cases\n")
