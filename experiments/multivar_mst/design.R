simulate_multivar_llqr <- function(n, case, seed, p = 4L) {
  stopifnot(n > p, p == 4L, case %in% 1:2)
  set.seed(seed)
  if (case == 1L) {
    x <- matrix(stats::rnorm(n * p), n, p)
    signal <- 1 + 1.2 * x[, 1L] - 0.9 * x[, 2L]
    signal <- signal + 0.5 * sin(x[, 1L] + x[, 3L])
    signal <- signal + 0.7 * x[, 3L] * x[, 4L]
    error <- 0.5 * stats::rt(n, df = 4)
  } else {
    sigma <- outer(seq_len(p), seq_len(p), function(i, j) 0.5^abs(i - j))
    x <- matrix(stats::rnorm(n * p), n, p) %*% chol(sigma)
    signal <- 0.5 + 0.8 * x[, 1L] + 0.6 * x[, 2L]^2
    signal <- signal - 0.4 * x[, 3L] + 0.5 * x[, 1L] * x[, 4L]
    scale <- 0.4 + 0.25 * abs(x[, 1L]) + 0.1 * abs(x[, 2L])
    error <- scale * stats::rnorm(n)
  }
  colnames(x) <- paste0("x", seq_len(p))
  list(x = x, y = as.numeric(signal + error), signal = signal, error = error,
       z = x, case = case, n = n, p = p, seed = seed)
}
