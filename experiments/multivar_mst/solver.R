compute_multivar_bandwidth <- function(x, tau = 0.5, h = NULL) {
  x <- as.matrix(x)
  if (!is.null(h)) {
    h <- as.numeric(h)
    if (length(h) == 1L) h <- rep(h, ncol(x))
  } else {
    factor <- 1.25 * (tau * (1 - tau) / stats::dnorm(stats::qnorm(tau))^2)^0.2
    h <- factor * apply(x, 2L, stats::sd) * nrow(x)^(-1 / (ncol(x) + 4))
  }
  if (length(h) != ncol(x) || any(!is.finite(h)) || any(h <= 0))
    stop("h must contain one positive finite bandwidth per coordinate.")
  h
}

gaussian_multivar_weights <- function(x, z, h) {
  u <- sweep(sweep(x, 2L, z, "-"), 2L, h, "/")
  exp(-0.5 * rowSums(u^2))
}

validate_multivar_inputs <- function(x, y, tau, z) {
  x <- as.matrix(x); y <- as.numeric(y)
  if (is.null(z)) z <- x else z <- as.matrix(z)
  if (nrow(x) != length(y) || ncol(x) != ncol(z) || nrow(z) < 1L)
    stop("x, y and z have incompatible dimensions.")
  if (any(!is.finite(x)) || any(!is.finite(y)) || any(!is.finite(z)))
    stop("x, y and z must be finite.")
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1)
    stop("tau must lie strictly between zero and one.")
  list(x = x, y = y, z = z, tau = tau)
}

llqr_direct_multivar <- function(x, y, tau = 0.5, z = NULL, h = NULL) {
  input <- validate_multivar_inputs(x, y, tau, z)
  x <- input$x; y <- input$y; z <- input$z
  h <- compute_multivar_bandwidth(x, tau, h)
  q <- ncol(x) + 1L; m <- nrow(z)
  beta <- matrix(NA_real_, m, q)
  underflow_count <- integer(m)
  for (j in seq_len(m)) {
    w <- gaussian_multivar_weights(x, z[j, ], h)
    underflow_count[j] <- sum(w == 0)
    fit <- quantreg::rq(y ~ x, tau = tau, weights = w, method = "br")
    coefficient <- as.numeric(stats::coef(fit))
    if (length(coefficient) != q || any(!is.finite(coefficient)))
      stop("quantreg returned invalid coefficients at evaluation point ", j, ".")
    beta[j, ] <- coefficient
  }
  list(
    ll_est = beta[, 1L] + rowSums(z * beta[, -1L, drop = FALSE]),
    gradient = beta[, -1L, drop = FALSE], beta = beta, z = z, h = h,
    solver_info = list(backend = "quantreg_rq_br", method = "br",
                       underflow_count = underflow_count)
  )
}

load_multivar_mst_kernel <- function(path = "experiments/multivar_mst/build/multivar_mst_optimized.so") {
  path <- normalizePath(path, mustWork = TRUE)
  dll <- dyn.load(path)
  list(dll = dll, symbol = getNativeSymbolInfo("ssqr_gaussian_tree", dll))
}

llqr_seq_screen_mst <- function(x, y, tau = 0.5, z = NULL, h = NULL,
                                threshold_factor = 0.1, diagnostics = FALSE,
                                kernel = NULL, full_active = FALSE) {
  input <- validate_multivar_inputs(x, y, tau, z)
  x <- input$x; y <- input$y; z <- input$z
  if (!is.numeric(threshold_factor) || length(threshold_factor) != 1L ||
      !is.finite(threshold_factor) || threshold_factor <= 0)
    stop("threshold_factor must be positive and finite.")
  h <- compute_multivar_bandwidth(x, tau, h)
  n <- nrow(x); d <- ncol(x); q <- d + 1L; m <- nrow(z)
  gamma <- if (full_active) 1e100 else
    threshold_factor * log(log(n))^2 / sqrt(log(n))
  if (is.null(kernel)) kernel <- load_multivar_mst_kernel()
  out <- .Fortran(kernel$symbol,
    a = as.double(cbind(1, x)), y = as.double(y), x = as.double(x),
    n = as.integer(n), q = as.integer(q), d = as.integer(d),
    z = as.double(z), ne = as.integer(m), h = as.double(h),
    tau = as.double(tau), tol = 1e-14, maxit = 1000000L,
    threshold = rep(as.double(gamma), m), min_keep = 1L, cache_flags = 27L,
    beta = double(m * q), hseq = integer(m * q), diagnostics = integer(m * 18L),
    visit_order = integer(m), parent = integer(m), edge_weight = double(m),
    mst_seconds = 0, ierr = 0L, failed_eval = 0L)
  if (out$ierr != 0L)
    stop("seq-screen-MST failed at evaluation point ", out$failed_eval,
         " (ierr=", out$ierr, ").")
  beta <- matrix(out$beta, m, q)
  hseq <- matrix(out$hseq, m, q)
  if (any(!is.finite(beta)) || any(hseq < 1L | hseq > n) ||
      any(apply(hseq, 1L, anyDuplicated) != 0L))
    stop("seq-screen-MST returned invalid coefficients or H indices.")
  path <- matrix(out$diagnostics, m, 18L, dimnames = list(NULL, c(
    "n_active", "first_tableau_rows", "final_tableau_rows", "iterations",
    "repairs", "init_mode", "init_trigger", "full_recovery", "independent",
    "certificate_hits", "residual_rows", "first_full_m_recovery",
    "initial_threshold_hits", "threshold_expansion_steps", "effective_threshold_hits",
    "basis_forced_rows", "first_aggregate_rows", "first_retained_size")))
  root <- out$visit_order[1L]
  path[root, 13:18] <- NA_integer_
  underflow_count <- n - path[, "n_active"]
  info <- list(
    backend = if (full_active) "tree_full_active" else "seq_screen_mst",
    fallback_triggered = FALSE, returned_backend = "weighted_qr_tree_core",
    threshold_initial = gamma, repair_count = sum(path[, "repairs"]),
    full_active_recovery_count = sum(path[, "full_recovery"]),
    independent_init_count = sum(path[, "independent"]),
    underflow_count = underflow_count, mst_seconds = out$mst_seconds)
  result <- list(
    ll_est = beta[, 1L] + rowSums(z * beta[, -1L, drop = FALSE]),
    gradient = beta[, -1L, drop = FALSE], beta = beta, z = z, h = h,
    eval_order = out$visit_order, parent = out$parent, root = root,
    edge_weight = out$edge_weight, H_seq = hseq, solver_info = info)
  if (diagnostics) result$diagnostics <- path
  result
}

weighted_qr_objective <- function(x, y, z, h, tau, beta) {
  vapply(seq_len(nrow(z)), function(j) {
    w <- gaussian_multivar_weights(x, z[j, ], h)
    r <- y - cbind(1, x) %*% beta[j, ]
    sum(w * r * (tau - (r < 0)))
  }, numeric(1L))
}

weighted_qr_kkt_residual <- function(x, y, z, h, tau, beta, hseq) {
  design <- cbind(1, x)
  vapply(seq_len(nrow(z)), function(j) {
    w <- gaussian_multivar_weights(x, z[j, ], h)
    r <- as.numeric(y - design %*% beta[j, ])
    H <- hseq[j, ]; nonzero <- setdiff(seq_along(y), H)
    score <- tau - (r[nonzero] < 0)
    target <- -colSums(design[nonzero, , drop = FALSE] * (w[nonzero] * score))
    mat <- t(design[H, , drop = FALSE] * w[H])
    g <- tryCatch(solve(mat, target), error = function(e) rep(Inf, ncol(design)))
    stationarity <- max(abs(as.numeric(mat %*% g - target)))
    bounds <- max(c(0, g - tau, tau - 1 - g))
    max(stationarity / max(1, sum(w)), bounds)
  }, numeric(1L))
}
