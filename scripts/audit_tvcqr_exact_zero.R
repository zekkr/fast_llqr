#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

tvcqr_audit_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(name, " must be a finite scalar.")
  }
  x
}

tvcqr_audit_metadata_value <- function(metadata, fit, name, default) {
  value <- metadata[[name]] %||% fit[[name]] %||% default
  if (length(value) == 0L) default else value[[1L]]
}

tvcqr_audit_design <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  time_index <- seq_len(n) / n
  A_base <- cbind(1, x)
  cbind(A_base, A_base * time_index)
}

tvcqr_audit_active_matrix <- function(time_index, evaluation_points, h) {
  abs(outer(time_index, evaluation_points, "-") / h) < 1
}

# Compact, blockwise audit used by the replication pipeline. Raw residuals are
# reconstructed from beta_full_est; residual_est is deliberately never read.
audit_tvcqr_exact_zero_compact <- function(x, y, h, fit, metadata = list(),
                                           block_size = 128L) {
  started <- proc.time()[["elapsed"]]
  x <- as.matrix(x)
  y <- as.numeric(y)
  h <- tvcqr_audit_scalar(as.numeric(h), "h")
  block_size <- as.integer(tvcqr_audit_scalar(as.numeric(block_size), "block_size"))
  if (h <= 0) stop("h must be positive.")
  if (block_size <= 0L) stop("block_size must be positive.")
  if (nrow(x) != length(y) || length(y) < 3L || ncol(x) < 1L) {
    stop("x and y must have compatible dimensions and at least three observations.")
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("x and y must be finite.")
  }
  if (!is.list(fit) || is.null(fit$beta_full_est)) {
    stop("fit must contain beta_full_est.")
  }

  A <- tvcqr_audit_design(x)
  beta <- as.matrix(fit$beta_full_est)
  n <- length(y)
  n_eval <- nrow(beta)
  q <- ncol(A)
  if (!identical(ncol(beta), q) || n_eval < 1L || any(!is.finite(beta))) {
    stop("beta_full_est must be finite with 2 * (ncol(x) + 1) columns.")
  }

  evaluation_points <- seq_len(n_eval) / n_eval
  time_index <- seq_len(n) / n
  multipliers <- c(10, 100, 1000)
  residual_scale <- max(abs(y)) + max(abs(A)) * apply(abs(beta), 1L, max) + 1
  eps_matrix <- outer(multipliers * .Machine$double.eps, residual_scale)

  exact_all <- integer(n_eval)
  exact_current <- integer(n_eval)
  exact_next <- rep.int(NA_integer_, n_eval)
  near_next <- matrix(
    NA_integer_, nrow = n_eval, ncol = length(multipliers),
    dimnames = list(NULL, paste0("c", multipliers))
  )

  for (lo in seq.int(1L, n_eval, by = block_size)) {
    hi <- min(n_eval, lo + block_size - 1L)
    jj <- lo:hi
    residual <- -A %*% t(beta[jj, , drop = FALSE])
    residual <- sweep(residual, 1L, y, "+")
    literal <- residual == 0
    current_active <- tvcqr_audit_active_matrix(
      time_index, evaluation_points[jj], h
    )

    exact_all[jj] <- colSums(literal)
    exact_current[jj] <- colSums(literal & current_active)

    has_next <- jj < n_eval
    if (any(has_next)) {
      jj_next <- jj[has_next]
      next_active <- tvcqr_audit_active_matrix(
        time_index, evaluation_points[jj_next + 1L], h
      )
      next_literal <- literal[, has_next, drop = FALSE]
      exact_next[jj_next] <- colSums(next_literal & next_active)
      abs_residual <- abs(residual[, has_next, drop = FALSE])
      for (cc in seq_along(multipliers)) {
        threshold <- matrix(
          eps_matrix[cc, jj_next], nrow = n, ncol = length(jj_next), byrow = TRUE
        )
        near_next[jj_next, cc] <- colSums(next_active & abs_residual <= threshold)
      }
    }
  }

  next_rows <- if (n_eval > 1L) seq_len(n_eval - 1L) else integer()
  next_values <- exact_next[next_rows]
  violation_idx <- next_rows[next_values > q]
  max_next_idx <- if (length(next_rows)) next_rows[which.max(next_values)] else NA_integer_

  witness <- if (length(violation_idx)) {
    do.call(rbind, lapply(violation_idx, function(j) {
      residual_j <- as.numeric(y - A %*% beta[j, ])
      next_active <- abs((time_index - evaluation_points[j + 1L]) / h) < 1
      zero_indices <- which(next_active & residual_j == 0)
      data.frame(
        evaluation_index = j,
        evaluation_point = evaluation_points[j],
        next_evaluation_point = evaluation_points[j + 1L],
        n_zero_exact_all = exact_all[j],
        n_zero_exact_current = exact_current[j],
        n_zero_exact_next_window = exact_next[j],
        zero_observation_indices_next = paste(zero_indices, collapse = "|"),
        n_zero_c10_next_window = near_next[j, 1L],
        n_zero_c100_next_window = near_next[j, 2L],
        n_zero_c1000_next_window = near_next[j, 3L],
        epsilon_c10 = eps_matrix[1L, j],
        epsilon_c100 = eps_matrix[2L, j],
        epsilon_c1000 = eps_matrix[3L, j],
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(
      evaluation_index = integer(), evaluation_point = numeric(),
      next_evaluation_point = numeric(), n_zero_exact_all = integer(),
      n_zero_exact_current = integer(), n_zero_exact_next_window = integer(),
      zero_observation_indices_next = character(),
      n_zero_c10_next_window = integer(), n_zero_c100_next_window = integer(),
      n_zero_c1000_next_window = integer(), epsilon_c10 = numeric(),
      epsilon_c100 = numeric(), epsilon_c1000 = numeric()
    )
  }

  next_distribution <- tabulate(pmin(next_values, q + 1L) + 1L, nbins = q + 2L)
  names(next_distribution) <- c(paste0("n_next_zero_", 0:q), "n_next_zero_gt_q")
  summary <- data.frame(
    case = as.integer(tvcqr_audit_metadata_value(metadata, fit, "case", NA_integer_)),
    n = n,
    tau = as.numeric(tvcqr_audit_metadata_value(metadata, fit, "tau", NA_real_)),
    seed = as.integer(tvcqr_audit_metadata_value(metadata, fit, "seed", NA_integer_)),
    rep_id = as.integer(tvcqr_audit_metadata_value(metadata, fit, "rep_id", NA_integer_)),
    method = as.character(tvcqr_audit_metadata_value(
      metadata, fit, "method", "tvcqr_seq_ppro_fortran_1"
    )),
    q = q,
    n_evaluation_points = n_eval,
    n_next_window_points = length(next_rows),
    max_n_zero_exact_all = max(exact_all),
    max_n_zero_exact_current = max(exact_current),
    max_n_zero_exact_next_window = if (length(next_rows)) max(next_values) else NA_integer_,
    as.list(next_distribution),
    n_literal_violation_next_window = length(violation_idx),
    replication_literal_violation_next_window = length(violation_idx) > 0L,
    first_literal_violation_eval = if (length(violation_idx)) violation_idx[1L] else NA_integer_,
    max_literal_next_eval = max_next_idx,
    max_n_zero_c10_next_window = if (length(next_rows)) max(near_next[next_rows, 1L]) else NA_integer_,
    max_n_zero_c100_next_window = if (length(next_rows)) max(near_next[next_rows, 2L]) else NA_integer_,
    max_n_zero_c1000_next_window = if (length(next_rows)) max(near_next[next_rows, 3L]) else NA_integer_,
    n_nearzero_exceeds_q_next_c10 = if (length(next_rows)) sum(near_next[next_rows, 1L] > q) else 0L,
    n_nearzero_exceeds_q_next_c100 = if (length(next_rows)) sum(near_next[next_rows, 2L] > q) else 0L,
    n_nearzero_exceeds_q_next_c1000 = if (length(next_rows)) sum(near_next[next_rows, 3L] > q) else 0L,
    audit_elapsed_sec = proc.time()[["elapsed"]] - started,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(witness)) {
    witness <- cbind(
      summary[rep(1L, nrow(witness)), c("case", "n", "tau", "seed", "rep_id", "method")],
      witness, row.names = NULL
    )
  }
  list(summary = summary, literal_violation_witnesses = witness)
}

# Pointwise audit retained for deterministic fixtures and witness reproduction.
audit_tvcqr_exact_zero <- function(x, y, h, fit, metadata = list()) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  h <- tvcqr_audit_scalar(as.numeric(h), "h")
  if (h <= 0) stop("h must be positive.")
  A <- tvcqr_audit_design(x)
  beta <- as.matrix(fit$beta_full_est)
  n <- length(y)
  n_eval <- nrow(beta)
  q <- ncol(A)
  if (nrow(x) != n || ncol(beta) != q || any(!is.finite(beta))) {
    stop("Inputs and beta_full_est have incompatible dimensions or non-finite values.")
  }
  time_index <- seq_len(n) / n
  evaluation_points <- seq_len(n_eval) / n_eval
  c_values <- c(10, 100, 1000)
  A_max <- max(abs(A))
  y_max <- max(abs(y))

  do.call(rbind, lapply(seq_len(n_eval), function(j) {
    residual_raw <- as.numeric(y - A %*% beta[j, ])
    current_active <- abs((time_index - evaluation_points[j]) / h) < 1
    has_next <- j < n_eval
    next_active <- if (has_next) {
      abs((time_index - evaluation_points[j + 1L]) / h) < 1
    } else rep(FALSE, n)
    scale_j <- y_max + A_max * max(abs(beta[j, ])) + 1
    eps_values <- c_values * .Machine$double.eps * scale_j
    next_counts <- if (has_next) {
      vapply(eps_values, function(eps_j) {
        sum(next_active & abs(residual_raw) <= eps_j)
      }, integer(1L))
    } else rep(NA_integer_, length(eps_values))
    n_exact_next <- if (has_next) sum(next_active & residual_raw == 0) else NA_integer_

    data.frame(
      case = tvcqr_audit_metadata_value(metadata, fit, "case", NA_integer_),
      n = n,
      tau = tvcqr_audit_metadata_value(metadata, fit, "tau", NA_real_),
      seed = tvcqr_audit_metadata_value(metadata, fit, "seed", NA_integer_),
      method = tvcqr_audit_metadata_value(
        metadata, fit, "method", "tvcqr_seq_ppro_fortran_1"
      ),
      evaluation_index = j,
      evaluation_point = evaluation_points[j],
      next_evaluation_point = if (has_next) evaluation_points[j + 1L] else NA_real_,
      q = q,
      n_zero_exact_all = sum(residual_raw == 0),
      n_zero_exact_current = sum(current_active & residual_raw == 0),
      n_zero_exact_next_window = n_exact_next,
      literal_violation_next_window = if (has_next) n_exact_next > q else NA,
      n_zero_c10_next_window = next_counts[1L],
      n_zero_c100_next_window = next_counts[2L],
      n_zero_c1000_next_window = next_counts[3L],
      nearzero_exceeds_q_next_c10 = if (has_next) next_counts[1L] > q else NA,
      nearzero_exceeds_q_next_c100 = if (has_next) next_counts[2L] > q else NA,
      nearzero_exceeds_q_next_c1000 = if (has_next) next_counts[3L] > q else NA,
      residual_scale = scale_j,
      epsilon_c10 = eps_values[1L],
      epsilon_c100 = eps_values[2L],
      epsilon_c1000 = eps_values[3L],
      stringsAsFactors = FALSE
    )
  }))
}

if (sys.nframe() == 0L) {
  stop("Source this file and call audit_tvcqr_exact_zero_compact() or audit_tvcqr_exact_zero().")
}
