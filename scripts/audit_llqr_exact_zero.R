#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

llqr_audit_kernel_weights <- function(x, s, h, case) {
  u <- (s - x) / h
  if (case == 1L) {
    return(exp(-0.5 * u^2) / sqrt(2 * pi))
  }
  if (case == 2L) {
    return(ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0))
  }
  stop("case must be 1 or 2.")
}

llqr_audit_scalar <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(name, " must be a finite scalar.")
  }
  x
}

llqr_audit_metadata_value <- function(metadata, fit, name, default) {
  value <- metadata[[name]] %||% fit[[name]] %||% default
  if (length(value) == 0L) default else value[[1L]]
}

# Compact, blockwise audit for large simulation grids. The fit vectors may be
# aligned with an unsorted evaluation vector; all three are reordered together
# before the next-window quantities are formed.
audit_llqr_exact_zero_compact <- function(x, y, z, h, fit, case = 1L,
                                          metadata = list(), block_size = 128L) {
  started <- proc.time()[["elapsed"]]
  x <- as.numeric(x)
  y <- as.numeric(y)
  z <- as.numeric(z)
  h <- llqr_audit_scalar(as.numeric(h), "h")
  case <- as.integer(llqr_audit_scalar(as.numeric(case), "case"))
  block_size <- as.integer(llqr_audit_scalar(as.numeric(block_size), "block_size"))
  if (h <= 0) stop("h must be positive.")
  if (block_size <= 0L) stop("block_size must be positive.")
  if (!(case %in% c(1L, 2L))) stop("case must be 1 or 2.")
  if (length(x) != length(y) || length(x) < 3L) {
    stop("x and y must have the same length, with at least three observations.")
  }
  if (length(z) < 1L || any(!is.finite(x)) || any(!is.finite(y)) ||
      any(!is.finite(z))) {
    stop("x, y, and z must be finite, and z must be nonempty.")
  }
  if (!is.list(fit) || is.null(fit$ll_est) || is.null(fit$d_ll_est)) {
    stop("fit must contain ll_est and d_ll_est.")
  }

  ll_est <- as.numeric(fit$ll_est)
  d_ll_est <- as.numeric(fit$d_ll_est)
  if (length(ll_est) != length(z) || length(d_ll_est) != length(z) ||
      any(!is.finite(ll_est)) || any(!is.finite(d_ll_est))) {
    stop("ll_est and d_ll_est must be finite with one entry per evaluation point.")
  }

  evaluation_order <- order(z, method = "radix")
  z <- z[evaluation_order]
  ll_est <- ll_est[evaluation_order]
  d_ll_est <- d_ll_est[evaluation_order]

  n <- length(y)
  n_eval <- length(z)
  q <- 2L
  alpha0 <- ll_est - z * d_ll_est
  theta_max <- pmax(abs(alpha0), abs(d_ll_est))
  residual_scale <- max(abs(y)) + max(1, max(abs(x))) * theta_max + 1
  multipliers <- c(10, 100, 1000)
  eps_matrix <- outer(multipliers * .Machine$double.eps, residual_scale)

  exact_all <- integer(n_eval)
  exact_current <- integer(n_eval)
  exact_next <- rep.int(NA_integer_, n_eval)
  near_next <- matrix(NA_integer_, nrow = n_eval, ncol = 3L,
                      dimnames = list(NULL, paste0("c", multipliers)))

  active_matrix <- function(eval_points) {
    u <- outer(x, eval_points, "-") / h
    if (case == 1L) {
      return(exp(-0.5 * u^2) > 0)
    }
    abs(u) < 1
  }

  for (lo in seq.int(1L, n_eval, by = block_size)) {
    hi <- min(n_eval, lo + block_size - 1L)
    jj <- lo:hi
    residual <- -tcrossprod(x, d_ll_est[jj])
    residual <- sweep(residual, 2L, alpha0[jj], "-")
    residual <- sweep(residual, 1L, y, "+")
    literal <- residual == 0
    current_active <- active_matrix(z[jj])

    exact_all[jj] <- colSums(literal)
    exact_current[jj] <- colSums(literal & current_active)

    has_next <- jj < n_eval
    if (any(has_next)) {
      jj_next <- jj[has_next]
      next_active <- active_matrix(z[jj_next + 1L])
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

  next_rows <- seq_len(max(0L, n_eval - 1L))
  next_values <- exact_next[next_rows]
  violation_idx <- next_rows[next_values > q]
  max_next_idx <- if (length(next_rows)) next_rows[which.max(next_values)] else NA_integer_

  witness <- if (length(violation_idx)) {
    data.frame(
      evaluation_index = violation_idx,
      evaluation_point = z[violation_idx],
      next_evaluation_point = z[violation_idx + 1L],
      n_zero_exact_all = exact_all[violation_idx],
      n_zero_exact_current = exact_current[violation_idx],
      n_zero_exact_next_window = exact_next[violation_idx],
      n_zero_c10_next_window = near_next[violation_idx, 1L],
      n_zero_c100_next_window = near_next[violation_idx, 2L],
      n_zero_c1000_next_window = near_next[violation_idx, 3L],
      epsilon_c10 = eps_matrix[1L, violation_idx],
      epsilon_c100 = eps_matrix[2L, violation_idx],
      epsilon_c1000 = eps_matrix[3L, violation_idx],
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      evaluation_index = integer(), evaluation_point = numeric(),
      next_evaluation_point = numeric(), n_zero_exact_all = integer(),
      n_zero_exact_current = integer(), n_zero_exact_next_window = integer(),
      n_zero_c10_next_window = integer(), n_zero_c100_next_window = integer(),
      n_zero_c1000_next_window = integer(), epsilon_c10 = numeric(),
      epsilon_c100 = numeric(), epsilon_c1000 = numeric()
    )
  }

  summary <- data.frame(
    case = case,
    n = n,
    tau = as.numeric(llqr_audit_metadata_value(metadata, fit, "tau", NA_real_)),
    seed = as.integer(llqr_audit_metadata_value(metadata, fit, "seed", NA_integer_)),
    rep_id = as.integer(llqr_audit_metadata_value(metadata, fit, "rep_id", NA_integer_)),
    method = as.character(llqr_audit_metadata_value(metadata, fit, "method", "llqr")),
    n_evaluation_points = n_eval,
    n_next_window_points = length(next_rows),
    max_n_zero_exact_all = max(exact_all),
    max_n_zero_exact_current = max(exact_current),
    max_n_zero_exact_next_window = if (length(next_rows)) max(next_values) else NA_integer_,
    n_next_zero_0 = sum(next_values == 0L),
    n_next_zero_1 = sum(next_values == 1L),
    n_next_zero_2 = sum(next_values == 2L),
    n_next_zero_gt_q = sum(next_values > q),
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
    stringsAsFactors = FALSE
  )

  if (nrow(witness)) {
    witness <- cbind(summary[rep(1L, nrow(witness)), c("case", "n", "tau", "seed", "rep_id", "method")],
                     witness, row.names = NULL)
  }
  list(summary = summary, literal_violation_witnesses = witness)
}

audit_llqr_exact_zero <- function(x, y, z, h, fit, case = 1L, metadata = list()) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  z <- as.numeric(z)
  h <- llqr_audit_scalar(as.numeric(h), "h")
  case <- as.integer(llqr_audit_scalar(as.numeric(case), "case"))
  if (h <= 0) stop("h must be positive.")
  if (length(x) != length(y) || length(x) < 3L) {
    stop("x and y must have the same length, with at least three observations.")
  }
  if (length(z) < 1L || is.unsorted(z, strictly = FALSE)) {
    stop("z must be nonempty and ordered in the same evaluation order as fit.")
  }
  if (any(!is.finite(x)) || any(!is.finite(y)) || any(!is.finite(z))) {
    stop("x, y, and z must be finite.")
  }
  if (!is.list(fit) || is.null(fit$ll_est) || is.null(fit$d_ll_est)) {
    stop("fit must contain ll_est and d_ll_est.")
  }

  ll_est <- as.numeric(fit$ll_est)
  d_ll_est <- as.numeric(fit$d_ll_est)
  if (length(ll_est) != length(z) || length(d_ll_est) != length(z)) {
    stop("ll_est and d_ll_est must have one entry per evaluation point.")
  }
  if (any(!is.finite(ll_est)) || any(!is.finite(d_ll_est))) {
    stop("ll_est and d_ll_est must be finite.")
  }

  n <- length(y)
  n_eval <- length(z)
  q <- 2L
  A <- cbind(1, x)
  A_max <- max(abs(A))
  y_max <- max(abs(y))
  machine_eps <- .Machine$double.eps
  c_values <- c(10, 100, 1000)
  H_seq <- fit$H_seq
  if (!is.null(H_seq)) {
    H_seq <- as.matrix(H_seq)
    if (!identical(dim(H_seq), c(n_eval, q))) {
      stop("When supplied, H_seq must have dimension length(z) by 2.")
    }
  }

  point_rows <- vector("list", n_eval)
  for (j in seq_len(n_eval)) {
    alpha1 <- d_ll_est[j]
    alpha0 <- ll_est[j] - z[j] * alpha1
    theta <- c(alpha0, alpha1)
    residual_raw <- as.numeric(y - A %*% theta)
    scale_j <- y_max + A_max * max(abs(theta)) + 1
    eps_values <- c_values * machine_eps * scale_j

    current_active <- llqr_audit_kernel_weights(x, z[j], h, case) > 0
    has_next <- j < n_eval
    next_active <- if (has_next) {
      llqr_audit_kernel_weights(x, z[j + 1L], h, case) > 0
    } else {
      rep(FALSE, n)
    }

    exact_all <- residual_raw == 0
    exact_current <- exact_all & current_active
    exact_next <- exact_all & next_active
    current_abs <- sort(abs(residual_raw[current_active]))
    qplus1_gap <- if (length(current_abs) >= q + 1L) current_abs[q + 1L] else NA_real_

    H_idx <- if (is.null(H_seq)) integer() else as.integer(H_seq[j, ])
    H_valid <- length(H_idx) == q && !anyNA(H_idx) &&
      all(H_idx >= 1L & H_idx <= n) && !anyDuplicated(H_idx)
    max_abs_H <- if (H_valid) max(abs(residual_raw[H_idx])) else NA_real_
    exact_outside_H <- exact_all
    if (H_valid) exact_outside_H[H_idx] <- FALSE

    current_counts <- vapply(
      eps_values,
      function(eps_j) sum(current_active & abs(residual_raw) <= eps_j),
      integer(1)
    )
    next_counts <- if (has_next) {
      vapply(
        eps_values,
        function(eps_j) sum(next_active & abs(residual_raw) <= eps_j),
        integer(1)
      )
    } else {
      rep(NA_integer_, length(eps_values))
    }

    point_rows[[j]] <- data.frame(
      case = case,
      n = n,
      tau = metadata$tau %||% fit$tau %||% NA_real_,
      seed = metadata$seed %||% fit$seed %||% NA_integer_,
      method = metadata$method %||% fit$method %||% "llqr",
      evaluation_index = j,
      evaluation_point = z[j],
      next_evaluation_point = if (has_next) z[j + 1L] else NA_real_,
      q = q,
      n_active_current = sum(current_active),
      n_active_next_window = if (has_next) sum(next_active) else NA_integer_,
      active_design_rank = if (any(current_active)) {
        qr(A[current_active, , drop = FALSE])$rank
      } else 0L,
      n_zero_exact_all = sum(exact_all),
      n_zero_exact_current = sum(exact_current),
      n_zero_exact_next_window = if (has_next) sum(exact_next) else NA_integer_,
      literal_violation_next_window = if (has_next) sum(exact_next) > q else NA,
      n_zero_c10_current = current_counts[1],
      n_zero_c100_current = current_counts[2],
      n_zero_c1000_current = current_counts[3],
      n_zero_c10_next_window = next_counts[1],
      n_zero_c100_next_window = next_counts[2],
      n_zero_c1000_next_window = next_counts[3],
      nearzero_exceeds_q_next_c10 = if (has_next) next_counts[1] > q else NA,
      nearzero_exceeds_q_next_c100 = if (has_next) next_counts[2] > q else NA,
      nearzero_exceeds_q_next_c1000 = if (has_next) next_counts[3] > q else NA,
      n_zero_exact_outside_H = sum(exact_outside_H),
      max_abs_raw_residual_H = max_abs_H,
      qplus1_abs_residual_current = qplus1_gap,
      residual_scale = scale_j,
      epsilon_c10 = eps_values[1],
      epsilon_c100 = eps_values[2],
      epsilon_c1000 = eps_values[3],
      returned_backend = fit$returned_backend %||% NA_character_,
      init_mode = if (!is.null(fit$init_mode)) fit$init_mode[j] else NA_character_,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, point_rows)
}

summarize_llqr_exact_zero <- function(pointwise) {
  required <- c(
    "case", "n", "tau", "seed", "method", "n_zero_exact_next_window",
    "literal_violation_next_window", "n_zero_c10_next_window",
    "n_zero_c100_next_window", "n_zero_c1000_next_window"
  )
  if (!all(required %in% names(pointwise))) {
    stop("pointwise audit data are missing required columns.")
  }

  group_names <- c("case", "n", "tau", "seed", "method")
  group_key <- do.call(
    paste,
    c(
      lapply(pointwise[group_names], function(x) ifelse(is.na(x), "<NA>", as.character(x))),
      sep = "\r"
    )
  )
  pieces <- split(pointwise, group_key)
  rows <- lapply(pieces, function(d) {
    next_rows <- !is.na(d$n_zero_exact_next_window)
    data.frame(
      case = d$case[1],
      n = d$n[1],
      tau = d$tau[1],
      seed = d$seed[1],
      method = d$method[1],
      n_evaluation_points = nrow(d),
      n_next_window_points = sum(next_rows),
      max_n_zero_exact_next_window = if (any(next_rows)) {
        max(d$n_zero_exact_next_window[next_rows])
      } else NA_integer_,
      n_literal_violation_next_window = sum(d$literal_violation_next_window, na.rm = TRUE),
      replication_literal_violation_next_window = any(
        d$literal_violation_next_window,
        na.rm = TRUE
      ),
      max_n_zero_c10_next_window = if (any(next_rows)) {
        max(d$n_zero_c10_next_window[next_rows])
      } else NA_integer_,
      max_n_zero_c100_next_window = if (any(next_rows)) {
        max(d$n_zero_c100_next_window[next_rows])
      } else NA_integer_,
      max_n_zero_c1000_next_window = if (any(next_rows)) {
        max(d$n_zero_c1000_next_window[next_rows])
      } else NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  rownames_out <- do.call(rbind, rows)
  rownames(rownames_out) <- NULL
  rownames_out
}

write_llqr_exact_zero_report <- function(summary, path) {
  table_text <- paste(capture.output(print(summary, row.names = FALSE)), collapse = "\n")
  lines <- c(
    "# LLQR exact-zero audit",
    "",
    "The primary Theorem 2 diagnostic is the next-window literal count",
    "",
    "$$",
    "N_{j\\to j+1,0}=\\sum_{i:w_i(s_{j+1})>0}\\mathbf{1}\\{r_{ij}^{\\rm raw}=0\\}.",
    "$$",
    "",
    "A strict diagnostic flag is recorded when $N_{j\\to j+1,0}>q$. Counts based on",
    "$c\\,\\epsilon_{\\mathrm{machine}}$ for $c\\in\\{10,100,1000\\}$ are numerical",
    "sensitivity diagnostics and are not labelled as theorem violations.",
    "",
    "```text",
    table_text,
    "```",
    ""
  )
  writeLines(lines, path, useBytes = TRUE)
}

llqr_audit_parse_args <- function(args) {
  values <- list()
  for (arg in args) {
    if (grepl("^--input=", arg)) values$input <- sub("^--input=", "", arg)
    if (grepl("^--output-prefix=", arg)) {
      values$output_prefix <- sub("^--output-prefix=", "", arg)
    }
  }
  if (is.null(values$input) || is.null(values$output_prefix)) {
    stop("Usage: audit_llqr_exact_zero.R --input=INPUT.rds --output-prefix=PATH")
  }
  values
}

run_llqr_exact_zero_audit_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- llqr_audit_parse_args(args)
  input <- readRDS(opts$input)
  records <- input$records %||% if (
    is.list(input) && all(c("x", "y", "z", "h", "fit") %in% names(input))
  ) list(input) else input
  if (!is.list(records) || length(records) < 1L) {
    stop("Input RDS must contain one audit record or a list named records.")
  }

  pointwise <- do.call(rbind, lapply(records, function(record) {
    audit_llqr_exact_zero(
      x = record$x,
      y = record$y,
      z = record$z,
      h = record$h,
      fit = record$fit,
      case = record$case %||% 1L,
      metadata = record$metadata %||% list()
    )
  }))
  rownames(pointwise) <- NULL
  summary <- summarize_llqr_exact_zero(pointwise)

  output_dir <- dirname(opts$output_prefix)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write.csv(pointwise, paste0(opts$output_prefix, "_pointwise.csv"), row.names = FALSE)
  write.csv(summary, paste0(opts$output_prefix, "_summary.csv"), row.names = FALSE)
  write_llqr_exact_zero_report(summary, paste0(opts$output_prefix, "_report.md"))
  invisible(list(pointwise = pointwise, summary = summary))
}

if (sys.nframe() == 0L) {
  run_llqr_exact_zero_audit_cli()
}
