# ============================================================================ #
# Multivariate LLQR utilities and simplex-based solvers
# ============================================================================ #

#' Compute a stable LLQR bandwidth for multivariate evaluation.
#'
#' For one-dimensional inputs, this reuses the existing dpill-based rule.
#' For multivariate inputs, it returns a coordinate-wise bandwidth vector based
#' on marginal scales and the usual n^(-1/(d+4)) rate.
#'
#' @param x Covariate matrix.
#' @param y Response vector.
#' @param tau Quantile level in (0, 1).
#' @param h Optional user-supplied bandwidth. Returned unchanged if not NULL.
#' @return A positive scalar bandwidth for d = 1, or a positive vector of
#'   length ncol(x) for d > 1.
compute_llqr_multivar_bandwidth <- function(x, y, tau, h = NULL) {
  x <- as.matrix(x)
  y <- as.numeric(as.matrix(y))

  if (!is.null(h)) {
    return(as.numeric(h))
  }

  m <- nrow(x)
  nvar <- ncol(x)
  tau_factor <- 1.25 * (tau * (1 - tau) / (dnorm(qnorm(tau)))^2)^0.2

  if (nvar == 1L) {
    red_dim <- floor(0.2 * m)
    ord_y <- order(y)
    lo <- red_dim + 1L
    hi <- m - red_dim
    index_y <- if (lo > hi) ord_y else ord_y[lo:hi]

    h0 <- tryCatch(
      KernSmooth::dpill(x[index_y, , drop = FALSE], y[index_y]),
      error = function(e) NaN,
      warning = function(w) suppressWarnings(NaN)
    )

    h1 <- tau_factor * h0
    if (any(is.nan(h1)) || any(!is.finite(h1)) || any(h1 <= 0)) {
      h1 <- 1.25 * max(m^(-1 / (nvar + 4)), min(2, stats::sd(y)) * m^(-1 / (nvar + 4)))
    }
    return(as.numeric(h1))
  }

  x_scale <- apply(x, 2L, stats::sd)
  x_scale[!is.finite(x_scale) | x_scale <= 0] <- 1
  rate <- m^(-1 / (nvar + 4))
  as.numeric(tau_factor * x_scale * rate)
}

.llqr_multivar_validate_inputs <- function(x, y, tau, z = NULL) {
  x <- as.matrix(x)
  y <- as.numeric(as.matrix(y))

  if (is.null(z)) {
    z <- x
  } else {
    z <- as.matrix(z)
  }

  if (length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("tau must be a single finite number strictly between 0 and 1.")
  }
  if (nrow(x) != length(y)) {
    stop("length(y) must equal nrow(x).")
  }
  if (ncol(x) != ncol(z)) {
    stop("x and z must have the same number of columns.")
  }
  if (anyNA(x) || anyNA(y) || anyNA(z)) {
    stop("Missing values are not supported.")
  }
  if (nrow(x) <= ncol(x)) {
    stop("nrow(x) must be strictly larger than ncol(x).")
  }

  list(x = x, y = y, z = z)
}

.llqr_multivar_scale_vector <- function(h, nvar) {
  h <- as.numeric(h)
  if (length(h) == 1L) {
    scale_vec <- rep(h, nvar)
  } else if (length(h) == nvar) {
    scale_vec <- h
  } else {
    stop("For multivariate LLQR, h must be a scalar or a vector of length ncol(x).")
  }

  if (any(!is.finite(scale_vec)) || any(scale_vec <= 0)) {
    stop("Bandwidth values must be finite and strictly positive.")
  }

  scale_vec
}

.llqr_multivar_compute_weights <- function(x, z_eval, h) {
  nvar <- ncol(x)
  scale_vec <- .llqr_multivar_scale_vector(h, nvar)
  centered <- sweep(x, 2L, z_eval, FUN = "-")
  centered <- sweep(centered, 2L, scale_vec, FUN = "/")

  if (nvar == 1L) {
    return(as.numeric(stats::dnorm(centered[, 1L])))
  }

  # A common positive multiplicative constant does not affect the optimizer.
  as.numeric(exp(-0.5 * rowSums(centered^2)))
}

.llqr_multivar_basis_weight_term <- function(r1, w, nvar) {
  idx <- r1 - 1L - nvar
  out <- numeric(length(r1))
  pos <- idx > 0L
  if (any(pos)) {
    out[pos] <- w[idx[pos]]
  }
  out
}

.llqr_multivar_extract_estimate <- function(b, IB, nvar) {
  tmp <- IB %in% 1:(nvar + 1L)
  estimate <- b[tmp][order(IB[tmp])]

  if (length(estimate) != (1L + nvar)) {
    estimate <- numeric(1L + nvar)
    for (ii in 1:(nvar + 1L)) {
      hit <- which(IB == ii)
      if (length(hit) > 0L) {
        estimate[ii] <- b[hit[1L]]
      }
    }
  }

  estimate
}

.llqr_multivar_update_last_row <- function(gammax, cc, IB, r1, w, tau, m, nvar) {
  tau * .llqr_multivar_basis_weight_term(r1, w, nvar) -
    as.numeric(matrix(cc[IB], nrow = 1L) %*% gammax[1:m, , drop = FALSE])
}

.llqr_multivar_init_state <- function(x, y) {
  m <- nrow(x)
  nvar <- ncol(x)

  gammax <- cbind(matrix(1, nrow = m, ncol = 1L), x)
  b <- c(y, 0)
  neg_y <- y < 0
  gammax[neg_y, ] <- -gammax[neg_y, , drop = FALSE]
  b[b < 0] <- -b[b < 0]

  list(
    gammax = rbind(gammax, rep(0, nvar + 1L)),
    b = b,
    IB = ifelse(y >= 0, (1:m) + 1L + nvar, (1:m) + 1L + nvar + m),
    r1 = 1:(nvar + 1L),
    r2 = numeric(nvar + 1L),
    freevarrow = c(rep(FALSE, m), TRUE)
  )
}

.llqr_multivar_solve_from_state <- function(state, w, tau, tol, maxit, bland, m, nvar) {
  gammax <- state$gammax
  b <- state$b
  IB <- state$IB
  r1 <- state$r1
  r2 <- state$r2
  freevarrow <- state$freevarrow

  rr <- matrix(0, nrow = 2L, ncol = 1L + nvar)
  cc <- c(rep(0, 1L + nvar), tau * w, (1 - tau) * w)
  gammax[m + 1L, ] <- .llqr_multivar_update_last_row(
    gammax = gammax, cc = cc, IB = IB, r1 = r1, w = w, tau = tau, m = m, nvar = nvar
  )

  j <- 0L
  repeat {
    if (j >= maxit) {
      warning("Simplex iteration did not converge within maxit.")
      break
    }

    rr[1L, ] <- gammax[m + 1L, ]
    rr[2L, ] <- (.llqr_multivar_basis_weight_term(r1, w, nvar) - rr[1L, ]) * (r2 != 0)
    rr[1L, r2 == 0] <- -abs(rr[1L, r2 == 0])

    rrl <- min(rr)
    if (rrl >= -tol) {
      break
    }

    if (bland) {
      if (any(rr[1L, ] < -tol)) {
        tsep <- which(rr[1L, ] < -tol)
        tmp <- r1[tsep]
        t <- min(tmp)
        t_rr <- tsep[which(tmp == t)[1L]]
        tsep <- 1L
      } else {
        tsep <- which(rr[2L, ] < -tol)
        tmp <- r2[tsep]
        t <- min(tmp)
        t_rr <- tsep[which(tmp == t)[1L]]
        tsep <- 2L
      }
    } else {
      hit <- which(rr == rrl, arr.ind = TRUE)[1L, ]
      tsep <- hit[1L]
      t_rr <- hit[2L]
      t <- if (tsep == 1L) r1[t_rr] else r2[t_rr]
    }

    if (r2[t_rr] != 0L) {
      yy <- if (tsep == 1L) gammax[, t_rr] else -gammax[, t_rr]
      k_ratio <- b / yy
      if (bland) {
        k <- which((min(k_ratio[yy > 0 & !freevarrow]) == k_ratio) & (yy > 0))
        if (length(k) != 1L) {
          tmp <- IB[k]
          k <- k[which(tmp == min(tmp))[1L]]
        }
      } else {
        k <- which((min(k_ratio[yy > 0 & !freevarrow]) == k_ratio) & (yy > 0))[1L]
      }

      if (tsep != 1L) {
        idxw <- r1[t_rr] - 1L - nvar
        if (idxw > 0L) {
          yy[m + 1L] <- yy[m + 1L] + w[idxw]
        }
      }
    } else {
      yy <- gammax[, t_rr]
      if (yy[m + 1L] < 0) {
        k_ratio <- b / yy
        if (bland) {
          k <- which((min(k_ratio[yy > 0 & !freevarrow]) == k_ratio) & (yy > 0))
          if (length(k) != 1L) {
            tmp <- IB[k]
            k <- k[which(tmp == min(tmp))[1L]]
          }
        } else {
          k <- which((min(k_ratio[yy > 0 & !freevarrow]) == k_ratio) & (yy > 0))[1L]
        }
      } else {
        k_ratio <- -b / yy
        if (bland) {
          k <- which((min(k_ratio[yy < 0 & !freevarrow]) == k_ratio) & (yy < 0))
          if (length(k) != 1L) {
            tmp <- IB[k]
            k <- k[which(tmp == min(tmp))[1L]]
          }
        } else {
          k <- which((min(k_ratio[yy < 0 & !freevarrow]) == k_ratio) & (yy < 0))[1L]
        }
      }
      freevarrow[k] <- TRUE
    }

    ee <- yy / yy[k]
    ee[k] <- 1 - 1 / yy[k]

    if (IB[k] <= (nvar + m + 1L)) {
      gammax[, t_rr] <- 0
      gammax[k, t_rr] <- 1
      r1[t_rr] <- IB[k]
      r2[t_rr] <- IB[k] + m
    } else {
      gammax[, t_rr] <- 0
      gammax[k, t_rr] <- -1
      idxw <- IB[k] - m - nvar - 1L
      gammax[m + 1L, t_rr] <- if (idxw > 0L) w[idxw] else 0
      r1[t_rr] <- IB[k] - m
      r2[t_rr] <- IB[k]
    }

    gammax <- gammax - tcrossprod(ee, gammax[k, ])
    b <- b - ee * b[k]
    IB[k] <- t
    j <- j + 1L
  }

  list(
    state = list(
      gammax = gammax,
      b = b,
      IB = IB,
      r1 = r1,
      r2 = r2,
      freevarrow = freevarrow
    ),
    it_num = j,
    H = r1 - 1L - nvar,
    estimate = .llqr_multivar_extract_estimate(b, IB, nvar)
  )
}

#' Build an evaluation traversal for multivariate LLQR.
#'
#' @param z Evaluation-point matrix.
#' @param order_method One of "mst", "input", "random", or "sorted".
#' @param root_method For MST only, one of "center" or "first".
#' @param distance_scale Optional positive scalar or vector used to scale
#'   coordinates before building the MST.
#' @param seed Optional seed used when order_method = "random".
#' @return A list with traversal order, parent vector, root index, and edge
#'   weights.
build_llqr_multivar_traversal <- function(z,
                                          order_method = c("mst", "input", "random", "sorted"),
                                          root_method = c("center", "first"),
                                          distance_scale = NULL,
                                          seed = NULL) {
  z <- as.matrix(z)
  rounds <- nrow(z)
  nvar <- ncol(z)
  order_method <- match.arg(order_method)
  root_method <- match.arg(root_method)

  if (rounds == 0L) {
    stop("z must contain at least one evaluation point.")
  }

  if (!is.null(distance_scale)) {
    scale_vec <- .llqr_multivar_scale_vector(distance_scale, nvar)
    z_dist <- sweep(z, 2L, scale_vec, FUN = "/")
  } else {
    z_dist <- z
  }

  if (order_method == "sorted") {
    if (nvar != 1L) {
      stop("order_method = 'sorted' is only defined for one-dimensional z.")
    }
    ord <- order(z[, 1L])
    parent <- integer(rounds)
    edge_weight <- numeric(rounds)
    if (rounds >= 2L) {
      parent[ord[-1L]] <- ord[-rounds]
      edge_weight[ord[-1L]] <- abs(z_dist[ord[-1L], 1L] - z_dist[ord[-rounds], 1L])
    }
    return(list(order = ord, parent = parent, root = ord[1L], edge_weight = edge_weight))
  }

  if (order_method == "input") {
    ord <- seq_len(rounds)
    parent <- integer(rounds)
    edge_weight <- numeric(rounds)
    if (rounds >= 2L) {
      parent[ord[-1L]] <- ord[-rounds]
      diffs <- z_dist[ord[-1L], , drop = FALSE] - z_dist[ord[-rounds], , drop = FALSE]
      edge_weight[ord[-1L]] <- sqrt(rowSums(diffs^2))
    }
    return(list(order = ord, parent = parent, root = ord[1L], edge_weight = edge_weight))
  }

  if (order_method == "random") {
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      } else {
        NULL
      }
      on.exit({
        if (is.null(old_seed)) {
          rm(".Random.seed", envir = .GlobalEnv)
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(seed)
    }

    ord <- sample.int(rounds)
    parent <- integer(rounds)
    edge_weight <- numeric(rounds)
    if (rounds >= 2L) {
      parent[ord[-1L]] <- ord[-rounds]
      diffs <- z_dist[ord[-1L], , drop = FALSE] - z_dist[ord[-rounds], , drop = FALSE]
      edge_weight[ord[-1L]] <- sqrt(rowSums(diffs^2))
    }
    return(list(order = ord, parent = parent, root = ord[1L], edge_weight = edge_weight))
  }

  if (root_method == "first") {
    root <- 1L
  } else {
    z_mean <- colMeans(z_dist)
    root <- which.min(rowSums((z_dist - matrix(z_mean, nrow = rounds, ncol = nvar, byrow = TRUE))^2))
  }

  in_tree <- rep(FALSE, rounds)
  parent <- integer(rounds)
  edge_weight <- rep(Inf, rounds)
  nearest_dist2 <- rep(Inf, rounds)
  ord <- integer(rounds)

  in_tree[root] <- TRUE
  ord[1L] <- root
  parent[root] <- 0L
  edge_weight[root] <- 0
  nearest_dist2[root] <- 0

  if (rounds >= 2L) {
    diff0 <- sweep(z_dist, 2L, z_dist[root, ], FUN = "-")
    nearest_dist2 <- rowSums(diff0^2)
    nearest_dist2[root] <- 0
    parent[] <- root
    parent[root] <- 0L

    for (tt in 2L:rounds) {
      cand <- which(!in_tree)
      v <- cand[which.min(nearest_dist2[cand])]

      in_tree[v] <- TRUE
      ord[tt] <- v
      edge_weight[v] <- sqrt(nearest_dist2[v])

      remaining <- which(!in_tree)
      if (length(remaining) > 0L) {
        diffv <- sweep(z_dist[remaining, , drop = FALSE], 2L, z_dist[v, ], FUN = "-")
        d2 <- rowSums(diffv^2)
        improve <- d2 < nearest_dist2[remaining]
        if (any(improve)) {
          idx <- remaining[improve]
          nearest_dist2[idx] <- d2[improve]
          parent[idx] <- v
        }
      }
    }
  }

  list(order = ord, parent = parent, root = root, edge_weight = edge_weight)
}

#' Cold-start multivariate LLQR solved by the custom simplex routine.
#'
#' @param x Covariate matrix.
#' @param y Response vector.
#' @param tau Quantile level.
#' @param z Evaluation-point matrix. Defaults to x.
#' @param h Optional bandwidth.
#' @param tol Convergence tolerance.
#' @param maxit Maximum simplex iterations per evaluation point.
#' @param bland Whether to use Bland's rule.
#' @param track_order If TRUE, return results in the input order of z.
#' @return A list containing estimates, iteration counts, bandwidth, and
#'   traversal metadata.
llqr_tau_multivar <- function(x, y, tau, z = NULL, h = NULL,
                              tol = 1e-14, maxit = 10000,
                              bland = FALSE, track_order = FALSE) {
  validated <- .llqr_multivar_validate_inputs(x = x, y = y, tau = tau, z = z)
  x <- validated$x
  y <- validated$y
  z <- validated$z

  m <- nrow(x)
  nvar <- ncol(x)
  rounds <- nrow(z)

  h <- compute_llqr_multivar_bandwidth(x = x, y = y, tau = tau, h = h)
  order_method <- if (nvar == 1L) "sorted" else "input"
  traversal <- build_llqr_multivar_traversal(
    z = z,
    order_method = order_method,
    distance_scale = if (nvar > 1L) h else NULL
  )

  eval_order <- traversal$order
  ll_est_work <- numeric(rounds)
  it_num_work <- integer(rounds)
  H_work <- matrix(0, nrow = rounds, ncol = 1L + nvar)
  d_ll_est_work <- if (nvar == 1L) numeric(rounds) else matrix(0, nrow = rounds, ncol = nvar)
  if (nvar > 1L) {
    colnames(d_ll_est_work) <- colnames(x)
  }

  init_state <- .llqr_multivar_init_state(x = x, y = y)

  for (tt in seq_len(rounds)) {
    idx <- eval_order[tt]
    z_eval <- z[idx, , drop = FALSE]
    w <- .llqr_multivar_compute_weights(x = x, z_eval = z_eval[1, ], h = h)
    res <- .llqr_multivar_solve_from_state(
      state = init_state, w = w, tau = tau, tol = tol, maxit = maxit,
      bland = bland, m = m, nvar = nvar
    )
    estimate <- res$estimate

    ll_est_work[tt] <- sum(c(1, z_eval[1, ]) * estimate)
    if (nvar == 1L) {
      d_ll_est_work[tt] <- estimate[2L]
    } else {
      d_ll_est_work[tt, ] <- estimate[2:(nvar + 1L)]
    }
    it_num_work[tt] <- res$it_num
    H_work[tt, ] <- res$H
  }

  if (track_order) {
    inv_perm <- order(eval_order)
    ll_est <- ll_est_work[inv_perm]
    it_num <- it_num_work[inv_perm]
    H <- H_work[inv_perm, , drop = FALSE]
    d_ll_est <- if (nvar == 1L) d_ll_est_work[inv_perm] else d_ll_est_work[inv_perm, , drop = FALSE]
  } else {
    ll_est <- ll_est_work
    it_num <- it_num_work
    H <- H_work
    d_ll_est <- d_ll_est_work
  }

  list(
    ll_est = ll_est,
    d_ll_est = d_ll_est,
    it_num = it_num,
    h = h,
    H = H,
    eval_order = eval_order,
    order_method = order_method
  )
}

#' Sequential multivariate LLQR with warm starts along a chosen traversal.
#'
#' @param x Covariate matrix.
#' @param y Response vector.
#' @param tau Quantile level.
#' @param z Evaluation-point matrix. Defaults to x.
#' @param h Optional bandwidth.
#' @param tol Convergence tolerance.
#' @param maxit Maximum simplex iterations per evaluation point.
#' @param bland Whether to use Bland's rule.
#' @param track_order If TRUE, return results in the input order of z.
#' @param order_method Traversal rule: "mst", "input", "random", or "sorted".
#' @param root_method Root selection for MST traversal.
#' @param distance_scale Either "bandwidth" or "none". The MST is built on
#'   bandwidth-scaled coordinates when "bandwidth" is used.
#' @param traversal_seed Optional seed used for random traversal.
#' @return A list containing estimates, iteration counts, bandwidth, and
#'   traversal metadata.
llqr_tau_seq_multivar <- function(x, y, tau, z = NULL, h = NULL,
                                  tol = 1e-14, maxit = 10000,
                                  bland = FALSE, track_order = FALSE,
                                  order_method = NULL,
                                  root_method = c("center", "first"),
                                  distance_scale = c("bandwidth", "none"),
                                  traversal_seed = NULL) {
  validated <- .llqr_multivar_validate_inputs(x = x, y = y, tau = tau, z = z)
  x <- validated$x
  y <- validated$y
  z <- validated$z

  m <- nrow(x)
  nvar <- ncol(x)
  rounds <- nrow(z)
  root_method <- match.arg(root_method)
  distance_scale <- match.arg(distance_scale)

  h <- compute_llqr_multivar_bandwidth(x = x, y = y, tau = tau, h = h)
  if (is.null(order_method)) {
    order_method <- if (nvar == 1L) "sorted" else "mst"
  }

  traversal <- build_llqr_multivar_traversal(
    z = z,
    order_method = order_method,
    root_method = root_method,
    distance_scale = if (distance_scale == "bandwidth") h else NULL,
    seed = traversal_seed
  )

  eval_order <- traversal$order
  parent <- traversal$parent
  root <- traversal$root
  edge_weight <- traversal$edge_weight

  children_remaining <- integer(rounds)
  if (rounds > 1L) {
    for (ii in seq_len(rounds)) {
      if (parent[ii] > 0L) {
        children_remaining[parent[ii]] <- children_remaining[parent[ii]] + 1L
      }
    }
  }

  init_state <- .llqr_multivar_init_state(x = x, y = y)

  ll_est_work <- numeric(rounds)
  d_ll_est_work <- if (nvar == 1L) numeric(rounds) else matrix(0, nrow = rounds, ncol = nvar)
  if (nvar > 1L) {
    colnames(d_ll_est_work) <- colnames(x)
  }
  it_num_work <- integer(rounds)
  H_work <- matrix(0, nrow = rounds, ncol = 1L + nvar)
  state_cache <- vector("list", rounds)

  for (tt in seq_len(rounds)) {
    idx <- eval_order[tt]
    z_eval <- z[idx, , drop = FALSE]
    w <- .llqr_multivar_compute_weights(x = x, z_eval = z_eval[1, ], h = h)

    donor_state <- if (tt == 1L) init_state else state_cache[[parent[idx]]]
    if (is.null(donor_state)) {
      stop(sprintf("Warm-start state for node %d (parent %d) is missing.", idx, parent[idx]))
    }

    res <- .llqr_multivar_solve_from_state(
      state = donor_state, w = w, tau = tau, tol = tol, maxit = maxit,
      bland = bland, m = m, nvar = nvar
    )
    estimate <- res$estimate

    ll_est_work[tt] <- sum(c(1, z_eval[1, ]) * estimate)
    if (nvar == 1L) {
      d_ll_est_work[tt] <- estimate[2L]
    } else {
      d_ll_est_work[tt, ] <- estimate[2:(nvar + 1L)]
    }
    it_num_work[tt] <- res$it_num
    H_work[tt, ] <- res$H

    if (children_remaining[idx] > 0L) {
      state_cache[idx] <- list(res$state)
    }

    if (tt > 1L) {
      pidx <- parent[idx]
      children_remaining[pidx] <- children_remaining[pidx] - 1L
      if (children_remaining[pidx] == 0L) {
        state_cache[pidx] <- list(NULL)
      }
    }
  }

  if (track_order) {
    inv_perm <- order(eval_order)
    ll_est <- ll_est_work[inv_perm]
    d_ll_est <- if (nvar == 1L) d_ll_est_work[inv_perm] else d_ll_est_work[inv_perm, , drop = FALSE]
    it_num <- it_num_work[inv_perm]
    H <- H_work[inv_perm, , drop = FALSE]
  } else {
    ll_est <- ll_est_work
    d_ll_est <- d_ll_est_work
    it_num <- it_num_work
    H <- H_work
  }

  list(
    ll_est = ll_est,
    d_ll_est = d_ll_est,
    it_num = it_num,
    h = h,
    H = H,
    eval_order = eval_order,
    parent = parent,
    root = root,
    edge_weight = edge_weight,
    order_method = order_method,
    root_method = root_method,
    distance_scale = distance_scale
  )
}
