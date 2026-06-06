# ============================================================================ #
# Generates simulation data for one-dim Local linear quantile regression
# ============================================================================ #
generate_data <- function(n, case = 1, seed = NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }

  case <- llqr_validate_case(case)
  if(case == 1L){
    x <- rnorm(n)
  } else if(case == 2L){
    x <- runif(n)
  }
  
  error <- rnorm(n)
  y <- 1 + 2*x^2 + error
  return(list(x = x,y = y))
}

# ============================================================================ #
# 2x2 Matrix Inverse Solver
# ============================================================================ #
inv22 <- function(mat) {
  # Check if input is a 2x2 matrix
  if (!all(dim(mat) == c(2, 2))) {
    stop("Input must be a 2x2 matrix")
  }
  
  # Calculate determinant
  det <- mat[1, 1] * mat[2, 2] - mat[1, 2] * mat[2, 1]
  
  # Check if matrix is invertible (determinant != 0)
  if (det == 0) {
    stop("Matrix is not invertible (determinant is zero)")
  }
  
  # Compute inverse using the formula for 2x2 matrix:
  # [ d, -b ]
  # [ -c, a ] divided by determinant (ad - bc)
  inv_mat <- matrix(c(mat[2, 2], -mat[2, 1], -mat[1, 2], mat[1, 1]), nrow = 2) / det
  
  return(inv_mat)
}

# ============================================================================ #
# LLQR case, kernel, and bandwidth helpers
# ============================================================================ #
llqr_validate_case <- function(case) {
  case_num <- suppressWarnings(as.numeric(case))
  if (length(case_num) != 1L || is.na(case_num) ||
      !is.finite(case_num) || case_num != floor(case_num) ||
      !(case_num %in% 1:2)) {
    stop("Invalid LLQR case specification. Use 1 or 2.")
  }
  as.integer(case_num)
}

llqr_validate_h_factor <- function(h.factor) {
  h.factor <- as.numeric(h.factor)
  if (length(h.factor) != 1L || is.na(h.factor) ||
      !is.finite(h.factor) || h.factor <= 0) {
    stop("h.factor must be a positive finite scalar.")
  }
  h.factor
}

llqr_uses_epanechnikov <- function(case) {
  llqr_validate_case(case) == 2L
}

llqr_kernel_weights <- function(u, case = 1) {
  case <- llqr_validate_case(case)
  u <- as.numeric(u)
  if (case == 2L) {
    w <- numeric(length(u))
    inside <- abs(u) <= 1
    w[inside] <- 0.75 * (1 - u[inside]^2)
    return(w)
  }
  dnorm(u)
}

llqr_default_bandwidth <- function(x, y, tau, h = NULL, case = 1, h.factor = 1) {
  case <- llqr_validate_case(case)
  h.factor <- llqr_validate_h_factor(h.factor)
  x <- as.matrix(x)
  y <- as.matrix(y)
  m <- nrow(x)
  nvar <- ncol(x)

  if (!is.null(h)) {
    h <- as.numeric(h)
    if (length(h) != 1L || is.na(h) || !is.finite(h) || h <= 0) {
      stop("h must be a positive finite scalar when provided.")
    }
    return(h)
  }

  if (case == 2L) {
    return(as.numeric(h.factor * m^(-0.2)))
  }

  red_dim <- floor(0.2 * m)
  index_y <- order(y)[red_dim:(m - red_dim)]
  h_val <- KernSmooth::dpill(x[index_y, , drop = FALSE], y[index_y])
  h_val <- 1.25 * h_val * (tau * (1 - tau)/(dnorm(qnorm(tau)))^2)^0.2
  if (is.nan(h_val)) {
    h_val <- 1.25 * max(m^(-1/(nvar + 4)), min(2, sd(y))*m^(-1/(nvar + 4)))
  }
  as.numeric(h_val)
}

llqr_threshold_scale <- function(m, case = 1) {
  case <- llqr_validate_case(case)
  if (case == 1L) {
    return(log(log(m)) / sqrt(log(m)))
  }
  log(m)^(1/2) * m^(-2/5)
}

# ============================================================================ #
# LLQR Candidate Certification Helpers
# ============================================================================ #
llqr_quantile_loss <- function(r, tau) {
  r <- as.numeric(r)
  r * (tau - (r < 0))
}

llqr_full_objective <- function(estimate, A, y, w, tau) {
  estimate <- as.numeric(estimate)
  A <- as.matrix(A)
  y <- as.numeric(y)
  w <- as.numeric(w)
  r <- as.numeric(y - A %*% estimate)
  sum(w * llqr_quantile_loss(r, tau))
}

llqr_set_equal_int <- function(x, y) {
  x <- sort(unique(as.integer(x)))
  y <- sort(unique(as.integer(y)))
  identical(x, y)
}

llqr_estimate_from_local_linear <- function(ll_value, d_ll_value, z_value) {
  c(as.numeric(ll_value) - as.numeric(z_value) * as.numeric(d_ll_value), as.numeric(d_ll_value))
}

llqr_classify_candidate <- function(cert) {
  if (!isTRUE(cert$candidate_H_valid)) {
    return("invalid_H")
  }
  if (!isTRUE(cert$candidate_in_zero_set)) {
    return("H_not_in_zero_set")
  }
  "certified"
}

build_llqr_full_cert_record <- function(round, backend,
                                        estimate_candidate, H_candidate,
                                        A, y, w, tau,
                                        baseline_estimate = NULL,
                                        baseline_H = NULL,
                                        residual_tol = 1e-8,
                                        rank_tol = 1e-10,
                                        obj_tol = 1e-10) {
  cert <- certify_llqr_candidate_full(
    estimate_candidate = estimate_candidate,
    H_candidate = H_candidate,
    A = A,
    y = y,
    w = w,
    tau = tau,
    baseline_estimate = baseline_estimate,
    baseline_H = baseline_H,
    residual_tol = residual_tol,
    rank_tol = rank_tol,
    obj_tol = obj_tol
  )

  list(
    round = as.integer(round),
    backend = backend,
    estimate_candidate = as.numeric(estimate_candidate),
    H_candidate = as.integer(H_candidate),
    seq_H = if (!is.null(baseline_H)) as.integer(baseline_H) else NULL,
    certification = cert,
    candidate_class = llqr_classify_candidate(cert),
    fallback_triggered = !isTRUE(cert$cert_ok),
    obj_full = cert$obj_full,
    obj_full_seq = cert$obj_full_seq,
    obj_gap_vs_seq = cert$obj_gap_vs_seq,
    H_recovered_from_full = cert$H_recovered_from_full,
    h_set_match_vs_seq = cert$h_set_match_vs_seq
  )
}

ensure_llqr_fortran_library_loaded <- function(lib_filename) {
  lib_path <- file.path("src", "fortran", lib_filename)
  loaded_paths <- vapply(getLoadedDLLs(), function(x) x[["path"]], character(1), USE.NAMES = FALSE)
  if (!normalizePath(lib_path, winslash = "/", mustWork = FALSE) %in%
      normalizePath(loaded_paths, winslash = "/", mustWork = FALSE)) {
    if (!file.exists(lib_path)) {
      stop(sprintf("Required Fortran library not found: %s", lib_path))
    }
    dyn.load(lib_path)
  }
}

certify_llqr_candidate_full <- function(estimate_candidate, H_candidate, A, y, w, tau,
                                        baseline_estimate = NULL, baseline_H = NULL,
                                        residual_tol = 1e-8, rank_tol = 1e-10,
                                        obj_tol = 1e-10) {
  estimate_candidate <- as.numeric(estimate_candidate)
  H_candidate <- as.integer(H_candidate)
  A <- as.matrix(A)
  y <- as.numeric(y)
  w <- as.numeric(w)

  p <- ncol(A)
  r <- as.numeric(y - A %*% estimate_candidate)
  zero_idx <- which(abs(r) <= residual_tol)

  candidate_H_in_range <- length(H_candidate) == p &&
    !anyNA(H_candidate) &&
    !any(H_candidate < 1L | H_candidate > nrow(A)) &&
    !anyDuplicated(H_candidate)
  candidate_rank_ok <- candidate_H_in_range &&
    (qr(A[H_candidate, , drop = FALSE], tol = rank_tol)$rank == p)
  candidate_zero_ok <- candidate_H_in_range &&
    all(abs(r[H_candidate]) <= residual_tol)
  candidate_H_valid <- candidate_rank_ok && candidate_zero_ok
  candidate_in_zero_set <- candidate_H_in_range && all(H_candidate %in% zero_idx)

  obj_full <- llqr_full_objective(
    estimate = estimate_candidate,
    A = A,
    y = y,
    w = w,
    tau = tau
  )
  obj_full_seq <- if (!is.null(baseline_estimate)) {
    llqr_full_objective(
      estimate = baseline_estimate,
      A = A,
      y = y,
      w = w,
      tau = tau
    )
  } else {
    NA_real_
  }
  obj_gap_vs_seq <- if (is.finite(obj_full_seq)) obj_full - obj_full_seq else NA_real_
  obj_ok_vs_seq <- if (is.finite(obj_gap_vs_seq)) obj_gap_vs_seq <= obj_tol else NA

  h_set_match_vs_seq <- if (!is.null(baseline_H)) {
    llqr_set_equal_int(H_candidate, baseline_H)
  } else {
    NA
  }

  cert_ok <- candidate_H_valid && candidate_in_zero_set

  list(
    cert_ok = cert_ok,
    obj_full = obj_full,
    obj_full_seq = obj_full_seq,
    obj_gap_vs_seq = obj_gap_vs_seq,
    obj_ok_vs_seq = obj_ok_vs_seq,
    H_recovered_from_full = zero_idx,
    candidate_H_valid = candidate_H_valid,
    candidate_in_zero_set = candidate_in_zero_set,
    h_set_match_vs_seq = h_set_match_vs_seq,
    residual = r
  )
}

# ============================================================================ #
# Sequential algorithm for one-dim Local linear quantile regression
# ============================================================================ #
llqr_seq <- function(x, y, tau = 0.5, z = NULL, h = NULL, tol = 1e-14, maxit = 1e6, 
                         bland = F, track_order = F, case = 1, h.factor = 1){
  case <- llqr_validate_case(case)
  h.factor <- llqr_validate_h_factor(h.factor)
  # x must be one-dimensional
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  # If z is not provided, use x; otherwise, convert z to a matrix
  if (is.null(z)) {
    z <- as.matrix(x)
  } else {
    z <- as.matrix(z)
  }
  
  # Always sort z and store the original order if track_order is TRUE
  original_order <- order(z)
  z <- z[original_order, , drop = FALSE]
  
  m <- nrow(x) # subject number
  nvar <- ncol(x) # number of var
  rounds <- nrow(z) # number of evaluation points
  
  h <- llqr_default_bandwidth(x = x, y = y, tau = tau, h = h, case = case, h.factor = h.factor)
  
  
  eva_z <- z[1] - x
  w <- llqr_kernel_weights(eva_z / h, case = case)
  cc <- c(rep(0, 1+nvar), tau*w, (1-tau)*w)
  
  # A matrix, b matrix
  gammax <- cbind(matrix(1, nrow = m, ncol = 1), x)
  b <- c(y, 0)
  gammax[y<0,] <- -gammax[y<0,]
  b[b<0] <- -b[b<0]
  
  # Ib
  IB <- (y>=0)*((1:m)+1+nvar)+(y<0)*((1:m)+1+nvar+m)
  
  gammax <- rbind(gammax, -cc[IB] %*% gammax)
  freevarrow <- c(vector(mode = 'logical', m), TRUE) # these variables cannot be non-basic variables
  r1 <- 1:(nvar+1)
  r2 <- vector('numeric', nvar+1)
  rr <- matrix(0, nrow = 2, ncol = 1+nvar) # last row in the table
  
  # estbeta_all <- matrix(0, nrow = nvar+1, ncol = rounds)
  it_num <- rep(0,rounds)
  ll_est <- rep(0,rounds)
  d_ll_est <- rep(0,rounds)
  H <- matrix(0, nrow = m, ncol = 1+nvar)
  
  for (rd in 1:rounds){
    if (rd>=2){
      eva_z <- z[rd] - x
      w <- llqr_kernel_weights(eva_z / h, case = case)
      
      cc[(2+nvar):(1+nvar+m)] <- tau*w
      cc[(2+nvar+m):(1+nvar+2*m)] <- (1-tau)*w
      gammax[(m+1),] <- tau*w[r1-1-nvar] - cc[IB] %*% gammax[1:m,] 
      # freevarrow <- c(vector(mode = 'logical', m), T) # these variables cannot be non-basic variables
      # gammax[(m+1),] <- tau - c(cc[IB],0) %*% gammax
    }
    
    j <- 0
    while (j<maxit) {
      # print(j)
      # start iteration
      # step 2
      rr[1,] <- gammax[(1+m),]
      rr[2,] <- (w[(r1-1-nvar) * (r1-1-nvar>0) + (r1-1-nvar<=0)] - rr[1,]) * (r2!=0) 
      rr[1, r2==0] <- -abs(rr[1, r2==0])
      
      # terminate 
      rrl <- min(rr)
      if (rrl>= -tol){
        break
      }
      
      # step 3
      # choose t as entering variable
      if (bland){
        if (any(rr[1,]< -tol)){
          tsep <- which(rr[1,]< -tol)
          tmp <- r1[tsep]
          t <- min(tmp)
          t_rr <- tsep[which(tmp==t)]
          tsep <- 1
        }else{
          tsep <- which(rr[2,]< -tol)
          tmp <- r2[tsep]
          t <- min(tmp)
          t_rr <- tsep[which(tmp==t)]
          tsep <- 2
        }
      }else{
        tsep <- which(rr==rrl, arr.ind = TRUE)[1,] #which(rr==rrl, T)[1,]
        t_rr <- tsep[2]
        tsep <- tsep[1]
        if (tsep==1){
          t <- r1[t_rr]
        }else{
          t <- r2[t_rr]
        }
      }
      
      # print(t)
      if (r2[t_rr]!=0){ # the entering variable is one of last 2n variables
        # choose k as leaving variable
        # step 4
        if (tsep==1){
          yy <- gammax[, t_rr]
        }else{
          yy <- -gammax[, t_rr]
        }
        
        # step 5: choose k as leaving variable
        k <- b/yy
        if (bland){
          k <- which((min(k[yy>0 & !freevarrow])==k)&(yy>0))
          if (length(k)!=1){
            tmp <- IB[k]
            k <- k[which(tmp==min(tmp))]
          }
        }else{
          k <- which((min(k[yy>0 & !freevarrow])==k)&(yy>0))[1]
        }
        
        # pivoting step 6'
        if (tsep!=1){ # the entering variable is one of p+n+2~p+1+2n variables
          yy[m+1] <- yy[m+1]+w[r1[t_rr]-1-nvar] 
        }
      }
      else{ # the entering variable is one of first (p+1) variables, which are free variables
        yy <- gammax[, t_rr]
        if (yy[m+1]<0){   # r_p < 0, which means beta_p can change from 0 to a positive value
          # step 5
          k <- b/yy
          if (bland){
            k <- which((min(k[yy>0 & !freevarrow])==k)&(yy>0))
            if (length(k)!=1){
              tmp <- IB[k]
              k <- k[which(tmp==min(tmp))]
            }
          }else{
            k <- which((min(k[yy>0 & !freevarrow])==k)&(yy>0))[1]
          }
          
        }
        else{ # r_p >= 0, which means beta_p can change from 0 to a negative value
          k <- -b/yy
          if (bland){
            k <- which((min(k[yy<0 & !freevarrow])==k)&(yy<0))
            if (length(k)!=1){
              tmp <- IB[k]
              k <- k[which(tmp==min(tmp))]
            }
          }else{
            k <- which((min(k[yy<0 & !freevarrow])==k)&(yy<0))[1]  
          }
          
        }
        # pivoting step 6'
        freevarrow[k] <- T
      }
      # print(IB[k])
      
      ee <- yy/yy[k]
      ee[k] <- 1 - 1/yy[k]
      
      if (IB[k]<=(nvar+m+1)){ 
        gammax[, t_rr] <- 0 
        gammax[k,t_rr] <- 1
        r1[t_rr] <- IB[k]
        r2[t_rr] <- IB[k]+m
      }else{                  
        gammax[, t_rr] <- 0
        gammax[k,t_rr] <- -1
        gammax[(m+1),t_rr] <- w[IB[k] - m - nvar - 1] 
        r1[t_rr] <- IB[k]-m
        r2[t_rr] <- IB[k]
      }
      
      # pivoting step 6
      gammax <- gammax - tcrossprod(ee, gammax[k,])
      b <- b - ee * b[k]
      IB[k] <- t
      
      j <- j+1
    }
    
    if (j==maxit){
      warning('Not converge')
    }
    it_num[rd] <- j
    H[rd,] <- r1 - 1 - nvar
    
    tmp <- IB %in% 1:(nvar+1)
    estimate <- b[tmp][order(IB[tmp])]
    
    if (length(estimate)!=(1+nvar)){
      estimate <- vector('numeric', nvar+1)
      for (i in 1:(nvar+1)){
        try(estimate[i] <- b[IB==i], silent = T) # some b[IB==i] may not exist, which means this beta is 0.
      }
    }
    
    # estbeta_all[, rd] <- estimate
    ll_est[rd] <- crossprod(c(1,z[rd]),estimate)
    d_ll_est[rd] <- estimate[2]
    
  }
  
  # Reorder results back to the original order if track_order is TRUE
  if (track_order) {
    ll_est <- ll_est[order(original_order)]
    d_ll_est <- d_ll_est[order(original_order)]
    # it_num <- it_num[order(original_order)]
  }
  
  return(list(ll_est = ll_est, d_ll_est = d_ll_est, it_num = it_num, h = h, H_seq = H))
}

# ============================================================================ #
# Preprocessing algorithm for one-dim LLQR
# ============================================================================ #
#Local linear quantile regression for one-dim predictors 
llqr_ppro <- function(x, y, tau = 0.5, z = NULL, h = NULL, Mm.factor = 1e-3, 
                      track_order = F, case = 1, h.factor = 1, pmethod = NULL){
  case <- llqr_validate_case(case)
  h.factor <- llqr_validate_h_factor(h.factor)
  # x must be one-dimensional
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  # If z is not provided, use x; otherwise, convert z to a matrix
  if (is.null(z)) {
    z <- as.matrix(x)
  } else {
    z <- as.matrix(z)
  }
  
  # Always sort z and store the original order if track_order is TRUE
  original_order <- order(z)
  z <- z[original_order, , drop = FALSE]
  
  m <- nrow(x) # subject number
  nvar <- ncol(x) # number of var
  rounds <- nrow(z) # number of evaluation points
  
  h <- llqr_default_bandwidth(x = x, y = y, tau = tau, h = h, case = case, h.factor = h.factor)
  
  x.norms <- apply(x, 1, function(row) sqrt(sum(row^2)))
  mm <- llqr_threshold_scale(m = m, case = case)
  
  # initialize the output
  ll_est <- rep(0,rounds)
  d_ll_est <- rep(0,rounds)
  residual_est <- matrix(0, nrow = m, ncol = m)
  
  # Initial estimation at the first evaluation point
  eva_z <- z[1] - x
  xx <- cbind(matrix(1, nrow = m, ncol = 1), eva_z) # n*(p+1)
  w <- llqr_kernel_weights(eva_z / h, case = case)
  wxx <- apply(xx, 2, function(x) x * w)
  wy <- y * w
  row_all_zero_x <- apply(wxx, 1, function(z) all(z == 0))
  zero_rows <- row_all_zero_x & (wy == 0)
  wxx <- wxx[!zero_rows, , drop = FALSE] # remove zero rows
  wy <- wy[!zero_rows]
  if (!is.null(pmethod)) {
    fit <- quantreg::rq.fit(x = wxx, y = wy, tau = tau, method = pmethod) # rq.fit is faster than rq
  } else {
    fit <- quantreg::rq.fit(x = wxx, y = wy, tau = tau)
  }
  b <- fit$coef
  # r <- z$resid
  r <- y - xx %*% b
  
  ll_est[1] <- b[1]
  d_ll_est[1] <- b[2]
  residual_est[1, ] <- r
  
  
  for (rd in 2:rounds){
    # print(rd)
    not_optimal <- TRUE
    not_new_sl_sh <- TRUE
    mmm <- mm
    eva_z <- z[rd] - x
    w <- llqr_kernel_weights(eva_z / h, case = case)
    xx <- cbind(matrix(1, nrow = m, ncol = 1), eva_z) # n*(p+1)
    wxx <- apply(xx, 2, function(x) x * w)
    wy <- y * w
    
    while (not_optimal) {
      r <- residual_est[rd-1,]
      if (not_new_sl_sh) {
        M <- Mm.factor * mmm * log(log(m))
        sl <- r < - M
        sh <- r > M
      }
      
      wxxs <- wxx[!sh & !sl, ]
      wys <- wy[!sh & !sl]
      if (any(sl)) {
        glob.wx <- colSums(wxx[sl, , drop = FALSE])
        glob.wy <- sum(wy[sl])
        wxxs <- rbind(wxxs, glob.wx)
        wys <- c(wys, glob.wy)
      }
      if (any(sh)) {
        ghib.wx <- colSums(wxx[sh, , drop = FALSE])
        ghib.wy <- sum(wy[sh])
        wxxs <- rbind(wxxs, ghib.wx)
        wys <- c(wys, ghib.wy)  
      }
      rows_all_zero_x <- apply(wxxs, 1, function(z) all(z == 0))
      zero_rows <- rows_all_zero_x & (wys == 0)
      wxxs <- wxxs[!zero_rows, , drop = FALSE] # remove zero rows
      wys <- wys[!zero_rows]
      ms <- nrow(wxxs) # subsample size
      # print(ns)
      if (!is.null(pmethod)) {
        fit <- quantreg::rq.fit(x = wxxs, y = wys, tau = tau, method = pmethod) # rq.fit is faster than rq
      } else {
        fit <- quantreg::rq.fit(x = wxxs, y = wys, tau = tau)
      }  
      # z <- quantreg::rq.fit(x = xxs, y = ys, tau = tau, weights = ws, method = method)
      b <- fit$coef
      # r <- y - crossprod(t(xx), b)
      r <- y - xx %*% b
      sh.bad <- (r < 0) & sh
      sl.bad <- (r > 0) & sl
      bad.signs <- sum(sh.bad | sl.bad)
      if (bad.signs > 0) {
        if (bad.signs > 0.1 * ms) {
          mmm <- 2 * mmm
          not_new_sl_sh <- TRUE
          # cat("Too many fixups:  doubling m at evaluation point ", rd, "\n")
        } else {
          sh <- sh & !sh.bad
          sl <- sl & !sl.bad
          not_new_sl_sh <- FALSE
          # cat("Fixed some signs at evaluation point ", rd, "\n")
        }
      }
      else not_optimal <- FALSE
    }
    
    ll_est[rd] <- b[1]
    d_ll_est[rd] <- b[2]
    residual_est[rd, ] <- r
  }
  
  # Reorder results back to the original order if track_order is TRUE
  if (track_order) {
    ll_est <- ll_est[order(original_order)]
    d_ll_est <- d_ll_est[order(original_order)]
    # it_num <- it_num[order(original_order)]
  }
  
  return(list(ll_est = ll_est, d_ll_est = d_ll_est, residual_est = residual_est, h = h))
}



# ============================================================================ #
# Sequential plus preprocessing algorithm for one-dim LLQR
# ============================================================================ #
llqr_seq_ppro <- function(x, y, tau = 0.5, z = NULL, h = NULL, tol = 1e-14, maxit = 1e6, 
                              Mm.factor = 1e-3, bland = F, track_order = F,
                              case = 1, h.factor = 1, min_subsample_size = NULL,
                              store_residual = FALSE,
                              fallback = FALSE,
                              debug_trace = FALSE, debug_rounds = NULL){
  case <- llqr_validate_case(case)
  h.factor <- llqr_validate_h_factor(h.factor)
  # x must be one-dimensional
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  # If z is not provided, use x; otherwise, convert z to a matrix
  if (is.null(z)) {
    z <- as.matrix(x)
  } else {
    z <- as.matrix(z)
  }
  
  # Always sort z and store the original order if track_order is TRUE
  original_order <- order(z)
  z <- z[original_order, , drop = FALSE]
  
  m <- nrow(x) # subject number
  nvar <- ncol(x) # number of var
  rounds <- nrow(z) # number of evaluation points
  
  h <- llqr_default_bandwidth(x = x, y = y, tau = tau, h = h, case = case, h.factor = h.factor)
  
  if (is.null(min_subsample_size)){
    min_subsample_size <- max(5 * (nvar + 1), ceiling(0.2 * m))
  }
  residual_tol <- 1e-6
  rank_tol <- 1e-10
  max_empty_pivot_retries <- 3L

  is_valid_llqr_H <- function(H_idx, r_vec) {
    H_idx <- as.integer(H_idx)
    if (length(H_idx) != (nvar + 1) ||
        anyNA(H_idx) ||
        any(H_idx < 1L | H_idx > m) ||
        anyDuplicated(H_idx)) {
      return(FALSE)
    }
    rank_ok <- if (nvar == 1L) {
      abs(A[H_idx[2], 2] - A[H_idx[1], 2]) > rank_tol
    } else {
      qr(A[H_idx, , drop = FALSE], tol = rank_tol)$rank == (nvar + 1)
    }
    rank_ok && (max(abs(r_vec[H_idx])) <= residual_tol)
  }
  
  x.norms <- apply(x, 1, function(row) sqrt(sum(row^2)))
  #mm <- sqrt(log(m)) * (1 / sqrt(m * h) + h^2) * max(x.norms)
  mm <- llqr_threshold_scale(m = m, case = case)
  
  eva_z <- z[1] - x
  w <- llqr_kernel_weights(eva_z / h, case = case)
  cc <- c(rep(0, 1+nvar), tau*w, (1-tau)*w)
  
  # A matrix, b matrix
  A <- cbind(matrix(1, nrow = m, ncol = 1), x)
  gammax <- A
  gammax[y<0,] <- -gammax[y<0,]
  b <- c(y,0)
  b[b<0] <- -b[b<0]
  
  # Ib
  IB <- (y>=0)*((1:m)+1+nvar)+(y<0)*((1:m)+1+nvar+m)
  
  gammax <- rbind(gammax, -cc[IB] %*% gammax)
  freevarrow <- c(vector(mode = 'logical', m), TRUE) # these variables cannot be non-basic variables
  r1 <- 1:(nvar+1)
  r2 <- vector('numeric', nvar+1)
  rr <- matrix(0, nrow = 2, ncol = 1+nvar) # last row in the table
  
  # initialize the output
  # estbeta_all <- matrix(0, nrow = nvar+1, ncol = rounds)
  it_num <- rep(0,rounds)
  ll_est <- rep(0,rounds)
  d_ll_est <- rep(0,rounds)
  # Optional residual-history storage: keeping all m x m residuals is O(m^2) memory.
  residual_est <- if (store_residual) matrix(0, nrow = rounds, ncol = m) else NULL
  acceptance_diagnostics <- vector("list", rounds)
  first_failed_eval <- NA_integer_
  failure_reason <- NULL
  failure_info <- NULL
  debug_round_mask <- rep(FALSE, rounds)
  if (debug_trace) {
    if (is.null(debug_rounds)) {
      debug_round_mask[] <- TRUE
    } else {
      debug_rounds <- unique(as.integer(debug_rounds))
      debug_rounds <- debug_rounds[!is.na(debug_rounds)]
      debug_rounds <- debug_rounds[debug_rounds >= 1L & debug_rounds <= rounds]
      debug_round_mask[debug_rounds] <- TRUE
    }
  }
  round_debug <- if (debug_trace) vector("list", rounds) else NULL
  init_round_debug <- function(rd) {
    if (!debug_trace || !debug_round_mask[rd]) {
      return(invisible(NULL))
    }
    round_debug[[rd]] <<- list(
      round = rd,
      attempts = list(),
      final_status = "pending"
    )
  }
  append_round_attempt <- function(rd, attempt_info) {
    if (!debug_trace || !debug_round_mask[rd]) {
      return(invisible(NULL))
    }
    if (is.null(round_debug[[rd]])) {
      init_round_debug(rd)
    }
    round_debug[[rd]]$attempts[[length(round_debug[[rd]]$attempts) + 1L]] <<- attempt_info
  }
  finalize_round_debug <- function(rd, status, info = list()) {
    if (!debug_trace || !debug_round_mask[rd]) {
      return(invisible(NULL))
    }
    if (is.null(round_debug[[rd]])) {
      init_round_debug(rd)
    }
    round_debug[[rd]]$final_status <<- status
    if (length(info) > 0L) {
      for (nm in names(info)) {
        round_debug[[rd]][[nm]] <<- info[[nm]]
      }
    }
  }
  r_prev <- rep(0, m)
  n_sub <- rep(0, rounds)
  H_seq <- matrix(0L, nrow = rounds, ncol = nvar+1)
  n_sub[1] <- m

  make_return <- function(returned_backend,
                          fallback_triggered = FALSE,
                          fallback_reason = NULL,
                          M_value = NA_real_) {
    ll_est_out <- ll_est
    d_ll_est_out <- d_ll_est
    if (track_order) {
      ll_est_out <- ll_est_out[order(original_order)]
      d_ll_est_out <- d_ll_est_out[order(original_order)]
    }

    list(
      ll_est = ll_est_out,
      d_ll_est = d_ll_est_out,
      it_num = it_num,
      residual_est = residual_est,
      h = h,
      M = M_value,
      n_sub = n_sub,
      H_seq = H_seq,
      acceptance_diagnostics = acceptance_diagnostics,
      round_debug = round_debug,
      first_failed_eval = first_failed_eval,
      failure_reason = failure_reason,
      failure_info = failure_info,
      fallback_triggered = fallback_triggered,
      fallback_reason = fallback_reason,
      returned_backend = returned_backend
    )
  }

  fail_ppro_suffix <- function(rd_start, reason, info = NULL, current_M = NA_real_) {
    first_failed_eval <<- as.integer(rd_start)
    failure_reason <<- reason
    failure_info <<- info
    idx <- rd_start:rounds
    ll_est[idx] <<- NA_real_
    d_ll_est[idx] <<- NA_real_
    it_num[idx] <<- NA_integer_
    n_sub[idx] <<- NA_integer_
    H_seq[idx, ] <<- NA_integer_
    if (store_residual) {
      residual_est[idx, ] <<- NA_real_
    }
    if (is.null(acceptance_diagnostics[[rd_start]])) {
      acceptance_diagnostics[[rd_start]] <<- list(
        round = as.integer(rd_start),
        reason = reason,
        M = current_M,
        info = info
      )
    }
    finalize_round_debug(
      rd = rd_start,
      status = paste0("ppro_failed_", reason),
      info = list(failure_info = info)
    )
    make_return(returned_backend = "ppro_failed", M_value = current_M)
  }

  finish_with_seq_fallback <- function(rd_start, failure_reason_arg, failure_info_arg = NULL) {
    first_failed_eval <<- as.integer(rd_start)
    failure_reason <<- failure_reason_arg
    failure_info <<- failure_info_arg
    fallback_maxit <- max(as.numeric(maxit), 1e6)
    seq_fit <- llqr_seq(
      x = x,
      y = y,
      tau = tau,
      z = z,
      h = h,
      tol = tol,
      maxit = fallback_maxit,
      bland = bland,
      case = case,
      h.factor = h.factor,
      track_order = FALSE
    )

    ll_est <<- as.numeric(seq_fit$ll_est[seq_len(rounds)])
    d_ll_est <<- as.numeric(seq_fit$d_ll_est[seq_len(rounds)])
    it_num <<- as.integer(seq_fit$it_num[seq_len(rounds)])
    n_sub <<- rep.int(m, rounds)
    H_seq <<- matrix(as.integer(seq_fit$H_seq[seq_len(rounds), , drop = FALSE]),
                    nrow = rounds, ncol = nvar + 1)
    if (store_residual) {
      residual_est <<- matrix(NA_real_, nrow = rounds, ncol = m)
      z_values <- as.numeric(z)
      for (i in seq_len(rounds)) {
        beta <- llqr_estimate_from_local_linear(
          ll_value = ll_est[i],
          d_ll_value = d_ll_est[i],
          z_value = z_values[i]
        )
        residual_est[i, ] <<- as.numeric(y) - as.numeric(cbind(1, x) %*% beta)
      }
    }
    make_return(
      returned_backend = "seq_fallback",
      fallback_triggered = TRUE,
      fallback_reason = failure_reason_arg,
      M_value = NA_real_
    )
  }
  
  ## rd =1 
  {
    j <- 0
    rd <- 1
    while (j<maxit) {
      # print(j)
      # start iteration
      # step 2
      rr[1,] <- gammax[(1+m),]
      rr[2,] <- (w[(r1-1-nvar) * (r1-1-nvar>0) + (r1-1-nvar<=0)] - rr[1,]) * (r2!=0) ################################
      rr[1, r2==0] <- -abs(rr[1, r2==0])
      
      # terminate 
      rrl <- min(rr)
      if (rrl>= -tol){
        break
      }
      
      # step 3
      # choose t as entering variable
      if (bland){
        if (any(rr[1,]< -tol)){
          tsep <- which(rr[1,]< -tol)
          tmp <- r1[tsep]
          t <- min(tmp)
          t_rr <- tsep[which(tmp==t)]
          tsep <- 1
        }else{
          tsep <- which(rr[2,]< -tol)
          tmp <- r2[tsep]
          t <- min(tmp)
          t_rr <- tsep[which(tmp==t)]
          tsep <- 2
        }
      }else{
        tsep <- which(rr==rrl, arr.ind = TRUE)[1,] #which(rr==rrl, T)[1,]
        t_rr <- tsep[2]
        tsep <- tsep[1]
        if (tsep==1){
          t <- r1[t_rr]
        }else{
          t <- r2[t_rr]
        }
      }
      
      # print(t)
      if (r2[t_rr]!=0){ # the entering variable is one of last 2n variables
        # choose k as leaving variable
        # step 4
        if (tsep==1){
          yy <- gammax[, t_rr]
        }else{
          yy <- -gammax[, t_rr]
        }
        
        # step 5: choose k as leaving variable
        k <- b/yy
        if (bland){
          k <- which((min(k[yy>0 & !freevarrow])==k)&(yy>0))
          if (length(k)!=1){
            tmp <- IB[k]
            k <- k[which(tmp==min(tmp))]
          }
        }else{
          k <- which((min(k[yy>0 & !freevarrow])==k)&(yy>0))[1]
        }
        
        # pivoting step 6'
        if (tsep!=1){ # the entering variable is one of p+n+2~p+1+2n variables
          yy[m+1] <- yy[m+1]+w[r1[t_rr]-1-nvar] #####################################
        }
      }
      else{ # the entering variable is one of first (p+1) variables, which are free variables
        yy <- gammax[, t_rr]
        if (yy[m+1]<0){   # r_p < 0, which means beta_p can change from 0 to a positive value
          # step 5
          k <- b/yy
          if (bland){
            k <- which((min(k[yy>0 & !freevarrow])==k)&(yy>0))
            if (length(k)!=1){
              tmp <- IB[k]
              k <- k[which(tmp==min(tmp))]
            }
          }else{
            k <- which((min(k[yy>0 & !freevarrow])==k)&(yy>0))[1]
          }
          
        }
        else{ # r_p >= 0, which means beta_p can change from 0 to a negative value
          k <- -b/yy
          if (bland){
            k <- which((min(k[yy<0 & !freevarrow])==k)&(yy<0))
            if (length(k)!=1){
              tmp <- IB[k]
              k <- k[which(tmp==min(tmp))]
            }
          }else{
            k <- which((min(k[yy<0 & !freevarrow])==k)&(yy<0))[1]  
          }
          
        }
        # pivoting step 6'
        freevarrow[k] <- T
      }
      # print(IB[k])
      
      ee <- yy/yy[k]
      ee[k] <- 1 - 1/yy[k]
      
      if (IB[k]<=(nvar+m+1)){ 
        gammax[, t_rr] <- 0 
        gammax[k,t_rr] <- 1
        r1[t_rr] <- IB[k]
        r2[t_rr] <- IB[k]+m
      }else{                  
        gammax[, t_rr] <- 0
        gammax[k,t_rr] <- -1
        gammax[(m+1),t_rr] <- w[IB[k] - m - nvar - 1] ################################
        r1[t_rr] <- IB[k]-m
        r2[t_rr] <- IB[k]
      }
      
      # pivoting step 6
      gammax <- gammax - tcrossprod(ee, gammax[k,])
      b <- b - ee * b[k]
      IB[k] <- t
      
      j <- j+1
    }
    
    if (j==maxit){
      info <- list(round = rd, maxit = maxit, stage = "initial_full_sample")
      if (isTRUE(fallback)) {
        return(finish_with_seq_fallback(rd, "maxit", info))
      }
      return(fail_ppro_suffix(rd, "maxit", info, current_M = NA_real_))
    }
    it_num[rd] <- j
    
    u <- numeric(m)
    v <- numeric(m)
    
    u_tmp <- IB %in% (nvar+2):(nvar+m+1)
    v_tmp <- IB %in% (nvar+m+2):(nvar+2*m+1)
    tmp <- IB %in% (1:(nvar+1))
    estimate <- b[tmp][order(IB[tmp])]
    u[IB[u_tmp]-(nvar+1)] <- b[1:m][u_tmp] 
    v[IB[v_tmp]-(nvar+1+m)] <- b[1:m][v_tmp]
    
    if (length(estimate)!=(1+nvar)){
      estimate <- vector('numeric', nvar+1)
      for (i in 1:(nvar+1)){
        try(estimate[i] <- b[IB==i], silent = T) # some b[IB==i] may not exist, which means this beta is 0.
      }
    }
    
    # estbeta_all[, rd] <- estimate
    ll_est[rd] <- crossprod(c(1,z[rd]),estimate)
    d_ll_est[rd] <- estimate[1+nvar]
    r_prev <- u - v
    if (store_residual) {
      residual_est[rd, ] <- r_prev
    }
    H <- r1 - 1 - nvar
    H_seq[1, ] <- H
    finalize_round_debug(
      rd = 1L,
      status = "initial_full_sample",
      info = list(
        H_prev = NA_integer_,
        H_final = H,
        estimate = as.numeric(estimate),
        residual = as.numeric(r_prev),
        ll_est = ll_est[1],
        d_ll_est = d_ll_est[1],
        n_sub = m
      )
    )
  }
  
  # we use a big (n+2) rows gammax to store the gammaxs
  gammaxs.temp <- matrix(NA, nrow = m + 2, ncol = nvar + 1)
  bs.temp <- rep(NA, m + 2)

  for (rd in 2:rounds){
    # print(rd)
    not_optimal <- TRUE
    not_new_sl_sh <- TRUE
    force_full_sample <- FALSE
    empty_pivot_count <- 0L
    attempt_counter <- 0L
    eva_z <- z[rd] - x
    w <- llqr_kernel_weights(eva_z / h, case = case)
    mmm <- mm
    M <- NA_real_
    init_round_debug(rd)
    
    j <- 0
    while (not_optimal) {
      attempt_counter <- attempt_counter + 1L
      
      # Only previous residuals are needed to build the current screening sets.
      r <- r_prev
      if (force_full_sample) {
        sl <- rep(FALSE, m)
        sh <- rep(FALSE, m)
      } else if (not_new_sl_sh){
        # Fix 4: Scale threshold to residual magnitude
        residual_scale <- median(abs(r))
        M <- max(Mm.factor * mmm * log(log(m)), 0.1 * residual_scale) #Mm.factor * mmm * log(log(m))
        
        # If too few observations would remain, increase M
        n_potential_S <- sum(abs(r) <= M)
        while (n_potential_S < min_subsample_size && M < max(abs(r))) {
          M <- M * 1.5
          n_potential_S <- sum(abs(r) <= M)
        }
        
        sl <- r < -M
        sh <- r > M
      }
      
      not_jl_or_jh <- !(sl | sh) # not_jl_or_jh is TRUE if r is not in the interval [-M, M]
      
      # ISSUE 2 FIX: Force H observations into subsample
      fix2 = TRUE
      if(fix2){      
        H_prev <- H_seq[rd-1,]
        for (h_idx in H_prev) {
          if (h_idx > 0 && h_idx <= m) {
            if (sl[h_idx] || sh[h_idx]) {
              
              sl[h_idx] <- FALSE
              sh[h_idx] <- FALSE
              not_jl_or_jh[h_idx] <- TRUE
            }
          }
        }
      }
      
      idx_not_jl_or_jh <- which(not_jl_or_jh) # the index of r that are not in the interval [-M, M]
      ms <- sum(not_jl_or_jh) # ms is the number of elements in r that are not in the interval [-M, M]
      ws <- w[not_jl_or_jh]
      
      if(rd == 2) {
        gammaxs.temp[1:ms,] <- A[not_jl_or_jh, ,drop = FALSE]
        bs.temp[1:ms] <- y[not_jl_or_jh]
      }
      
      if (any(sl)) {
        gammaxsl <- A[sl, , drop = FALSE]
        wsl <- w[sl]
        glob.wx <- colSums(gammaxsl * wsl)
        glob.wy <- sum(y[sl] * wsl)
        gammaxs.temp[m + 1,] <- glob.wx
        bs.temp[m + 1] <- glob.wy
        ws <- c(ws, 1) 
        ms <- ms + 1
      }
      if (any(sh)) {
        gammaxsh <- A[sh, , drop = FALSE]
        wsh <- w[sh]
        ghib.wx <- colSums(gammaxsh * wsh)
        ghib.wy <- sum(y[sh] * wsh)
        gammaxs.temp[m + 2,] <- ghib.wx
        bs.temp[m + 2] <- ghib.wy
        ws <- c(ws, 1)
        ms <- ms + 1 # num of observations in the subsample
      }
      # gammaxs_temp <- gammaxs
      
      idpos <- which((r[not_jl_or_jh] > 0))
      idneg <- which((r[not_jl_or_jh] < 0))
      H <- H_seq[rd-1,]
      H <- match(H, idx_not_jl_or_jh) # the 2 observations that are interpolated, index in the subsample
      idpos <- idpos[!idpos %in% H] # avoid numerical error
      idneg <- idneg[!idneg %in% H]
      u.in.IBs <- idpos + (nvar + 1) # new index in the subsample
      v.in.IBs <- idneg + (nvar + 1) + ms # new index in the subsample
      r1 <- H + nvar + 1 
      r2 <- r1 + ms
      
      if (any(sl) && any(sh)) {
        Hbar <- c(idpos,idneg,ms-1,ms)
        IBs <- c(1:(nvar + 1), u.in.IBs, v.in.IBs, nvar+1+2*ms-1, nvar+1+ms) # length(IBs) = ms
        # P <- diag(c(rep(1,length(u.in.IBs)),rep(-1,length(v.in.IBs)),-1,1))
        P <- c(rep(1,length(u.in.IBs)),rep(-1,length(v.in.IBs)),-1,1)
        freevarrow <- c(rep(TRUE,nvar + 1), rep(FALSE,length(u.in.IBs)), 
                        rep(FALSE,length(v.in.IBs)), TRUE, TRUE, TRUE) # 1~2, v_L and u_H not into the nonbasic set
      } else if (any(sl)) {
        Hbar <- c(idpos,idneg,ms)
        IBs <- c(1:(nvar + 1), u.in.IBs, v.in.IBs, nvar+1+2*ms) # length(IBs) = ms
        # P <- diag(c(rep(1,length(u.in.IBs)),rep(-1,length(v.in.IBs)+1)))
        P <- c(rep(1,length(u.in.IBs)),rep(-1,length(v.in.IBs)),-1)
        freevarrow <- c(rep(TRUE,nvar + 1), rep(FALSE,length(u.in.IBs)), 
                        rep(FALSE,length(v.in.IBs)),TRUE, TRUE) # 1~2, v_L not into the nonbasic set
      } else if (any(sh)) {
        Hbar <- c(idpos,idneg,ms)
        IBs <- c(1:(nvar + 1), u.in.IBs, v.in.IBs, nvar+1+ms) # length(IBs) = ms
        # P <- diag(c(rep(1,length(u.in.IBs)),rep(-1,length(v.in.IBs)),1))
        P <- c(rep(1,length(u.in.IBs)),rep(-1,length(v.in.IBs)),1)
        freevarrow <- c(rep(TRUE,nvar + 1), rep(FALSE,length(u.in.IBs)), 
                        rep(FALSE,length(v.in.IBs)), TRUE, TRUE) # 1~2, u_H not into the nonbasic set
      } else {
        Hbar <- c(idpos,idneg)
        IBs <- c(1:(nvar + 1), u.in.IBs, v.in.IBs)
        # P <- diag(c(rep(1,length(u.in.IBs)),rep(-1,length(v.in.IBs))))
        P <- c(rep(1,length(u.in.IBs)),rep(-1,length(v.in.IBs)))
        freevarrow <- c(rep(TRUE,nvar + 1), rep(FALSE,length(u.in.IBs)), 
                        rep(FALSE,length(v.in.IBs)), TRUE) # 1~2 not into the nonbasic set
      }
      
      if (rd == 2){
        if (any(sl) & any(sh)) {
          gammaxs <- gammaxs.temp[c(1:(ms-2),m+1,m+2),]
          bs <- bs.temp[c(1:(ms-2),m+1,m+2)]
          lambda <- c(tau * ws[idpos], (1 - tau) * ws[idneg], 1 - tau, tau)
        } else if (any(sl)) {
          gammaxs <- gammaxs.temp[c(1:(ms-1),m+1),]
          bs <- bs.temp[c(1:(ms-1),m+1)]
          lambda <- c(tau * ws[idpos], (1 - tau) * ws[idneg], 1 - tau)
        } else if (any(sh)) {
          gammaxs <- gammaxs.temp[c(1:(ms-1),m+2),]
          bs <- bs.temp[c(1:(ms-1),m+2)]
          lambda <- c(tau * ws[idpos], (1 - tau) * ws[idneg], tau)
        } else {
          gammaxs <- gammaxs.temp[1:ms,]
          bs <- bs.temp[1:ms]
          lambda <- c(tau * ws[idpos], (1 - tau) * ws[idneg])
        }
        
        xhinv <- inv22(gammaxs[H, ])  
        # Pxhbarxhinv <- P %*% gammaxs[Hbar, ] %*% xhinv  
        Pxhbar <- gammaxs[Hbar, ] * P
        Pxhbarxhinv <- Pxhbar %*% xhinv
        gammaxs <- rbind(xhinv, - Pxhbarxhinv)
        gammaxs <- rbind(gammaxs, tau * ws[H] + t(lambda) %*% Pxhbarxhinv)
        # bs <- c(xhinv %*% bs[H], - Pxhbarxhinv %*% bs[H] + P %*% bs[Hbar])
        bs <- c(xhinv %*% bs[H], - Pxhbarxhinv %*% bs[H] + bs[Hbar] * P)
        bs <- c(bs, 0)
      }
      else {
        
        xhinv <- gammaxs.temp[1:(nvar + 1),]
        
        idx_Hbar_pos2 <- idx_not_jl_or_jh[idpos]
        matched_rows_pos <- idx_Hbar_pos2 %in% idx_Hbar_pos
        matched_idpos <-  idx_Hbar_pos2[matched_rows_pos]
        rows.pos <- match(matched_idpos, idx_Hbar_pos) # The row indices in gammaxs.pos of last time to be transfered to new gammaxs 
        rows.pos2 <- which(matched_rows_pos)
        indices <- (nvar+2):(nvar+1+length(u.in.IBs))
        gammaxs.temp[indices[rows.pos2], ] <- gammaxs.pos[rows.pos, ]
        bs.temp[indices[rows.pos2]] <- bs.pos[rows.pos]
        
        if(any(!matched_rows_pos)) {
          rows.pos2.nm <- idx_Hbar_pos2[!matched_rows_pos] # The rows that not in last gammaxs in original data
          Pxhbarxhinv.pos <- A[rows.pos2.nm,] %*% xhinv # P X(hbar) X(h)^{-1}
          gammaxs.temp[indices[!matched_rows_pos],] <- - Pxhbarxhinv.pos
          bs.temp[indices[!matched_rows_pos]] <- - Pxhbarxhinv.pos %*% y[idx_not_jl_or_jh[H]] + y[rows.pos2.nm] # -P X(hbar) X(h)^{-1} y(h) + P y(hbar)
        }
        
        idx_Hbar_neg2 <- idx_not_jl_or_jh[idneg]
        matched_rows_neg <- idx_Hbar_neg2 %in% idx_Hbar_neg
        matched_idneg <-  idx_Hbar_neg2[matched_rows_neg]
        rows.neg <- match(matched_idneg, idx_Hbar_neg) # The row indices in gammaxs.neg of last time to be transfered to new gammaxs
        rows.neg2 <- which(matched_rows_neg)
        indices <- (nvar+1+length(u.in.IBs)+1):(nvar+1+length(u.in.IBs)+length(v.in.IBs))
        gammaxs.temp[indices[rows.neg2], ] <- gammaxs.neg[rows.neg,]
        bs.temp[indices[rows.neg2]] <- bs.neg[rows.neg]
        
        if(any(!matched_rows_neg)) {
          rows.neg2.nm <- idx_Hbar_neg2[!matched_rows_neg] # The rows that not in last gammaxs in original data
          Pxhbarxhinv.neg <- - A[rows.neg2.nm,] %*% xhinv # P X(hbar) X(h)^{-1}
          gammaxs.temp[indices[!matched_rows_neg], ]<- - Pxhbarxhinv.neg
          bs.temp[indices[!matched_rows_neg]] <- - Pxhbarxhinv.neg %*% y[idx_not_jl_or_jh[H]] - y[rows.neg2.nm] # -P X(hbar) X(h)^{-1} y(h) + P y(hbar)
        }
        
        if (any(sl) & any(sh)) {
          gammaxs.temp[m+1,] <- gammaxs.temp[m+1,] %*% xhinv # -P X_L^T X(h)^{-1} = X_L^T X(h)^{-1}
          gammaxs.temp[m+2,] <- - gammaxs.temp[m+2,] %*% xhinv # -P X_H^T X(h)^{-1} = -X_H^T X(h)^{-1}
          bs.temp[m+1] <- gammaxs.temp[m+1,] %*% y[idx_not_jl_or_jh[H]] - bs.temp[m+1] # -P X_L^T X(h)^{-1} y(h) + P y_L
          bs.temp[m+2] <- gammaxs.temp[m+2,] %*% y[idx_not_jl_or_jh[H]] + bs.temp[m+2] # -P X_H^T X(h)^{-1} y(h) + P y_H
          gammaxs <- gammaxs.temp[c(1:(ms-2),m+1,m+2),]
          bs <- bs.temp[c(1:(ms-2),m+1,m+2)]
          lambda <- c(tau * ws[idpos], (1 - tau) * ws[idneg], 1 - tau, tau)
        } else if (any(sl)) {
          gammaxs.temp[m+1,] <- gammaxs.temp[m+1,] %*% xhinv # -P X_L^T X(h)^{-1}
          bs.temp[m+1] <- gammaxs.temp[m+1,] %*% y[idx_not_jl_or_jh[H]] - bs.temp[m+1] # -P X_L^T X(h)^{-1} y(h) + P y_L
          gammaxs <- gammaxs.temp[c(1:(ms-1),m+1),]
          bs <- bs.temp[c(1:(ms-1),m+1)]
          lambda <- c(tau * ws[idpos], (1 - tau) * ws[idneg], 1 - tau)
        } else if (any(sh)) {
          gammaxs.temp[m+2,] <- - gammaxs.temp[m+2,] %*% xhinv # -P X_H^T X(h)^{-1}
          bs.temp[m+2] <- gammaxs.temp[m+2,] %*% y[idx_not_jl_or_jh[H]] + bs.temp[m+2] # -P X_H^T X(h)^{-1} y(h) + P y_H
          gammaxs <- gammaxs.temp[c(1:(ms-1),m+2),]
          bs <- bs.temp[c(1:(ms-1),m+2)]
          lambda <- c(tau * ws[idpos], (1 - tau) * ws[idneg], tau)
        } else {
          gammaxs <- gammaxs.temp[1:ms,]
          bs <- bs.temp[1:ms]
          lambda <- c(tau * ws[idpos], (1 - tau) * ws[idneg])
        } 
        
        Pxhbarxhinv <- - gammaxs[(nvar+2):ms,] 
        if (length(lambda) != nrow(Pxhbarxhinv)) {
          print("Error: Length of lambda does not match number of columns in Pxhbarxhinv")
        }
        gammaxs <- rbind(gammaxs, tau * ws[H] + t(lambda) %*% Pxhbarxhinv)
        bs <- c(bs,0)
      }
      
      no_pivot_attempt <- FALSE
      while (j < maxit){
        # print(j)
        # start iteration
        # step 2
        rr[1, ] <- gammaxs[ms+1, ]
        rr[2, ] <- (ws[r1 - 1 - nvar] - rr[1, ]) 
        
        # terminate 
        rrl <- min(rr)
        if (rrl >= -tol){
          break
        }
        
        # step 3
        # choose t as the entering variable
        # selected from the set of ui and vi of (p+1) points that are interpolated by the regression quantiles
        if (bland){
          if (any(rr[1,] < -tol)){
            tsep <- which(rr[1,] < -tol)
            tmp <- r1[tsep]
            t <- min(tmp)
            t_rr <- tsep[which(tmp == t)]
            tsep <- 1
          } else {
            tsep <- which(rr[2,] < -tol)
            tmp <- r2[tsep]
            t <- min(tmp)
            t_rr <- tsep[which(tmp == t)]
            tsep <- 2
          }
        } else {
          tsep <- which(rr == rrl, arr.ind = TRUE)[1, ]
          t_rr <- tsep[2]
          tsep <- tsep[1]
          if (tsep == 1){
            t <- r1[t_rr]
          } else {
            t <- r2[t_rr]
          }
        }
        # print the entering variable
        # print(t)
        
        # choose k as the leaving variable
        # step 4
        if (tsep == 1){
          yy <- gammaxs[, t_rr]
        } else {
          yy <- -gammaxs[, t_rr]
        }
        
        # step 5
        k <- bs / yy
        
        if (length(k[yy > 0 & !freevarrow]) == 0) {
          no_pivot_attempt <- TRUE
          break
        }
        
        if (bland){
          k <- which((min(k[yy > 0 & !freevarrow]) == k) & (yy > 0))
          if (length(k) != 1){
            tmp <- IBs[k]
            k <- k[which(tmp == min(tmp))]
          }
        } else {
          k <- which((min(k[yy > 0 & !freevarrow]) == k) & (yy > 0))[1]
        }
        # print the leaving variable
        # print(k)
        
        # pivoting step 6'
        if (tsep != 1){ # the entering variable is one of v_i
          yy[ms + 1] <- yy[ms + 1] + ws[r1[t_rr] - 1 - nvar] #####################################
        }
        
        ee <- yy / yy[k]
        ee[k] <- 1 - 1 / yy[k]
        
        if (IBs[k] <= (ms + nvar + 1)){ 
          gammaxs[, t_rr] <- 0 
          gammaxs[k, t_rr] <- 1
          r1[t_rr] <- IBs[k]
          r2[t_rr] <- IBs[k] + ms
        } else {                  
          gammaxs[, t_rr] <- 0
          gammaxs[k, t_rr] <- -1
          gammaxs[(ms + 1), t_rr] <- ws[IBs[k] - ms - nvar - 1] 
          r1[t_rr] <- IBs[k] - ms
          r2[t_rr] <- IBs[k]
        }
        
        # pivoting step 6
        gammaxs <- gammaxs - tcrossprod(ee, gammaxs[k, ])
        bs <- bs - ee * bs[k]
        IBs[k] <- t
        
        j <- j + 1
      }

      if (isTRUE(no_pivot_attempt)) {
        no_pivot_info <- list(
          round = rd,
          attempt_id = attempt_counter,
          force_full_sample = force_full_sample,
          empty_pivot_count = empty_pivot_count,
          M = M,
          H_prev = H_seq[rd - 1, ],
          r_prev = as.numeric(r_prev),
          sl_idx = which(sl),
          sh_idx = which(sh),
          idx_not_jl_or_jh = idx_not_jl_or_jh,
          ms = ms
        )
        append_round_attempt(
          rd = rd,
          attempt_info = c(list(action = "retry_no_pivot", no_pivot_attempt = TRUE), no_pivot_info)
        )
        if (isTRUE(force_full_sample)) {
          acceptance_diagnostics[[rd]] <- c(list(reason = "no_pivot_full_sample"), no_pivot_info)
          if (isTRUE(fallback)) {
            return(finish_with_seq_fallback(rd, "no_pivot_full_sample", no_pivot_info))
          }
          return(fail_ppro_suffix(rd, "no_pivot_full_sample", no_pivot_info, current_M = M))
        }
        empty_pivot_count <- empty_pivot_count + 1L
        mmm <- 2 * mmm
        not_new_sl_sh <- TRUE
        if (empty_pivot_count >= max_empty_pivot_retries) {
          force_full_sample <- TRUE
        }
        next
      }

      if (j >= maxit) {
        maxit_info <- list(
          round = rd,
          attempt_id = attempt_counter,
          force_full_sample = force_full_sample,
          empty_pivot_count = empty_pivot_count,
          M = M,
          H_prev = H_seq[rd - 1, ],
          r_prev = as.numeric(r_prev),
          sl_idx = which(sl),
          sh_idx = which(sh),
          idx_not_jl_or_jh = idx_not_jl_or_jh,
          ms = ms,
          maxit = maxit
        )
        append_round_attempt(
          rd = rd,
          attempt_info = c(list(action = "retry_maxit", no_pivot_attempt = FALSE), maxit_info)
        )
        if (isTRUE(force_full_sample)) {
          acceptance_diagnostics[[rd]] <- c(list(reason = "maxit_full_sample"), maxit_info)
          if (isTRUE(fallback)) {
            return(finish_with_seq_fallback(rd, "maxit_full_sample", maxit_info))
          }
          return(fail_ppro_suffix(rd, "maxit_full_sample", maxit_info, current_M = M))
        }
        empty_pivot_count <- empty_pivot_count + 1L
        mmm <- 2 * mmm
        not_new_sl_sh <- TRUE
        if (empty_pivot_count >= max_empty_pivot_retries) {
          force_full_sample <- TRUE
        }
        next
      }
      
      estimate <- bs[1:(nvar + 1)]
      
      # Check signs of residuals
      r <- y - A %*% estimate
      sh.bad <- (r < 0) & sh
      sl.bad <- (r > 0) & sl
      bad.signs <- sum(sh.bad | sl.bad)
      H_candidate <- r1 - 1 - nvar
      H_candidate <- idx_not_jl_or_jh[H_candidate]
      cert <- certify_llqr_candidate_full(
        estimate_candidate = estimate,
        H_candidate = H_candidate,
        A = A,
        y = y,
        w = w,
        tau = tau,
        residual_tol = residual_tol,
        rank_tol = rank_tol
      )
      cert_initial <- cert
      cert_refit <- NULL
      refit_attempted <- FALSE
      refit_recovered <- FALSE
      H_refit_eligible <- (bad.signs == 0) &&
        !isTRUE(cert$cert_ok) &&
        length(H_candidate) == (nvar + 1) &&
        !anyNA(H_candidate) &&
        !any(H_candidate < 1L | H_candidate > nrow(A)) &&
        !anyDuplicated(H_candidate)
      if (isTRUE(H_refit_eligible)) {
        H_refit_rank_ok <- qr(A[H_candidate, , drop = FALSE], tol = rank_tol)$rank == ncol(A)
        if (isTRUE(H_refit_rank_ok)) {
          refit_attempted <- TRUE
          estimate_refit <- tryCatch(
            as.numeric(solve(A[H_candidate, , drop = FALSE], y[H_candidate])),
            error = function(e) NULL
          )
          if (!is.null(estimate_refit) &&
              length(estimate_refit) == (nvar + 1) &&
              all(is.finite(estimate_refit))) {
            cert_refit <- certify_llqr_candidate_full(
              estimate_candidate = estimate_refit,
              H_candidate = H_candidate,
              A = A,
              y = y,
              w = w,
              tau = tau,
              residual_tol = residual_tol,
              rank_tol = rank_tol
            )
            if (isTRUE(cert_refit$cert_ok)) {
              estimate <- estimate_refit
              r <- cert_refit$residual
              bs[1:(nvar + 1)] <- estimate_refit
              sh.bad <- (r < 0) & sh
              sl.bad <- (r > 0) & sl
              bad.signs <- sum(sh.bad | sl.bad)
              cert <- cert_refit
              refit_recovered <- bad.signs == 0
            }
          }
        }
      }
      accept_subsample <- (bad.signs == 0) && isTRUE(cert$cert_ok)
      append_round_attempt(
        rd = rd,
        attempt_info = list(
          attempt_id = attempt_counter,
          action = if (accept_subsample) "accept" else "reject_retry",
          no_pivot_attempt = FALSE,
          force_full_sample = force_full_sample,
          empty_pivot_count = empty_pivot_count,
          M = M,
          H_prev = H_seq[rd - 1, ],
          r_prev = as.numeric(r_prev),
          sl_idx = which(sl),
          sh_idx = which(sh),
          idx_not_jl_or_jh = idx_not_jl_or_jh,
          ms = ms,
          H_candidate = H_candidate,
          estimate_candidate = as.numeric(estimate),
          bad_signs = bad.signs,
          local_accept = accept_subsample,
          certification = cert,
          certification_initial = cert_initial,
          certification_refit = cert_refit,
          refit_attempted = refit_attempted,
          refit_recovered = refit_recovered
        )
      )
      if (!accept_subsample) {
        if (bad.signs > 0.1 * ms) { 
          mmm <- 2 * mmm
          not_new_sl_sh <- TRUE
          # cat("Too many fixups:  doubling m at evaluation point", rd, "\n")
        } else if (bad.signs > 0) {
          sh <- sh & !sh.bad
          sl <- sl & !sl.bad
          not_new_sl_sh <- FALSE
          # cat("Some fixups: fixing ", rd, "\n")
        } else {
          cert_failure_reason <- if (!isTRUE(cert$candidate_H_valid)) {
            "candidate_H_invalid"
          } else if (!isTRUE(cert$candidate_in_zero_set)) {
            "candidate_not_in_zero_set"
          } else {
            "candidate_cert_failed"
          }
          cert_failure_info <- list(
            round = rd,
            M = M,
            ms = ms,
            bad_signs = bad.signs,
            H_prev = H_seq[rd - 1, ],
            H_candidate = H_candidate,
            estimate_candidate = as.numeric(estimate),
            sl_idx = which(sl),
            sh_idx = which(sh),
            idx_not_jl_or_jh = idx_not_jl_or_jh,
            certification = cert,
            certification_initial = cert_initial,
            certification_refit = cert_refit,
            refit_attempted = refit_attempted,
            refit_recovered = refit_recovered
          )
          if ((ms >= m) || (!any(sl) && !any(sh)) || isTRUE(force_full_sample)) {
            acceptance_diagnostics[[rd]] <- c(
              list(reason = cert_failure_reason),
              cert_failure_info
            )
            finalize_round_debug(
              rd = rd,
              status = paste0("ppro_failed_", cert_failure_reason),
              info = list(
                H_prev = H_seq[rd - 1, ],
                H_final = H_candidate,
                estimate = as.numeric(estimate),
                residual = as.numeric(r),
                ll_est = as.numeric(crossprod(c(1, z[rd]), estimate)),
                d_ll_est = estimate[1 + nvar],
                n_sub = ms
              )
            )
            if (isTRUE(fallback)) {
              return(finish_with_seq_fallback(rd, cert_failure_reason, cert_failure_info))
            }
            return(fail_ppro_suffix(rd, cert_failure_reason, cert_failure_info, current_M = M))
          }
          mmm <- 2 * mmm
          not_new_sl_sh <- TRUE
        }
      }
      if (accept_subsample) {
        not_optimal <- FALSE
        H <- H_candidate
        r[H] <- 0 # set the residuals of H to 0, avoid numerical issues
        if (rd < m){
          gammaxs.temp[1:(nvar+1),] <- gammaxs[1:(nvar+1),]
          bs.temp[1:(nvar+1)] <- bs[1:(nvar+1)]
          ms.org <- ms - any(sl) - any(sh) # The subsample size without x_L and x_H
          id_gammaxs_Hbar <- IBs[(nvar+2):ms.org] 
          p_Hbar <- ifelse(id_gammaxs_Hbar > nvar+1+ms, -1, 1) # 1 for u and -1 for v
          id_gammaxs_Hbar <- ifelse(id_gammaxs_Hbar > nvar + 1 + ms,
                                    id_gammaxs_Hbar - nvar - 1 - ms,
                                    id_gammaxs_Hbar - nvar - 1) # The (2p+1)~ms rows of gammaxs; indices of Hbar in the subsample
          idx_Hbar <- idx_not_jl_or_jh[id_gammaxs_Hbar] # The indices of Hbar in the original data
          idx_Hbar_pos <- idx_Hbar[p_Hbar==1] # the indices of Hbar of last time in  original data
          idx_Hbar_neg <- idx_Hbar[p_Hbar==-1]
          gammaxs.pos <- gammaxs[(nvar+2):ms.org,, drop=FALSE][p_Hbar==1,,drop=FALSE]
          gammaxs.neg <- gammaxs[(nvar+2):ms.org,, drop=FALSE][p_Hbar==-1,,drop=FALSE]
          bs.pos <- bs[(nvar+2):ms.org][p_Hbar==1]
          bs.neg <- bs[(nvar+2):ms.org][p_Hbar==-1]
        }
      }
    }
    
    
    if (j == maxit){
      info <- list(round = rd, maxit = maxit, stage = "post_acceptance_guard")
      if (isTRUE(fallback)) {
        return(finish_with_seq_fallback(rd, "maxit", info))
      }
      return(fail_ppro_suffix(rd, "maxit", info, current_M = M))
    }
    it_num[rd] <- j
    ll_est[rd] <- crossprod(c(1,z[rd]),estimate)
    d_ll_est[rd] <- estimate[1+nvar]
    r_prev <- r
    if (store_residual) {
      residual_est[rd, ] <- r
    }
    n_sub[rd] <- ms
    H_seq[rd,] <- H
    finalize_round_debug(
      rd = rd,
      status = "accepted",
      info = list(
        H_prev = H_seq[rd - 1, ],
        H_final = H,
        estimate = as.numeric(estimate),
        residual = as.numeric(r),
        ll_est = ll_est[rd],
        d_ll_est = d_ll_est[rd],
        n_sub = ms
      )
    )
  }
  
  return(make_return(
    returned_backend = "ppro",
    M_value = if (exists("M", inherits = FALSE)) M else NA_real_
  ))
}

# ============================================================================ #
# Sequential algorithm for one-dim LLQR (Fortran version)
# ============================================================================ #
llqr_seq_fortran_wrapper <- function(x, y, tau = 0.5, z = NULL, h = NULL, tol = 1e-14,
                                     maxit = 1e6, bland = FALSE, track_order = FALSE,
                                     case = 1, h.factor = 1) {
  ensure_llqr_fortran_library_loaded("llqr_seq.so")
  case <- llqr_validate_case(case)
  h.factor <- llqr_validate_h_factor(h.factor)
  
  # Convert to vectors
  x <- as.vector(x)
  y <- as.vector(y)
  if (is.null(z)) {
    z <- x
  }
  z <- as.vector(z)
  original_order <- order(z)
  z <- z[original_order]
  
  # Get dimensions
  m <- length(x)
  nvar <- 1  # Always 1 for univariate LLQR
  rounds <- length(z)
  
  h <- llqr_default_bandwidth(x = x, y = y, tau = tau, h = h, case = case, h.factor = h.factor)
  
  # Prepare output arrays
  ll_est <- numeric(rounds)
  d_ll_est <- numeric(rounds)
  it_num <- integer(rounds)
  residual_est <- matrix(0.0, nrow = rounds, ncol = m)
  H_mat <- matrix(0L, nrow = rounds, ncol = nvar + 1)
  
  # Call Fortran subroutine
  result <- .Fortran("llqr_seq_fortran",
                     x = as.double(x),
                     y = as.double(y),
                     z = as.double(z),
                     m = as.integer(m),
                     nvar = as.integer(nvar),
                     rounds = as.integer(rounds),
                     tau = as.double(tau),
                     h = as.double(h),
                     tol = as.double(tol),
                     maxit = as.integer(maxit),
                     case_int = as.integer(case),
                     bland_int = as.integer(bland),
                     ll_est = as.double(ll_est),
                     d_ll_est = as.double(d_ll_est),
                     it_num = as.integer(it_num),
                     residual_est = as.double(residual_est),
                     H_mat = as.integer(H_mat))
  
  # Return results matching R function structure
  ll_est <- result$ll_est
  d_ll_est <- result$d_ll_est
  if (track_order) {
    ll_est <- ll_est[order(original_order)]
    d_ll_est <- d_ll_est[order(original_order)]
  }

  list(
    ll_est = ll_est,
    d_ll_est = d_ll_est,
    it_num = result$it_num,
    residual_est = matrix(result$residual_est, nrow = rounds, ncol = m),
    h = result$h,
    H_seq = matrix(result$H_mat, nrow = rounds, ncol = nvar + 1)
  )
}

# ============================================================================ #
# Sequential plus preprocessing algorithm for one-dim LLQR (Fortran version)
# ============================================================================ #
llqr_seq_ppro_fortran_wrapper <- function(x, y, tau = 0.5, z = NULL, h = NULL,
                                      Mm.factor = 1e-3, case = 1, h.factor = 1, tol = 1e-14,
                                      maxit = 1e6, bland = FALSE,
                                      track_order = FALSE,
                                      fallback = FALSE,
                                      return_raw_backend = FALSE,
                                      debug_trace = FALSE,
                                      debug_rounds = NULL,
                                      min_subsample_size = NULL,
                                      store_residual = FALSE,
                                      always_same_h_refit = TRUE) {
  ensure_llqr_fortran_library_loaded("llqr_ppro.so")
  case <- llqr_validate_case(case)
  h.factor <- llqr_validate_h_factor(h.factor)
  x <- as.vector(x)
  y <- as.vector(y)
  
  # Auto-load library if not already loaded
  # if (!.llqr_ppro_loaded) {
  #   load_llqr_ppro()
  # }
  if (is.null(z)) {
    z <- x
  }
  z <- as.vector(z)
  original_order <- order(z)
  z <- z[original_order]
  
  # Setup
  m <- length(y)
  nvar <- 1  # univariate (can be extended for multivariate)
  rounds <- length(z)
  acceptance_diagnostics <- vector("list", rounds)
  certification_log <- vector("list", rounds)
  first_bad_round <- NA_integer_
  fallback_triggered <- FALSE
  cert_fail_detected <- FALSE
  fallback_reason <- NULL
  debug_round_mask <- rep(FALSE, rounds)
  if (debug_trace) {
    if (is.null(debug_rounds)) {
      debug_round_mask[] <- TRUE
    } else {
      debug_rounds <- unique(as.integer(debug_rounds))
      debug_rounds <- debug_rounds[!is.na(debug_rounds)]
      debug_rounds <- debug_rounds[debug_rounds >= 1L & debug_rounds <= rounds]
      debug_round_mask[debug_rounds] <- TRUE
    }
  }
  round_debug <- if (debug_trace) vector("list", rounds) else NULL
  init_round_debug <- function(rd) {
    if (!debug_trace || !debug_round_mask[rd]) {
      return(invisible(NULL))
    }
    round_debug[[rd]] <<- list(
      round = rd,
      attempts = list(),
      final_status = "pending"
    )
  }
  append_round_attempt <- function(rd, attempt_info) {
    if (!debug_trace || !debug_round_mask[rd]) {
      return(invisible(NULL))
    }
    if (is.null(round_debug[[rd]])) {
      init_round_debug(rd)
    }
    round_debug[[rd]]$attempts[[length(round_debug[[rd]]$attempts) + 1L]] <<- attempt_info
  }
  finalize_round_debug <- function(rd, status, info = list()) {
    if (!debug_trace || !debug_round_mask[rd]) {
      return(invisible(NULL))
    }
    if (is.null(round_debug[[rd]])) {
      init_round_debug(rd)
    }
    round_debug[[rd]]$final_status <<- status
    if (length(info) > 0L) {
      for (nm in names(info)) {
        round_debug[[rd]][[nm]] <<- info[[nm]]
      }
    }
  }

  h <- llqr_default_bandwidth(x = x, y = y, tau = tau, h = h, case = case, h.factor = h.factor)
  residual_tol <- 1e-6

  if (is.null(min_subsample_size)) {
    min_subsample_size_in <- max(5L * (nvar + 1L), ceiling(0.2 * m))
  } else {
    min_subsample_size_in <- as.integer(ceiling(as.numeric(min_subsample_size)))
    if (length(min_subsample_size_in) != 1L || is.na(min_subsample_size_in) ||
        !is.finite(min_subsample_size_in) || min_subsample_size_in <= 0L) {
      stop("min_subsample_size must be a positive finite scalar or NULL.")
    }
  }

  bland_int <- if (bland) 1L else 0L
  
  # Call Fortran
  result <- .Fortran("llqr_ppro_fortran",
                     x = as.double(x),
                     y = as.double(y),
                     z = as.double(z),
                     m = as.integer(m),
                     nvar = as.integer(nvar),
                     rounds = as.integer(rounds),
                     tau = as.double(tau),
                     h = as.double(h),
                     tol = as.double(tol),
                     maxit = as.integer(maxit),
                     Mm_factor = as.double(Mm.factor),
                     case_int = as.integer(case),
                     bland_int = as.integer(bland_int),
                     min_subsample_size_in = as.integer(min_subsample_size_in),
                     ll_est = double(rounds),
                     d_ll_est = double(rounds),
                     it_num = integer(rounds),
                     residual_est = matrix(0.0, nrow = rounds, ncol = m),
                     H_mat = matrix(0L, nrow = rounds, ncol = nvar + 1),
                     n_sub_out = integer(rounds),
                     ierr = integer(1),
                     failed_eval = integer(1),
                     always_same_h_refit_int = as.integer(isTRUE(always_same_h_refit)))

  raw_backend <- list(
    ll_est = as.numeric(result$ll_est),
    d_ll_est = as.numeric(result$d_ll_est),
    it_num = as.integer(result$it_num),
    residual_est = result$residual_est,
    H_seq = matrix(result$H_mat, nrow = rounds, ncol = nvar + 1),
    M = NA_real_,
    n_sub = as.integer(result$n_sub_out),
    z = z,
    ierr = as.integer(result$ierr),
    failed_eval = as.integer(result$failed_eval)
  )

  make_return <- function(ll_est_value, d_ll_est_value, it_num_value, residual_est_value,
                          H_seq_value, n_sub_value, returned_backend,
                          first_failed_eval = NA_integer_, failure_reason = NULL,
                          failure_info = NULL, fallback_triggered = FALSE,
                          fallback_reason = NULL) {
    ll_est_out <- as.numeric(ll_est_value)
    d_ll_est_out <- as.numeric(d_ll_est_value)
    if (track_order) {
      ll_est_out <- ll_est_out[order(original_order)]
      d_ll_est_out <- d_ll_est_out[order(original_order)]
    }
    residual_est_out <- if (isTRUE(store_residual)) residual_est_value else NULL

    list(
      ll_est = ll_est_out,
      d_ll_est = d_ll_est_out,
      it_num = as.integer(it_num_value),
      residual_est = residual_est_out,
      H_seq = H_seq_value,
      h = h,
      M = NA_real_,
      n_sub = n_sub_value,
      acceptance_diagnostics = acceptance_diagnostics,
      certification_log = certification_log,
      round_debug = round_debug,
      first_failed_eval = first_failed_eval,
      first_bad_round = first_failed_eval,
      failure_reason = failure_reason,
      failure_info = failure_info,
      fallback_triggered = fallback_triggered,
      cert_fail_detected = !is.na(first_failed_eval),
      fallback_reason = fallback_reason,
      backend_ierr = as.integer(result$ierr),
      returned_backend = returned_backend,
      raw_backend = if (isTRUE(return_raw_backend)) raw_backend else NULL
    )
  }

  fail_ppro_suffix <- function(rd_start, reason, info = NULL) {
    ll_est_failed <- raw_backend$ll_est
    d_ll_est_failed <- raw_backend$d_ll_est
    it_num_failed <- raw_backend$it_num
    residual_failed <- raw_backend$residual_est
    H_failed <- raw_backend$H_seq
    n_sub_failed <- raw_backend$n_sub
    idx <- rd_start:rounds
    ll_est_failed[idx] <- NA_real_
    d_ll_est_failed[idx] <- NA_real_
    it_num_failed[idx] <- NA_integer_
    residual_failed[idx, ] <- NA_real_
    H_failed[idx, ] <- NA_integer_
    n_sub_failed[idx] <- NA_integer_
    if (is.null(acceptance_diagnostics[[rd_start]])) {
      acceptance_diagnostics[[rd_start]] <<- list(
        round = as.integer(rd_start),
        reason = reason,
        info = info
      )
    }
    finalize_round_debug(
      rd = rd_start,
      status = paste0("ppro_failed_", reason),
      info = list(failure_info = info)
    )
    make_return(
      ll_est_value = ll_est_failed,
      d_ll_est_value = d_ll_est_failed,
      it_num_value = it_num_failed,
      residual_est_value = residual_failed,
      H_seq_value = H_failed,
      n_sub_value = n_sub_failed,
      returned_backend = "ppro_failed",
      first_failed_eval = as.integer(rd_start),
      failure_reason = reason,
      failure_info = info
    )
  }

  finish_with_seq_fallback <- function(rd_start, failure_reason_arg = NULL, failure_info_arg = NULL) {
    fallback_maxit <- max(as.numeric(maxit), 1e6)
    fallback_fit <- llqr_seq_fortran_wrapper(
      x = x,
      y = y,
      tau = tau,
      z = z,
      h = h,
      tol = tol,
      maxit = fallback_maxit,
      bland = bland,
      track_order = FALSE,
      case = case,
      h.factor = h.factor
    )
    make_return(
      ll_est_value = fallback_fit$ll_est,
      d_ll_est_value = fallback_fit$d_ll_est,
      it_num_value = fallback_fit$it_num,
      residual_est_value = fallback_fit$residual_est,
      H_seq_value = fallback_fit$H_seq,
      n_sub_value = rep.int(m, rounds),
      returned_backend = "seq_fallback",
      first_failed_eval = as.integer(rd_start),
      failure_reason = failure_reason_arg,
      failure_info = failure_info_arg,
      fallback_triggered = TRUE,
      fallback_reason = failure_reason_arg
    )
  }

  infer_failed_round <- function(H_seq_value) {
    invalid <- which(apply(H_seq_value, 1L, function(h_row) {
      any(is.na(h_row)) || any(h_row < 1L | h_row > m) || anyDuplicated(as.integer(h_row)) > 0L
    }))
    if (length(invalid) > 0L) {
      return(as.integer(invalid[1L]))
    }
    1L
  }

  backend_failed_round <- function() {
    rd_backend <- as.integer(result$failed_eval)
    if (length(rd_backend) == 1L && !is.na(rd_backend) &&
        rd_backend >= 1L && rd_backend <= rounds) {
      return(rd_backend)
    }
    infer_failed_round(raw_backend$H_seq)
  }

  if (!identical(as.integer(result$ierr), 0L)) {
    ierr_value <- as.integer(result$ierr)
    reason <- switch(
      as.character(ierr_value),
      "1" = "ppro_cert_failed",
      "2" = "ppro_simplex_nonconverged",
      "3" = "candidate_H_invalid",
      "4" = "candidate_H_singular",
      "5" = "fortran_invariant_violation",
      "fortran_ierr"
    )
    rd_failed <- backend_failed_round()
    info <- list(ierr = ierr_value, failed_eval = as.integer(result$failed_eval))
    acceptance_diagnostics[[rd_failed]] <- list(
      round = rd_failed,
      reason = reason,
      ierr = ierr_value
    )
    if (isTRUE(fallback)) {
      return(finish_with_seq_fallback(rd_failed, reason, info))
    }
    return(fail_ppro_suffix(rd_failed, reason, info))
  }

  ll_est_raw <- raw_backend$ll_est
  d_ll_est_raw <- raw_backend$d_ll_est
  H_seq_raw <- raw_backend$H_seq
  residual_est_raw <- raw_backend$residual_est
  n_sub <- raw_backend$n_sub
  A <- cbind(1, x)

  for (rd in seq_len(rounds)) {
    init_round_debug(rd)
    estimate_candidate <- llqr_estimate_from_local_linear(
      ll_value = ll_est_raw[rd],
      d_ll_value = d_ll_est_raw[rd],
      z_value = z[rd]
    )
    H_candidate <- as.integer(H_seq_raw[rd, ])
    w <- llqr_kernel_weights((z[rd] - x) / h, case = case)
    cert_record <- build_llqr_full_cert_record(
      round = rd,
      backend = "fortran_ppro",
      estimate_candidate = estimate_candidate,
      H_candidate = H_candidate,
      A = A,
      y = y,
      w = w,
      tau = tau,
      residual_tol = residual_tol
    )
    certification_log[[rd]] <- cert_record
    append_round_attempt(
      rd = rd,
      attempt_info = list(
        attempt_id = 1L,
        action = if (isTRUE(cert_record$certification$cert_ok)) "accept" else "reject_cert",
        H_candidate = H_candidate,
        estimate_candidate = estimate_candidate,
        certification = cert_record$certification,
        candidate_class = cert_record$candidate_class
      )
    )
    if (!isTRUE(cert_record$certification$cert_ok)) {
      reason <- if (!isTRUE(cert_record$certification$candidate_H_valid)) {
        "candidate_H_invalid"
      } else if (!isTRUE(cert_record$certification$candidate_in_zero_set)) {
        "candidate_not_in_zero_set"
      } else {
        "candidate_cert_failed"
      }
      acceptance_diagnostics[[rd]] <- cert_record
      finalize_round_debug(
        rd = rd,
        status = paste0("ppro_failed_", reason),
        info = list(
          H_prev = if (rd > 1L) H_seq_raw[rd - 1L, ] else NA_integer_,
          H_final = H_candidate,
          estimate = estimate_candidate,
          residual = as.numeric(y - A %*% estimate_candidate),
          ll_est = ll_est_raw[rd],
          d_ll_est = d_ll_est_raw[rd],
          n_sub = n_sub[rd]
        )
      )
      if (isTRUE(fallback)) {
        return(finish_with_seq_fallback(rd, reason, cert_record))
      }
      return(fail_ppro_suffix(rd, reason, cert_record))
    }
    finalize_round_debug(
      rd = rd,
      status = "accepted",
      info = list(
        H_prev = if (rd > 1L) H_seq_raw[rd - 1L, ] else NA_integer_,
        H_final = H_candidate,
        estimate = estimate_candidate,
        residual = as.numeric(y - A %*% estimate_candidate),
        ll_est = ll_est_raw[rd],
        d_ll_est = d_ll_est_raw[rd],
        n_sub = n_sub[rd]
      )
    )
  }

  make_return(
    ll_est_value = ll_est_raw,
    d_ll_est_value = d_ll_est_raw,
    it_num_value = raw_backend$it_num,
    residual_est_value = residual_est_raw,
    H_seq_value = H_seq_raw,
    n_sub_value = n_sub,
    returned_backend = "ppro"
  )
}
