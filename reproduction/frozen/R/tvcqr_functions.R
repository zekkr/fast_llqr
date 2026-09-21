# ============================================================================ #
# Generates simulation data for time-varying coefficient quantile regression
# ============================================================================ #
resolve_tvcqr_case <- function(case) {
  if (is.numeric(case)) {
    if (length(case) != 1L || is.na(case) || !case %in% c(1, 2)) {
      stop("Invalid case specification. Use 1 or 2.")
    }
    case_id <- as.integer(case)
  } else if (is.character(case)) {
    case_key <- tolower(trimws(case))
    if (case_key %in% c("1", "case1", "case1_iid", "iid")) {
      case_id <- 1L
    } else if (case_key %in% c("2", "case2", "case2_dependent", "dependent")) {
      case_id <- 2L
    } else {
      stop("Invalid case specification. Use 1, 2, or 'case2_dependent'.")
    }
  } else {
    stop("Invalid case specification. Use 1 or 2.")
  }
  
  case_key <- if (case_id == 1L) "case1_iid" else "case2_dependent"
  case_label <- if (case_id == 1L) {
    "Case 1"
  } else {
    "Case 2: locally stationary dependent covariates and errors"
  }
  
  list(
    case_id = case_id,
    case_key = case_key,
    case_label = case_label
  )
}

tvcqr_theta_paths <- function(time_index) {
  list(
    theta0 = sin(2 * pi * time_index),
    theta1 = rep(0.5, length(time_index)),
    theta2 = 2 * log(1 + 2 * time_index),
    theta3 = exp(-(time_index - 0.5)^2)
  )
}

build_truncated_locally_stationary_series <- function(innovations,
                                                      coeff_path,
                                                      J = 100L,
                                                      burn_in = 500L,
                                                      scale = 1) {
  J <- as.integer(J)
  burn_in <- as.integer(burn_in)
  
  lag_matrix <- stats::embed(innovations, J + 1L)
  lag_matrix <- lag_matrix[(burn_in + 1L):(burn_in + length(coeff_path)), , drop = FALSE]
  coeff_matrix <- outer(coeff_path, 0:J, "^")
  
  scale * rowSums(coeff_matrix * lag_matrix)
}

generate_ts_case2_dependent <- function(n, J = 100L, burn_in = 500L) {
  J <- as.integer(J)
  burn_in <- as.integer(burn_in)
  
  if (is.na(J) || J < 0L) {
    stop("J must be a non-negative integer.")
  }
  if (is.na(burn_in) || burn_in < 0L) {
    stop("burn_in must be a non-negative integer.")
  }
  
  time_index <- (1:n) / n
  theta <- tvcqr_theta_paths(time_index)
  
  a_t <- 1 / 2 - (time_index - 1 / 2)^2
  b_t <- 1 / 2 - time_index / 2
  c_t <- 1 / 4 + time_index / 2
  
  total_length <- n + burn_in + J
  zeta <- rnorm(total_length)
  eta <- rnorm(total_length)
  eps <- rnorm(total_length)
  xi_aux <- (eta + eps) / sqrt(2)
  
  # Case 2 replaces the iid design by truncated locally stationary linear processes.
  error <- build_truncated_locally_stationary_series(
    innovations = zeta,
    coeff_path = a_t,
    J = J,
    burn_in = burn_in,
    scale = 1 / 4
  )
  x1 <- build_truncated_locally_stationary_series(
    innovations = xi_aux,
    coeff_path = b_t,
    J = J,
    burn_in = burn_in
  )
  x2 <- build_truncated_locally_stationary_series(
    innovations = eta,
    coeff_path = c_t,
    J = J,
    burn_in = burn_in
  )
  x3 <- stats::rchisq(n, df = 3) / 3
  
  x <- cbind(x1, x2, x3)
  y <- theta$theta0 + theta$theta1 * x[, 1] + theta$theta2 * x[, 2] +
    theta$theta3 * x[, 3] + error
  
  list(
    x = x,
    y = y,
    error = error,
    time_index = time_index,
    theta0 = theta$theta0,
    theta1 = theta$theta1,
    theta2 = theta$theta2,
    theta3 = theta$theta3,
    case = 2L,
    case_key = "case2_dependent",
    case_label = "Case 2: locally stationary dependent covariates and errors",
    J = J,
    burn_in = burn_in
  )
}

generate_ts <- function(n, case = 1, seed = NULL, J = 100L, burn_in = 500L) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  case_info <- resolve_tvcqr_case(case)
  time_index <- (1:n) / n
  theta <- tvcqr_theta_paths(time_index)
  
  if (case_info$case_id == 2L) {
    return(generate_ts_case2_dependent(n = n, J = J, burn_in = burn_in))
  }
  
  x <- matrix(rnorm(n * 3), n, 3)
  error <- rnorm(n, mean = 0, sd = 1)
  y <- theta$theta0 + theta$theta1 * x[, 1] + theta$theta2 * x[, 2] +
    theta$theta3 * x[, 3] + error
  
  return(list(
    x = x,
    y = y,
    error = error,
    time_index = time_index,
    theta0 = theta$theta0,
    theta1 = theta$theta1,
    theta2 = theta$theta2,
    theta3 = theta$theta3,
    case = case_info$case_id,
    case_key = case_info$case_key,
    case_label = case_info$case_label,
    J = as.integer(J),
    burn_in = as.integer(burn_in)
  ))
}

# ============================================================================ #
# TVCQR Candidate Certification Helpers
# ============================================================================ #
tvcqr_quantile_loss <- function(r, tau) {
  r <- as.numeric(r)
  r * (tau - (r < 0))
}

tvcqr_full_objective <- function(estimate, A, y, w, tau) {
  estimate <- as.numeric(estimate)
  A <- as.matrix(A)
  y <- as.numeric(y)
  w <- as.numeric(w)
  r <- as.numeric(y - A %*% estimate)
  sum(w * tvcqr_quantile_loss(r, tau))
}

tvcqr_set_equal_int <- function(x, y) {
  x <- sort(unique(as.integer(x)))
  y <- sort(unique(as.integer(y)))
  identical(x, y)
}

tvcqr_classify_candidate <- function(cert) {
  if (!isTRUE(cert$candidate_H_in_range)) {
    return("invalid_H")
  }
  if (!isTRUE(cert$candidate_rank_ok)) {
    return("singular_H")
  }
  if (!isTRUE(cert$candidate_H_valid)) {
    return("H_not_zero")
  }
  if (!isTRUE(cert$candidate_in_zero_set)) {
    return("H_not_in_zero_set")
  }
  if (!isTRUE(cert$objective_finite)) {
    return("objective_not_finite")
  }
  "certified"
}

certify_tvcqr_candidate_full <- function(estimate_candidate, H_candidate, A, y, w, tau,
                                         baseline_estimate = NULL, baseline_H = NULL,
                                         residual_tol = 1e-8, rank_tol = 1e-10,
                                         obj_tol = 1e-10) {
  estimate_candidate <- as.numeric(estimate_candidate)
  H_candidate <- as.integer(H_candidate)
  A <- as.matrix(A)
  y <- as.numeric(y)
  w <- as.numeric(w)

  p <- ncol(A)
  beta_finite <- length(estimate_candidate) == p && all(is.finite(estimate_candidate))
  r <- if (beta_finite) {
    as.numeric(y - A %*% estimate_candidate)
  } else {
    rep(NA_real_, length(y))
  }
  zero_idx <- which(abs(r) <= residual_tol)

  candidate_H_in_range <- length(H_candidate) == p &&
    !anyNA(H_candidate) &&
    !any(H_candidate < 1L | H_candidate > nrow(A)) &&
    !anyDuplicated(H_candidate)
  candidate_rank_ok <- candidate_H_in_range &&
    (qr(A[H_candidate, , drop = FALSE], tol = rank_tol)$rank == p)
  candidate_zero_ok <- beta_finite && candidate_H_in_range &&
    all(abs(r[H_candidate]) <= residual_tol)
  candidate_H_valid <- candidate_rank_ok && candidate_zero_ok
  candidate_in_zero_set <- candidate_H_in_range && all(H_candidate %in% zero_idx)

  obj_full <- if (beta_finite) {
    tvcqr_full_objective(
      estimate = estimate_candidate,
      A = A,
      y = y,
      w = w,
      tau = tau
    )
  } else {
    NA_real_
  }
  objective_finite <- is.finite(obj_full)
  obj_full_seq <- if (!is.null(baseline_estimate)) {
    tvcqr_full_objective(
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
    tvcqr_set_equal_int(H_candidate, baseline_H)
  } else {
    NA
  }

  psi_mid <- if (beta_finite) tau - as.numeric(r < 0) else rep(NA_real_, length(y))
  kkt_score_mid <- if (beta_finite) {
    as.numeric(crossprod(A, w * psi_mid))
  } else {
    rep(NA_real_, p)
  }
  kkt_score_mid_max_abs <- suppressWarnings(max(abs(kkt_score_mid), na.rm = TRUE))
  if (!is.finite(kkt_score_mid_max_abs)) {
    kkt_score_mid_max_abs <- NA_real_
  }

  cert_ok <- beta_finite && candidate_H_valid && candidate_in_zero_set && objective_finite

  list(
    cert_ok = cert_ok,
    beta_finite = beta_finite,
    obj_full = obj_full,
    obj_full_seq = obj_full_seq,
    obj_gap_vs_seq = obj_gap_vs_seq,
    obj_ok_vs_seq = obj_ok_vs_seq,
    objective_finite = objective_finite,
    H_recovered_from_full = zero_idx,
    candidate_H_in_range = candidate_H_in_range,
    candidate_rank_ok = candidate_rank_ok,
    candidate_zero_ok = candidate_zero_ok,
    candidate_H_valid = candidate_H_valid,
    candidate_in_zero_set = candidate_in_zero_set,
    h_set_match_vs_seq = h_set_match_vs_seq,
    kkt_score_mid = kkt_score_mid,
    kkt_score_mid_max_abs = kkt_score_mid_max_abs,
    residual = r
  )
}

build_tvcqr_full_cert_record <- function(round, backend,
                                         estimate_candidate, H_candidate,
                                         A, y, w, tau,
                                         baseline_estimate = NULL,
                                         baseline_H = NULL,
                                         residual_tol = 1e-8,
                                         rank_tol = 1e-10,
                                         obj_tol = 1e-10) {
  cert <- certify_tvcqr_candidate_full(
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
    candidate_class = tvcqr_classify_candidate(cert),
    fallback_triggered = !isTRUE(cert$cert_ok),
    obj_full = cert$obj_full,
    obj_full_seq = cert$obj_full_seq,
    obj_gap_vs_seq = cert$obj_gap_vs_seq,
    H_recovered_from_full = cert$H_recovered_from_full,
    h_set_match_vs_seq = cert$h_set_match_vs_seq
  )
}

# ============================================================================ #
# Time-varying coefficient quantile regression by Local linear estimator (Using 'Quantreg')
# ============================================================================ #
tvc_rq <- function (x, y, tau = 0.5, h = NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  m <- nrow(x) # Number of time nodes
  nvar <- ncol(x) # number of pred var
  
  if (is.null(h)) {
    h <- m^{-0.2}
  }
  time_index <- (1:m)/m
  theta_ll_est <- matrix(0, nrow = m, ncol = nvar+1)
  x1 <- cbind(matrix(1, nrow = m, ncol = 1), x)
  
  for (t in 1:m) {
    x2 <- cbind(x, apply(x1, 2, function(x) x * (1:m-t) / m))
    w <- 0.75 * (1 - ((t/m-time_index)/h)^2) * (abs(t/m-time_index) <= h)
    q <- quantreg::rq(y ~ x2, tau = tau, weights = w)
    theta_ll_est[t,] <- coef(q)[1:(nvar+1)]
  }
  
  return(list(theta_ll_est = theta_ll_est, h = h))
}

# ============================================================================ #
# Sequential algorithm for TVCQR 
# ============================================================================ #
tvcqr_seq <- function(x, y, tau = 0.5, h = NULL, h.factor = 1, tol = 1e-14, maxit = 1e6, bland = F){
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  m <- nrow(x) # Number of time nodes
  nvar <- ncol(x) # number of pred var
  
  if (is.null(h)) {
    h <- m^{-0.2} * h.factor
  }
  
  time_index <- (1:m)/m
  w <- 0.75 * (1 - ((1/m-time_index)/h)^2) * (abs(1/m-time_index) <= h) 
  cc <- c(rep(0, 2*(1+nvar)), tau*w, (1-tau)*w)
  
  # A matrix, b matrix
  gammax <- cbind(matrix(1, nrow = m, ncol = 1), x)
  gammax <- cbind(gammax, apply(gammax, 2, function(x) x * (1:m)/m))
  b <- c(y, 0)
  gammax[y<0,] <- -gammax[y<0,]
  b[b<0] <- -b[b<0]
  
  # Ib
  IB <- (y>=0)*((1:m)+2*(1+nvar))+(y<0)*((1:m)+2*(1+nvar)+m)
  
  gammax <- rbind(gammax, -cc[IB] %*% gammax)
  freevarrow <- c(vector(mode = 'logical', m), TRUE) # these variables cannot be non-basic variables
  r1 <- 1:(2*(nvar+1))
  r2 <- vector('numeric', 2*(nvar+1))
  rr <- matrix(0, nrow = 2, ncol = 2*(1+nvar)) # last row in the table
  
  # Initialize the output
  it_num <- rep(0,m)
  theta_ll_est <- matrix(0, nrow = m, ncol = nvar+1)
  pivot_id_list <- vector("list", m)
  residual_est <- matrix(0, nrow = m, ncol = m)
  H_seq <- matrix(0,nrow = m, ncol = 2*(1+nvar))
  
  for (eva_t in 1:m){
    if (eva_t>=2){
      w <- 0.75 * (1 - ((eva_t/m-time_index)/h)^2) * (abs(eva_t/m-time_index) <= h) 
      cc[(3+2*nvar):(2+2*nvar+m)] <- tau*w
      cc[(3+2*nvar+m):(2+2*nvar+2*m)] <- (1-tau)*w
      gammax[(m+1),] <- tau*w[r1-2-2*nvar] - cc[IB] %*% gammax[1:m,] 
    }
    
    j <- 0
    while (j<maxit) {
      # print(j)
      # start iteration
      # step 2
      rr[1,] <- gammax[(1+m),]
      rr[2,] <- (w[(r1-2-2*nvar) * (r1-2-2*nvar>0) + (r1-2-2*nvar<=0)] - rr[1,]) * (r2!=0) 
      rr[1, r2==0] <- -abs(rr[1, r2==0])
      
      # terminate 
      rrl <- min(rr)
      if (rrl>= -tol){
        break
      }
      
      # step 3
      # choose t
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
        tsep <- which(rr==rrl, arr.ind = TRUE)[1,]
        t_rr <- tsep[2]
        tsep <- tsep[1]
        if (tsep==1){
          t <- r1[t_rr]
        }else{
          t <- r2[t_rr]
        }
      }
      # print the entering variable
      # print(t)
      
      if (r2[t_rr]!=0){
        # choose k
        # step 4
        if (tsep==1){
          yy <- gammax[, t_rr]
        }else{
          yy <- -gammax[, t_rr]
        }
        
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
        
        
        # pivoting step 6'
        if (tsep!=1){
          yy[m+1] <- yy[m+1]+w[r1[t_rr]-2-2*nvar] 
        }
      }else{
        yy <- gammax[, t_rr]
        if (yy[m+1]<0){   
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
        else{
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
        freevarrow[k] <- T # B_k cannot be non-basic
      }
      # print the leaving variable
      # print(IB[k])
      
      pivot_id_list[[eva_t]] <- rbind(pivot_id_list[[eva_t]], c(t, IB[k]))
      
      ee <- yy/yy[k]
      ee[k] <- 1 - 1/yy[k]
      
      if (IB[k]<=(m+2*nvar+2)){ 
        gammax[, t_rr] <- 0 
        gammax[k,t_rr] <- 1
        r1[t_rr] <- IB[k]
        r2[t_rr] <- IB[k]+m
      }else{                  
        gammax[, t_rr] <- 0
        gammax[k,t_rr] <- -1
        gammax[(m+1),t_rr] <- w[IB[k] - m - 2*nvar - 2]
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
    it_num[eva_t] <- j
    H_seq[eva_t,] <- r1 - 2*nvar - 2
    
    u <- numeric(m)
    v <- numeric(m)
    
    tmp <- IB %in% 1:(2*(nvar+1))
    u_tmp <- IB %in% (2*(nvar+1)+1):(2*(nvar+1)+m)
    v_tmp <- IB %in% (2*(nvar+1)+m+1):(2*(nvar+1)+2*m)
    estimate <- b[tmp][order(IB[tmp])] 
    # order(v) returns the position of the elements in the original vector after sorting in ascending order
    # so that the estimate elements are in the order of beta_0,...,beta_p, gamma_0,...,gamma_p.
    u[IB[u_tmp]-(2*(nvar+1))] <- b[1:m][u_tmp] # since b is of dim m+1, we need to subtract m+1
    v[IB[v_tmp]-(2*(nvar+1)+m)] <- b[1:m][v_tmp]
    
    if (length(estimate)!=2*(1+nvar)){
      estimate <- vector('numeric', 2*(nvar+1))
      for (i in 1:2*(nvar+1)){
        try(estimate[i] <- b[IB==i], silent = T) # some b[IB==i] may not exist, which means this beta is 0.
      }
    }
    
    theta_ll_est[eva_t, ] <- estimate[1:(nvar+1)] + (eva_t/m)*estimate[(nvar+2):(2*(nvar+1))]
    residual_est[eva_t, ] <- u - v
    
  }
  
  return(list(theta_ll_est = theta_ll_est, it_num = it_num, pivot_id_list = pivot_id_list, residual_est = residual_est, H_seq = H_seq))
}

# ============================================================================ #
# Proprecessing algorithm for TVCQR 
# ============================================================================ #
tvc_rq_ppro <- function (x, y, tau = 0.5, h = NULL, Mm.factor = 1e-4, pmethod = NULL) {
  n <- length(y)
  if (tau < 0 | tau > 1) 
    stop("tau outside (0,1)")
  
  if (nrow(x) != n) 
    stop("x and y don't match n")
  p <- ncol(x)
  
  if (is.null(h)) {
    h <- n^{-0.2}
  }
  x.norms <- apply(x, 1, function(row) sqrt(sum(row^2)))
  m <- log(n)^{4} * h^2 * max(x.norms)  # if assuming x sub-gaussian, directly use h^2 * log(m)^{9/2} rate
  
  theta_ll_est <- matrix(NA, nrow = n, ncol = p+1)
  d_theta_ll_est <- matrix(NA, nrow = n, ncol = p+1)
  residual_est <- matrix(NA, nrow = n, ncol = n)
  
  # Estimate the residuals at time 1
  xx <- cbind(matrix(1, nrow = n, ncol = 1), x)
  xx <- cbind(xx, apply(xx, 2, function(x) x * (1:n-1) / n)) # n*2(p+1)
  w <- 0.75 * (1 - ((1:n-1)/ (n*h))^2) * ((abs(1:n-1) / n) <= h) 
  wxx <- apply(xx, 2, function(x) x * w)
  wy <- y * w
  row_all_zero_x <- apply(wxx, 1, function(z) all(z == 0))
  zero_rows <- row_all_zero_x & (wy == 0)
  wxx <- wxx[!zero_rows, , drop = FALSE] # remove zero rows
  wy <- wy[!zero_rows]
  if (!is.null(pmethod)) {
    z <- quantreg::rq.fit(x = wxx, y = wy, tau = tau, method = pmethod) # rq.fit is faster than rq
  } else {
    z <- quantreg::rq.fit(x = wxx, y = wy, tau = tau)
  }
  b <- z$coef
  # r <- z$resid
  r <- y - xx %*% b
  
  theta_ll_est[1,] <- b[1:(p+1)]
  d_theta_ll_est[1,] <- b[(p+2):(2*(p+1))]
  residual_est[1,] <- r
  
  
  for (t in 2:n) {
    not_optimal <- TRUE
    not_new_sl_sh <- TRUE
    mm <- m
    xx <- cbind(matrix(1, nrow = n, ncol = 1), x)
    xx <- cbind(xx, apply(xx, 2, function(x) x * (1:n-t) / n)) # n*2(p+1)
    w <- 0.75 * (1 - ((1:n-t)/ (n*h))^2) * ((abs(1:n-t) / n) <= h) 
    wxx <- apply(xx, 2, function(x) x * w)
    wy <- y * w
    while (not_optimal) {
      r <- residual_est[t-1,]
      if (not_new_sl_sh) {
        M <- Mm.factor * mm * log(log(n))
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
      ns <- nrow(wxxs) # subsample size
      # print(ns)
      if (!is.null(pmethod)) {
        z <- quantreg::rq.fit(x = wxxs, y = wys, tau = tau, method = pmethod) # rq.fit is faster than rq
      } else {
        z <- quantreg::rq.fit(x = wxxs, y = wys, tau = tau)
      }
      # z <- quantreg::rq.fit(x = wxxs, y = wys, tau = tau, method = method)  
      # z <- quantreg::rq.fit(x = xxs, y = ys, tau = tau, weights = ws, method = method)
      b <- z$coef
      # r <- y - crossprod(t(xx), b)
      r <- y - xx %*% b
      sh.bad <- (r <= 0) & sh
      sl.bad <- (r >= 0) & sl
      bad.signs <- sum(sh.bad | sl.bad)
      if (bad.signs > 0) {
        if (bad.signs > 0.1 * ns) {
          mm <- 2 * mm
          not_new_sl_sh <- TRUE
          cat("Too many fixups:  doubling m at time ", t, "\n")
        } else {
          sh <- sh & !sh.bad
          sl <- sl & !sl.bad
          not_new_sl_sh <- FALSE
          cat("Fixed some signs at time ", t, "\n")
        }
      }
      else not_optimal <- FALSE
    }
    theta_ll_est[t,] <- b[1:(p+1)]
    d_theta_ll_est[t,] <- b[(p+2):(2*(p+1))]
    # residual_est[t, ] <- y - crossprod(t(xx), b)
    # residual_est[t, ] <- y - xx %*% b
    residual_est[t, ] <- r
  }
  
  list(theta_ll_est = theta_ll_est, residual_est = residual_est, h = h)
}

# ============================================================================ #
# Sequential plus preprocessing algorithm for TVCQR # eps seems not used
# ============================================================================ #
tvcqr_seq_ppro <- function(x, y, tau = 0.5, h = NULL, h.factor = 1, tol = 1e-14, maxit = 1e6,
                           bland = FALSE, Mm.factor = 1e-4, eps = 1e-06, cpp_helper = FALSE,
                           store_residual = FALSE, fallback = FALSE,
                           debug_trace = FALSE, debug_rounds = NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  m <- nrow(x) # Number of time nodes
  nvar <- ncol(x) # number of pred var
  
  if (is.null(h)) {
    h <- m^{-0.2} * h.factor
  }
  
  x.norms <- apply(x, 1, function(row) sqrt(sum(row^2)))
  # mm <- sqrt(log(m)) * (1 / sqrt(m * h) + h^2) * max(x.norms)
  mm <- log(m)^{4} * h^2 *max(x.norms) # if assuming x sub-gaussian, directly use h^2 * log(m)^{9/2} rate
  
  time_index <- (1:m) / m
  w <- 0.75 * (1 - ((1/m - time_index) / h)^2) * (abs(1/m - time_index) <= h) 
  cc <- c(rep(0, 2 * (1 + nvar)), tau * w, (1 - tau) * w)
  
  # A matrix, b matrix
  A <- cbind(matrix(1, nrow = m, ncol = 1), x)
  A <- cbind(A, apply(A, 2, function(x) x * (1:m) / m))
  gammax <- A
  gammax[y < 0, ] <- -gammax[y < 0, ]
  b <- c(y, 0)
  b[b<0] <- -b[b<0]
  
  # Ib 
  IB <- (y >= 0) * ((1:m) + 2 * (1 + nvar)) + (y < 0) * ((1:m) + 2 * (1 + nvar) + m) 
  
  gammax <- rbind(gammax, -cc[IB] %*% gammax)
  freevarrow <- c(vector(mode = 'logical', m), TRUE) # these variables cannot be non-basic variables, ensure that beta variables keep in the basis
  r1 <- 1:(2 * (nvar + 1))
  r2 <- vector('numeric', 2 * (nvar + 1))
  rr <- matrix(0, nrow = 2, ncol = 2 * (1 + nvar)) # last row in the table
  
  # Initialize the output
  it_num <- rep(0, m)
  theta_ll_est <- matrix(0, nrow = m, ncol = nvar + 1)
  # Optional residual-history storage: keeping all m x m residuals is O(m^2) memory.
  residual_est <- if (store_residual) matrix(0, nrow = m, ncol = m) else NULL
  r_prev <- rep(0, m)
  n_sub <- rep(0, m)
  H_seq <- matrix(0, nrow = m, ncol = 2*(nvar+1))
  acceptance_diagnostics <- vector("list", m)
  first_failed_eval <- NA_integer_
  failure_reason <- NULL
  failure_info <- NULL
  round_debug <- if (debug_trace) vector("list", m) else NULL
  # test_sl_sh <- rep(0, m)
  n_sub[1] <- m
  min_subsample_size <- max(5 * (2 * (nvar + 1)), ceiling(0.2 * m))
  residual_tol <- 1e-6
  rank_tol <- 1e-10
  max_empty_pivot_retries <- 3L

  is_valid_tvcqr_H <- function(H_idx, r_vec) {
    p <- 2 * (nvar + 1)
    H_idx <- as.integer(H_idx)
    if (length(H_idx) != p ||
        anyNA(H_idx) ||
        any(H_idx < 1L | H_idx > m) ||
        anyDuplicated(H_idx)) {
      return(FALSE)
    }
    AH <- A[H_idx, , drop = FALSE]
    rank_ok <- qr(AH, tol = rank_tol)$rank == p
    rank_ok && (max(abs(r_vec[H_idx])) <= residual_tol)
  }

  certify_tvcqr_candidate_local <- function(estimate_candidate, H_candidate, r_vec) {
    H_candidate <- as.integer(H_candidate)
    r_vec <- as.numeric(r_vec)
    p <- 2 * (nvar + 1)
    zero_idx <- which(abs(r_vec) <= residual_tol)
    candidate_H_in_range <- length(H_candidate) == p &&
      !anyNA(H_candidate) &&
      !any(H_candidate < 1L | H_candidate > m) &&
      !anyDuplicated(H_candidate)
    candidate_rank_ok <- candidate_H_in_range &&
      (qr(A[H_candidate, , drop = FALSE], tol = rank_tol)$rank == p)
    candidate_zero_ok <- candidate_H_in_range &&
      all(abs(r_vec[H_candidate]) <= residual_tol)
    candidate_H_valid <- candidate_rank_ok && candidate_zero_ok
    candidate_in_zero_set <- candidate_H_in_range && all(H_candidate %in% zero_idx)
    list(
      cert_ok = candidate_H_valid && candidate_in_zero_set,
      candidate_H_valid = candidate_H_valid,
      candidate_in_zero_set = candidate_in_zero_set,
      H_recovered_from_full = zero_idx,
      residual = r_vec
    )
  }

  make_return <- function(returned_backend,
                          fallback_triggered = FALSE,
                          fallback_reason = NULL,
                          M_value = NA_real_) {
    list(
      theta_ll_est = theta_ll_est,
      it_num = it_num,
      residual_est = residual_est,
      M = M_value,
      n_sub = n_sub,
      H_seq = H_seq,
      h = h,
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

  fail_ppro_suffix <- function(t_start, reason, info = NULL, current_M = NA_real_) {
    first_failed_eval <<- as.integer(t_start)
    failure_reason <<- reason
    failure_info <<- info
    idx <- t_start:m
    theta_ll_est[idx, ] <<- NA_real_
    it_num[idx] <<- NA_integer_
    n_sub[idx] <<- NA_integer_
    H_seq[idx, ] <<- NA_integer_
    if (store_residual) {
      residual_est[idx, ] <<- NA_real_
    }
    if (is.null(acceptance_diagnostics[[t_start]])) {
      acceptance_diagnostics[[t_start]] <<- list(
        time = as.integer(t_start),
        reason = reason,
        M = current_M,
        info = info
      )
    }
    make_return(returned_backend = "ppro_failed", M_value = current_M)
  }

  finish_with_seq_fallback <- function(t_start, failure_reason_arg, failure_info_arg = NULL) {
    first_failed_eval <<- as.integer(t_start)
    failure_reason <<- failure_reason_arg
    failure_info <<- failure_info_arg
    fallback_maxit <- max(as.numeric(maxit), 1e6)
    seq_fit <- tvcqr_seq(
      x = x,
      y = y,
      tau = tau,
      h = h,
      h.factor = h.factor,
      tol = tol,
      maxit = fallback_maxit,
      bland = bland
    )

    theta_ll_est <<- seq_fit$theta_ll_est
    it_num <<- as.integer(seq_fit$it_num)
    n_sub <<- rep.int(m, m)
    H_seq <<- seq_fit$H_seq
    if (store_residual) {
      residual_est <<- seq_fit$residual_est
    }
    make_return(
      returned_backend = "seq_fallback",
      fallback_triggered = TRUE,
      fallback_reason = failure_reason_arg,
      M_value = NA_real_
    )
  }
  
  # t = 1
  {  
    j <- 0
    while (j < maxit){
      # print(j)
      # start iteration
      # step 2
      rr[1, ] <- gammax[(1 + m), ]
      rr[2, ] <- (w[(r1 - 2 - 2 * nvar) * (r1 - 2 - 2 * nvar > 0) + (r1 - 2 - 2 * nvar <= 0)] - rr[1, ]) * (r2 != 0)
      rr[1, r2 == 0] <- -abs(rr[1, r2 == 0])
      
      # terminate 
      rrl <- min(rr)
      if (rrl >= -tol){
        break
      }
      
      # step 3
      # choose t
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
      
      if (r2[t_rr] != 0){
        # choose k
        # step 4
        if (tsep == 1){
          yy <- gammax[, t_rr]
        } else {
          yy <- -gammax[, t_rr]
        }
        
        # step 5
        k <- b / yy
        if (bland){
          k <- which((min(k[yy > 0 & !freevarrow]) == k) & (yy > 0))
          if (length(k) != 1){
            tmp <- IB[k]
            k <- k[which(tmp == min(tmp))]
          }
        } else {
          k <- which((min(k[yy > 0 & !freevarrow]) == k) & (yy > 0))[1]
        }
        
        # pivoting step 6'
        if (tsep != 1){
          yy[m + 1] <- yy[m + 1] + w[r1[t_rr] - 2 - 2 * nvar] #####################################
        }
      } else {
        yy <- gammax[, t_rr]
        if (yy[m + 1] < 0){   ################################
          # step 5
          k <- b / yy
          if (bland){
            k <- which((min(k[yy > 0 & !freevarrow]) == k) & (yy > 0))
            if (length(k) != 1){
              tmp <- IB[k]
              k <- k[which(tmp == min(tmp))]
            }
          } else {
            k <- which((min(k[yy > 0 & !freevarrow]) == k) & (yy > 0))[1]
          }
        }
        else{
          k <- -b / yy
          if (bland){
            k <- which((min(k[yy < 0 & !freevarrow]) == k) & (yy < 0))
            if (length(k) != 1){
              tmp <- IB[k]
              k <- k[which(tmp == min(tmp))]
            }
          } else {
            k <- which((min(k[yy < 0 & !freevarrow]) == k) & (yy < 0))[1]  
          }
          
        }
        # pivoting step 6'
        freevarrow[k] <- TRUE # B_k cannot be non-basic
      }
      # print the leaving variable
      # print(IB[k])
      
      ee <- yy / yy[k]
      ee[k] <- 1 - 1 / yy[k]
      
      if (IB[k] <= (m + 2 * nvar + 2)){ 
        gammax[, t_rr] <- 0 
        gammax[k, t_rr] <- 1
        r1[t_rr] <- IB[k]
        r2[t_rr] <- IB[k] + m
      } else {                  
        gammax[, t_rr] <- 0
        gammax[k, t_rr] <- -1
        gammax[(m + 1), t_rr] <- w[IB[k] - m - 2 * nvar - 2] ################################
        r1[t_rr] <- IB[k] - m
        r2[t_rr] <- IB[k]
      }
      
      # pivoting step 6
      gammax <- gammax - tcrossprod(ee, gammax[k, ])
      b <- b - ee * b[k]
      IB[k] <- t
      
      j <- j + 1
    }
    
    if (j == maxit){
      info <- list(time = 1L, maxit = maxit, stage = "initial_full_sample")
      if (isTRUE(fallback)) {
        return(finish_with_seq_fallback(1L, "maxit", info))
      }
      return(fail_ppro_suffix(1L, "maxit", info, current_M = NA_real_))
    }
    it_num[1] <- j
    
    u <- numeric(m)
    v <- numeric(m)
    
    tmp <- IB %in% 1:(2*(nvar+1))
    u_tmp <- IB %in% (2*(nvar+1)+1):(2*(nvar+1)+m)
    v_tmp <- IB %in% (2*(nvar+1)+m+1):(2*(nvar+1)+2*m)
    estimate <- b[tmp][order(IB[tmp])] 
    # order(v) returns the position of the elements in the original vector after sorting in ascending order
    # so that the estimate elements are in the order of beta_0,...,beta_p, gamma_0,...,gamma_p.
    u[IB[u_tmp]-(2*(nvar+1))] <- b[1:m][u_tmp] # since b is of dim m+1, we need to subtract m+1
    v[IB[v_tmp]-(2*(nvar+1)+m)] <- b[1:m][v_tmp]
    
    if (length(estimate) != 2 * (1 + nvar)){
      estimate <- vector('numeric', 2 * (nvar + 1))
      for (i in 1:2 * (nvar + 1)){
        try(estimate[i] <- b[IB == i], silent = TRUE) # some b[IB==i] may not exist, which means this beta is 0.
      }
    }
    
    theta_ll_est[1, ] <- estimate[1:(nvar + 1)] + (1 / m) * estimate[(nvar + 2):(2 * (nvar + 1))]
    # residual_est[1, ] <- y - x %*% theta_ll_est[1, ]
    r_prev <- u - v
    if (store_residual) {
      residual_est[1, ] <- r_prev
    }
    H <- r1 - 2 - 2 * nvar
    H_seq[1, ] <- H
  }
  
  # we use a big (n+2) rows gammax to store the gammaxs
  gammaxs.temp <- matrix(NA, nrow = m + 2, ncol = 2*(nvar + 1))
  bs.temp <- rep(NA, m + 2)

  for (eva_t in 2:m){
    # print(eva_t)
    not_optimal <- TRUE
    not_new_sl_sh <- TRUE
    force_full_sample <- FALSE
    empty_pivot_count <- 0L
    w <- 0.75 * (1 - ((eva_t / m - time_index) / h)^2) * (abs(eva_t / m - time_index) <= h) 
    mmm <- mm 
    j <- 0
    while (not_optimal) {
      
      # Only previous residuals are needed to construct the next screening set.
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
      
      not_jl_or_jh <- !(sl | sh)
      
      # ISSUE 2 FIX: Force H observations into subsample
      fix2 = TRUE
      if(fix2){      
        H_prev <- H_seq[eva_t-1,]
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
      
      
      
      idx_not_jl_or_jh <- which(not_jl_or_jh) # The index of subsample in the original data
      ms <- sum(not_jl_or_jh)  
      ws <- w[not_jl_or_jh]
      
      if (eva_t == 2) {
        gammaxs.temp[1:ms,] <- A[not_jl_or_jh, , drop = FALSE]
        bs.temp[1:ms] <- y[not_jl_or_jh]
      }
      
      if (any(sl)) { 
        gammaxsl <- A[sl, , drop = FALSE]
        wsl <- w[sl]
        if (cpp_helper) {
          glob.wx <- fast_weighted_colsums_v2(gammaxsl, wsl) 
        } else {
          glob.wx <- colSums(gammaxsl * wsl)  
        } 
        glob.wy <- sum(y[sl] * wsl) 
        gammaxs.temp[m + 1,] <- glob.wx
        bs.temp[m + 1] <- glob.wy
        ws <- c(ws, 1)
        ms <- ms + 1
      }
      
      if (any(sh)) {
        gammaxsh <- A[sh, , drop = FALSE]
        wsh <- w[sh]
        if (cpp_helper) {
          ghib.wx <- fast_weighted_colsums_v2(gammaxsh, wsh) 
        } else {
          ghib.wx <- colSums(gammaxsh * wsh)  
        }  
        ghib.wy <- sum(y[sh] * wsh)  
        gammaxs.temp[m + 2,] <- ghib.wx
        bs.temp[m + 2] <- ghib.wy
        ws <- c(ws, 1)
        ms <- ms + 1
      }
      
      idpos <- which((r[not_jl_or_jh] > 0)) # the observations whose residuals are positive, index in the subsample
      idneg <- which((r[not_jl_or_jh] < 0)) # the observations whose residuals are negative, index the subsample
      H <- H_seq[eva_t-1,] # the 2p observations that are interpolated, index in the original sample
      H <- match(H, idx_not_jl_or_jh) # the 2p observations that are interpolated, index in the subsample
      idpos <- idpos[!idpos %in% H] # avoid numerical error
      idneg <- idneg[!idneg %in% H]
      u.in.IBs <- idpos + 2*(nvar + 1) 
      v.in.IBs <- idneg + 2*(nvar + 1) + ms 
      r1 <- H + 2*(nvar + 1)
      r2 <- r1 + ms
      
      if (any(sl) && any(sh)) {
        Hbar <- c(idpos, idneg, ms - 1, ms)
        IBs <- c(1:(2*(nvar+1)), u.in.IBs, v.in.IBs, 2*(nvar+1) + 2*ms - 1,  2*(nvar+1) + ms)
        P <- diag(c(rep(1, length(u.in.IBs)), rep(-1, length(v.in.IBs)),-1,1))
        freevarrow <- c(rep(TRUE, 2*(nvar+1)), rep(FALSE, length(u.in.IBs)),  
                        rep(FALSE, length(v.in.IBs)), TRUE, TRUE, TRUE)
      } else if (any(sl)) {
        Hbar <- c(idpos, idneg, ms)
        IBs <- c(1:(2*(nvar+1)), u.in.IBs, v.in.IBs, 2*(nvar+1) + 2*ms)
        P <- diag(c(rep(1, length(u.in.IBs)), rep(-1, length(v.in.IBs) + 1)))
        freevarrow <- c(rep(TRUE, 2*(nvar+1)), rep(FALSE, length(u.in.IBs)), 
                        rep(FALSE, length(v.in.IBs)), TRUE, TRUE)
      } else if (any(sh)) {
        Hbar <- c(idpos, idneg, ms)
        IBs <- c(1:(2*(nvar+1)), u.in.IBs, v.in.IBs, 2*(nvar+1) + ms)
        P <- diag(c(rep(1, length(u.in.IBs)), rep(-1, length(v.in.IBs)),1))
        freevarrow <- c(rep(TRUE, 2*(nvar+1)), rep(FALSE, length(u.in.IBs)), 
                        rep(FALSE, length(v.in.IBs)), TRUE, TRUE)
      } else {
        Hbar <- c(idpos, idneg)
        IBs <- c(1:(2*(nvar+1)), u.in.IBs, v.in.IBs)
        P <- diag(c(rep(1, length(u.in.IBs)), rep(-1, length(v.in.IBs))))
        freevarrow <- c(rep(TRUE, 2*(nvar+1)), rep(FALSE, length(u.in.IBs)), 
                        rep(FALSE, length(v.in.IBs)), TRUE)
      }
      
      if (eva_t == 2){
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
        
        # gammaxs <- na.omit(gammaxs.temp)
        # bs <- na.omit(bs.temp)
        
        xhinv <- solve(gammaxs[H, ])  
        Pxhbarxhinv <- P %*% gammaxs[Hbar, ] %*% xhinv  
        gammaxs <- rbind(xhinv, - Pxhbarxhinv)
        gammaxs <- rbind(gammaxs, tau * ws[H] + t(lambda) %*% Pxhbarxhinv)
        bs <- c(xhinv %*% bs[H], - Pxhbarxhinv %*% bs[H] + P %*% bs[Hbar])
        bs <- c(bs, 0)
        
      } 
      else {
        
        xhinv <- gammaxs.temp[1:(2*(nvar + 1)),]
        
        idx_Hbar_pos2 <- idx_not_jl_or_jh[idpos]
        matched_rows_pos <- idx_Hbar_pos2 %in% idx_Hbar_pos
        matched_idpos <-  idx_Hbar_pos2[matched_rows_pos]
        rows.pos <- match(matched_idpos, idx_Hbar_pos) # The row indices in gammaxs.pos of last time to be transfered to new gammaxs 
        rows.pos2 <- which(matched_rows_pos)
        # gammaxs.temp[(2*(nvar+1)+1):(2*(nvar+1)+length(u.in.IBs)), ][rows.pos2,] <- gammaxs.pos[rows.pos,]
        indices <- (2*(nvar+1)+1):(2*(nvar+1)+length(u.in.IBs))
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
        indices <- (2*(nvar+1)+length(u.in.IBs)+1):(2*(nvar+1)+length(u.in.IBs)+length(v.in.IBs))
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
        
        Pxhbarxhinv <- - gammaxs[(2*(nvar+1)+1):ms,] 
        # if (length(lambda) != nrow(Pxhbarxhinv)) {
        #   print("Error: Length of lambda does not match number of columns in Pxhbarxhinv")
        # }
        gammaxs <- rbind(gammaxs, tau * ws[H] + t(lambda) %*% Pxhbarxhinv)
        bs <- c(bs,0)
      }
      
      # Simplex iteration
      no_pivot_attempt <- FALSE
      while (j < maxit){
        # print(j)
        # start iteration
        # step 2
        rr[1, ] <- gammaxs[(1 + ms), ]
        
        
        rr[2, ] <- (ws[r1 - 2 - 2 * nvar] - rr[1, ])
        
        
        
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
          yy[ms + 1] <- yy[ms + 1] + ws[r1[t_rr] - 2 - 2 * nvar] 
        }
        
        ee <- yy / yy[k]
        ee[k] <- 1 - 1 / yy[k]
        
        if (IBs[k] <= (ms + 2 * nvar + 2)){ 
          gammaxs[, t_rr] <- 0 
          gammaxs[k, t_rr] <- 1
          r1[t_rr] <- IBs[k]
          r2[t_rr] <- IBs[k] + ms
        } else {                  
          gammaxs[, t_rr] <- 0
          gammaxs[k, t_rr] <- -1
          gammaxs[(ms + 1), t_rr] <- ws[IBs[k] - ms - 2 * nvar - 2] 
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
          time = eva_t,
          force_full_sample = force_full_sample,
          empty_pivot_count = empty_pivot_count,
          M = M,
          H_prev = H_seq[eva_t - 1, ],
          r_prev = as.numeric(r_prev),
          sl_idx = which(sl),
          sh_idx = which(sh),
          idx_not_jl_or_jh = idx_not_jl_or_jh,
          ms = ms
        )
        if (isTRUE(force_full_sample)) {
          acceptance_diagnostics[[eva_t]] <- c(list(reason = "no_pivot_full_sample"), no_pivot_info)
          if (isTRUE(fallback)) {
            return(finish_with_seq_fallback(eva_t, "no_pivot_full_sample", no_pivot_info))
          }
          return(fail_ppro_suffix(eva_t, "no_pivot_full_sample", no_pivot_info, current_M = M))
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
          time = eva_t,
          force_full_sample = force_full_sample,
          empty_pivot_count = empty_pivot_count,
          M = M,
          H_prev = H_seq[eva_t - 1, ],
          r_prev = as.numeric(r_prev),
          sl_idx = which(sl),
          sh_idx = which(sh),
          idx_not_jl_or_jh = idx_not_jl_or_jh,
          ms = ms,
          maxit = maxit
        )
        if (isTRUE(force_full_sample)) {
          acceptance_diagnostics[[eva_t]] <- c(list(reason = "maxit_full_sample"), maxit_info)
          if (isTRUE(fallback)) {
            return(finish_with_seq_fallback(eva_t, "maxit_full_sample", maxit_info))
          }
          return(fail_ppro_suffix(eva_t, "maxit_full_sample", maxit_info, current_M = M))
        }
        empty_pivot_count <- empty_pivot_count + 1L
        mmm <- 2 * mmm
        not_new_sl_sh <- TRUE
        if (empty_pivot_count >= max_empty_pivot_retries) {
          force_full_sample <- TRUE
        }
        next
      }
      
      estimate <- bs[1:(2*(nvar + 1))]
      # Check signs of residuals
      if (cpp_helper) {
        r <- y - fast_mat_mult(A, estimate)
      } else {
        r <- y - A %*% estimate
      }
      sh.bad <- (r <= 0) & sh
      sl.bad <- (r >= 0) & sl
      bad.signs <- sum(sh.bad | sl.bad)
      H_candidate <- r1 - 2 - 2 * nvar
      H_candidate <- idx_not_jl_or_jh[H_candidate]
      cert <- certify_tvcqr_candidate_local(
        estimate_candidate = estimate,
        H_candidate = H_candidate,
        r_vec = r
      )
      accept_subsample <- (bad.signs == 0) && isTRUE(cert$cert_ok)
      if (!accept_subsample) {
        if (bad.signs > 0.1 * ms) { 
          mmm <- 2 * mmm
          not_new_sl_sh <- TRUE
          #cat("Too many fixups:  doubling m at time ", eva_t, "\n")
        } else if (bad.signs > 0) {
          sh <- sh & !sh.bad
          sl <- sl & !sl.bad
          not_new_sl_sh <- FALSE
          #cat("Some fixups: fixing ", eva_t, "\n")
        } else {
          cert_failure_reason <- if (!isTRUE(cert$candidate_H_valid)) {
            "candidate_H_invalid"
          } else if (!isTRUE(cert$candidate_in_zero_set)) {
            "candidate_not_in_zero_set"
          } else {
            "candidate_cert_failed"
          }
          cert_failure_info <- list(
            time = eva_t,
            M = M,
            ms = ms,
            bad_signs = bad.signs,
            H_prev = H_seq[eva_t - 1, ],
            H_candidate = H_candidate,
            estimate_candidate = as.numeric(estimate),
            sl_idx = which(sl),
            sh_idx = which(sh),
            idx_not_jl_or_jh = idx_not_jl_or_jh,
            certification = cert
          )
          if ((ms >= m) || (!any(sl) && !any(sh)) || isTRUE(force_full_sample)) {
            acceptance_diagnostics[[eva_t]] <- c(list(reason = cert_failure_reason), cert_failure_info)
            if (isTRUE(fallback)) {
              return(finish_with_seq_fallback(eva_t, cert_failure_reason, cert_failure_info))
            }
            return(fail_ppro_suffix(eva_t, cert_failure_reason, cert_failure_info, current_M = M))
          }
          mmm <- 2 * mmm
          not_new_sl_sh <- TRUE
        }
      }
      if (accept_subsample) { # reach optimality
        not_optimal <- FALSE
        # the index corresponding to r1 in the 1~m raw data
        H <- H_candidate
        r[H] <- 0 # set the residuals of H to 0, avoid numerical issues
        if (eva_t < m){
          gammaxs.temp[1:(2*(nvar+1)),] <- gammaxs[1:(2*(nvar+1)),]
          bs.temp[1:(2*(nvar+1))] <- bs[1:(2*(nvar+1))]
          ms.org <- ms - any(sl) - any(sh) # The subsample size without x_L and x_H
          id_gammaxs_Hbar <- IBs[(2*(nvar+1)+1):ms.org] 
          p_Hbar <- ifelse(id_gammaxs_Hbar > 2*(nvar+1)+ms, -1, 1) # 1 for u and -1 for v
          id_gammaxs_Hbar <- ifelse(id_gammaxs_Hbar > 2*(nvar+1) + ms,
                                    id_gammaxs_Hbar - 2*(nvar+1) - ms,
                                    id_gammaxs_Hbar - 2*(nvar+1)) # The (2p+1)~ms rows of gammaxs; indices of Hbar in the subsample
          idx_Hbar <- idx_not_jl_or_jh[id_gammaxs_Hbar] # The indices of Hbar in the original data
          idx_Hbar_pos <- idx_Hbar[p_Hbar==1] # the indices of Hbar of last time in  original data
          idx_Hbar_neg <- idx_Hbar[p_Hbar==-1]
          gammaxs.pos <- gammaxs[(2*(nvar+1)+1):ms.org,, drop=FALSE][p_Hbar==1,,drop=FALSE]
          gammaxs.neg <- gammaxs[(2*(nvar+1)+1):ms.org,, drop=FALSE][p_Hbar==-1,,drop=FALSE]
          bs.pos <- bs[(2*(nvar+1)+1):ms.org][p_Hbar==1]
          bs.neg <- bs[(2*(nvar+1)+1):ms.org][p_Hbar==-1]
        }
      }
    }
    
    if (j == maxit){
      info <- list(time = eva_t, maxit = maxit, stage = "post_acceptance_guard")
      if (isTRUE(fallback)) {
        return(finish_with_seq_fallback(eva_t, "maxit", info))
      }
      return(fail_ppro_suffix(eva_t, "maxit", info, current_M = M))
    }
    it_num[eva_t] <- j 
    theta_ll_est[eva_t, ] <- estimate[1:(nvar + 1)] + (eva_t / m) * estimate[(nvar + 2):(2 * (nvar + 1))]
    r_prev <- r
    if (store_residual) {
      residual_est[eva_t, ] <- r
    }
    n_sub[eva_t] <- ms
    H_seq[eva_t,] <- H
    # test_sl_sh[eva_t] <- any(sl) + any(sh)
    
  }
  
  return(make_return(
    returned_backend = "ppro",
    M_value = if (exists("M", inherits = FALSE)) M else NA_real_
  ))
}

# ============================================================================ #
# Sequential algorithm for TVCQR (Based on Fortran)
# ============================================================================ #
tvcqr_seq_fortran_wrapper <- function(x, y, tau = 0.5, h = NULL, tol = 1e-14, 
                                      maxit = 1e6, bland = FALSE) {
  
  # Convert inputs to proper format
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  # Get dimensions
  m <- nrow(x)      # Number of time points
  nvar <- ncol(x)   # Number of predictor variables
  
  # Handle default bandwidth
  if (is.null(h)) {
    h <- 0.0  # Fortran will set default h = m^(-0.2)
  }
  
  # Prepare output arrays
  theta_ll_est <- matrix(0.0, nrow = m, ncol = nvar + 1)
  it_num <- integer(m)
  residual_est <- matrix(0.0, nrow = m, ncol = m)
  H_mat <- matrix(0L, nrow = m, ncol = 2 * (nvar + 1))
  
  # Call Fortran subroutine
  # Note: Fortran expects column-major order, which R uses by default
  result <- .Fortran("tvcqr_seq_fortran",
                     x = as.double(x),                    # m x nvar matrix
                     y = as.double(y),                    # m vector
                     m = as.integer(m),                   # number of observations
                     nvar = as.integer(nvar),             # number of predictors
                     tau = as.double(tau),                # quantile level
                     h = as.double(h),                    # bandwidth (in/out)
                     tol = as.double(tol),                # convergence tolerance
                     maxit = as.integer(maxit),           # max iterations
                     bland_int = as.integer(bland),       # 0 or 1 for Bland's rule
                     theta_ll_est = as.double(theta_ll_est),  # output: m x (nvar+1)
                     it_num = as.integer(it_num),             # output: m vector
                     residual_est = as.double(residual_est),  # output: m x m
                     H_mat = as.integer(H_mat)                # output: m x 2(nvar+1)
  )
  
  # Extract results and convert back to matrices
  theta_ll_est <- matrix(result$theta_ll_est, nrow = m, ncol = nvar + 1)
  residual_est <- matrix(result$residual_est, nrow = m, ncol = m)
  H_mat <- matrix(result$H_mat, nrow = m, ncol = 2 * (nvar + 1))
  
  # Note: The original R function returns pivot_id_list, but this is mainly
  # for debugging and not essential for the algorithm's output.
  # We're omitting it here to match the paper's implementation.
  
  # Return results matching original function's structure
  return(list(
    theta_ll_est = theta_ll_est,
    it_num = result$it_num,
    residual_est = residual_est,
    H_seq = H_mat
  ))
}

# ============================================================================ #
# Sequential plus preprocessing algorithm for TVCQR (Based on Fortran)
# ============================================================================ #
tvcqr_seq_ppro_fortran_wrapper <- function(x, y, tau = 0.5, h = NULL, h.factor = 1, 
                                           tol = 1e-14, maxit = 1e6, bland = FALSE, 
                                           Mm.factor = 1e-4, eps = 1e-06,
                                           store_residual = FALSE, fallback = FALSE,
                                           debug_trace = FALSE,
                                           min_subsample_size = NULL,
                                           always_same_h_refit = TRUE,
                                           threshold_lower_bound = TRUE,
                                           threshold_scale_mode = c("loglog", "log")) {
  
  # First, let's check if the Fortran function is properly loaded
  # This helps users identify if they need to compile and load the shared library
  if (!is.loaded("tvcqr_seq_ppro_fortran")) {
    stop("Fortran function 'tvcqr_seq_ppro_fortran' is not loaded.\n",
         "Please compile and load the shared library first using:\n",
         "  1. Compile: gfortran -shared -fPIC -o tvcqr_seq_ppro_fortran.so ",
         "tvcqr_seq_ppro_fortran.f90 -llapack -lblas\n",
         "  2. Load in R: dyn.load('tvcqr_seq_ppro_fortran.so')")
  }
  
  # Convert inputs to matrices to ensure proper dimensions
  # This handles both vector and matrix inputs gracefully
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  # Extract dimensions for clarity and validation
  m <- nrow(x)      # Number of time points
  nvar <- ncol(x)   # Number of predictor variables
  
  # Validate input dimensions
  if (length(y) != m) {
    stop("Length of y (", length(y), ") must equal number of rows in x (", m, ")")
  }
  
  # Validate tau parameter (quantile level)
  if (tau <= 0 || tau >= 1) {
    stop("tau must be between 0 and 1 (exclusive). Current value: ", tau)
  }

  if (is.null(min_subsample_size)) {
    min_subsample_size_in <- max(5L * (2L * (nvar + 1L)), ceiling(0.2 * m))
  } else {
    min_subsample_size_numeric <- as.numeric(min_subsample_size)
    if (length(min_subsample_size_numeric) != 1L ||
        is.na(min_subsample_size_numeric) ||
        !is.finite(min_subsample_size_numeric) ||
        min_subsample_size_numeric < 0) {
      stop("min_subsample_size must be a non-negative finite scalar or NULL.")
    }
    min_subsample_size_in <- as.integer(ceiling(min_subsample_size_numeric))
  }
  if (!is.logical(threshold_lower_bound) ||
      length(threshold_lower_bound) != 1L ||
      is.na(threshold_lower_bound)) {
    stop("threshold_lower_bound must be a non-missing logical scalar.")
  }
  threshold_scale_mode <- match.arg(threshold_scale_mode)
  threshold_scale_mode_int <- switch(
    threshold_scale_mode,
    loglog = 1L,
    log = 2L
  )
  if (!threshold_lower_bound) {
    Mm.factor_numeric <- as.numeric(Mm.factor)
    if (length(Mm.factor_numeric) != 1L ||
        is.na(Mm.factor_numeric) ||
        !is.finite(Mm.factor_numeric) ||
        Mm.factor_numeric <= 0) {
      stop("Mm.factor must be a positive finite scalar when threshold_lower_bound is FALSE.")
    }
  }
  
  # Handle bandwidth parameter
  # If h is NULL or 0, the Fortran code will calculate it using h.factor
  if (is.null(h)) {
    h_value <- 0.0  # Signal to Fortran to calculate bandwidth
  } else {
    h_value <- as.double(h)
    if (h_value <= 0) {
      warning("Provided h <= 0. Will calculate bandwidth using h.factor.")
    }
  }
  
  # Prepare output arrays with proper dimensions
  # These arrays will be filled by the Fortran subroutine
  theta_ll_est <- matrix(0.0, nrow = m, ncol = nvar + 1)
  beta_full_est <- matrix(0.0, nrow = m, ncol = 2 * (nvar + 1))
  it_num <- integer(m)
  residual_buffer <- if (isTRUE(store_residual)) {
    matrix(0.0, nrow = m, ncol = m)
  } else {
    double(1L)
  }
  M_out <- 0.0
  first_n_sub <- integer(m)
  repair_count <- integer(m)
  final_n_sub <- integer(m)
  init_mode <- integer(m)
  init_trigger <- integer(m)
  H_seq <- matrix(0L, nrow = m, ncol = 2 * (nvar + 1))
  same_h_refit_attempted <- integer(m)
  same_h_refit_recovered <- integer(m)
  acceptance_diagnostics <- vector("list", m)
  certification_log <- vector("list", m)
  residual_tol <- 1e-6
  rank_tol <- 1e-10
  
  # Call the Fortran subroutine using .Fortran interface
  # Note: .Fortran always passes by value and returns modified copies
  result <- .Fortran("tvcqr_seq_ppro_fortran",
                     # Input arguments - converted to appropriate types
                     x = as.double(x),                    # Flatten matrix to vector
                     y = as.double(y),                    
                     m = as.integer(m),                   
                     nvar = as.integer(nvar),             
                     tau = as.double(tau),                
                     h = as.double(h_value),              # May be modified by Fortran
                     h_factor = as.double(h.factor),      
                     tol = as.double(tol),                
                     maxit = as.integer(maxit),           
                     bland_int = as.integer(bland),       # Convert logical to integer
                     Mm_factor = as.double(Mm.factor),    
                     eps = as.double(eps),                
                     store_residual_int = as.integer(isTRUE(store_residual)),
                     debug_int = as.integer(isTRUE(debug_trace)),
                     # Output arguments - pre-allocated arrays
                     theta_ll_est = as.double(theta_ll_est),
                     beta_full_est = as.double(beta_full_est),
                     it_num = as.integer(it_num),
                     residual_est = as.double(residual_buffer),
                     M_out = as.double(M_out),
                     first_n_sub = as.integer(first_n_sub),
                     repair_count = as.integer(repair_count),
                     final_n_sub = as.integer(final_n_sub),
                     init_mode = as.integer(init_mode),
                     init_trigger = as.integer(init_trigger),
                     H_seq = as.integer(H_seq),
                     same_h_refit_attempted = as.integer(same_h_refit_attempted),
                     same_h_refit_recovered = as.integer(same_h_refit_recovered),
                     ierr = as.integer(0),
                     failed_eval = as.integer(0),
                     min_subsample_size_in = as.integer(min_subsample_size_in),
                     always_same_h_refit_int = as.integer(isTRUE(always_same_h_refit)),
                     threshold_lower_bound_int = as.integer(threshold_lower_bound),
                     threshold_scale_mode_int = as.integer(threshold_scale_mode_int),
                     # Don't duplicate arrays (more efficient)
                     DUP = FALSE)

  make_return <- function(theta_value, beta_full_value, it_value, residual_value, M_value,
                          first_n_sub_value, repair_count_value, final_n_sub_value,
                          init_mode_value, init_trigger_value, H_value, returned_backend,
                          first_failed_eval = NA_integer_,
                          failure_reason = NULL, failure_info = NULL,
                          fallback_triggered = FALSE, fallback_reason = NULL) {
    failed_eval_value <- if (!is.na(first_failed_eval)) {
      as.integer(first_failed_eval)
    } else {
      as.integer(result$failed_eval)
    }
    init_mode_value <- as.integer(init_mode_value)
    init_trigger_value <- as.integer(init_trigger_value)
    init_mode_labels <- c(
      `0` = "first_point_full_cold",
      `1` = "warm_H",
      `2` = "shifted_independent",
      `3` = "full_active_recovery"
    )
    init_trigger_labels <- c(
      `0` = "none",
      `1` = "invalid_previous_H",
      `2` = "unmappable_previous_H",
      `3` = "singular_previous_H"
    )
    init_mode_out <- unname(init_mode_labels[as.character(init_mode_value)])
    init_trigger_out <- unname(init_trigger_labels[as.character(init_trigger_value)])
    list(
      theta_ll_est = theta_value,
      beta_full_est = beta_full_value,
      it_num = as.integer(it_value),
      residual_est = residual_value,
      M = M_value,
      first_n_sub = as.integer(first_n_sub_value),
      repair_count = as.integer(repair_count_value),
      final_n_sub = as.integer(final_n_sub_value),
      n_sub = as.integer(final_n_sub_value),
      init_mode = init_mode_out,
      init_trigger = init_trigger_out,
      independent_init_count = sum(init_mode_value == 2L, na.rm = TRUE),
      full_active_recovery_count = sum(init_mode_value == 3L, na.rm = TRUE),
      H_seq = H_value,
      h = result$h,
      acceptance_diagnostics = acceptance_diagnostics,
      certification_log = certification_log,
      same_h_refit_attempted = as.integer(result$same_h_refit_attempted),
      same_h_refit_recovered = as.integer(result$same_h_refit_recovered),
      first_failed_eval = first_failed_eval,
      first_bad_round = first_failed_eval,
      failed_eval = failed_eval_value,
      failure_reason = failure_reason,
      failure_info = failure_info,
      fallback_triggered = fallback_triggered,
      cert_fail_detected = !is.na(first_failed_eval),
      fallback_reason = fallback_reason,
      backend_ierr = as.integer(result$ierr),
      backend_failed_eval = as.integer(result$failed_eval),
      returned_backend = returned_backend
    )
  }

  finish_with_seq_fallback <- function(t_start, failure_reason_arg, failure_info_arg = NULL) {
    fallback_maxit <- max(as.numeric(maxit), 1e6)
    fallback_fit <- tvcqr_seq_fortran_wrapper(
      x = x,
      y = y,
      tau = tau,
      h = if (h_value > 0) h_value else NULL,
      tol = tol,
      maxit = fallback_maxit,
      bland = bland
    )
    make_return(
      theta_value = fallback_fit$theta_ll_est,
      beta_full_value = matrix(NA_real_, nrow = m, ncol = 2 * (nvar + 1)),
      it_value = fallback_fit$it_num,
      residual_value = if (isTRUE(store_residual)) fallback_fit$residual_est else NULL,
      M_value = NA_real_,
      first_n_sub_value = rep.int(m, m),
      repair_count_value = rep.int(NA_integer_, m),
      final_n_sub_value = rep.int(m, m),
      init_mode_value = rep.int(NA_integer_, m),
      init_trigger_value = rep.int(NA_integer_, m),
      H_value = fallback_fit$H_seq,
      returned_backend = "seq_fallback",
      first_failed_eval = as.integer(t_start),
      failure_reason = failure_reason_arg,
      failure_info = failure_info_arg,
      fallback_triggered = TRUE,
      fallback_reason = failure_reason_arg
    )
  }

  infer_failed_time <- function(H_value) {
    invalid <- which(apply(H_value, 1L, function(h_row) {
      any(is.na(h_row)) || any(h_row < 1L | h_row > m) || anyDuplicated(as.integer(h_row)) > 0L
    }))
    if (length(invalid) > 0L) {
      return(as.integer(invalid[1L]))
    }
    1L
  }

  backend_failed_time <- function() {
    failed_eval_value <- as.integer(result$failed_eval)
    if (length(failed_eval_value) == 1L && !is.na(failed_eval_value) &&
        failed_eval_value >= 1L && failed_eval_value <= m) {
      return(failed_eval_value)
    }
    infer_failed_time(H_seq)
  }

  failure_reason_from_ierr <- function(ierr_value) {
    switch(
      as.character(ierr_value),
      "1" = "ppro_cert_failed",
      "2" = "ppro_simplex_nonconverged",
      "3" = "candidate_H_invalid",
      "4" = "candidate_H_singular",
      "5" = "fortran_invariant_violation",
      "fortran_ierr"
    )
  }

  fail_ppro_suffix <- function(t_start, reason, info = NULL,
                               theta_value, beta_full_value, it_value, residual_value,
                               M_value, first_n_sub_value, repair_count_value,
                               final_n_sub_value, init_mode_value, init_trigger_value,
                               H_value) {
    idx <- t_start:m
    theta_value[idx, ] <- NA_real_
    beta_full_value[idx, ] <- NA_real_
    it_value[idx] <- NA_integer_
    first_n_sub_value[idx] <- NA_integer_
    repair_count_value[idx] <- NA_integer_
    final_n_sub_value[idx] <- NA_integer_
    init_mode_value[idx] <- NA_integer_
    init_trigger_value[idx] <- NA_integer_
    H_value[idx, ] <- NA_integer_
    if (!is.null(residual_value)) {
      residual_value[idx, ] <- NA_real_
    }
    if (is.null(acceptance_diagnostics[[t_start]])) {
      acceptance_diagnostics[[t_start]] <<- list(
        time = as.integer(t_start),
        reason = reason,
        info = info
      )
    }
    make_return(
      theta_value = theta_value,
      beta_full_value = beta_full_value,
      it_value = it_value,
      residual_value = residual_value,
      M_value = M_value,
      first_n_sub_value = first_n_sub_value,
      repair_count_value = repair_count_value,
      final_n_sub_value = final_n_sub_value,
      init_mode_value = init_mode_value,
      init_trigger_value = init_trigger_value,
      H_value = H_value,
      returned_backend = "ppro_failed",
      first_failed_eval = as.integer(t_start),
      failure_reason = reason,
      failure_info = info
    )
  }
  
  # Reshape the flattened arrays back to matrices
  # Fortran stores matrices in column-major order, same as R
  theta_ll_est <- matrix(result$theta_ll_est, nrow = m, ncol = nvar + 1, byrow = FALSE)
  beta_full_est <- matrix(result$beta_full_est, nrow = m, ncol = 2 * (nvar + 1), byrow = FALSE)
  residual_est <- if (isTRUE(store_residual)) {
    matrix(result$residual_est, nrow = m, ncol = m, byrow = FALSE)
  } else {
    NULL
  }
  H_seq <- matrix(result$H_seq, nrow = m, ncol = 2 * (nvar + 1), byrow = FALSE)

  if (!identical(as.integer(result$ierr), 0L)) {
    ierr_value <- as.integer(result$ierr)
    reason <- failure_reason_from_ierr(ierr_value)
    failed_eval_value <- as.integer(result$failed_eval)
    t_failed <- backend_failed_time()
    info <- list(
      ierr = ierr_value,
      failed_eval = failed_eval_value,
      backend_ierr = ierr_value,
      backend_failed_eval = failed_eval_value,
      failure_source = "fortran"
    )
    acceptance_diagnostics[[t_failed]] <- list(
      time = t_failed,
      reason = reason,
      ierr = ierr_value,
      failed_eval = failed_eval_value
    )
    if (isTRUE(fallback)) {
      return(finish_with_seq_fallback(t_failed, reason, info))
    }
    return(fail_ppro_suffix(
      t_start = t_failed,
      reason = reason,
      info = info,
      theta_value = theta_ll_est,
      beta_full_value = beta_full_est,
      it_value = result$it_num,
      residual_value = residual_est,
      M_value = result$M_out,
      first_n_sub_value = result$first_n_sub,
      repair_count_value = result$repair_count,
      final_n_sub_value = result$final_n_sub,
      init_mode_value = result$init_mode,
      init_trigger_value = result$init_trigger,
      H_value = H_seq
    ))
  }

  A_base <- cbind(1, x)
  time_index <- seq_len(m) / m
  A_full <- cbind(A_base, A_base * time_index)
  h_backend <- as.numeric(result$h)
  for (eva_t in seq_len(m)) {
    H_candidate <- as.integer(H_seq[eva_t, ])
    w_t <- 0.75 * (1 - (((eva_t / m) - time_index) / h_backend)^2) *
      (abs((eva_t / m) - time_index) <= h_backend)
    certification_log[[eva_t]] <- build_tvcqr_full_cert_record(
      round = eva_t,
      backend = "ppro",
      estimate_candidate = beta_full_est[eva_t, ],
      H_candidate = H_candidate,
      A = A_full,
      y = y,
      w = w_t,
      tau = tau,
      residual_tol = residual_tol,
      rank_tol = rank_tol
    )
    certification_log[[eva_t]]$certification$backend_ierr_ok <- TRUE
    certification_log[[eva_t]]$certification$backend_bad_sign_verification_passed <- TRUE
    certification_log[[eva_t]]$same_h_refit_attempted <- as.integer(result$same_h_refit_attempted[eva_t])
    certification_log[[eva_t]]$same_h_refit_recovered <- as.integer(result$same_h_refit_recovered[eva_t])

    cert <- certification_log[[eva_t]]$certification
    if (!isTRUE(cert$cert_ok)) {
      reason <- if (!isTRUE(cert$candidate_H_in_range)) {
        "candidate_H_invalid"
      } else if (!isTRUE(cert$candidate_rank_ok)) {
        "candidate_H_singular"
      } else if (!isTRUE(cert$candidate_H_valid)) {
        "candidate_H_invalid"
      } else {
        "candidate_not_in_zero_set"
      }
      info <- list(
        ierr = 0L,
        failed_eval = as.integer(eva_t),
        backend_ierr = as.integer(result$ierr),
        backend_failed_eval = as.integer(result$failed_eval),
        failure_source = "wrapper_certification",
        time = eva_t,
        H_candidate = H_candidate,
        candidate_H_in_range = cert$candidate_H_in_range,
        candidate_rank_ok = cert$candidate_rank_ok,
        candidate_zero_ok = cert$candidate_zero_ok,
        candidate_H_valid = cert$candidate_H_valid,
        candidate_in_zero_set = cert$candidate_in_zero_set,
        H_recovered_from_full = cert$H_recovered_from_full,
        obj_full = cert$obj_full,
        objective_finite = cert$objective_finite,
        kkt_score_mid_max_abs = cert$kkt_score_mid_max_abs
      )
      acceptance_diagnostics[[eva_t]] <- c(list(reason = reason), info)
      if (isTRUE(fallback)) {
        return(finish_with_seq_fallback(eva_t, reason, info))
      }
      return(fail_ppro_suffix(
        t_start = eva_t,
        reason = reason,
        info = info,
        theta_value = theta_ll_est,
        beta_full_value = beta_full_est,
        it_value = result$it_num,
        residual_value = residual_est,
        M_value = result$M_out,
        first_n_sub_value = result$first_n_sub,
        repair_count_value = result$repair_count,
        final_n_sub_value = result$final_n_sub,
        init_mode_value = result$init_mode,
        init_trigger_value = result$init_trigger,
        H_value = H_seq
      ))
    }
  }
  
  # Return a named list with all results
  # This structure makes it easy to access individual components
  return(make_return(
    theta_value = theta_ll_est,
    beta_full_value = beta_full_est,
    it_value = result$it_num,
    residual_value = residual_est,
    M_value = result$M_out,
    first_n_sub_value = result$first_n_sub,
    repair_count_value = result$repair_count,
    final_n_sub_value = result$final_n_sub,
    init_mode_value = result$init_mode,
    init_trigger_value = result$init_trigger,
    H_value = H_seq,
    returned_backend = "ppro"
  ))
}

compare_tvcqr_ppro_strict <- function(fit_r, fit_f, tol_theta = 1e-8, tol_resid = 1e-6) {
  same_int_vec <- function(a, b) {
    a <- as.integer(a)
    b <- as.integer(b)
    length(a) == length(b) && all((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))
  }

  h_row_set_equal <- function(a, b) {
    a <- as.integer(a)
    b <- as.integer(b)
    if (length(a) != length(b) || anyNA(a) || anyNA(b)) {
      return(FALSE)
    }
    identical(sort(a), sort(b))
  }

  null_or_na <- function(x) {
    is.null(x) || length(x) == 0L || all(is.na(x))
  }

  theta_diff <- suppressWarnings(max(abs(fit_r$theta_ll_est - fit_f$theta_ll_est), na.rm = TRUE))
  if (!is.finite(theta_diff)) {
    theta_diff <- NA_real_
  }
  theta_ok <- isTRUE(!is.na(theta_diff) && theta_diff <= tol_theta)

  H_r <- fit_r$H_seq
  H_f <- fit_f$H_seq
  H_dim_ok <- identical(dim(H_r), dim(H_f))
  ordered_H_ok <- H_dim_ok && isTRUE(all.equal(H_r, H_f, tolerance = 0, check.attributes = FALSE))
  H_ordered_mismatches <- if (H_dim_ok) {
    which(rowSums(H_r != H_f, na.rm = FALSE) > 0L)
  } else {
    seq_len(max(nrow(H_r), nrow(H_f)))
  }
  H_row_set_mismatches <- if (H_dim_ok) {
    which(!vapply(seq_len(nrow(H_r)), function(i) h_row_set_equal(H_r[i, ], H_f[i, ]), logical(1)))
  } else {
    seq_len(max(nrow(H_r), nrow(H_f)))
  }
  H_row_set_ok <- H_dim_ok && length(H_row_set_mismatches) == 0L

  n_sub_ok <- same_int_vec(fit_r$n_sub, fit_f$n_sub)
  n_sub_mismatches <- if (length(fit_r$n_sub) == length(fit_f$n_sub)) {
    which(!((is.na(fit_r$n_sub) & is.na(fit_f$n_sub)) |
              (!is.na(fit_r$n_sub) & !is.na(fit_f$n_sub) & fit_r$n_sub == fit_f$n_sub)))
  } else {
    seq_len(max(length(fit_r$n_sub), length(fit_f$n_sub)))
  }

  it_num_ok <- same_int_vec(fit_r$it_num, fit_f$it_num)
  it_num_mismatches <- if (length(fit_r$it_num) == length(fit_f$it_num)) {
    which(!((is.na(fit_r$it_num) & is.na(fit_f$it_num)) |
              (!is.na(fit_r$it_num) & !is.na(fit_f$it_num) & fit_r$it_num == fit_f$it_num)))
  } else {
    seq_len(max(length(fit_r$it_num), length(fit_f$it_num)))
  }

  residual_compared <- !is.null(fit_r$residual_est) && !is.null(fit_f$residual_est)
  residual_max_diff <- NA_real_
  residual_ok <- NA
  if (residual_compared) {
    residual_max_diff <- suppressWarnings(max(abs(fit_r$residual_est - fit_f$residual_est), na.rm = TRUE))
    if (!is.finite(residual_max_diff)) {
      residual_max_diff <- NA_real_
    }
    residual_ok <- isTRUE(!is.na(residual_max_diff) && residual_max_diff <= tol_resid)
  }

  backend_ok <- identical(fit_r$returned_backend, "ppro") && identical(fit_f$returned_backend, "ppro")
  backend_ierr_ok <- is.null(fit_f$backend_ierr) || identical(as.integer(fit_f$backend_ierr), 0L)
  failure_reason_ok <- null_or_na(fit_r$failure_reason) && null_or_na(fit_f$failure_reason)

  blocking_pass <- theta_ok &&
    ordered_H_ok &&
    H_row_set_ok &&
    n_sub_ok &&
    backend_ok &&
    backend_ierr_ok &&
    failure_reason_ok &&
    (!isTRUE(residual_compared) || isTRUE(residual_ok))

  list(
    blocking_pass = blocking_pass,
    diagnostics_pass = blocking_pass && it_num_ok,
    theta_ok = theta_ok,
    max_theta_diff = theta_diff,
    ordered_H_ok = ordered_H_ok,
    ordered_H_mismatches = H_ordered_mismatches,
    H_row_set_ok = H_row_set_ok,
    H_row_set_mismatches = H_row_set_mismatches,
    n_sub_ok = n_sub_ok,
    n_sub_mismatches = n_sub_mismatches,
    it_num_ok = it_num_ok,
    it_num_mismatches = it_num_mismatches,
    residual_compared = residual_compared,
    residual_ok = residual_ok,
    residual_max_diff = residual_max_diff,
    backend_ok = backend_ok,
    backend_ierr_ok = backend_ierr_ok,
    failure_reason_ok = failure_reason_ok
  )
}
