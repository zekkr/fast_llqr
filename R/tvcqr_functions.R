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
      sh.bad <- (r < 0) & sh
      sl.bad <- (r > 0) & sl
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
                           store_residual = FALSE) {
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
  # test_sl_sh <- rep(0, m)
  n_sub[1] <- m
  min_subsample_size <- max(5 * (nvar + 1), ceiling(0.2 * m))
  residual_tol <- 1e-8
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
      warning('Not converge')
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

  finish_with_seq_fallback <- function(t_start, current_M) {
    seq_fit <- tvcqr_seq(
      x = x,
      y = y,
      tau = tau,
      h = h,
      h.factor = h.factor,
      tol = tol,
      maxit = maxit,
      bland = bland
    )

    theta_ll_est[t_start:m, ] <<- seq_fit$theta_ll_est[t_start:m, , drop = FALSE]
    it_num[t_start:m] <<- NA_real_
    if (store_residual) {
      residual_est[t_start:m, ] <<- seq_fit$residual_est[t_start:m, , drop = FALSE]
    }
    n_sub[t_start:m] <<- m
    H_seq[t_start:m, ] <<- seq_fit$H_seq[t_start:m, , drop = FALSE]

    return(list(
      theta_ll_est = theta_ll_est,
      it_num = it_num,
      residual_est = residual_est,
      M = current_M,
      n_sub = n_sub,
      H_seq = H_seq
    ))
  }
  
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
      sh.bad <- (r < 0) & sh
      sl.bad <- (r > 0) & sl
      bad.signs <- sum(sh.bad | sl.bad)
      H_candidate <- r1 - 2 - 2 * nvar
      H_candidate <- idx_not_jl_or_jh[H_candidate]
      accept_subsample <- (bad.signs == 0) && is_valid_tvcqr_H(H_candidate, r)
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
          if ((ms >= m) && !any(sl) && !any(sh) && !force_full_sample) {
            return(finish_with_seq_fallback(eva_t, M))
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
      warning('Not converge')
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
  
  return(list(theta_ll_est = theta_ll_est, it_num = it_num, residual_est = residual_est, M = M, n_sub = n_sub, H_seq = H_seq))
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
                                           store_residual = TRUE) {
  
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
  it_num <- integer(m)
  residual_est <- if (isTRUE(store_residual)) matrix(0.0, nrow = m, ncol = m) else 0.0
  M_out <- 0.0
  n_sub <- integer(m)
  H_seq <- matrix(0L, nrow = m, ncol = 2 * (nvar + 1))
  
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
                     # Output arguments - pre-allocated arrays
                     theta_ll_est = as.double(theta_ll_est),
                     it_num = as.integer(it_num),
                     residual_est = as.double(residual_est),
                     M_out = as.double(M_out),
                     n_sub = as.integer(n_sub),
                     H_seq = as.integer(H_seq),
                     ierr = as.integer(0),
                     # Don't duplicate arrays (more efficient)
                     DUP = FALSE)

  if (!identical(as.integer(result$ierr), 0L)) {
    fallback <- tvcqr_seq_fortran_wrapper(
      x = x,
      y = y,
      tau = tau,
      h = if (h_value > 0) h_value else NULL,
      tol = tol,
      maxit = maxit,
      bland = bland
    )
    fallback$M <- NA_real_
    fallback$n_sub <- rep(m, m)
    if (!isTRUE(store_residual)) {
      fallback$residual_est <- NULL
    }
    return(fallback)
  }
  
  # Reshape the flattened arrays back to matrices
  # Fortran stores matrices in column-major order, same as R
  theta_ll_est <- matrix(result$theta_ll_est, nrow = m, ncol = nvar + 1, byrow = FALSE)
  residual_est <- if (isTRUE(store_residual)) {
    matrix(result$residual_est, nrow = m, ncol = m, byrow = FALSE)
  } else {
    NULL
  }
  H_seq <- matrix(result$H_seq, nrow = m, ncol = 2 * (nvar + 1), byrow = FALSE)
  
  # Return a named list with all results
  # This structure makes it easy to access individual components
  return(list(
    theta_ll_est = theta_ll_est,    # Time-varying coefficient estimates
    it_num = result$it_num,          # Number of iterations at each time point
    residual_est = residual_est,     # Residuals at each time point
    M = result$M_out,                # Final threshold value used
    n_sub = result$n_sub,            # Subsample size at each time point
    H_seq = H_seq,                   # Interpolation indices at each time point
    h = result$h                     # Bandwidth used (useful if it was calculated)
  ))
}
