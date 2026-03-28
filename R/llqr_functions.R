# ============================================================================ #
# Generates simulation data for one-dim Local linear quantile regression
# ============================================================================ #
generate_data <- function(n, case = 1, seed = NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if(case == 1){
    x <- rnorm(n)
  } else if(case == 2){
    x <- runif(n)
  } else {
    stop("Invalid case specification. Use 1 or 2.")
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
# Sequential algorithm for one-dim Local linear quantile regression
# ============================================================================ #
llqr_seq <- function(x, y, tau = 0.5, z = NULL, h = NULL, tol = 1e-14, maxit = 1e6, 
                         bland = F, track_order = F){
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
  
  if (is.null(h)) {
    # choosing bandwidth using rule of thumb in Yu and Jones 1998.
    red_dim <- floor(0.2 * m) # Rounding of Numbers
    index_y <- order(y)[red_dim:(m - red_dim)] # increasing = TRUE
    h <- KernSmooth::dpill(x[index_y, ], y[index_y]) # bugs exist when multivariate case
    h <- 1.25 * h * (tau * (1 - tau)/(dnorm(qnorm(tau)))^2)^0.2
    if (h == "NaN") {
      h <- 1.25 * max(m^(-1/(nvar + 4)), min(2, sd(y))*m^(-1/(nvar + 4)))
    }
  }
  
  
  eva_z <- z[1] - x
  w <- dnorm(eva_z/h)
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
      w <- dnorm(eva_z/h)
      
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
                      track_order = F, case = 1, pmethod = NULL){
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
  
  if (is.null(h)) {
    # choosing bandwidth using rule of thumb in Yu and Jones 1998.
    red_dim <- floor(0.2 * m) # Rounding of Numbers
    index_y <- order(y)[red_dim:(m - red_dim)] # increasing = TRUE
    h <- KernSmooth::dpill(x[index_y, ], y[index_y]) 
    h <- 1.25 * h * (tau * (1 - tau)/(dnorm(qnorm(tau)))^2)^0.2
    if (h == "NaN") {
      h <- 1.25 * max(m^(-1/(nvar + 4)), min(2, sd(y))*m^(-1/(nvar + 4)))
    }
  }
  
  x.norms <- apply(x, 1, function(row) sqrt(sum(row^2)))
  if (case == 1){
    mm <- log(log(m))/sqrt(log(m))
  } else if (case == 2) {
    mm <-  log(m)^(1/2) * m^{-2/5}
  }
  
  # initialize the output
  ll_est <- rep(0,rounds)
  d_ll_est <- rep(0,rounds)
  residual_est <- matrix(0, nrow = m, ncol = m)
  
  # Initial estimation at the first evaluation point
  eva_z <- z[1] - x
  xx <- cbind(matrix(1, nrow = m, ncol = 1), eva_z) # n*(p+1)
  w <- dnorm(eva_z/h)
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
    w <- dnorm(eva_z/h)
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
                              Mm.factor = 1, bland = F, track_order = F, 
                              case = 1, min_subsample_size = NULL,
                              store_residual = FALSE){
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
  
  if (is.null(h)) {
    # choosing bandwidth using rule of thumb in Yu and Jones 1998.
    red_dim <- floor(0.2 * m) # Rounding of Numbers
    index_y <- order(y)[red_dim:(m - red_dim)] # increasing = TRUE
    h <- KernSmooth::dpill(x[index_y, ], y[index_y]) 
    h <- 1.25 * h * (tau * (1 - tau)/(dnorm(qnorm(tau)))^2)^0.2
    if (h == "NaN") {
      h <- 1.25 * max(m^(-1/(nvar + 4)), min(2, sd(y))*m^(-1/(nvar + 4)))
    }
  }
  
  if (is.null(min_subsample_size)){
    min_subsample_size <- max(5 * (nvar + 1), ceiling(0.2 * m))
  }
  residual_tol <- 1e-8
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
  if (case == 1){
    mm <- log(log(m))/sqrt(log(m))
  } else if (case == 2) {
    mm <-  log(m)^(1/2) * m^{-2/5}
  }
  
  eva_z <- z[1] - x
  w <- dnorm(eva_z/h)
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
  residual_est <- if (store_residual) matrix(0, nrow = m, ncol = m) else NULL
  r_prev <- rep(0, m)
  n_sub <- rep(0, m)
  H_seq <- matrix(0, nrow = m, ncol = nvar+1)
  n_sub[1] <- m
  
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
      warning('Not converge')
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
  }
  
  # we use a big (n+2) rows gammax to store the gammaxs
  gammaxs.temp <- matrix(NA, nrow = m + 2, ncol = nvar + 1)
  bs.temp <- rep(NA, m + 2)

  finish_with_seq_fallback <- function(rd_start, current_M) {
    seq_fit <- llqr_seq(
      x = x,
      y = y,
      tau = tau,
      z = z[rd_start:rounds, , drop = FALSE],
      h = h,
      tol = tol,
      maxit = maxit,
      track_order = FALSE
    )

    ll_est[rd_start:rounds] <<- seq_fit$ll_est
    d_ll_est[rd_start:rounds] <<- seq_fit$d_ll_est
    it_num[rd_start:rounds] <<- NA_real_
    if (store_residual) {
      residual_est[rd_start:rounds, ] <<- seq_fit$residual_est
    }
    n_sub[rd_start:rounds] <<- m
    H_seq[rd_start:rounds, ] <<- NA_real_

    if (track_order) {
      ll_est <<- ll_est[order(original_order)]
      d_ll_est <<- d_ll_est[order(original_order)]
    }

    return(list(
      ll_est = ll_est,
      d_ll_est = d_ll_est,
      it_num = it_num,
      residual_est = residual_est,
      h = h,
      M = current_M,
      n_sub = n_sub,
      H_seq = H_seq
    ))
  }
  
  for (rd in 2:rounds){
    # print(rd)
    not_optimal <- TRUE
    not_new_sl_sh <- TRUE
    force_full_sample <- FALSE
    empty_pivot_count <- 0L
    eva_z <- z[rd] - x
    w <- dnorm(eva_z/h)
    mmm <- mm
    
    j <- 0
    while (not_optimal) {
      
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
        P <- c(rep(1,length(u.in.IBs)),rep(-1,length(v.in.IBs)),1)
        freevarrow <- c(rep(TRUE,nvar + 1), rep(FALSE,length(u.in.IBs)), 
                        rep(FALSE,length(v.in.IBs)),TRUE, TRUE) # 1~2, v_L not into the nonbasic set
      } else if (any(sh)) {
        Hbar <- c(idpos,idneg,ms)
        IBs <- c(1:(nvar + 1), u.in.IBs, v.in.IBs, nvar+1+ms) # length(IBs) = ms
        # P <- diag(c(rep(1,length(u.in.IBs)),rep(-1,length(v.in.IBs)),1))
        P <- c(rep(1,length(u.in.IBs)),rep(-1,length(v.in.IBs)),-1)
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
      accept_subsample <- (bad.signs == 0) && is_valid_llqr_H(H_candidate, r)
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
          # Once screening has already collapsed to the full sample, repeated
          # invalid_H retries add no new information. Fall back to the stable
          # sequential solver for the remaining evaluation points instead.
          if ((ms >= m) && !any(sl) && !any(sh) && !force_full_sample) {
            return(finish_with_seq_fallback(rd, M))
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
      warning('Not converge')
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
  }
  
  # Reorder results back to the original order if track_order is TRUE
  if (track_order) {
    ll_est <- ll_est[order(original_order)]
    d_ll_est <- d_ll_est[order(original_order)]
    # it_num <- it_num[order(original_order)]
  }
  
  return(list(ll_est = ll_est, d_ll_est = d_ll_est, it_num = it_num, residual_est = residual_est, h = h, M = M, n_sub = n_sub, H_seq = H_seq))
}

# ============================================================================ #
# Sequential algorithm for one-dim LLQR (Fortran version)
# ============================================================================ #
llqr_seq_fortran_wrapper <- function(x, y, tau = 0.5, z = NULL, h = NULL, tol = 1e-14,
                                     maxit = 1e6, bland = FALSE) {
  
  # Convert to vectors
  x <- as.vector(x)
  y <- as.vector(y)
  if (is.null(z)) {
    z <- x
  }
  z <- as.vector(z)
  # Keep output row order consistent with R llqr_seq/llqr_seq_ppro (sorted-z order).
  z <- z[order(z)]
  
  # Get dimensions
  m <- length(x)
  nvar <- 1  # Always 1 for univariate LLQR
  rounds <- length(z)
  
  # Handle bandwidth - calculate same way as R function
  if (is.null(h)) {
    # Using rule of thumb from Yu and Jones 1998 (same as llqr_tau_seq)
    red_dim <- floor(0.2 * m)
    index_y <- order(y)[red_dim:(m - red_dim)]
    h <- KernSmooth::dpill(x[index_y], y[index_y])
    h <- 1.25 * h * (tau * (1 - tau)/(dnorm(qnorm(tau)))^2)^0.2
    if (is.nan(h)) {
      h <- 1.25 * max(m^(-1/(nvar + 4)), min(2, sd(y))*m^(-1/(nvar + 4)))
    }
  }
  
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
                     bland_int = as.integer(bland),
                     ll_est = as.double(ll_est),
                     d_ll_est = as.double(d_ll_est),
                     it_num = as.integer(it_num),
                     residual_est = as.double(residual_est),
                     H_mat = as.integer(H_mat))
  
  # Return results matching R function structure
  list(
    ll_est = result$ll_est,
    d_ll_est = result$d_ll_est,
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
                                      Mm.factor = 1e-3, case = 1, tol = 1e-14,
                                      maxit = 1e6, bland = TRUE) {
  
  # Auto-load library if not already loaded
  # if (!.llqr_ppro_loaded) {
  #   load_llqr_ppro()
  # }
  if (is.null(z)) {
    z <- x
  }
  z <- as.vector(z)
  # Keep output row order consistent with R llqr_seq/llqr_seq_ppro (sorted-z order).
  z <- z[order(z)]
  if (!(case %in% c(1, 2))) {
    stop("Invalid case specification. Use 1 (normal) or 2 (uniform).")
  }
  
  # Setup
  m <- length(y)
  nvar <- 1  # univariate (can be extended for multivariate)
  rounds <- length(z)
  
  # Handle bandwidth - calculate same way as R function
  if (is.null(h)) {
    # Using rule of thumb from Yu and Jones 1998 (same as llqr_tau_seq)
    red_dim <- floor(0.2 * m)
    index_y <- order(y)[red_dim:(m - red_dim)]
    h <- KernSmooth::dpill(x[index_y], y[index_y])
    h <- 1.25 * h * (tau * (1 - tau)/(dnorm(qnorm(tau)))^2)^0.2
    if (is.nan(h)) {
      h <- 1.25 * max(m^(-1/(nvar + 4)), min(2, sd(y))*m^(-1/(nvar + 4)))
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
                     ll_est = double(rounds),
                     d_ll_est = double(rounds),
                     it_num = integer(rounds),
                     residual_est = matrix(0.0, nrow = rounds, ncol = m),
                     H_mat = matrix(0L, nrow = rounds, ncol = nvar + 1),
                     ierr = integer(1))

  if (!identical(as.integer(result$ierr), 0L)) {
    fallback <- llqr_seq_fortran_wrapper(
      x = x,
      y = y,
      tau = tau,
      z = z,
      h = h,
      tol = tol,
      maxit = maxit,
      bland = bland
    )
    fallback$M <- NA_real_
    fallback$n_sub <- rep(length(y), length(fallback$ll_est))
    return(fallback)
  }
  
  # Return results matching R's output format
  return(list(
    ll_est = result$ll_est,
    d_ll_est = result$d_ll_est,
    it_num = result$it_num,
    residual_est = result$residual_est,
    H_seq = result$H_mat  # Note: Named H_seq to match R output
  ))
}
