#!/usr/bin/env Rscript
# Independent small-grid release gate. No changes to stable seq or test fixtures.
root <- normalizePath(getwd()); args <- commandArgs(TRUE)
outdir <- if(length(args)) args[1L] else "tmp/u11-validation"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
source("R/u11_functions.R"); research <- load_u11()
source("reproduction/frozen/adapters.R")
source("reproduction/frozen/model_methods.R")
source("reproduction/frozen/R/llqr_functions.R")
source("reproduction/frozen/R/tvcqr_functions.R")
source("reproduction/frozen/design.R")
library(fastllqr)
dll <- research$dll
Sys.setenv(SSQR_CACHE_FLAGS=27,SSQR_PROVIDER_FLAGS=1)
lean <- list(llqr=dyn.load("tmp/u11-build/llqr_seq_lean_sortskip.so"),
             tvcqr=dyn.load("tmp/u11-build/tvcqr_seq_lean_nohistory.so"))
h_equal <- function(a,b) identical(dim(a),dim(b)) && all(vapply(seq_len(nrow(a)),function(i) setequal(a[i,],b[i,]),logical(1)))
kkt <- function(A,y,b,H,w,tau) {
  r <- y-drop(A%*%b); u <- w*ifelse(r>=0,tau,tau-1); u[H] <- 0
  u[H] <- solve(t(A[H,,drop=FALSE]),-drop(crossprod(A,u)))
  max(0,abs(drop(crossprod(A,u))), (tau-1)*w-u,u-tau*w)
}
rows <- list()
for (case in 1:5) for(n in c(100L,200L)) for(tau in c(.2,.5,.8)) for(seed in 2026:2027) {
  ll <- case %in% c(1,2,5); kernel <- if(case==5) "gaussian" else "epanechnikov"
  d <- if(case==1) generate_logistic_case2(n,seed) else if(ll) generate_data(n,if(case==5) 1 else 2,seed) else generate_ts(n,case-2,seed)
  x<-d$x;y<-d$y; z<-if(case==1)d$z else if(ll)sort(x) else seq_len(n)/n
  h<-if(case==5) llqr_default_bandwidth(x,y,tau,case=1) else n^(-.2)
  if(ll) {
    p<-fastllqr::llqr_seq_ppro(x,y,tau,z,h,kernel,diagnostics=TRUE)
    r<-research$llqr_seq_ppro(x,y,tau,z,h,kernel,diagnostics=TRUE)
    f<-ssqr_llqr(dll,x,y,tau,if(case==5)1L else 2L,h,z)
    s<-run_llqr_lean(lean$llqr,x,y,z,tau,h,if(case==5)1L else 2L)
    Hs<-matrix(s$H_mat,ncol=2);A<-cbind(1,x)
    beta<-cbind(p$ll_est-z*p$d_ll_est,p$d_ll_est)
    direct<-llqr_local_fit(x,y,tau,z,h,case=if(case==5)1 else 2)$ll_est
    est<-p$ll_est
  } else {
    p<-fastllqr::tvcqr_seq_ppro(x,y,tau,h,diagnostics=TRUE)
    r<-research$tvcqr_seq_ppro(x,y,tau,h,diagnostics=TRUE)
    f<-ssqr_tvcqr(dll,x,y,tau,h)
    s<-run_tvcqr_lean(lean$tvcqr,x,y,tau)
    Hs<-matrix(s$H_mat,ncol=2*(ncol(x)+1));X<-cbind(1,x);A<-cbind(X,X*z)
    beta<-p$beta_full_est;direct<-tvc_rq(x,y,tau,h)$theta_ll_est;est<-p$theta_ll_est
  }
  stopifnot(f$ierr==0,identical(p$diagnostics$path,f$diagnostics),h_equal(p$diagnostics$H_seq,f$H),h_equal(p$diagnostics$H_seq,Hs),
    h_equal(p$diagnostics$H_seq,r$diagnostics$H_seq),identical(p,r),
    max(abs(est-f$estimate))<1e-8,max(abs(est-direct))<1e-8)
  gaps<-viol<-numeric(length(z))
  for(j in seq_along(z)) {
    u<-(z[j]-(if(ll)x else z))/h
    w<-if(kernel=="gaussian") dnorm(u) else .75*pmax(1-u*u,0)
    rr<-y-drop(A%*%beta[j,]);obj<-sum(w*rr*(tau-(rr<0)))
    b<-quantreg::rq.wfit(A[w>0,,drop=FALSE],y[w>0],tau,w[w>0],method="br")$coef
    dr<-y-drop(A%*%b);objref<-sum(w*dr*(tau-(dr<0)))
    gaps[j]<-abs(obj-objref)/(1+abs(objref))
    viol[j]<-kkt(A,y,beta[j,],p$diagnostics$H_seq[j,],w,tau)
  }
  stopifnot(max(gaps)<1e-8,max(viol)<1e-6,!attr(p,"solver_info")$fallback_triggered)
  rows[[length(rows)+1L]]<-data.frame(case,n,tau,seed,points=length(z),H_match=TRUE,
    max_estimate_error=max(abs(est-direct)),max_objective_gap=max(gaps),max_KKT_violation=max(viol),
    full_active_recovery=attr(p,"solver_info")$full_active_recovery_count,seq_fallback=FALSE)
}
res<-do.call(rbind,rows);write.csv(res,file.path(outdir,"validation.csv"),row.names=FALSE)
cat("PASS datasets=",nrow(res)," points=",sum(res$points)," H_mismatches=0 seq_fallback=0 full_active_recovery=",sum(res$full_active_recovery),"\n",sep="")
print(sapply(res[c("max_estimate_error","max_objective_gap","max_KKT_violation")],max))
