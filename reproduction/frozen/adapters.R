# Model adapters: design, weights/thresholds and output interpretation only.
ssqr_experiment_dir <- function() {
  Sys.getenv("SSQR_EXPERIMENT_DIR", "experiments/hpc_u11_rep500")
}
ssqr_load <- function(mode="checked") {
  dyn.load(file.path(ssqr_experiment_dir(),"build",paste0("ssqr_",mode,".so")))
}
ssqr_call <- function(dll,A,y,grid,coordinate,kernel,h,tau,threshold,min_keep=1L,
                      provider_flags=as.integer(Sys.getenv("SSQR_PROVIDER_FLAGS","1"))) {
  n<-length(y); q<-ncol(A); ne<-length(grid)
  stopifnot(nrow(A)==n,length(coordinate)==n,all(is.finite(c(A,y,grid,coordinate,h,tau,threshold))),
            h>0,tau>0,tau<1,length(threshold)%in%c(1L,ne))
  out<-.Fortran("ssqr_kernel_path",PACKAGE=dll[["name"]],a=as.double(A),y=as.double(y),
    n=as.integer(n),q=as.integer(q),grid=as.double(grid),ne=as.integer(ne),coordinate=as.double(coordinate),
    kernel=as.integer(kernel),h=as.double(h),tau=as.double(tau),tol=1e-14,maxit=1000000L,
    threshold=rep_len(as.double(threshold),ne),min_keep=as.integer(min_keep),
    cache_flags=as.integer(Sys.getenv("SSQR_CACHE_FLAGS","27")),
    provider_flags=as.integer(provider_flags),
    beta=double(ne*q),hseq=integer(ne*q),diagnostics=integer(ne*18L),ierr=0L,failed_eval=0L)
  diagnostics<-structure(matrix(out$diagnostics,ne,18L),dimnames=list(NULL,c(
    "n_active","first_tableau_rows","final_tableau_rows","iterations","repairs","init_mode",
    "init_trigger","full_recovery","independent","certificate_hits","residual_rows",
    "first_full_m_recovery","initial_threshold_hits","threshold_expansion_steps",
    "effective_threshold_hits","basis_forced_rows","first_aggregate_rows","first_retained_size")))
  diagnostics[1L,13:18]<-NA_integer_
  list(beta=matrix(out$beta,ne,q),H=matrix(out$hseq,ne,q),
       diagnostics=diagnostics,
       ierr=out$ierr,failed_eval=out$failed_eval,seq_fallback=FALSE)
}
ssqr_llqr <- function(dll,x,y,tau,case,h,z=sort(x),Mm.factor=.1,provider_flags=NULL) {
  n<-length(y); ordered<-!is.unsorted(z)
  if(ordered){grid<-z}else{ord<-order(z);grid<-z[ord]}
  mm<-if(case==1L)log(log(n))/sqrt(log(n)) else sqrt(log(n))*n^(-.4)
  if(is.null(provider_flags))provider_flags<-as.integer(Sys.getenv("SSQR_PROVIDER_FLAGS","1"))
  threshold<-Mm.factor*mm*log(log(n))
  out<-ssqr_call(dll,cbind(1,x),y,grid,x,if(case==1L)1L else 2L,h,tau,
                 threshold,provider_flags=provider_flags)
  if(!ordered){
    undo<-order(ord);out$beta<-out$beta[undo,,drop=FALSE];out$H<-out$H[undo,,drop=FALSE]
    out$diagnostics<-out$diagnostics[undo,,drop=FALSE]
  }
  out$estimate<-out$beta[,1L]+z*out$beta[,2L];out$derivative<-out$beta[,2L]
  out$threshold_initial<-threshold
  out
}
ssqr_tvcqr <- function(dll,x,y,tau,h=length(y)^(-.2),Mm.factor=1e-5,provider_flags=NULL) {
  n<-length(y); time<-seq_len(n)/n; X<-cbind(1,x); p<-ncol(X)
  mm<-log(n)^4*h^2*max(sqrt(rowSums(x^2)))
  if(is.null(provider_flags))provider_flags<-as.integer(Sys.getenv("SSQR_PROVIDER_FLAGS","1"))
  threshold<-Mm.factor*mm*log(log(n))
  out<-ssqr_call(dll,cbind(X,X*time),y,time,time,2L,h,tau,threshold,
                 provider_flags=provider_flags)
  out$estimate<-out$beta[,seq_len(p),drop=FALSE]+out$beta[,p+seq_len(p),drop=FALSE]*time
  out$derivative<-out$beta[,p+seq_len(p),drop=FALSE]
  out$threshold_initial<-threshold
  out
}
