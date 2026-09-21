args<-commandArgs(TRUE);build<-if(length(args))args[1] else 'tmp/u11-build'
source('R/llqr_functions.R');source('reproduction/frozen/design.R')
source('reproduction/frozen/adapters.R');source('reproduction/frozen/model_methods.R')
core<-dyn.load(file.path(build,'ssqr_optimized.so'));lean<-dyn.load(file.path(build,'llqr_seq_lean_sortskip.so'))
Sys.setenv(SSQR_CACHE_FLAGS=27,SSQR_PROVIDER_FLAGS=1)
rows<-list()
for(case in 1:2)for(n in c(100,200))for(tau in c(.2,.5,.8))for(seed in 2026:2027){
 d<-if(case==1)generate_logistic_case2(n,seed) else generate_data(n,2,seed)
 x<-d$x;y<-d$y;z<-if(case==1)d$z else sort(x);h<-n^(-.2);A<-cbind(1,x)
 f<-ssqr_llqr(core,x,y,tau,2,h,z);s<-run_llqr_lean(lean,x,y,z,tau,h,2)
 b<-llqr_local_fit(x,y,tau,z,h,2);H<-matrix(s$H_mat,ncol=2)
 stopifnot(f$ierr==0,f$failed_eval==0,!f$seq_fallback,all(vapply(seq_along(z),function(j)setequal(f$H[j,],H[j,]),logical(1))))
 gap<-viol<-numeric(length(z))
 for(j in seq_along(z)){
  w<-llqr_kernel_weights((z[j]-x)/h,2);r<-y-drop(A%*%f$beta[j,]);bb<-c(b$ll_est[j]-z[j]*b$d_ll_est[j],b$d_ll_est[j]);rr<-y-drop(A%*%bb)
  ob<-sum(w*rr*(tau-(rr<0)));gap[j]<-abs(sum(w*r*(tau-(r<0)))-ob)/(1+abs(ob))
  u<-w*ifelse(r>=0,tau,tau-1);u[H[j,]]<-0;u[H[j,]]<-solve(t(A[H[j,],,drop=FALSE]),-drop(crossprod(A,u)))
  viol[j]<-max(0,abs(drop(crossprod(A,u))),(tau-1)*w-u,u-tau*w)
 }
 err<-max(abs(f$estimate-b$ll_est));stopifnot(err<1e-8,max(gap)<1e-8,max(viol)<1e-6)
 rows[[length(rows)+1]]<-data.frame(case,n,tau,seed,points=length(z),estimate_error=err,objective_gap=max(gap),KKT=max(viol),full_active_recovery=sum(f$diagnostics[-1,'full_recovery'])+f$diagnostics[1,'first_full_m_recovery'],seq_fallback=0,H_mismatches=0)
}
out<-do.call(rbind,rows)
if(length(args)>1)write.csv(out,args[2],row.names=FALSE)
cat('PASS datasets=',nrow(out),' points=',sum(out$points),' recovery=',sum(out$full_active_recovery),' H_mismatches=0 seq_fallback=0\n',sep='')
print(sapply(out[c('estimate_error','objective_gap','KKT')],max))
