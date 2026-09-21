options(warn=1)
source('R/llqr_functions.R')
old <- new.env(); sys.source('reproduction/frozen/R/llqr_functions.R', old)
# Observe the public rq boundary rather than infer the full-row contract from output.
.rq_seen <- list()
trace('rq', where=asNamespace('quantreg'), print=FALSE, tracer=quote({
  .GlobalEnv$.rq_seen[[length(.GlobalEnv$.rq_seen)+1L]] <- list(
    n=length(eval(formula[[2]], environment(formula))), w=weights,
    method=method, dots=list(...))
}))
rows <- list()
for (case in 1:2) for(tau in c(.2,.5,.8)) for(seed in 2026:2027) {
  d <- generate_data(200,case,seed);x<-d$x;y<-d$y
  if(seed==2027) x[2] <- x[1] # duplicate predictor
  z<-c(.7,.3,.7,.5);h<-if(case==1) .8 else .22
  before<-length(.rq_seen)
  a<-llqr_local_fit(x,y,tau,z,h,case); b<-old$llqr_local_fit(x,y,tau,z,h,case)
  calls<-.rq_seen[seq.int(before+1L,before+length(z))]
  stopifnot(length(.rq_seen)-before==length(z),all(vapply(calls,function(v)v$n==length(x)&&length(v$w)==length(x)&&v$method=='br'&&!length(v$dots),logical(1))))
  gap<-kkt<-numeric(length(z))
  A<-cbind(1,x)
  for(j in seq_along(z)) {
    w<-llqr_kernel_weights((z[j]-x)/h,case)
    stopifnot(identical(calls[[j]]$w,w))
    ba<-c(a$ll_est[j]-z[j]*a$d_ll_est[j],a$d_ll_est[j]);bb<-c(b$ll_est[j]-z[j]*b$d_ll_est[j],b$d_ll_est[j])
    r<-y-drop(A%*%ba);rr<-y-drop(A%*%bb)
    gap[j]<-abs(sum(w*r*(tau-(r<0)))-sum(w*rr*(tau-(rr<0))))/(1+sum(w*rr*(tau-(rr<0))))
    H<-head(order(ifelse(w>0,abs(r),Inf)),2)
    u<-w*ifelse(r>=0,tau,tau-1);u[H]<-0
    u[H]<-solve(t(A[H,,drop=FALSE]),-drop(crossprod(A,u)))
    kkt[j]<-max(0,abs(drop(crossprod(A,u))),(tau-1)*w-u,u-tau*w)
  }
  stopifnot(max(gap)<1e-8,max(kkt)<1e-6,max(abs(a$ll_est-b$ll_est))<1e-8)
  rows[[length(rows)+1]]<-data.frame(case,tau,seed,estimate=max(abs(a$ll_est-b$ll_est)),objective=max(gap),KKT=max(kkt))
}
# Exact support boundaries and insufficient support are passed through rq unchanged.
x<-c(-1,-.5,0,.5,1,2);y<-c(2,1,0,2,3,4)
before<-length(.rq_seen); a<-llqr_local_fit(x,y,.5,z=0,h=1,case=2)
stopifnot(identical(tail(.rq_seen,1)[[1]]$w,c(0,.5625,.75,.5625,0,0)))
for(z in c(0,10)) {
  e<-try(llqr_local_fit(c(0,1,2),c(1,3,2),z=z,h=.1,case=2),silent=TRUE)
  ref<-try(quantreg::rq(c(1,3,2)~c(0,1,2),weights=llqr_kernel_weights((z-c(0,1,2))/.1,2)),silent=TRUE)
  stopifnot(inherits(e,'try-error')==inherits(ref,'try-error'))
}
untrace('rq',where=asNamespace('quantreg'))
x<-do.call(rbind,rows);print(x)
cat('PASS 12 datasets / 48 comparison points; all rq calls retain full rows; boundary/degenerate cases checked; fallback N/A\n')
