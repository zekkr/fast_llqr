#!/usr/bin/env Rscript
source("reproduction/build.R")
build_u11("tmp/u11-paper-profile","paper-hpc",getwd())
file.copy("tmp/u11-paper-profile/ssqr_optimized.so","tmp/u11-paper-profile/ssqr_paper.so",overwrite=TRUE)
paper<-dyn.load("tmp/u11-paper-profile/ssqr_paper.so")
source("reproduction/frozen/adapters.R");source("reproduction/frozen/design.R")
source("reproduction/frozen/R/llqr_functions.R");source("reproduction/frozen/R/tvcqr_functions.R")
h_equal<-function(a,b)all(vapply(seq_len(nrow(a)),function(i)setequal(a[i,],b[i,]),logical(1)))
rows<-list()
for(case in 1:4){
 n<-1000L;h<-n^(-.2);tau<-.5
 d<-if(case==1)generate_logistic_case2(n,2026) else if(case==2)generate_data(n,2,2026) else generate_ts(n,case-2,2026)
 z<-if(case==1)d$z else if(case==2)sort(d$x) else seq_len(n)/n
 package_fit<-function()if(case<=2)fastllqr::llqr_seq_ppro(d$x,d$y,tau,z,h,"epanechnikov",diagnostics=TRUE) else fastllqr::tvcqr_seq_ppro(d$x,d$y,tau,h,diagnostics=TRUE)
 frozen_fit<-function()if(case<=2)ssqr_llqr(paper,d$x,d$y,tau,2L,h,z) else ssqr_tvcqr(paper,d$x,d$y,tau,h)
 a<-package_fit();b<-frozen_fit();stopifnot(b$ierr==0L,h_equal(a$diagnostics$H_seq,b$H))
 est<-if(case<=2)a$ll_est else a$theta_ll_est
 stopifnot(max(abs(est-b$estimate))<1e-8)
 funcs<-list(package_portable=package_fit,archived_core_paper_flags=frozen_fit)
 for(rep in 1:6)for(method in if(rep%%2)names(funcs) else rev(names(funcs))){
   gc(FALSE);t<-proc.time()[["elapsed"]];ans<-funcs[[method]]();elapsed<-proc.time()[["elapsed"]]-t
   rows[[length(rows)+1L]]<-data.frame(case,n,rep,method,seconds=elapsed,H_match=TRUE,max_estimate_difference=max(abs(est-b$estimate)))
 }
}
write.csv(do.call(rbind,rows),"reproduction/reports/profile_benchmark.csv",row.names=FALSE)
print(aggregate(seconds~case+method,do.call(rbind,rows),mean))
