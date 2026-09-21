#!/usr/bin/env Rscript
experiment_dir <- Sys.getenv("SSQR_EXPERIMENT_DIR", "experiments/hpc_u11_rep500")
source(file.path(experiment_dir,"design.R"))
source(file.path(experiment_dir,"metrics.R"))
source(file.path(experiment_dir,"model_methods.R"))
source(file.path(experiment_dir,"adapters.R"))
source("R/llqr_functions.R")
check <- function(x,label) {stopifnot(isTRUE(x));cat("PASS",label,"\n")}
g <- llqr_interior_grid(c(0, .05, .95, 1))
check(identical(g$z,.5) && g$placeholder && g$n_interior == 0L, "empty grid uses theorem placeholder")
g <- llqr_interior_grid(c(0,.1,1))
check(identical(g$z,.1) && !g$placeholder && g$n_interior == 1L, "one point and included boundary")
check(identical(llqr_interior_grid(c(.9,.1,.95))$z,c(.1,.9)), "grid sorts both included endpoints")
h <- matrix(c(7L,8L),1,2)
check(valid_h_path(h,10L,2L,1L), "H observation IDs may exceed number of evaluations")
check(!valid_h_path(matrix(c(7L,11L),1,2),10L,2L,1L), "out-of-sample H rejected")
check(audit_h_path("unified_u11",TRUE,h,TRUE,h[,2:1,drop=FALSE],10L,2L,1L)$accepted, "single-point H row-set match")
# Read diagnostic summarizer definitions without executing the driver.
expr <- as.list(parse(file.path(experiment_dir,"driver_array.R")))
for (name in c("safe_max_int","safe_sum_int","empty_compact_diagnostics","u11_diagnostics")) {
 e <- Filter(function(e)is.call(e)&&identical(e[[1]],as.name("<-"))&&identical(e[[2]],as.name(name)),expr)
 stopifnot(length(e)==1L);eval(e[[1]],.GlobalEnv)
}
dll <- ssqr_load("checked"); lean <- load_lean_kernels(experiment_dir)$llqr
Sys.setenv(SSQR_CACHE_FLAGS="27",SSQR_PROVIDER_FLAGS="1")
dat <- generate_logistic_case2(200L,2026L)
check(length(dat$x)==200L && length(dat$z)<200L && all(dat$z>=.1 & dat$z<=.9), "full sample retained while grid restricted")
for (z in list(.5,dat$z)) {
 n <- length(dat$x);ne <- length(z);bw <- n^(-.2)
 a <- ssqr_llqr(dll,dat$x,dat$y,.5,2L,bw,z)
 b <- run_llqr_lean(lean,dat$x,dat$y,z,.5,bw,2L)
 ref <- llqr_local_fit(dat$x,dat$y,.5,z,bw,2L)
 check(a$ierr==0L && length(a$estimate)==ne && length(b$ll_est)==ne &&
       audit_h_path("unified_u11",TRUE,a$H,TRUE,b$H_mat,n,2L,ne)$accepted,
       paste("real fit dimensions and H on",ne,"points"))
 check(paper_discrepancy(a$estimate,ref$ll_est,ne)$value<1e-8 &&
       paper_discrepancy(b$ll_est,ref$ll_est,ne)$value<1e-8,"direct numerical agreement")
 ds <- u11_diagnostics(a$diagnostics)
 check(ds$first_pass_total==ne-1L && ds$retained_decomposition_ok, "diagnostic transition count")
 if (ne==1L) check(is.na(ds$first_pass_rate) && ds$repair_total==0L && ds$first_retained_size_max==0L,"single-point diagnostics are defined")
}
