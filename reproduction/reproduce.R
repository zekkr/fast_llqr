#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
option<-function(key,default=NULL){i<-match(paste0("--",key),args);if(is.na(i))default else {if(i==length(args))stop("Missing value: ",key);args[i+1L]}}
if(system2(Sys.getenv("PYTHON","python3"),"reproduction/check_frozen.py")!=0L)stop("Frozen input check failed")
mode<-option("mode","smoke");stopifnot(mode%in%c("archived","smoke","full"))
root<-normalizePath(getwd());stopifnot(file.exists("reproduction/FROZEN_SHA256.json"))
out<-option("output",file.path(root,"output",paste0("reproduction_",mode,"_",format(Sys.time(),"%Y%m%d_%H%M%S"))))
if(file.exists(out))stop("Output already exists: ",out)
if(mode=="archived") {
  status<-system2(Sys.getenv("PYTHON","python3"),c("reproduction/rebuild_archived.py",shQuote(out)))
  quit(save="no",status=status)
}
reps<-as.integer(option("reps",if(mode=="full")"500" else "2"));seed<-as.integer(option("seed-base","2025"))
if(mode=="full")stopifnot(reps==500L,seed==2025L)
if(mode=="smoke")stopifnot(reps>0L,reps<=10L)
baseline<-option("llqr-baseline","formula");stopifnot(baseline%in%c("formula","historical"))
profile<-option("profile","portable")
grid<-if(mode=="full")expand.grid(paper_case=1:4,tau=c(.2,.5,.8),n=c(1000L,2000L,5000L,10000L)) else data.frame(paper_case=1:4,tau=.5,n=100L)
index<-option("config");if(!is.null(index)){index<-as.integer(index);stopifnot(index>=1L,index<=nrow(grid));grid<-grid[index,,drop=FALSE]}
dir.create(out,recursive=TRUE);out<-normalizePath(out)
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,BLIS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
# All worker subprocesses inherit the fixed single-thread settings.
work<-file.path(out,"workspace");dir.create(work)
file.copy("reproduction/frozen/R",work,recursive=TRUE)
if(baseline=="formula")stopifnot(file.copy("R/llqr_functions.R",file.path(work,"R","llqr_functions.R"),overwrite=TRUE))
expdir<-file.path(work,"experiment");dir.create(expdir)
files<-list.files("reproduction/frozen",full.names=TRUE,no..=TRUE,all.files=TRUE)
file.copy(files,expdir,recursive=TRUE)
source("reproduction/build.R");build_u11(file.path(expdir,"build"),profile,root)
cpu<-if(file.exists("/proc/cpuinfo"))grep("^model name",readLines("/proc/cpuinfo"),value=TRUE)[1L] else if(Sys.info()[["sysname"]]=="Darwin")system2("sysctl",c("-n","machdep.cpu.brand_string"),stdout=TRUE) else Sys.info()[["machine"]]
os_release<-if(file.exists("/etc/os-release"))readLines("/etc/os-release") else character()
meta<-list(llqr_baseline=baseline,cpu=cpu,os_release=os_release,external_libraries=extSoftVersion(),session=capture.output(sessionInfo()),system=Sys.info(),profile=profile,
  packages=lapply(c("quantreg","KernSmooth"),function(p)list(package=p,version=as.character(packageVersion(p)))),
  source_sha=system2("git","rev-parse HEAD",stdout=TRUE),threads=Sys.getenv(c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS")))
saveRDS(meta,file.path(out,"environment.rds"));writeLines(capture.output(str(meta)),file.path(out,"environment.txt"))
write.csv(grid,file.path(out,"configuration_map.csv"),row.names=FALSE)
Sys.setenv(SSQR_PROJECT_ROOT=work,SSQR_EXPERIMENT_DIR=expdir,SSQR_OUTPUT_ROOT=out,
 SSQR_NUM_REP=reps,SSQR_SEED_BASE=seed,SSQR_CHUNK_SIZE=reps,SLURM_ARRAY_TASK_ID=1,
 SLURM_CPUS_PER_TASK=1,SSQR_ALLOW_SMOKE=if(mode=="smoke")"1" else "0",
 SSQR_CACHE_FLAGS=27,SSQR_PROVIDER_FLAGS=1,SSQR_BUILD_MODE="optimized")
rscript<-file.path(R.home("bin"),"Rscript")
for(i in seq_len(nrow(grid))){
 g<-grid[i,]; logistic<-g$paper_case==1L
 Sys.setenv(SSQR_RUN_TAG=paste0("paper_case",g$paper_case),SSQR_MODEL=if(g$paper_case<=2)"llqr" else "tvcqr",
   SSQR_CASE=if(g$paper_case<=2)2 else g$paper_case-2,SSQR_N=g$n,SSQR_TAU=g$tau)
 for(script in c("driver_array.R","merge_config.R")){
   name<-if(logistic)script else paste0("original_",script)
   status<-system2(rscript,c("--vanilla",shQuote(file.path(expdir,name))))
   if(status!=0L)stop("Experiment failed: ",name,"; output retained at ",out)
 }
}
# Check scientific acceptance as well as file completeness.
paths<-list.files(out,pattern="^replication_metrics[.]rds$",recursive=TRUE,full.names=TRUE)
rows<-lapply(paths,readRDS)
stopifnot(length(rows)==nrow(grid),all(vapply(rows,function(x)nrow(x)==3L*reps && all(x$accepted_ok) && all(x$solver_ok) && !any(x$fallback_triggered),logical(1))))
summary<-do.call(rbind,lapply(seq_along(rows),function(i){x<-rows[[i]];do.call(rbind,lapply(split(x,x$method),function(z)data.frame(model=z$model[1],raw_case=z$case[1],paper_case=as.integer(sub("paper_case","",z$run_tag[1])),n=z$n[1],tau=z$tau[1],method=z$method[1],replications=nrow(z),time_mean_sec=mean(z$elapsed_sec),time_sd_sec=sd(z$elapsed_sec),max_discrepancy=max(z$discrepancy))))}))
write.csv(summary,file.path(out,"new_run_summary.csv"),row.names=FALSE)
cat("PASS: ",nrow(grid)," configurations; ",reps," paired replications each. Results: ",out,"\n",sep="")
