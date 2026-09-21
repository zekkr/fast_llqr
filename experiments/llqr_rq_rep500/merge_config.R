#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
geti<-function(k,d=NA_integer_){x<-Sys.getenv(k,"");if(nzchar(x))as.integer(x)else as.integer(d)}
getn<-function(k,d=NA_real_){x<-Sys.getenv(k,"");if(nzchar(x))as.numeric(x)else as.numeric(d)}
project<-normalizePath(Sys.getenv("SSQR_PROJECT_ROOT",getwd()),mustWork=TRUE);setwd(project)
experiment_dir<-Sys.getenv("SSQR_EXPERIMENT_DIR","experiments/llqr_rq_rep500")
source("reproduction/frozen/retry.R")
model<-tolower(Sys.getenv("SSQR_MODEL",""));case_id<-geti("SSQR_CASE");tau<-getn("SSQR_TAU");n<-geti("SSQR_N")
num_rep<-geti("SSQR_NUM_REP",500L);seed_base<-geti("SSQR_SEED_BASE",2025L)
run_tag<-Sys.getenv("SSQR_RUN_TAG","");output_root<-Sys.getenv("SSQR_OUTPUT_ROOT","")
smoke<-Sys.getenv("SSQR_ALLOW_SMOKE","0")=="1"
stopifnot(model=="llqr",case_id%in%1:2,(smoke||num_rep==500L),(smoke||seed_base==2025L),nzchar(run_tag),nzchar(output_root))
tag<-sprintf("case%d_tau%02d_n%d",case_id,round(100*tau),n);dir<-file.path(output_root,run_tag,model,tag);partial<-file.path(dir,"partials")
observed<-list.files(partial,pattern="^rep[0-9]+[.]rds$");expected<-sprintf("rep%04d.rds",seq_len(num_rep))
rows<-list();attempts<-list();missing<-integer();malformed<-integer()
for(rep_id in seq_len(num_rep)){
  path<-file.path(partial,sprintf("rep%04d.rds",rep_id));if(!file.exists(path)){missing<-c(missing,rep_id);next}
  value<-tryCatch(readRDS(path),error=function(e)e);identity<-list(run_tag=run_tag,model=model,case=case_id,tau=tau,n=n,rep_id=rep_id)
  valid<-!inherits(value,"error")&&isTRUE(tryCatch(validate_attempt_chain(value,identity,seed_base),error=function(e)FALSE))
  if(!valid){malformed<-c(malformed,rep_id);next}
  rows[[length(rows)+1L]]<-value$final[match(experiment_methods,value$final$method),,drop=FALSE]
  attempts[[length(attempts)+1L]]<-value$attempts
}
merged<-if(length(rows))do.call(rbind,rows)else data.frame();all_attempts<-if(length(attempts))do.call(rbind,attempts)else data.frame()
if(nrow(merged)){merged<-merged[order(merged$rep_id,match(merged$method,experiment_methods)),,drop=FALSE];rownames(merged)<-NULL}
dir.create(dir,recursive=TRUE,showWarnings=FALSE);saveRDS(merged,file.path(dir,"replication_metrics.rds"),compress="xz")
write.csv(merged,file.path(dir,"replication_metrics.csv"),row.names=FALSE,na="")
saveRDS(all_attempts,file.path(dir,"attempt_metrics.rds"),compress="xz");write.csv(all_attempts,file.path(dir,"attempt_metrics.csv"),row.names=FALSE,na="")
extra<-setdiff(observed,expected);integrity<-data.frame(run_tag,model,case=case_id,tau,n,num_rep,
  n_complete=length(rows),n_missing=length(missing),missing_ids=paste(missing,collapse=","),
  n_malformed=length(malformed),malformed_ids=paste(malformed,collapse=","),n_extra=length(extra),extra_files=paste(extra,collapse=","),
  complete=length(rows)==num_rep&&!length(missing)&&!length(malformed)&&!length(extra),stringsAsFactors=FALSE)
write.csv(integrity,file.path(dir,"integrity.csv"),row.names=FALSE,na="")
cat(sprintf("merged=%d/%d missing=%d malformed=%d extra=%d\n",length(rows),num_rep,length(missing),length(malformed),length(extra)))
if(!integrity$complete)quit(save="no",status=2L)

if (!nrow(merged) || !all(merged$accepted_ok %in% TRUE) || !all(merged$discrepancy_status == "ok")) stop("Unaccepted fits; do not publish")
