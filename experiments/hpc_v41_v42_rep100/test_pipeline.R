#!/usr/bin/env Rscript
# Small synthetic orchestration fixture; no scientific results are generated.
source("experiments/hpc_v41_v42_rep100/retry_helpers.R")
source_path <- "data/v41_v42_rep100/llqr_regen_smoke_rep2/llqr/case2_tau50_n1000/partials/rep0001.rds"
if (!file.exists(source_path)) stop("Run the documented rep2 local smoke first")
template <- readRDS(source_path)$final
settings <- list(policy="paper_baseline_error", max_attempts=20L, stride=1000000)
root <- tempfile("llqr_regen_pipeline_")
dir.create(root)
Sys.setenv(FASTQR_PROJECT_DIR=getwd(), FASTQR_V4142_BASE_DIR=root,
  FASTQR_RUN_TAG="pipeline", FASTQR_MODEL="llqr", FASTQR_MODELS="llqr",
  FASTQR_NUM_REP=4L, FASTQR_SEED_BASE=2026L,
  FASTQR_BASELINE_RETRY_POLICY="paper_baseline_error", FASTQR_MAX_ATTEMPTS_PER_REP=20L,
  FASTQR_RETRY_STRIDE=1000000, FASTQR_PREVIOUS_RUN_TAG="")
for (case_id in 1:2) for (tau in c(.2,.5,.8)) for (n in c(1000L,2000L,5000L,10000L)) {
  tag <- sprintf("case%d_tau%02d_n%d", case_id, as.integer(100*tau), n)
  folder <- file.path(root,"pipeline","llqr",tag,"partials")
  dir.create(folder, recursive=TRUE)
  for (rep_id in 1:4) {
    result <- run_baseline_retries(rep_id,2026L,settings,function(r,seed) {
      x <- template
      x$run_tag<-"pipeline"; x$case<-case_id; x$tau<-tau; x$n<-n
      x$rep_id<-r; x$seed<-seed
      x$method_position<-((seq_len(4L)-1L-(r-1L)%%4L)%%4L)+1L
      x$elapsed_sec<-c(1,.1,.09,.095); x$data_generation_sec<-.001
      # Each config has one replacement; candidates fail on that discarded attempt.
      x$method_ok<-TRUE; x$threw_error<-FALSE; x$error_stage<-NA_character_
      if (r==1L && seed==2027L) {
        x$method_ok[c(1,3,4)]<-FALSE; x$threw_error[1L]<-TRUE; x$error_stage[1L]<-"fit"
      }
      x
    })
    saveRDS(result,file.path(folder,sprintf("rep%04d.rds",rep_id)))
  }
  Sys.setenv(FASTQR_CASE=case_id,FASTQR_TAU=tau,FASTQR_N=n)
  status <- system2(file.path(R.home("bin"),"Rscript"),"experiments/hpc_v41_v42_rep100/merge_config.R", stdout=FALSE)
  stopifnot(status==0L)
}
status <- system2(file.path(R.home("bin"),"Rscript"),"experiments/hpc_v41_v42_rep100/summarize_run.R")
stopifnot(status==0L)
table_dir <- file.path(root,"pipeline","tables")
summary <- read.csv(file.path(table_dir,"method_summary.csv"))
ratios <- read.csv(file.path(table_dir,"candidate_vs_lean_seq.csv"))
attempts <- read.csv(file.path(table_dir,"attempt_metrics_all.csv"))
retries <- read.csv(file.path(table_dir,"retry_summary.csv"))
stopifnot(nrow(summary)==96L, nrow(ratios)==48L, nrow(attempts)==480L,
          all(summary$n_success==4L), all(summary$timing_order_balanced),
          sum(retries$n_replaced_replications)==24L, sum(!attempts$method_ok)==72L,
          all(ratios$candidate_over_lean_seq < 1))
cat("Passed complete LLQR-only 24-config synthetic merge/summary pipeline\n")
cat("Fixture path:",root,"\n")
