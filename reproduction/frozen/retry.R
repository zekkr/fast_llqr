experiment_methods <- c("direct_baseline","lean_seq","unified_u11")

baseline_retryable <- function(rows) {
  x<-rows[rows$method=="direct_baseline",,drop=FALSE]
  nrow(x)==1L&&isTRUE(x$threw_error)&&identical(x$error_stage,"fit")
}

run_attempt_chain <- function(model,rep_id,seed_base,attempt_fun,max_attempts=20L,stride=1000000L) {
  limit<-as.integer(max_attempts)
  original<-as.double(seed_base)+rep_id
  if(original+(limit-1)*stride>.Machine$integer.max)stop("seed overflow")
  saved<-list()
  for(attempt in seq_len(limit)) {
    seed<-as.integer(original+(attempt-1)*stride)
    rows<-attempt_fun(rep_id,seed)
    if(!is.data.frame(rows)||nrow(rows)!=length(experiment_methods)||!setequal(rows$method,experiment_methods))
      stop("attempt must return exactly the registered experiment methods")
    rows$original_seed<-as.integer(original);rows$attempt<-attempt
    rows$max_attempts<-limit;rows$retry_stride<-as.integer(stride)
    eligible<-baseline_retryable(rows)
    rows$baseline_retryable<-eligible
    rows$retry_exhausted<-eligible&&attempt==limit
    saved[[attempt]]<-rows
    if(!eligible)break
  }
  all<-do.call(rbind,saved);all$is_final_attempt<-all$attempt==length(saved)
  rownames(all)<-NULL;final<-all[all$is_final_attempt,,drop=FALSE];rownames(final)<-NULL
  list(schema_version=2L,final=final,attempts=all)
}

validate_attempt_chain <- function(value,identity,seed_base) {
  if(!is.list(value)||!identical(value$schema_version,2L)||
     !is.data.frame(value$attempts)||!is.data.frame(value$final))return(FALSE)
  a<-value$attempts
  required<-c(names(identity),"rep_id","seed","method","method_position","original_seed",
    "attempt","max_attempts","retry_stride","baseline_retryable","retry_exhausted",
    "is_final_attempt","solver_ok","accepted_ok","finite_estimate","elapsed_sec",
    "discrepancy","discrepancy_status","ierr","failed_eval","fallback_triggered",
    "h_check_status","h_match_vs_lean","h_mismatch_eval_count","first_h_mismatch_eval",
    "threw_error","error_stage")
  if(!nrow(a)||!all(required%in%names(a))||anyNA(a[c("rep_id","seed","method_position","attempt")]))return(FALSE)
  for(k in names(identity))if(!isTRUE(all(a[[k]]==identity[[k]])))return(FALSE)
  rep_id<-identity$rep_id;last<-max(a$attempt);limit<-20L
  if(last<1L||last>limit||!identical(sort(unique(a$attempt)),seq_len(last)))return(FALSE)
  if(!isTRUE(all(a$original_seed==seed_base+rep_id))||!isTRUE(all(a$max_attempts==limit))||
     !isTRUE(all(a$retry_stride==1000000L)))return(FALSE)
  for(k in seq_len(last)) {
    x<-a[a$attempt==k,,drop=FALSE];shift<-(rep_id-1L)%%length(experiment_methods)
    expected_pos<-((match(x$method,experiment_methods)-1L-shift)%%length(experiment_methods))+1L
    eligible<-baseline_retryable(x)
    if(nrow(x)!=length(experiment_methods)||!setequal(x$method,experiment_methods)||
       !isTRUE(all(x$seed==seed_base+rep_id+(k-1L)*1000000L))||
       !isTRUE(all(x$method_position==expected_pos))||
       !isTRUE(all(x$is_final_attempt==(k==last)))||
       !isTRUE(all(x$baseline_retryable==eligible))||
       !isTRUE(all(x$retry_exhausted==(eligible&&k==limit))))return(FALSE)
    if(k<last&&!eligible)return(FALSE)
    if(k==last&&eligible&&k<limit)return(FALSE)
  }
  final<-a[a$is_final_attempt,,drop=FALSE];rownames(final)<-NULL
  identical(final,value$final)
}
