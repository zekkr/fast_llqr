# Compact experiment metrics. Full estimates and H paths are transient only.
paper_discrepancy <- function(candidate, baseline, expected_k) {
  a <- as.numeric(candidate); b <- as.numeric(baseline)
  if (length(a) != expected_k || length(b) != expected_k) {
    return(list(value=NA_real_,status="dimension_mismatch"))
  }
  if (any(!is.finite(a)) || any(!is.finite(b))) {
    return(list(value=NA_real_,status="nonfinite_estimate"))
  }
  value <- mean(abs(a-b)/pmax(abs(b),1e-10))
  list(value=value,status=if(is.finite(value))"ok" else "nonfinite_discrepancy")
}

valid_h_path <- function(h,n,q,ne=n) {
  if (!is.numeric(h) || length(h)!=ne*q || any(!is.finite(h)) ||
      any(h!=floor(h)) || any(h<1L | h>n)) return(FALSE)
  z <- matrix(as.integer(h),nrow=ne,ncol=q)
  all(vapply(seq_len(ne),function(i)length(unique(z[i,]))==q,logical(1L)))
}

audit_h_path <- function(method,solver_ok,candidate_h,lean_ok,lean_h,n,q,ne=n) {
  out <- function(status,match=FALSE,count=NA_integer_,first=NA_integer_) {
    list(status=status,match=match,count=count,first=first,
         accepted=isTRUE(solver_ok)&&(method=="direct_baseline"||method=="lean_seq"||isTRUE(match)))
  }
  if(method=="direct_baseline")return(out("not_applicable_direct",NA))
  if(!isTRUE(solver_ok))return(out("method_failed"))
  if(!isTRUE(lean_ok))return(out("lean_failed"))
  if(!valid_h_path(lean_h,n,q,ne))return(out("invalid_lean_h"))
  if(!valid_h_path(candidate_h,n,q,ne))return(out("invalid_candidate_h"))
  a<-matrix(as.integer(candidate_h),ne,q);b<-matrix(as.integer(lean_h),ne,q)
  bad<-which(!vapply(seq_len(ne),function(i)identical(sort.int(a[i,]),sort.int(b[i,])),logical(1L)))
  if(length(bad))out("row_set_mismatch",FALSE,length(bad),bad[[1L]]) else out("match",TRUE,0L)
}

empty_fit <- function(method,position,elapsed=NA_real_,message=NA_character_,stage="fit") {
  list(method=method,position=as.integer(position),elapsed_sec=as.numeric(elapsed),
       estimate=NULL,H=NULL,solver_ok=FALSE,finite_estimate=FALSE,iteration_limit_hit=NA,ierr=NA_integer_,
       failed_eval=NA_integer_,first_tableau_rows_max=NA_integer_,
       first_retained_size_max=NA_integer_,initial_threshold_hits_max=NA_integer_,
       effective_threshold_hits_max=NA_integer_,threshold_expansion_total=NA_integer_,
       threshold_expansion_points=NA_integer_,threshold_expansion_steps_max=NA_integer_,
       threshold_initial=NA_real_,threshold_effective_max=NA_real_,basis_forced_rows_max=NA_integer_,
       first_aggregate_rows_max=NA_integer_,retained_decomposition_ok=NA,
       first_pass_count=NA_integer_,first_pass_total=NA_integer_,first_pass_rate=NA_real_,
       repaired_points=NA_integer_,repair_total=NA_integer_,
       internal_recovery_total=NA_integer_,first_full_m_recovery=NA_integer_,
       certificate_hits=NA_integer_,residual_rows_checked=NA_integer_,
       same_h_refit_attempted=NA_integer_,same_h_refit_recovered=NA_integer_,
       error_message=as.character(message),threw_error=TRUE,error_stage=stage)
}
