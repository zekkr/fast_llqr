#!/usr/bin/env Rscript
experiment_dir <- Sys.getenv("SSQR_EXPERIMENT_DIR", "experiments/hpc_u11_rep500")
source(file.path(experiment_dir, "metrics.R"))
source(file.path(experiment_dir, "retry.R"))
check <- function(value, label) {
  if (!isTRUE(value)) stop(label)
  cat("PASS", label, "\n")
}

d <- paper_discrepancy(c(1, 2), c(0, 2), 2L)
check(d$status == "ok" && d$value == 5e9, "paper discrepancy denominator floor")
check(paper_discrepancy(1, c(1, 2), 2L)$status == "dimension_mismatch", "dimension mismatch invalid")
h <- matrix(c(1, 2, 3, 2, 3, 1), 3, 2)
z <- audit_h_path("unified_u11", TRUE, h, TRUE, h[, 2:1], 3L, 2L)
check(z$accepted && z$count == 0L, "H row-set ignores column order")
bad <- h[, 2:1]
bad[1L, 1L] <- 3L
z <- audit_h_path("unified_u11", TRUE, h, TRUE, bad, 3L, 2L)
check(!z$accepted && z$count == 1L, "H mismatch rejected")

mock <- function(fail_first = FALSE, candidate_fail = FALSE) function(rep_id, seed) {
  positions <- ((seq_along(experiment_methods) - 1L -
    (rep_id - 1L) %% length(experiment_methods)) %% length(experiment_methods)) + 1L
  x <- data.frame(
    run_tag = "unit", model = "llqr", case = 1L, tau = 0.2, n = 1000L,
    rep_id, seed, method = experiment_methods, method_position = positions,
    elapsed_sec = 0.1, solver_ok = TRUE, accepted_ok = TRUE,
    finite_estimate = TRUE, discrepancy = 0, discrepancy_status = "ok",
    ierr = c(NA, NA, 0L), failed_eval = c(NA, NA, 0L),
    fallback_triggered = FALSE,
    h_check_status = c("not_applicable_direct", "match", "match"),
    h_match_vs_lean = c(NA, TRUE, TRUE), h_mismatch_eval_count = c(NA, 0L, 0L),
    first_h_mismatch_eval = NA_integer_, threw_error = FALSE,
    error_stage = NA_character_
  )
  if (fail_first && seed == 2042L) {
    x$solver_ok[1L] <- x$accepted_ok[1L] <- FALSE
    x$threw_error[1L] <- TRUE
    x$error_stage[1L] <- "fit"
  }
  if (candidate_fail) x$solver_ok[3L] <- x$accepted_ok[3L] <- FALSE
  x
}

identity <- list(run_tag = "unit", model = "llqr", case = 1L, tau = 0.2, n = 1000L, rep_id = 17L)
x <- run_attempt_chain("llqr", 17L, 2025L, mock(TRUE))
check(validate_attempt_chain(x, identity, 2025L) && nrow(x$attempts) == 6L &&
  all(x$final$seed == 1002042L), "baseline-only regeneration")
x <- run_attempt_chain("llqr", 17L, 2025L, mock(FALSE, TRUE))
check(nrow(x$attempts) == 3L, "candidate failure does not regenerate")
x <- run_attempt_chain("tvcqr", 17L, 2025L, mock(TRUE))
check(nrow(x$attempts) == 6L, "TVCQR baseline failure follows the same regeneration rule")

orders <- do.call(rbind, lapply(1:500, function(rep_id) mock()(rep_id, 2025L + rep_id)))
tab <- table(orders$method, orders$method_position)
check(all(apply(tab, 1L, function(v) diff(range(v)) <= 1L)), "rep500 three-method order is balanced")

no_retry_seq <- function(rep_id, seed) {
  x <- mock()(rep_id,seed);x$solver_ok[2L] <- x$accepted_ok[2L] <- FALSE
  x$threw_error[2L] <- TRUE;x$error_stage[2L] <- "fit";x
}
check(nrow(run_attempt_chain("llqr",17L,2025L,no_retry_seq)$attempts)==3L,
      "seq exception does not regenerate")
always_fail <- function(rep_id, seed) {
  x <- mock()(rep_id,seed);x$solver_ok[1L] <- x$accepted_ok[1L] <- FALSE
  x$threw_error[1L] <- TRUE;x$error_stage[1L] <- "fit";x
}
x <- run_attempt_chain("llqr",17L,2025L,always_fail)
check(nrow(x$attempts)==60L && all(x$final$retry_exhausted) && max(x$attempts$attempt)==20L,
      "baseline retry cap retains all twenty failed attempts")
