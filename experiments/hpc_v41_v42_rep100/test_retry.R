#!/usr/bin/env Rscript
source("experiments/hpc_v41_v42_rep100/retry_helpers.R")
settings <- list(policy = "paper_baseline_error", max_attempts = 20L, stride = 1000000)
identity <- list(run_tag = "unit", model = "llqr", case = 1L, tau = .2, n = 2000L)
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop(label)
  checks <<- checks + 1L
  cat("ok:", label, "\n")
}
mock_attempt <- function(baseline_error = function(seed) FALSE, candidate_error = function(seed) FALSE,
                         baseline_nonfinite = FALSE, generation_error = FALSE) {
  function(rep_id, seed) {
    x <- data.frame(run_tag = "unit", model = "llqr", case = 1L, tau = .2, n = 2000L,
      rep_id = rep_id, seed = seed, method = experiment_methods,
      method_position = ((seq_len(4L)-1L-(rep_id-1L)%%4L)%%4L)+1L,
      method_ok = TRUE, threw_error = FALSE, error_stage = NA_character_,
      error_message = NA_character_, elapsed_sec = c(1,.1,.12,.11))
    if (baseline_error(seed)) {
      x$method_ok[1L] <- FALSE; x$threw_error[1L] <- TRUE; x$error_stage[1L] <- "fit"
      x$error_message[1L] <- "Singular design matrix"
    }
    if (candidate_error(seed)) x$method_ok[3:4] <- FALSE
    if (baseline_nonfinite) x$method_ok[1L] <- FALSE
    if (generation_error) {
      x$method_ok <- FALSE; x$threw_error <- TRUE; x$error_stage <- "data_generation"
    }
    x
  }
}
run <- function(fun, st = settings) run_baseline_retries(17L, 2026L, st, fun)
valid <- function(x, st = settings) validate_retry_result(x, 17L, 2026L, identity, st)
x <- run(mock_attempt())
check(valid(x) && nrow(x$attempts)==4L && all(x$final$seed==2043L), "first attempt success")
x <- run(mock_attempt(function(seed) seed==2043L))
check(valid(x) && nrow(x$attempts)==8L && all(x$final$seed==1002043L), "baseline error regenerates all four methods")
check(all(x$attempts$method_position == rep(1:4,2L)), "rotation unchanged between attempts")
x <- run(mock_attempt(function(seed) TRUE))
check(valid(x) && nrow(x$attempts)==80L && all(x$final$retry_exhausted) && all(x$final$seed==19002043L), "twenty attempts exhausted")
x <- run(mock_attempt(candidate_error=function(seed) TRUE))
check(valid(x) && nrow(x$attempts)==4L && all(!x$final$method_ok[3:4]), "candidate failure never regenerates data")
x <- run(mock_attempt(function(seed) seed==2043L, function(seed) seed==2043L))
check(valid(x) && nrow(x$attempts)==8L && all(!x$attempts$method_ok[3:4]) && all(x$final$method_ok), "simultaneous failures retained in discarded attempt")
check(nrow(run(mock_attempt(baseline_nonfinite=TRUE))$attempts)==4L, "baseline nonfinite result is terminal without exception")
check(nrow(run(mock_attempt(generation_error=TRUE))$attempts)==4L, "data-generation error never triggers baseline retry")
off <- modifyList(settings, list(policy="none"))
xoff <- run(mock_attempt(function(seed) TRUE),off)
check(valid(xoff,off) && nrow(xoff$attempts)==4L, "legacy no-retry policy remains available")
corrupt <- x; corrupt$attempts$seed[5L] <- 7L
check(!valid(corrupt), "reject mismatched paired seed")
corrupt <- x; corrupt$attempts$threw_error[1L] <- FALSE
check(!valid(corrupt), "reject replacement not justified by baseline exception")
corrupt <- x; corrupt$final$seed <- 2043L
check(!valid(corrupt), "reject final output inconsistent with attempt log")
corrupt <- x; corrupt$attempts <- corrupt$attempts[-2L,]
check(!valid(corrupt), "reject missing method in attempt")
corrupt <- x; corrupt$attempts$case <- 2L
check(!valid(corrupt), "reject wrong configuration identity")
orders <- do.call(rbind, lapply(1:100, function(r) mock_attempt()(r,2026L+r)))
check(all(table(orders$method,orders$method_position)==25L), "100-rep rotation exactly balanced")
cat(sprintf("Passed %d retry/merge checks\n",checks))
