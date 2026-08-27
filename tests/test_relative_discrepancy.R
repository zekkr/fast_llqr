#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 0L) {
  stop("Run this test via Rscript.")
}

test_path <- normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
repo_dir <- normalizePath(file.path(dirname(test_path), ".."), mustWork = TRUE)

source(file.path(repo_dir, "R", "performance_measurement.R"))
source(file.path(repo_dir, "R", "llqr_simulation_helpers.R"))
source(file.path(repo_dir, "R", "tvcqr_simulation_helpers.R"))

assert_close <- function(actual, expected, label, tolerance = 1e-15) {
  if (length(actual) != length(expected) ||
      any(!is.finite(actual)) ||
      !isTRUE(all.equal(actual, expected, tolerance = tolerance))) {
    stop(sprintf(
      "%s failed: actual=%s expected=%s",
      label,
      paste(format(actual, digits = 17), collapse = ","),
      paste(format(expected, digits = 17), collapse = ",")
    ))
  }
}

llqr_reference <- c(0, 5e-11, 1e-10, -2, 0)
llqr_method <- c(1e-10, 1.5e-10, 2e-10, -1, 0)
llqr_expected <- mean(c(1, 1, 1, 0.5, 0))
llqr_results <- list(
  config = list(n = 5L, num_rep = 1L),
  method_names = c("llqr", "screen"),
  estimates_list = list(
    llqr = list(llqr_reference),
    screen = list(llqr_method)
  )
)

assert_close(
  calculate_average_relative_bias_llqr(llqr_results, "screen"),
  llqr_expected,
  "LLQR reporting function"
)
assert_close(
  unname(compute_llqr_max_average_relative_bias(llqr_results)["screen"]),
  llqr_expected,
  "LLQR simulation summary"
)
llqr_find_output <- capture.output(
  llqr_find <- find_llqr_max_bias_replication(llqr_results, "screen")
)
assert_close(llqr_find$max_bias, llqr_expected, "LLQR maximum-replication diagnostic")

tvcqr_reference <- matrix(c(0, 5e-11, 1e-10, -2), nrow = 1L)
tvcqr_method <- matrix(c(1e-10, 1.5e-10, 2e-10, -1), nrow = 1L)
tvcqr_expected <- mean(c(1, 1, 1, 0.5))
tvcqr_results <- list(
  config = list(n = 1L, num_rep = 1L),
  method_names = c("tvc_rq", "screen"),
  estimates_list = list(
    tvc_rq = list(tvcqr_reference),
    screen = list(tvcqr_method)
  )
)

assert_close(
  calculate_average_relative_bias_tvcqr(tvcqr_results, "screen"),
  tvcqr_expected,
  "TVCQR reporting function"
)
assert_close(
  unname(compute_tvcqr_max_average_relative_bias(tvcqr_results)["screen"]),
  tvcqr_expected,
  "TVCQR simulation summary"
)
tvcqr_find_output <- capture.output(
  tvcqr_find <- find_tvcqr_max_bias_replication(tvcqr_results, "screen")
)
assert_close(tvcqr_find$max_bias, tvcqr_expected, "TVCQR maximum-replication diagnostic")

invalid_llqr <- llqr_results
invalid_llqr$estimates_list$screen[[1L]][2L] <- NA_real_
invalid_output <- capture.output(
  invalid_value <- calculate_average_relative_bias_llqr(invalid_llqr, "screen")
)
if (!is.na(invalid_value)) {
  stop("Non-finite LLQR input should invalidate the full replication discrepancy.")
}
if (!is.na(suppressWarnings(compute_llqr_max_average_relative_bias(invalid_llqr)["screen"]))) {
  stop("LLQR maximum should be unavailable when any replication discrepancy is invalid.")
}

invalid_tvcqr <- tvcqr_results
invalid_tvcqr$estimates_list$screen[[1L]] <- matrix(1:3, nrow = 1L)
invalid_output <- capture.output(
  invalid_value <- calculate_average_relative_bias_tvcqr(invalid_tvcqr, "screen")
)
if (!is.na(invalid_value)) {
  stop("A TVCQR dimension mismatch should invalidate the full replication discrepancy.")
}
if (!is.na(suppressWarnings(compute_tvcqr_max_average_relative_bias(invalid_tvcqr)["screen"]))) {
  stop("TVCQR maximum should be unavailable when any replication discrepancy is invalid.")
}

invalid_llqr_baseline <- llqr_results
invalid_llqr_baseline$estimates_list$llqr[[1L]][2L] <- NA_real_
invalid_output <- capture.output(
  invalid_value <- calculate_average_relative_bias_llqr(invalid_llqr_baseline, "llqr")
)
if (!is.na(invalid_value)) {
  stop("A non-finite LLQR baseline must not be reported as zero discrepancy.")
}
if (!is.na(suppressWarnings(compute_llqr_max_average_relative_bias(invalid_llqr_baseline)["llqr"]))) {
  stop("The LLQR baseline maximum must be unavailable when its estimates are invalid.")
}

invalid_tvcqr_baseline <- tvcqr_results
invalid_tvcqr_baseline$estimates_list$tvc_rq[[1L]] <- matrix(1:3, nrow = 1L)
invalid_output <- capture.output(
  invalid_value <- calculate_average_relative_bias_tvcqr(invalid_tvcqr_baseline, "tvc_rq")
)
if (!is.na(invalid_value)) {
  stop("A malformed TVCQR baseline must not be reported as zero discrepancy.")
}
if (!is.na(suppressWarnings(compute_tvcqr_max_average_relative_bias(invalid_tvcqr_baseline)["tvc_rq"]))) {
  stop("The TVCQR baseline maximum must be unavailable when its estimates are invalid.")
}

cat("Relative discrepancy tests passed.\n")
