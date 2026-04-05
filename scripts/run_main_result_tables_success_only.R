#!/usr/bin/env Rscript

setup_file <- if (file.exists("scripts/setup_hpc.R")) {
  "scripts/setup_hpc.R"
} else {
  "scripts/setup.R"
}
source(setup_file)

tau_values <- c(0.2, 0.5, 0.8)
sample_sizes <- c(200, 500, 1000, 2000, 5000)
num_replications <- 1000L
output_dir <- "results/table"

cat("=== Main Result Table Export (success reps only) ===\n")
cat("This export relies on merged result objects.\n")
cat("Failed/missing replications are stored as NA/NULL and are skipped by summary statistics via na.rm=TRUE.\n\n")

jobs <- list(
  list(model = "llqr", case = 1L, filename_prefix = "llqr_result_case1"),
  list(model = "llqr", case = 2L, filename_prefix = "llqr_result_case2"),
  list(model = "tvcqr", case = 1L, filename_prefix = "tvcqr_result_case1"),
  list(model = "tvcqr", case = 2L, filename_prefix = "tvcqr_result_case2")
)

for (job in jobs) {
  cat(sprintf("Processing %s case %d ...\n", toupper(job$model), job$case))
  run_performance_analysis(
    case = job$case,
    tau = tau_values,
    n = sample_sizes,
    rep = num_replications,
    model = job$model,
    metrics = c("computation_time", "average_relative_bias"),
    probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
    output_format = "xlsx",
    output_dir = output_dir,
    filename_prefix = job$filename_prefix
  )
  cat("\n")
}

cat("Saved files:\n")
cat(file.path(output_dir, "llqr_result_case1.xlsx"), "\n")
cat(file.path(output_dir, "llqr_result_case2.xlsx"), "\n")
cat(file.path(output_dir, "tvcqr_result_case1.xlsx"), "\n")
cat(file.path(output_dir, "tvcqr_result_case2.xlsx"), "\n")
