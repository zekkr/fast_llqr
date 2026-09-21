#!/usr/bin/env Rscript
env_int <- function(name, default = NA_integer_) {
  value <- Sys.getenv(name, "")
  if (nzchar(value)) as.integer(value) else as.integer(default)
}
case_id <- env_int("MST_CASE"); n <- env_int("MST_N")
num_rep <- env_int("MST_NUM_REP", 100L)
run_tag <- Sys.getenv("MST_RUN_TAG", "")
output_root <- Sys.getenv("MST_OUTPUT_ROOT", "results/hpc/multivar_mst_runs")
config_dir <- file.path(output_root, run_tag, sprintf("case%d_n%d", case_id, n))
files <- list.files(file.path(config_dir, "partial"), "[.]csv$", full.names = TRUE)
if (!length(files)) stop("No partial files found for ", config_dir)
rows <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
rows <- rows[order(rows$rep_id), , drop = FALSE]
if (!identical(rows$rep_id, seq_len(num_rep)))
  stop("Missing, duplicate or unexpected replication IDs for ", config_dir)
if (any(!rows$direct_ok | !rows$screen_ok | !rows$audit_ok | !rows$h_set_match))
  stop("A solver or audit failure is present in ", config_dir)
write.csv(rows, file.path(config_dir, "replication_metrics.csv"), row.names = FALSE)
writeLines("PASS", file.path(config_dir, "CONFIG_PASS"))
