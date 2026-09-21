#!/usr/bin/env Rscript
env_int <- function(name, default = NA_integer_) {
  value <- Sys.getenv(name, ""); if (nzchar(value)) as.integer(value) else as.integer(default)
}
case_id <- env_int("MST_CASE"); n <- env_int("MST_N")
num_rep <- env_int("MST_NUM_REP", 5L); run_tag <- Sys.getenv("MST_RUN_TAG", "")
output_root <- Sys.getenv("MST_OUTPUT_ROOT", "results/hpc/multivar_mst_runs")
config_dir <- file.path(output_root, run_tag, sprintf("case%d_n%d", case_id, n))
files <- list.files(file.path(config_dir, "partial"), "[.]csv$", full.names = TRUE)
if (!length(files)) stop("No pilot partial files found for ", config_dir)
rows <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
rows <- rows[order(rows$rep_id, rows$threshold_factor), , drop = FALSE]
expected <- expand.grid(rep_id = seq_len(num_rep), threshold_factor = c(0.1, 0.2, 0.5, 1.0))
actual <- rows[, c("rep_id", "threshold_factor")]
actual_key <- paste(actual$rep_id, actual$threshold_factor)
expected_key <- paste(expected$rep_id, expected$threshold_factor)
if (nrow(rows) != nrow(expected) || length(unique(actual_key)) != nrow(expected) ||
    !setequal(actual_key, expected_key))
  stop("Missing, duplicate or unexpected pilot rows for ", config_dir)
if (any(!rows$direct_ok | !rows$screen_ok | !rows$audit_ok | !rows$h_set_match))
  stop("A solver or audit failure is present in ", config_dir)
write.csv(rows, file.path(config_dir, "pilot_metrics.csv"), row.names = FALSE)
writeLines("PASS", file.path(config_dir, "CONFIG_PASS"))
