#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

as_int <- function(name, default = NA_integer_) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  as.integer(value)
}

as_num <- function(name, default = NA_real_) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  as.numeric(value)
}

project_dir <- normalizePath(Sys.getenv("FASTQR_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
setwd(project_dir)
source(file.path(project_dir, "experiments/hpc_v41_v42_rep100/retry_helpers.R"))
settings <- retry_settings()

model <- tolower(Sys.getenv("FASTQR_MODEL", unset = ""))
case_id <- as_int("FASTQR_CASE")
tau <- as_num("FASTQR_TAU")
n <- as_int("FASTQR_N")
num_rep <- as_int("FASTQR_NUM_REP", 100L)
seed_base <- as_int("FASTQR_SEED_BASE", 2026L)
run_tag <- Sys.getenv("FASTQR_RUN_TAG", unset = "")
base_dir <- Sys.getenv("FASTQR_V4142_BASE_DIR", unset = "data/v41_v42_rep100")

if (!model %in% c("llqr", "tvcqr")) stop("FASTQR_MODEL must be llqr or tvcqr.")
if (!case_id %in% c(1L, 2L)) stop("FASTQR_CASE must be 1 or 2.")
if (!nzchar(run_tag)) stop("FASTQR_RUN_TAG is required.")

expected_methods <- c("direct_baseline", "lean_seq", "v41", "v42")
config_tag <- sprintf("case%d_tau%02d_n%d", case_id, as.integer(round(100 * tau)), n)
config_dir <- file.path(base_dir, run_tag, model, config_tag)
partial_dir <- file.path(config_dir, "partials")
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)

rows <- list()
attempt_rows <- list()
missing_ids <- integer()
malformed_ids <- integer()

for (rep_id in seq_len(num_rep)) {
  path <- file.path(partial_dir, sprintf("rep%04d.rds", rep_id))
  if (!file.exists(path)) {
    missing_ids <- c(missing_ids, rep_id)
    next
  }
  value <- tryCatch(readRDS(path), error = function(e) e)
  if (is.data.frame(value) && settings$policy == "none") {
    legacy <- value
    legacy$threw_error <- FALSE
    legacy$error_stage <- NA_character_
    legacy$data_generation_sec <- NA_real_
    value <- run_baseline_retries(rep_id, seed_base, settings, function(...) legacy)
  }
  valid <- !inherits(value, "error") && isTRUE(tryCatch(validate_retry_result(
    value, rep_id, seed_base,
    list(run_tag = run_tag, model = model, case = case_id, tau = tau, n = n), settings
  ), error = function(e) FALSE))
  if (!valid) {
    malformed_ids <- c(malformed_ids, rep_id)
    next
  }
  rows[[length(rows) + 1L]] <- value$final[match(expected_methods, value$final$method), , drop = FALSE]
  attempt_rows[[length(attempt_rows) + 1L]] <- value$attempts
}

merged <- if (length(rows)) do.call(rbind, rows) else data.frame()
if (nrow(merged)) {
  merged <- merged[order(merged$rep_id, match(merged$method, expected_methods)), , drop = FALSE]
  rownames(merged) <- NULL
}

csv_path <- file.path(config_dir, "replication_metrics.csv")
rds_path <- file.path(config_dir, "replication_metrics.rds")
write.csv(merged, csv_path, row.names = FALSE, na = "")
saveRDS(merged, rds_path, compress = "xz")
attempts <- if (length(attempt_rows)) do.call(rbind, attempt_rows) else data.frame()
write.csv(attempts, file.path(config_dir, "attempt_metrics.csv"), row.names = FALSE, na = "")
saveRDS(attempts, file.path(config_dir, "attempt_metrics.rds"), compress = "xz")

method_failures <- if (nrow(merged)) {
  aggregate(!merged$method_ok, list(method = merged$method), sum)
} else {
  data.frame(method = expected_methods, x = num_rep)
}
names(method_failures)[[2L]] <- "n_failed"

integrity <- data.frame(
  run_tag = run_tag,
  model = model,
  case = case_id,
  tau = tau,
  n = n,
  num_rep = num_rep,
  n_complete_rep_files = length(rows),
  n_missing = length(missing_ids),
  missing_rep_ids = paste(missing_ids, collapse = ","),
  n_malformed = length(malformed_ids),
  malformed_rep_ids = paste(malformed_ids, collapse = ","),
  all_replications_present = length(rows) == num_rep && !length(missing_ids) && !length(malformed_ids),
  stringsAsFactors = FALSE
)
for (method in expected_methods) {
  integrity[[paste0("failed_", method)]] <- method_failures$n_failed[match(method, method_failures$method)]
}
write.csv(integrity, file.path(config_dir, "integrity.csv"), row.names = FALSE, na = "")

cat(sprintf(
  "merged=%d/%d missing=%d malformed=%d output=%s\n",
  length(rows), num_rep, length(missing_ids), length(malformed_ids), csv_path
))
