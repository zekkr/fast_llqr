#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

preflight_dir <- Sys.getenv("SSQR_PREFLIGHT_DIR", "")
pushed_sha <- Sys.getenv("SSQR_PUSHED_SHA", "")
stopifnot(nzchar(preflight_dir), grepl("^[0-9a-f]{40}$", pushed_sha))

smoke_pass <- file.path(preflight_dir, "HPC_PIPELINE_SMOKE_PASS")
files <- sort(Sys.glob(file.path(preflight_dir, "cache_mode_audit", "config*.rds")))
stopifnot(file.exists(smoke_pass), identical(trimws(readLines(smoke_pass)), pushed_sha))
stopifnot(length(files) == 48L)
x <- do.call(rbind, lapply(files, readRDS))
stopifnot(nrow(x) == 144L, length(unique(x$task_id)) == 48L)
stopifnot(
  all(x$cache_flags %in% c(25L, 27L, 31L)),
  all(x$ierr == 0L), all(x$failed_eval == 0L),
  all(x$beta_bitwise_equal_25), all(x$beta_max_abs_diff_25 == 0),
  all(x$h_exact_equal_25), all(x$retained_diagnostics_equal_25),
  all(x$h_rowset_match_lean), all(x$h_mismatch_eval_count == 0L),
  all(x$retained_decomposition_ok)
)
write.csv(x, file.path(preflight_dir, "cache_mode_audit.csv"), row.names = FALSE)
writeLines(pushed_sha, file.path(preflight_dir, "PREFLIGHT_PASS"))
cat(sprintf("PASS HPC smoke and 48-config cache audit sha=%s\n", pushed_sha))
