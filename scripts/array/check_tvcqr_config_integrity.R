setup_file <- if (file.exists("scripts/setup_hpc.R")) {
  "scripts/setup_hpc.R"
} else {
  "scripts/setup.R"
}
source(setup_file)

as_int <- function(x, default) {
  x <- Sys.getenv(x, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  as.integer(x)
}

as_num <- function(x, default) {
  x <- Sys.getenv(x, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  as.numeric(x)
}

normalize_msg <- function(msg) {
  if (is.null(msg) || length(msg) == 0 || all(is.na(msg)) || !nzchar(paste(msg, collapse = ""))) {
    return("EMPTY_ERROR_MSG")
  }
  msg <- paste(msg, collapse = " | ")
  msg <- gsub("[[:space:]]+", " ", msg)
  trimws(msg)
}

case_id <- as_int("FASTQR_CASE", 1L)
tau <- as_num("FASTQR_TAU", 0.5)
n <- as_int("FASTQR_N", 200L)
num_rep <- as_int("FASTQR_NUM_REP", 1000L)
case_label <- resolve_tvcqr_case(case_id)$case_label
base_dir <- Sys.getenv("FASTQR_TVCQR_BASE_DIR", unset = "data/tvcqr_simu_results")
out_dir <- Sys.getenv("FASTQR_TVCQR_CHECK_OUT_DIR", unset = "results/table")

tau_str <- sprintf("tau%02d", as.integer(round(100 * tau)))
partial_dir <- file.path(base_dir, ".array_tmp", sprintf("case%d_%s_n%d_rep%d", case_id, tau_str, n, num_rep))
final_file <- file.path(base_dir, sprintf("case%d_%s_n%d_rep%d.RData", case_id, tau_str, n, num_rep))

fs <- if (dir.exists(partial_dir)) list.files(partial_dir, pattern = "^rep[0-9]+\\.RData$", full.names = TRUE) else character()
ids <- as.integer(sub("^rep([0-9]+)\\.RData$", "\\1", basename(fs)))
ids <- ids[!is.na(ids)]
missing <- setdiff(seq_len(num_rep), ids)
failed <- integer()
failed_msg <- character()

for (f in fs) {
  env <- new.env()
  ok <- try(load(f, envir = env), silent = TRUE)
  rep_id <- as.integer(sub("^rep([0-9]+)\\.RData$", "\\1", basename(f)))
  if (inherits(ok, "try-error") || !exists("partial_results", envir = env)) {
    failed <- c(failed, rep_id)
    failed_msg <- c(failed_msg, "load failed / partial_results missing")
    next
  }
  pr <- get("partial_results", envir = env)
  if (is.null(pr$status) || pr$status != "success") {
    failed <- c(failed, rep_id)
    failed_msg <- c(failed_msg, normalize_msg(pr$error_msg))
  }
}

failed <- sort(unique(failed))
summary_df <- data.frame(
  case = case_id,
  case_label = case_label,
  tau = tau,
  n = n,
  final_exists = file.exists(final_file),
  n_files = length(fs),
  n_missing = length(missing),
  missing_rate = length(missing) / num_rep,
  n_failed = length(failed),
  failed_rate = length(failed) / num_rep,
  first_missing = if (length(missing)) paste(head(missing, 20L), collapse = ",") else "",
  first_failed = if (length(failed)) paste(head(failed, 20L), collapse = ",") else "",
  first_failed_msg = if (length(failed_msg)) failed_msg[1L] else "",
  stringsAsFactors = FALSE
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("tvcqr_integrity_case%d_%s_n%d_rep%d.csv", case_id, tau_str, n, num_rep))
write.csv(summary_df, out_file, row.names = FALSE)
print(summary_df, row.names = FALSE)
cat("Saved:", out_file, "\n")
