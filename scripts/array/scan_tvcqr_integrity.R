setup_file <- if (file.exists("scripts/setup_hpc.R")) {
  "scripts/setup_hpc.R"
} else {
  "scripts/setup.R"
}
source(setup_file)

suppressPackageStartupMessages(library(parallel))

as_int <- function(x, default) {
  x <- Sys.getenv(x, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  as.integer(x)
}

parse_int_vec <- function(env_name, default = integer()) {
  x <- Sys.getenv(env_name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

parse_num_vec <- function(env_name, default = numeric()) {
  x <- Sys.getenv(env_name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

normalize_msg <- function(msg) {
  if (is.null(msg) || length(msg) == 0 || all(is.na(msg)) || !nzchar(paste(msg, collapse = ""))) {
    return("EMPTY_ERROR_MSG")
  }
  msg <- paste(msg, collapse = " | ")
  msg <- gsub("[[:space:]]+", " ", msg)
  trimws(msg)
}

base_dir <- Sys.getenv("FASTQR_TVCQR_BASE_DIR", unset = "data/tvcqr_simu_results")
partial_root <- file.path(base_dir, ".array_tmp")
cases <- parse_int_vec("FASTQR_CASES", c(1L))
ns <- parse_int_vec("FASTQR_NS", c(200L, 500L, 1000L, 2000L, 5000L))
taus <- parse_num_vec("FASTQR_TAUS", c(0.2, 0.5, 0.8))
num_rep <- as_int("FASTQR_NUM_REP", 1000L)
scan_cores <- as_int("FASTQR_SCAN_CORES", max(1L, min(8L, detectCores(logical = TRUE))))
out_prefix <- Sys.getenv("FASTQR_TVCQR_SCAN_PREFIX", unset = "results/table/tvcqr_current")
rep_id_dir <- sprintf("%s_rep_ids", out_prefix)

scan_one_config <- function(case_id, tau, n) {
  tau_str <- sprintf("tau%02d", as.integer(round(100 * tau)))
  partial_dir <- file.path(partial_root, sprintf("case%d_%s_n%d_rep%d", case_id, tau_str, n, num_rep))
  final_file <- file.path(base_dir, sprintf("case%d_%s_n%d_rep%d.RData", case_id, tau_str, n, num_rep))

  if (!dir.exists(partial_dir)) {
    summary_df <- data.frame(
      case = case_id,
      tau = tau,
      n = n,
      final_exists = file.exists(final_file),
      n_files = 0L,
      n_missing = num_rep,
      missing_rate = 1,
      n_failed = NA_integer_,
      failed_rate = NA_real_,
      first_missing = paste(head(seq_len(num_rep), 20L), collapse = ","),
      first_failed = NA_character_,
      first_failed_msg = NA_character_,
      n_rerun = num_rep,
      rep_ids_file = "",
      stringsAsFactors = FALSE
    )
    dir.create(rep_id_dir, recursive = TRUE, showWarnings = FALSE)
    rep_id_file <- file.path(rep_id_dir, sprintf("case%d_%s_n%d_rep%d.txt", case_id, tau_str, n, num_rep))
    writeLines(as.character(seq_len(num_rep)), rep_id_file)
    summary_df$rep_ids_file <- rep_id_file
    return(list(summary = summary_df, errors = data.frame()))
  }

  fs <- list.files(partial_dir, pattern = "^rep[0-9]+\\.RData$", full.names = TRUE)
  ids <- as.integer(sub("^rep([0-9]+)\\.RData$", "\\1", basename(fs)))
  ids <- ids[!is.na(ids)]
  missing <- setdiff(seq_len(num_rep), ids)

  check_file <- function(f) {
    env <- new.env()
    ok <- try(load(f, envir = env), silent = TRUE)
    rep_id <- as.integer(sub("^rep([0-9]+)\\.RData$", "\\1", basename(f)))

    if (inherits(ok, "try-error") || !exists("partial_results", envir = env)) {
      return(list(rep_id = rep_id, failed = TRUE, msg = "load failed / partial_results missing"))
    }

    pr <- get("partial_results", envir = env)
    if (is.null(pr$status) || pr$status != "success") {
      return(list(rep_id = rep_id, failed = TRUE, msg = normalize_msg(pr$error_msg)))
    }

    list(rep_id = rep_id, failed = FALSE, msg = "")
  }

  file_results <- if (length(fs) > 0L) {
    mclapply(fs, check_file, mc.cores = min(scan_cores, length(fs)))
  } else {
    list()
  }

  failed_ids <- sort(unique(vapply(file_results, function(x) if (isTRUE(x$failed)) x$rep_id else NA_integer_, integer(1L))))
  failed_ids <- failed_ids[!is.na(failed_ids)]
  failed_entries <- Filter(function(x) isTRUE(x$failed), file_results)
  failed_msgs <- vapply(failed_entries, function(x) x$msg, character(1L))
  rerun_ids <- sort(unique(c(missing, failed_ids)))

  err_df <- if (length(failed_entries) > 0L) {
    data.frame(
      case = case_id,
      tau = tau,
      n = n,
      rep_id = vapply(failed_entries, function(x) x$rep_id, integer(1L)),
      error_msg = failed_msgs,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }

  summary_df <- data.frame(
    case = case_id,
    tau = tau,
    n = n,
    final_exists = file.exists(final_file),
    n_files = length(fs),
    n_missing = length(missing),
    missing_rate = length(missing) / num_rep,
    n_failed = length(failed_ids),
    failed_rate = length(failed_ids) / num_rep,
    first_missing = if (length(missing)) paste(head(sort(missing), 20L), collapse = ",") else "",
    first_failed = if (length(failed_ids)) paste(head(failed_ids, 20L), collapse = ",") else "",
    first_failed_msg = if (length(failed_msgs)) failed_msgs[1L] else "",
    n_rerun = length(rerun_ids),
    rep_ids_file = "",
    stringsAsFactors = FALSE
  )

  if (length(rerun_ids) > 0L) {
    dir.create(rep_id_dir, recursive = TRUE, showWarnings = FALSE)
    rep_id_file <- file.path(rep_id_dir, sprintf("case%d_%s_n%d_rep%d.txt", case_id, tau_str, n, num_rep))
    writeLines(as.character(rerun_ids), rep_id_file)
    summary_df$rep_ids_file <- rep_id_file
  }

  list(summary = summary_df, errors = err_df)
}

grid <- expand.grid(case = cases, tau = taus, n = ns, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
results <- mcmapply(
  FUN = scan_one_config,
  case_id = grid$case,
  tau = grid$tau,
  n = grid$n,
  SIMPLIFY = FALSE,
  mc.cores = min(scan_cores, nrow(grid))
)

summary_df <- do.call(rbind, lapply(results, `[[`, "summary"))
summary_df <- summary_df[order(summary_df$case, summary_df$n, summary_df$tau), ]
error_list <- lapply(results, `[[`, "errors")
error_list <- error_list[vapply(error_list, nrow, integer(1L)) > 0L]
error_df <- if (length(error_list) > 0L) do.call(rbind, error_list) else data.frame()
problem_df <- subset(summary_df, n_missing > 0 | n_failed > 0 | !final_exists)

problem_tsv <- sprintf("%s_problem_configs.tsv", out_prefix)
summary_csv <- sprintf("%s_integrity_summary.csv", out_prefix)
errors_csv <- sprintf("%s_failed_reps.csv", out_prefix)
dir.create(dirname(summary_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, summary_csv, row.names = FALSE)
if (nrow(error_df) > 0L) {
  write.csv(error_df, errors_csv, row.names = FALSE)
}
write.table(
  problem_df[, c("case", "tau", "n", "n_missing", "n_failed", "missing_rate", "failed_rate", "n_rerun", "rep_ids_file")],
  file = problem_tsv,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("=== TVCQR integrity summary ===\n")
print(summary_df, row.names = FALSE)
cat("\n=== Problem configs ===\n")
print(problem_df, row.names = FALSE)
cat("\nSaved:\n")
cat(summary_csv, "\n")
if (nrow(error_df) > 0L) cat(errors_csv, "\n")
cat(problem_tsv, "\n")
