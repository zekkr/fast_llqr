#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

as_int <- function(name, default = NA_integer_) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  as.integer(value)
}

project_dir <- normalizePath(Sys.getenv("FASTQR_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
setwd(project_dir)

run_tag <- Sys.getenv("FASTQR_RUN_TAG", unset = "")
base_dir <- Sys.getenv("FASTQR_V4142_BASE_DIR", unset = "data/v41_v42_rep100")
num_rep <- as_int("FASTQR_NUM_REP", 100L)
seed_base <- as_int("FASTQR_SEED_BASE", 2026L)
if (!nzchar(run_tag)) stop("FASTQR_RUN_TAG is required.")

models <- c("llqr", "tvcqr")
cases <- c(1L, 2L)
taus <- c(0.2, 0.5, 0.8)
ns <- c(1000L, 2000L, 5000L, 10000L)
methods <- c("direct_baseline", "lean_seq", "v41", "v42")
run_dir <- file.path(base_dir, run_tag)
table_dir <- file.path(run_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

all_metrics <- list()
all_integrity <- list()
missing_configs <- character()

for (model in models) {
  for (case_id in cases) {
    for (tau in taus) {
      for (n in ns) {
        config_tag <- sprintf("case%d_tau%02d_n%d", case_id, as.integer(round(100 * tau)), n)
        config_dir <- file.path(run_dir, model, config_tag)
        metrics_path <- file.path(config_dir, "replication_metrics.rds")
        integrity_path <- file.path(config_dir, "integrity.csv")
        if (!file.exists(metrics_path) || !file.exists(integrity_path)) {
          missing_configs <- c(missing_configs, paste(model, config_tag, sep = "/"))
          next
        }
        value <- readRDS(metrics_path)
        if (nrow(value)) all_metrics[[length(all_metrics) + 1L]] <- value
        all_integrity[[length(all_integrity) + 1L]] <- read.csv(
          integrity_path, stringsAsFactors = FALSE, check.names = FALSE
        )
      }
    }
  }
}

metrics <- if (length(all_metrics)) do.call(rbind, all_metrics) else data.frame()
integrity <- if (length(all_integrity)) do.call(rbind, all_integrity) else data.frame()

summarize_method <- function(model, case_id, tau, n, method) {
  x <- metrics[
    metrics$model == model & metrics$case == case_id &
      metrics$tau == tau & metrics$n == n & metrics$method == method,
    , drop = FALSE
  ]
  x <- x[order(x$rep_id), , drop = FALSE]
  complete_ids <- nrow(x) == num_rep && identical(as.integer(x$rep_id), seq_len(num_rep))
  success <- complete_ids && all(x$method_ok) && all(x$finite_estimate)
  timing_ok <- success && all(is.finite(x$elapsed_sec))
  discrepancy_ok <- success && all(x$discrepancy_status == "ok") && all(is.finite(x$discrepancy))

  data.frame(
    run_tag = run_tag,
    model = model,
    case = case_id,
    tau = tau,
    n = n,
    num_rep = num_rep,
    seed_base = seed_base,
    method = method,
    version_status = if (method == "v41") {
      "research_only_active_first"
    } else if (method == "v42") {
      "paper_contract_aligned"
    } else {
      "reference"
    },
    n_rows = nrow(x),
    n_success = if (nrow(x)) sum(x$method_ok & x$finite_estimate) else 0L,
    n_failed = num_rep - if (nrow(x)) sum(x$method_ok & x$finite_estimate) else 0L,
    n_nonfinite = if (nrow(x)) sum(!x$finite_estimate) else num_rep,
    ierr_nonzero_count = if (nrow(x)) sum(!is.na(x$ierr) & x$ierr != 0L) else NA_integer_,
    failed_eval_nonzero_count = if (nrow(x)) sum(!is.na(x$failed_eval) & x$failed_eval != 0L) else NA_integer_,
    fallback_count = if (nrow(x)) sum(x$fallback_triggered) else 0L,
    first_point_full_m_count = if (nrow(x) && method %in% c("v41", "v42")) {
      sum(x$first_point_full_m %in% TRUE)
    } else {
      NA_integer_
    },
    independent_init_count = if (nrow(x) && method %in% c("v41", "v42") && all(is.finite(x$independent_init_count))) {
      sum(x$independent_init_count)
    } else {
      NA_integer_
    },
    full_active_recovery_count = if (nrow(x) && method %in% c("v41", "v42") && all(is.finite(x$full_active_recovery_count))) {
      sum(x$full_active_recovery_count)
    } else {
      NA_integer_
    },
    repair_total = if (nrow(x) && method %in% c("v41", "v42") && all(is.finite(x$repair_total))) {
      sum(x$repair_total)
    } else {
      NA_integer_
    },
    time_mean_sec = if (timing_ok) mean(x$elapsed_sec) else NA_real_,
    time_median_sec = if (timing_ok) median(x$elapsed_sec) else NA_real_,
    time_min_sec = if (timing_ok) min(x$elapsed_sec) else NA_real_,
    time_max_sec = if (timing_ok) max(x$elapsed_sec) else NA_real_,
    time_status = if (timing_ok) "ok" else "invalid_incomplete_or_failed",
    max_average_relative_bias = if (discrepancy_ok) max(x$discrepancy) else NA_real_,
    max_mean_relative_numerical_discrepancy = if (discrepancy_ok) max(x$discrepancy) else NA_real_,
    discrepancy_status = if (discrepancy_ok) "ok" else "invalid_incomplete_or_nonfinite",
    stringsAsFactors = FALSE
  )
}

summary_rows <- list()
for (model in models) {
  for (case_id in cases) {
    for (tau in taus) {
      for (n in ns) {
        for (method in methods) {
          summary_rows[[length(summary_rows) + 1L]] <- summarize_method(model, case_id, tau, n, method)
        }
      }
    }
  }
}
method_summary <- do.call(rbind, summary_rows)

ratio_rows <- list()
for (model in models) {
  for (case_id in cases) {
    for (tau in taus) {
      for (n in ns) {
        config <- method_summary[
          method_summary$model == model & method_summary$case == case_id &
            method_summary$tau == tau & method_summary$n == n,
          , drop = FALSE
        ]
        lean_time <- config$time_mean_sec[match("lean_seq", config$method)]
        for (candidate in c("v41", "v42")) {
          row <- config[config$method == candidate, , drop = FALSE]
          ratio <- if (length(lean_time) == 1L && is.finite(lean_time) &&
                       nrow(row) == 1L && is.finite(row$time_mean_sec)) {
            row$time_mean_sec / lean_time
          } else {
            NA_real_
          }
          ratio_rows[[length(ratio_rows) + 1L]] <- data.frame(
            model = model,
            case = case_id,
            tau = tau,
            n = n,
            candidate = candidate,
            candidate_time_mean_sec = row$time_mean_sec,
            lean_seq_time_mean_sec = lean_time,
            candidate_over_lean_seq = ratio,
            faster_than_lean_seq = if (is.finite(ratio)) ratio < 1 else NA,
            max_average_relative_bias_vs_direct = row$max_average_relative_bias,
            n_success = row$n_success,
            n_failed = row$n_failed,
            ierr_nonzero_count = row$ierr_nonzero_count,
            failed_eval_nonzero_count = row$failed_eval_nonzero_count,
            full_active_recovery_count = row$full_active_recovery_count,
            fallback_count = row$fallback_count,
            version_status = row$version_status,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}
candidate_summary <- do.call(rbind, ratio_rows)

write.csv(method_summary, file.path(table_dir, "method_summary.csv"), row.names = FALSE, na = "")
write.csv(candidate_summary, file.path(table_dir, "candidate_vs_lean_seq.csv"), row.names = FALSE, na = "")
write.csv(integrity, file.path(table_dir, "integrity_summary.csv"), row.names = FALSE, na = "")
if (nrow(metrics)) {
  write.csv(metrics, file.path(table_dir, "replication_metrics_all.csv"), row.names = FALSE, na = "")
}

all_complete <- !length(missing_configs) && nrow(integrity) == 48L &&
  all(integrity$all_replications_present %in% TRUE)
all_stable <- all(method_summary$n_success == num_rep) &&
  all(method_summary$n_failed == 0L) &&
  all(method_summary$n_nonfinite == 0L)
candidate_stable <- all(
  candidate_summary$n_success == num_rep &
    candidate_summary$n_failed == 0L &
    candidate_summary$ierr_nonzero_count == 0L &
    candidate_summary$failed_eval_nonzero_count == 0L &
    candidate_summary$fallback_count == 0L
)
faster_count <- sum(candidate_summary$faster_than_lean_seq %in% TRUE)
valid_ratio_count <- sum(is.finite(candidate_summary$candidate_over_lean_seq))
valid_bias <- candidate_summary[is.finite(candidate_summary$max_average_relative_bias_vs_direct), , drop = FALSE]
max_bias_row <- if (nrow(valid_bias)) valid_bias[which.max(valid_bias$max_average_relative_bias_vs_direct), , drop = FALSE] else NULL

summary_lines <- c(
  sprintf("run_tag: %s", run_tag),
  sprintf("grid_complete: %s", all_complete),
  sprintf("all_methods_stable: %s", all_stable),
  sprintf("v41_v42_stable: %s", candidate_stable),
  sprintf("missing_config_count: %d", length(missing_configs)),
  sprintf("candidate_faster_than_lean: %d/%d valid comparisons", faster_count, valid_ratio_count),
  if (is.null(max_bias_row)) {
    "max_candidate_bias: unavailable"
  } else {
    sprintf(
      "max_candidate_bias: %.17g at %s case=%d tau=%.1f n=%d method=%s",
      max_bias_row$max_average_relative_bias_vs_direct,
      max_bias_row$model, max_bias_row$case, max_bias_row$tau,
      max_bias_row$n, max_bias_row$candidate
    )
  },
  if (length(missing_configs)) {
    paste0("missing_configs: ", paste(missing_configs, collapse = ","))
  } else {
    "missing_configs: none"
  }
)
writeLines(summary_lines, file.path(table_dir, "run_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

if (!all_complete || !all_stable || !candidate_stable) quit(save = "no", status = 2L)
