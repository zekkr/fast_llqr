#!/usr/bin/env Rscript

suppressWarnings({
  options(stringsAsFactors = FALSE)
})

parse_args <- function(args) {
  cfg <- list(
    models = c("llqr", "tvcqr"),
    cases = NULL,
    taus = NULL,
    ns = NULL,
    rep = NA_integer_,
    seed_base = NA_integer_,
    out = NULL,
    llqr_base_dir = Sys.getenv("FASTQR_LLQR_BASE_DIR", unset = "data/llqr_simu_results"),
    tvcqr_base_dir = Sys.getenv("FASTQR_TVCQR_BASE_DIR", unset = "data/tvcqr_simu_results")
  )

  parse_csv_chr <- function(x) {
    vals <- strsplit(x, ",", fixed = TRUE)[[1]]
    trimws(vals[nchar(trimws(vals)) > 0])
  }
  parse_csv_int <- function(x) as.integer(parse_csv_chr(x))
  parse_csv_num <- function(x) as.numeric(parse_csv_chr(x))

  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- if (length(kv) > 1) kv[2] else ""

    if (key == "models") cfg$models <- parse_csv_chr(val)
    if (key == "cases") cfg$cases <- parse_csv_int(val)
    if (key == "taus") cfg$taus <- parse_csv_num(val)
    if (key == "ns") cfg$ns <- parse_csv_int(val)
    if (key == "rep") cfg$rep <- as.integer(val)
    if (key == "seed-base") cfg$seed_base <- as.integer(val)
    if (key == "out") cfg$out <- val
    if (key == "llqr-base-dir") cfg$llqr_base_dir <- val
    if (key == "tvcqr-base-dir") cfg$tvcqr_base_dir <- val
  }

  cfg$models <- unique(cfg$models)
  cfg$models <- cfg$models[cfg$models %in% c("llqr", "tvcqr")]
  if (length(cfg$models) == 0) stop("No valid --models. Use llqr,tvcqr.")
  if (is.null(cfg$cases) || anyNA(cfg$cases)) stop("--cases is required.")
  if (is.null(cfg$taus) || anyNA(cfg$taus)) stop("--taus is required.")
  if (is.null(cfg$ns) || anyNA(cfg$ns)) stop("--ns is required.")
  if (is.na(cfg$rep) || cfg$rep <= 0) stop("--rep must be a positive integer.")
  if (is.na(cfg$seed_base)) stop("--seed-base is required.")
  if (is.null(cfg$out) || !nzchar(cfg$out)) {
    cfg$out <- sprintf("results/table/rep%d_seed%d_method_summary.csv",
                       cfg$rep, cfg$seed_base)
  }
  cfg
}

source_if_exists <- function(path) {
  if (!file.exists(path)) stop(sprintf("Required file not found: %s", path))
  source(path)
}

tau_tag <- function(tau) sprintf("tau%02d", as.integer(round(tau * 100)))

result_path <- function(model, case, tau, n, rep, llqr_base_dir, tvcqr_base_dir) {
  data_dir <- if (model == "llqr") {
    llqr_base_dir
  } else {
    tvcqr_base_dir
  }
  file.path(data_dir, sprintf("case%d_%s_n%d_rep%d.RData",
                              case, tau_tag(tau), n, rep))
}

partial_dir_path <- function(model, case, tau, n, rep, llqr_base_dir, tvcqr_base_dir) {
  data_dir <- if (model == "llqr") {
    llqr_base_dir
  } else {
    tvcqr_base_dir
  }
  file.path(data_dir, ".array_tmp",
            sprintf("case%d_%s_n%d_rep%d", case, tau_tag(tau), n, rep))
}

load_results_obj <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (exists("results", envir = env)) return(get("results", envir = env))
  if (exists("simulation_results", envir = env)) return(get("simulation_results", envir = env))
  stop(sprintf("No results object found in %s", path))
}

all_method_names <- function(results_obj) {
  methods <- character()
  if (!is.null(results_obj$method_names)) methods <- c(methods, results_obj$method_names)
  if (!is.null(results_obj$timing_matrix)) methods <- c(methods, colnames(results_obj$timing_matrix))
  if (!is.null(results_obj$estimates_list)) methods <- c(methods, names(results_obj$estimates_list))
  if (!is.null(results_obj$H_seq_list)) methods <- c(methods, names(results_obj$H_seq_list))
  unique(methods[nzchar(methods)])
}

display_method <- function(model, raw_method) {
  if (model == "tvcqr" && raw_method == "tvc_rq") return("tvcqr baseline")
  raw_method
}

reference_method <- function(model) {
  if (model == "llqr") "llqr" else "tvc_rq"
}

seq_baseline_method <- function(model) {
  if (model == "llqr") "llqr_seq" else "tvcqr_seq"
}

infer_mm_factor <- function(results_obj, raw_method) {
  mapping <- results_obj$Mm.factor_mapping
  if (!is.null(mapping) && all(c("Method", "Mm.factor") %in% names(mapping))) {
    hit <- mapping$Mm.factor[match(raw_method, mapping$Method)]
    if (length(hit) == 1L && !is.na(hit) && nzchar(hit) && hit != "N/A") {
      return(as.character(hit))
    }
  }

  mm_vec <- results_obj$config$Mm.factor
  if (is.null(mm_vec) || !grepl("ppro.*_(\\d+)$", raw_method)) return("N/A")
  idx <- as.integer(sub(".*_(\\d+)$", "\\1", raw_method))
  if (is.na(idx) || idx < 1L || idx > length(mm_vec)) return("N/A")
  as.character(mm_vec[idx])
}

finite_or_na <- function(x, fun) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  fun(x)
}

estimate_present <- function(x) {
  !is.null(x) && length(x) > 0L && !all(is.na(as.numeric(x)))
}

method_success_count <- function(results_obj, method, rep_limit) {
  timing <- results_obj$timing_matrix
  if (!is.null(timing) && method %in% colnames(timing)) {
    return(sum(is.finite(as.numeric(timing[seq_len(min(nrow(timing), rep_limit)), method]))))
  }
  est <- results_obj$estimates_list[[method]]
  if (!is.null(est)) {
    return(sum(vapply(est[seq_len(min(length(est), rep_limit))], estimate_present, logical(1))))
  }
  0L
}

inspect_partial_status <- function(partial_dir, rep_limit) {
  if (!dir.exists(partial_dir)) {
    return(list(
      source = "merged RData only",
      n_missing = NA_integer_,
      n_failed = NA_integer_,
      failed_reps = integer()
    ))
  }

  missing <- integer()
  failed <- integer()
  for (rep_id in seq_len(rep_limit)) {
    f <- file.path(partial_dir, sprintf("rep%04d.RData", rep_id))
    if (!file.exists(f)) {
      missing <- c(missing, rep_id)
      next
    }
    env <- new.env(parent = emptyenv())
    ok <- try(load(f, envir = env), silent = TRUE)
    if (inherits(ok, "try-error") || !exists("partial_results", envir = env)) {
      failed <- c(failed, rep_id)
      next
    }
    pr <- get("partial_results", envir = env)
    if (is.null(pr$status) || pr$status != "success") {
      failed <- c(failed, rep_id)
    }
  }
  list(
    source = "partial tmp directory",
    n_missing = length(missing),
    n_failed = length(failed),
    failed_reps = failed
  )
}

row_set_equal <- function(a, b) {
  a <- as.integer(a)
  b <- as.integer(b)
  if (anyNA(a) || anyNA(b)) return(FALSE)
  setequal(a, b)
}

compare_hseq_one_rep <- function(base_h, target_h) {
  if (is.null(base_h) || is.null(target_h)) return(FALSE)
  base_m <- as.matrix(base_h)
  target_m <- as.matrix(target_h)
  if (nrow(base_m) != nrow(target_m) || ncol(base_m) != ncol(target_m)) return(FALSE)
  if (anyNA(base_m) || anyNA(target_m)) return(FALSE)
  for (i in seq_len(nrow(base_m))) {
    if (!row_set_equal(base_m[i, ], target_m[i, ])) return(FALSE)
  }
  TRUE
}

hseq_summary <- function(results_obj, model, raw_method, rep_limit) {
  if (!grepl("ppro", raw_method)) {
    return(list(failed_count = NA_integer_, first_failed = NA_integer_,
                status = "not applicable"))
  }
  if (is.null(results_obj$H_seq_list)) {
    return(list(failed_count = NA_integer_, first_failed = NA_integer_,
                status = "not run; H_seq_list missing"))
  }
  base_method <- seq_baseline_method(model)
  base_list <- results_obj$H_seq_list[[base_method]]
  target_list <- results_obj$H_seq_list[[raw_method]]
  if (is.null(base_list) || is.null(target_list)) {
    return(list(failed_count = NA_integer_, first_failed = NA_integer_,
                status = sprintf("not run; H_seq missing for %s or %s",
                                 base_method, raw_method)))
  }

  n_rep <- min(length(base_list), length(target_list), rep_limit)
  failed <- integer()
  for (rep_id in seq_len(n_rep)) {
    if (!compare_hseq_one_rep(base_list[[rep_id]], target_list[[rep_id]])) {
      failed <- c(failed, rep_id)
    }
  }
  list(
    failed_count = length(failed),
    first_failed = if (length(failed) > 0L) failed[1] else NA_integer_,
    status = "checked against seq baseline"
  )
}

fallback_summary <- function(results_obj, raw_method, rep_limit) {
  metadata <- results_obj$method_metadata
  if (is.null(metadata)) {
    return(list(count = NA_integer_, rate = NA_real_,
                status = "not reliable; method_metadata not persisted"))
  }
  method_meta <- metadata[[raw_method]]
  if (is.null(method_meta) || !is.list(method_meta)) {
    return(list(count = NA_integer_, rate = NA_real_,
                status = "not reliable; method metadata missing"))
  }

  flags <- rep(NA, rep_limit)
  n <- min(length(method_meta), rep_limit)
  for (i in seq_len(n)) {
    item <- method_meta[[i]]
    if (is.list(item) && !is.null(item$fallback_triggered)) {
      flags[i] <- isTRUE(item$fallback_triggered)
    }
  }
  if (all(is.na(flags))) {
    return(list(count = NA_integer_, rate = NA_real_,
                status = "not reliable; fallback_triggered missing"))
  }
  count <- sum(flags, na.rm = TRUE)
  list(count = count, rate = count / sum(!is.na(flags)), status = "best effort from method_metadata")
}

bias_vector <- function(results_obj, model, raw_method) {
  if (model == "llqr") {
    calculate_average_relative_bias_llqr(results_obj, raw_method)
  } else {
    calculate_average_relative_bias_tvcqr(results_obj, raw_method)
  }
}

baseline_status <- function(results_obj, model, rep_limit) {
  ref <- reference_method(model)
  if (is.null(results_obj$estimates_list) || !ref %in% names(results_obj$estimates_list)) {
    return(list(ok = FALSE, status = "baseline_missing"))
  }
  n_success <- method_success_count(results_obj, ref, rep_limit)
  if (n_success == 0L) return(list(ok = FALSE, status = "baseline_failed"))
  if (n_success < rep_limit) {
    return(list(ok = FALSE, status = sprintf(
      "baseline_incomplete:%d/%d_success", n_success, rep_limit
    )))
  }
  list(ok = TRUE, status = "ok")
}

summarize_one_result <- function(model, case, tau, n, num_rep, seed_base, path,
                                 llqr_base_dir, tvcqr_base_dir) {
  results_obj <- load_results_obj(path)
  methods <- all_method_names(results_obj)
  if (length(methods) == 0L) stop(sprintf("No methods found in %s", path))

  partial_status <- inspect_partial_status(
    partial_dir_path(model, case, tau, n, num_rep, llqr_base_dir, tvcqr_base_dir),
    num_rep
  )
  base_status <- baseline_status(results_obj, model, num_rep)

  rows <- vector("list", length(methods))
  for (i in seq_along(methods)) {
    method <- methods[i]
    timing_col <- if (!is.null(results_obj$timing_matrix) &&
                      method %in% colnames(results_obj$timing_matrix)) {
      as.numeric(results_obj$timing_matrix[, method])
    } else {
      rep(NA_real_, num_rep)
    }

    n_success <- method_success_count(results_obj, method, num_rep)
    inferred_missing <- max(num_rep - n_success, 0L)
    n_missing <- if (!is.na(partial_status$n_missing)) partial_status$n_missing else inferred_missing
    n_failed <- if (!is.na(partial_status$n_failed)) partial_status$n_failed else inferred_missing
    notes <- if (!is.na(partial_status$n_missing)) {
      "missing/failed from partial tmp directory"
    } else {
      "missing/failed inferred from merged timing/estimates"
    }

    if (base_status$ok) {
      bias <- bias_vector(results_obj, model, method)
      bias_complete <- !is.null(bias) && length(bias) == num_rep && all(is.finite(bias))
      if (bias_complete) {
        max_bias <- max(bias)
        bias_status <- "ok"
        invalid_bias <- FALSE
      } else {
        max_bias <- NA_real_
        n_finite_bias <- if (is.null(bias)) 0L else sum(is.finite(bias))
        bias_status <- sprintf("invalid_bias_vector:%d/%d_finite", n_finite_bias, num_rep)
        invalid_bias <- TRUE
      }
    } else {
      max_bias <- NA_real_
      bias_status <- base_status$status
      invalid_bias <- TRUE
    }

    fb <- fallback_summary(results_obj, method, num_rep)
    hs <- hseq_summary(results_obj, model, method, num_rep)

    rows[[i]] <- data.frame(
      model = model,
      case = case,
      tau = tau,
      n = n,
      num_rep = num_rep,
      seed_base = seed_base,
      raw_method = method,
      display_method = display_method(model, method),
      Mm.factor = infer_mm_factor(results_obj, method),
      n_success = n_success,
      n_failed = n_failed,
      n_missing = n_missing,
      time_mean_sec = finite_or_na(timing_col, mean),
      time_median_sec = finite_or_na(timing_col, median),
      time_min_sec = finite_or_na(timing_col, min),
      time_max_sec = finite_or_na(timing_col, max),
      max_average_relative_bias = max_bias,
      bias_status = bias_status,
      config_invalid_for_bias = invalid_bias,
      fallback_count = fb$count,
      fallback_rate = fb$rate,
      fallback_status = fb$status,
      hseq_failed_count = hs$failed_count,
      hseq_first_failed_rep = hs$first_failed,
      hseq_status = hs$status,
      result_file = path,
      notes = notes,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  source_if_exists("R/performance_measurement.R")

  all_rows <- list()
  missing_files <- character()

  for (model in cfg$models) {
    for (case in cfg$cases) {
      for (tau in cfg$taus) {
        for (n in cfg$ns) {
          path <- result_path(
            model, case, tau, n, cfg$rep, cfg$llqr_base_dir, cfg$tvcqr_base_dir
          )
          if (!file.exists(path)) {
            missing_files <- c(missing_files, path)
            next
          }
          all_rows[[length(all_rows) + 1L]] <- summarize_one_result(
            model = model,
            case = case,
            tau = tau,
            n = n,
            num_rep = cfg$rep,
            seed_base = cfg$seed_base,
            path = path,
            llqr_base_dir = cfg$llqr_base_dir,
            tvcqr_base_dir = cfg$tvcqr_base_dir
          )
        }
      }
    }
  }

  if (length(all_rows) == 0L) {
    stop("No result files found for requested grid.")
  }

  summary_df <- do.call(rbind, all_rows)
  dir.create(dirname(cfg$out), recursive = TRUE, showWarnings = FALSE)
  write.csv(summary_df, cfg$out, row.names = FALSE)

  failed_configs <- unique(summary_df[
    summary_df$n_failed > 0 | summary_df$n_missing > 0,
    c("model", "case", "tau", "n")
  ])
  bias_rows <- summary_df[is.finite(summary_df$max_average_relative_bias), ]
  max_bias_row <- if (nrow(bias_rows) > 0L) {
    bias_rows[which.max(bias_rows$max_average_relative_bias), ]
  } else {
    NULL
  }

  cat("\n=== Method Summary ===\n")
  cat(sprintf("Output: %s\n", cfg$out))
  cat(sprintf("Rows: %d\n", nrow(summary_df)))
  cat(sprintf("Missing result files: %d\n", length(missing_files)))
  if (length(missing_files) > 0L) {
    cat("First missing files:\n")
    cat(paste(head(missing_files, 5), collapse = "\n"), "\n")
  }
  cat(sprintf("Failed/missing config-method rows: %d\n",
              sum(summary_df$n_failed > 0 | summary_df$n_missing > 0, na.rm = TRUE)))
  if (nrow(failed_configs) > 0L) {
    cat("Failed/missing configs:\n")
    print(failed_configs, row.names = FALSE)
  } else {
    cat("Failed/missing configs: none detected\n")
  }
  if (!is.null(max_bias_row)) {
    cat("Max average relative bias row:\n")
    print(max_bias_row[, c("model", "case", "tau", "n", "raw_method",
                           "Mm.factor", "max_average_relative_bias")],
          row.names = FALSE)
  } else {
    cat("Max average relative bias row: none available\n")
  }
  cat("=== Method Summary Done ===\n")
}

main()
