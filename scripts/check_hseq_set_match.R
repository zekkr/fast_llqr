#!/usr/bin/env Rscript

# Batch check: whether H_seq rows of ppro methods match baseline seq methods
# in set sense (row-wise, order-insensitive).
#
# Usage examples:
#   Rscript scripts/check_hseq_set_match.R
#   Rscript scripts/check_hseq_set_match.R --models=llqr,tvcqr --rep=1000
#   Rscript scripts/check_hseq_set_match.R --models=llqr --cases=1,2 --taus=0.2,0.5,0.8 --ns=200,500,1000,2000,5000 --rep=1000

suppressWarnings({
  options(stringsAsFactors = FALSE)
})

parse_args <- function(args) {
  cfg <- list(
    models = c("llqr", "tvcqr"),
    cases = NULL,
    taus = NULL,
    ns = NULL,
    rep = 1000L,
    output_dir = "results/table",
    output_file = NULL
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
    if (key == "output_dir") cfg$output_dir <- val
    if (key == "output_file") cfg$output_file <- val
  }

  cfg$models <- unique(cfg$models)
  cfg$models <- cfg$models[cfg$models %in% c("llqr", "tvcqr")]
  if (length(cfg$models) == 0) stop("No valid models. Use llqr and/or tvcqr.")
  if (is.na(cfg$rep) || cfg$rep <= 0) stop("--rep must be a positive integer.")

  if (is.null(cfg$output_file)) {
    cfg$output_file <- sprintf("hseq_set_match_rep%d.csv", cfg$rep)
  }

  cfg
}

parse_file_meta <- function(path) {
  bn <- basename(path)
  m <- regexec("^case([0-9]+)_tau([0-9]+)_n([0-9]+)_rep([0-9]+)\\.RData$", bn)
  g <- regmatches(bn, m)[[1]]
  if (length(g) != 5) return(NULL)
  list(
    case = as.integer(g[2]),
    tau = as.integer(g[3]) / 100,
    n = as.integer(g[4]),
    rep = as.integer(g[5])
  )
}

row_set_equal <- function(a, b) {
  a <- as.integer(a)
  b <- as.integer(b)
  if (anyNA(a) || anyNA(b)) return(FALSE)
  setequal(a, b)
}

compare_hseq_rep <- function(base_h, target_h) {
  if (is.null(base_h)) return("missing_base_hseq")
  if (is.null(target_h)) return("missing_target_hseq")

  base_m <- as.matrix(base_h)
  target_m <- as.matrix(target_h)

  if (nrow(base_m) != nrow(target_m) || ncol(base_m) != ncol(target_m)) {
    return("shape_mismatch")
  }

  for (i in seq_len(nrow(base_m))) {
    if (!row_set_equal(base_m[i, ], target_m[i, ])) {
      return("row_set_mismatch")
    }
  }
  "ok"
}

check_one_method <- function(results_obj, base_method, target_method, rep_limit) {
  if (is.null(results_obj$H_seq_list)) {
    return(data.frame(
      target_method = target_method,
      baseline_method = base_method,
      n_rep_checked = 0L,
      n_failed = NA_integer_,
      n_missing_base_hseq = NA_integer_,
      n_missing_target_hseq = NA_integer_,
      n_shape_mismatch = NA_integer_,
      n_row_set_mismatch = NA_integer_,
      first_failed_rep = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  base_list <- results_obj$H_seq_list[[base_method]]
  target_list <- results_obj$H_seq_list[[target_method]]

  if (is.null(base_list) || is.null(target_list)) {
    return(data.frame(
      target_method = target_method,
      baseline_method = base_method,
      n_rep_checked = 0L,
      n_failed = NA_integer_,
      n_missing_base_hseq = NA_integer_,
      n_missing_target_hseq = NA_integer_,
      n_shape_mismatch = NA_integer_,
      n_row_set_mismatch = NA_integer_,
      first_failed_rep = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  n_rep <- min(length(base_list), length(target_list), rep_limit)
  if (n_rep <= 0) {
    return(data.frame(
      target_method = target_method,
      baseline_method = base_method,
      n_rep_checked = 0L,
      n_failed = 0L,
      n_missing_base_hseq = 0L,
      n_missing_target_hseq = 0L,
      n_shape_mismatch = 0L,
      n_row_set_mismatch = 0L,
      first_failed_rep = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  reasons <- character(n_rep)
  for (i in seq_len(n_rep)) {
    reasons[i] <- compare_hseq_rep(base_list[[i]], target_list[[i]])
  }

  failed_idx <- which(reasons != "ok")
  n_failed <- length(failed_idx)

  data.frame(
    target_method = target_method,
    baseline_method = base_method,
    n_rep_checked = n_rep,
    n_failed = n_failed,
    n_missing_base_hseq = sum(reasons == "missing_base_hseq"),
    n_missing_target_hseq = sum(reasons == "missing_target_hseq"),
    n_shape_mismatch = sum(reasons == "shape_mismatch"),
    n_row_set_mismatch = sum(reasons == "row_set_mismatch"),
    first_failed_rep = if (n_failed > 0) failed_idx[1] else NA_integer_,
    stringsAsFactors = FALSE
  )
}

select_methods <- function(model, method_names) {
  if (model == "llqr") {
    base_method <- "llqr_seq"
    target_methods <- grep("^llqr_seq_ppro(_fortran)?_[0-9]+$", method_names, value = TRUE)
  } else {
    base_method <- "tvcqr_seq"
    target_methods <- grep("^tvcqr_seq_ppro(_fortran)?_[0-9]+$", method_names, value = TRUE)
  }
  list(base = base_method, targets = target_methods)
}

check_model_files <- function(model, cfg) {
  data_dir <- if (model == "llqr") "data/llqr_simu_results" else "data/tvcqr_simu_results"
  files <- list.files(
    data_dir,
    pattern = sprintf("^case[0-9]+_tau[0-9]+_n[0-9]+_rep%d\\.RData$", cfg$rep),
    full.names = TRUE
  )
  if (length(files) == 0) {
    cat(sprintf("[WARN] No files found for model=%s in %s with rep=%d\n", model, data_dir, cfg$rep))
    return(data.frame())
  }

  rows <- list()
  idx <- 1L

  for (f in files) {
    meta <- parse_file_meta(f)
    if (is.null(meta)) next
    if (!is.null(cfg$cases) && !meta$case %in% cfg$cases) next
    if (!is.null(cfg$taus) && !meta$tau %in% cfg$taus) next
    if (!is.null(cfg$ns) && !meta$n %in% cfg$ns) next

    e <- new.env()
    ok <- try(load(f, envir = e), silent = TRUE)
    if (inherits(ok, "try-error") || !exists("results", envir = e)) {
      cat(sprintf("[WARN] Failed to load results from %s\n", f))
      next
    }
    results_obj <- get("results", envir = e)

    method_names <- results_obj$method_names
    if (is.null(method_names) && !is.null(results_obj$H_seq_list)) {
      method_names <- names(results_obj$H_seq_list)
    }
    if (is.null(method_names)) {
      cat(sprintf("[WARN] No method names in %s\n", f))
      next
    }

    sel <- select_methods(model, method_names)
    if (!sel$base %in% method_names) {
      cat(sprintf("[WARN] Baseline method %s missing in %s\n", sel$base, f))
      next
    }
    if (length(sel$targets) == 0) {
      cat(sprintf("[WARN] No ppro target methods found in %s\n", f))
      next
    }

    for (tm in sel$targets) {
      one <- check_one_method(results_obj, sel$base, tm, rep_limit = cfg$rep)
      one$model <- model
      one$case <- meta$case
      one$tau <- meta$tau
      one$n <- meta$n
      one$rep <- meta$rep
      one$file <- basename(f)
      one$failure_rate <- ifelse(one$n_rep_checked > 0, one$n_failed / one$n_rep_checked, NA_real_)

      one <- one[, c(
        "model", "case", "tau", "n", "rep", "target_method", "baseline_method",
        "n_rep_checked", "n_failed", "failure_rate",
        "n_missing_base_hseq", "n_missing_target_hseq", "n_shape_mismatch", "n_row_set_mismatch",
        "first_failed_rep", "file"
      )]

      rows[[idx]] <- one
      idx <- idx + 1L
    }
  }

  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

  out_rows <- list()
  pos <- 1L
  for (m in cfg$models) {
    cat(sprintf("\n[INFO] Checking model=%s\n", m))
    df <- check_model_files(m, cfg)
    if (nrow(df) > 0) {
      out_rows[[pos]] <- df
      pos <- pos + 1L
    }
  }

  if (length(out_rows) == 0) {
    cat("\n[INFO] No rows produced.\n")
    return(invisible(NULL))
  }

  out <- do.call(rbind, out_rows)
  out <- out[order(out$model, out$case, out$tau, out$n, out$target_method), ]
  rownames(out) <- NULL

  out_path <- file.path(cfg$output_dir, cfg$output_file)
  write.csv(out, out_path, row.names = FALSE)
  cat(sprintf("\n[INFO] Saved summary: %s\n", out_path))

  agg <- aggregate(
    cbind(n_failed, n_rep_checked) ~ model + target_method,
    data = out, FUN = sum, na.rm = TRUE
  )
  agg$failure_rate <- with(agg, ifelse(n_rep_checked > 0, n_failed / n_rep_checked, NA_real_))
  agg <- agg[order(-agg$n_failed, agg$model, agg$target_method), ]

  cat("\n[INFO] Failure counts by method (aggregated):\n")
  print(agg, row.names = FALSE)

  cat("\n[INFO] Rows with failures (n_failed > 0):\n")
  print(out[out$n_failed > 0, ], row.names = FALSE)
}

main()

