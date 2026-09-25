# Reproducible direct versus seq-MST experiment for Supplementary Table S1.
# Run from the repository root; this script uses only base R.
source("R/llqr_multivar_functions.R")

simulate_llqr_multivar_data <- function(n, p, case = 1L, seed = NULL) {
  if (p != 4L) stop("The Table S1 models require p = 4.")
  if (!is.null(seed)) set.seed(seed)
  if (case == 1L) {
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    f <- 1 + 1.2 * x[, 1L] - 0.9 * x[, 2L]
    eps <- 0.5 * rt(n, df = 4)
  } else if (case == 2L) {
    sigma <- outer(seq_len(p), seq_len(p), function(i, j) 0.5^abs(i - j))
    x <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol(sigma)
    hetero <- 0.4 + 0.25 * abs(x[, 1L]) + 0.1 * abs(x[, 2L])
    f <- 0.5 + 0.8 * x[, 1L] + 0.6 * x[, 2L]^2
    eps <- hetero * rnorm(n)
  } else stop("case must be 1 or 2.")
  colnames(x) <- paste0("x", seq_len(p))
  list(x = x, y = as.numeric(f + eps))
}

parse_options <- function(args) {
  out <- list(reps = 100L, output_dir = NULL, resume = FALSE,
              cases = c(1L, 2L), n_values = c(300L, 500L, 800L))
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--resume") out$resume <- TRUE else {
      if (i == length(args)) stop("Missing value after ", key)
      value <- args[i + 1L]
      if (key == "--reps") out$reps <- as.integer(value)
      else if (key == "--output-dir") out$output_dir <- value
      else if (key == "--cases") out$cases <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
      else if (key == "--n-values") out$n_values <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
      else stop("Unknown option: ", key)
      i <- i + 1L
    }
    i <- i + 1L
  }
  if (is.null(out$output_dir) || !nzchar(out$output_dir)) stop("--output-dir is required.")
  if (length(out$reps) != 1L || is.na(out$reps) || out$reps < 1L) stop("--reps must be positive.")
  if (length(out$cases) < 1L || anyNA(out$cases) ||
      any(!out$cases %in% 1:2) || anyDuplicated(out$cases)) stop("--cases must select 1, 2, or 1,2.")
  if (length(out$n_values) < 1L || anyNA(out$n_values) ||
      any(!out$n_values %in% c(300L, 500L, 800L)) ||
      anyDuplicated(out$n_values)) stop("--n-values must select distinct values from 300,500,800.")
  out
}

atomic_save_rds <- function(object, path) {
  tmp <- tempfile(pattern = ".incomplete_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) stop("Cannot move completed file into place: ", path)
}

fit_row <- function(fit, method, case_id, n, rep_id) {
  data.frame(case = case_id, n = n, p = 4L, replication = rep_id,
             method = method, total_iterations = sum(fit$it_num),
             mean_iterations = mean(fit$it_num),
             median_iterations = median(fit$it_num),
             max_iterations = max(fit$it_num), elapsed = unname(fit$elapsed),
             stringsAsFactors = FALSE)
}

run_pair <- function(config, case_id, n, rep_id, seed) {
  dat <- simulate_llqr_multivar_data(n, 4L, case_id, seed)
  x <- dat$x; y <- dat$y
  h <- compute_llqr_multivar_bandwidth(x, y, config$tau)
  cold_time <- system.time({
    cold <- llqr_tau_multivar(x, y, config$tau, z = x, h = h,
                              tol = config$tol, maxit = config$maxit,
                              bland = config$bland, track_order = TRUE)
  })
  cold$elapsed <- cold_time["elapsed"]
  mst_time <- system.time({
    mst <- llqr_tau_seq_multivar(x, y, config$tau, z = x, h = h,
                                 tol = config$tol, maxit = config$maxit,
                                 bland = config$bland, track_order = TRUE,
                                 order_method = "mst", root_method = "center",
                                 distance_scale = "bandwidth")
  })
  mst$elapsed <- mst_time["elapsed"]
  diff <- max(abs(as.numeric(mst$ll_est) - as.numeric(cold$ll_est)))
  diag <- data.frame(case = case_id, n = n, p = 4L, replication = rep_id,
                     seed = seed, h_summary = paste(signif(h, 4), collapse = ","),
                     mst_max_abs_diff = diff,
                     mst_mean_edge_weight = mean(mst$edge_weight[-mst$root]),
                     mst_max_edge_weight = max(mst$edge_weight[-mst$root]),
                     stringsAsFactors = FALSE)
  result <- list(case = case_id, n = n, p = 4L, replication = rep_id,
                 seed = seed, x = x, y = y, raw = rbind(
                   fit_row(cold, "cold", case_id, n, rep_id),
                   fit_row(mst, "seq_mst", case_id, n, rep_id)),
                 diagnostics = diag)
  failure <- if (any(!is.finite(c(cold$ll_est, mst$ll_est, diff)))) "nonfinite fit" else
    if (any(c(cold$it_num, mst$it_num) >= config$maxit)) "maxit reached" else
      if (diff > config$accuracy_tol) "accuracy tolerance exceeded" else NULL
  list(result = result, failure = failure)
}

validate_checkpoint <- function(item, case_id, n, rep_id, seed, config) {
  if (!is.list(item) || !identical(item$case, case_id) ||
      !identical(item$n, n) || !identical(item$p, 4L) ||
      !identical(item$replication, rep_id) || !identical(item$seed, seed) ||
      !is.matrix(item$x) || !identical(dim(item$x), c(n, 4L)) ||
      length(item$y) != n || !is.data.frame(item$raw) || nrow(item$raw) != 2L ||
      !identical(as.character(item$raw$method), c("cold", "seq_mst")) ||
      !is.data.frame(item$diagnostics) || nrow(item$diagnostics) != 1L ||
      any(item$raw$max_iterations >= config$maxit) ||
      !is.finite(item$diagnostics$mst_max_abs_diff) ||
      item$diagnostics$mst_max_abs_diff > config$accuracy_tol) {
    stop("Invalid or incomplete checkpoint: case=", case_id, " n=", n, " rep=", rep_id)
  }
  invisible(TRUE)
}

write_results <- function(items, config, output_dir) {
  raw <- do.call(rbind, lapply(items, `[[`, "raw"))
  diag <- do.call(rbind, lapply(items, `[[`, "diagnostics"))
  expected <- length(config$cases) * length(config$n_values) * config$reps
  keys <- paste(raw$case, raw$n, raw$replication, raw$method)
  if (length(items) != expected || nrow(raw) != 2L * expected ||
      nrow(diag) != expected || anyDuplicated(keys)) stop("Incomplete or duplicate run results.")
  cold <- raw[raw$method == "cold", c("case", "n", "replication", "total_iterations")]
  mst <- raw[raw$method == "seq_mst", c("case", "n", "replication", "total_iterations")]
  paired <- merge(cold, mst, by = c("case", "n", "replication"),
                  suffixes = c("_cold", "_mst"))
  paired$reduction <- 1 - paired$total_iterations_mst / paired$total_iterations_cold
  summary <- aggregate(cbind(total_iterations, elapsed) ~ case + n + method,
                       data = raw, FUN = mean)
  names(summary)[names(summary) == "total_iterations"] <- "mean_total_iterations"
  names(summary)[names(summary) == "elapsed"] <- "mean_elapsed"
  reduction <- aggregate(reduction ~ case + n, data = paired, FUN = mean)
  names(reduction)[3L] <- "mean_paired_reduction"
  table_s1 <- merge(merge(summary[summary$method == "cold", c("case", "n", "mean_total_iterations")],
                          summary[summary$method == "seq_mst", c("case", "n", "mean_total_iterations")],
                          by = c("case", "n"), suffixes = c("_direct", "_mst")),
                    reduction, by = c("case", "n"))
  table_s1 <- table_s1[order(table_s1$case, table_s1$n), ]
  table_s1$display_direct <- format(round(table_s1$mean_total_iterations_direct),
                                    big.mark = ",", scientific = FALSE, trim = TRUE)
  table_s1$display_mst <- format(round(table_s1$mean_total_iterations_mst),
                                 big.mark = ",", scientific = FALSE, trim = TRUE)
  table_s1$display_reduction <- sprintf("%.1f%%", 100 * table_s1$mean_paired_reduction)
  write.csv(raw, file.path(output_dir, "raw.csv"), row.names = FALSE)
  write.csv(diag, file.path(output_dir, "diagnostics.csv"), row.names = FALSE)
  write.csv(paired, file.path(output_dir, "paired_reductions.csv"), row.names = FALSE)
  write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
  write.csv(table_s1, file.path(output_dir, "table_s1.csv"), row.names = FALSE)
  atomic_save_rds(list(config = config, raw = raw, diagnostics = diag,
                       paired = paired, summary = summary, table_s1 = table_s1),
                  file.path(output_dir, "results.rds"))
  cat("Completed", expected, "paired datasets; max fit difference",
      format(max(diag$mst_max_abs_diff), digits = 7),
      "; mean parent-child edge", format(mean(diag$mst_mean_edge_weight), digits = 7), "\n")
  print(table_s1, row.names = FALSE)
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
config <- list(cases = options$cases, n_values = options$n_values, p = 4L,
               tau = 0.5, reps = options$reps, tol = 1e-14,
               maxit = 20000L, bland = FALSE, seed_base = 20260405L,
               accuracy_tol = 1e-8, methods = c("cold", "seq_mst"))
output_dir <- options$output_dir
source_paths <- c("scripts/run_llqr_multivar_iteration_simulation.R",
                  "R/llqr_multivar_functions.R")
source_hashes <- unname(tools::md5sum(source_paths))
names(source_hashes) <- source_paths
manifest_path <- file.path(output_dir, "manifest.rds")
if (options$resume) {
  if (!file.exists(manifest_path)) stop("No manifest to resume in: ", output_dir)
  manifest <- readRDS(manifest_path)
  if (!identical(manifest$config, config) ||
      !identical(manifest$source_hashes, source_hashes)) {
    stop("Resume rejected: configuration or source hashes changed.")
  }
} else {
  if (file.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory already contains files: ", output_dir)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "source"))
  for (path in source_paths) file.copy(path, file.path(output_dir, "source", basename(path)))
  manifest <- list(config = config, source_hashes = source_hashes,
                   rng_kind = RNGkind(), session_info = capture.output(sessionInfo()),
                   started_at = as.character(Sys.time()))
  atomic_save_rds(manifest, manifest_path)
  writeLines(manifest$session_info, file.path(output_dir, "session_info.txt"))
}
checkpoint_dir <- file.path(output_dir, "checkpoints")
failure_dir <- file.path(output_dir, "failures")
dir.create(checkpoint_dir, showWarnings = FALSE)
dir.create(failure_dir, showWarnings = FALSE)
items <- vector("list", length(config$cases) * length(config$n_values) * config$reps)
k <- 1L
for (case_id in config$cases) {
  for (n in config$n_values) {
    for (rep_id in seq_len(config$reps)) {
      seed <- config$seed_base + 100000L * case_id + 1000L * n + 10L * config$p + rep_id
      id <- sprintf("case%d_n%d_rep%03d", case_id, n, rep_id)
      checkpoint <- file.path(checkpoint_dir, paste0(id, ".rds"))
      if (file.exists(checkpoint)) {
        item <- readRDS(checkpoint)
        validate_checkpoint(item, case_id, n, rep_id, seed, config)
        cat("RESUME", id, "\n")
      } else {
        cat("RUN", id, "\n")
        pair <- run_pair(config, case_id, n, rep_id, seed)
        if (!is.null(pair$failure)) {
          atomic_save_rds(pair, file.path(failure_dir, paste0(id, ".rds")))
          stop("Pair failed at ", id, ": ", pair$failure)
        }
        item <- pair$result
        validate_checkpoint(item, case_id, n, rep_id, seed, config)
        atomic_save_rds(item, checkpoint)
        cat("DONE", id, " direct=", item$raw$total_iterations[1L],
            " mst=", item$raw$total_iterations[2L], "\n", sep = "")
      }
      items[[k]] <- item
      k <- k + 1L
    }
  }
}
write_results(items, config, output_dir)
