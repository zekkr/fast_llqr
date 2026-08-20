# scripts/array/driver_llqr_array.R
# Run LLQR simulation in SLURM job array (chunked), write partial results per replication.
# ---- project root ----
PROJECT_DIR <- normalizePath(Sys.getenv("FASTQR_PROJECT_DIR", unset = getwd()))
setwd(PROJECT_DIR)


suppressPackageStartupMessages({
  library(microbenchmark)
})

setup_file <- if (file.exists("scripts/setup_hpc.R")) {
  "scripts/setup_hpc.R"
} else {
  "scripts/setup.R"
}

source(setup_file)
source("scripts/audit_llqr_exact_zero.R")

as_int <- function(x, default) {
  x <- Sys.getenv(x, unset = NA_character_)
  if (is.na(x) || nchar(x) == 0) return(default)
  as.integer(x)
}
as_num <- function(x, default) {
  x <- Sys.getenv(x, unset = NA_character_)
  if (is.na(x) || nchar(x) == 0) return(default)
  as.numeric(x)
}
as_bool <- function(x, default = FALSE) {
  x <- Sys.getenv(x, unset = NA_character_)
  if (is.na(x) || nchar(x) == 0) return(default)
  tolower(x) %in% c("1","true","t","yes","y")
}
as_threshold_scale_mode <- function(x, default = "loglog") {
  value <- tolower(trimws(Sys.getenv(x, unset = default)))
  if (is.na(value) || nchar(value) == 0) return(default)
  if (!(value %in% c("loglog", "log"))) {
    stop(sprintf("%s must be loglog or log, got: %s", x, value))
  }
  value
}
as_timeout_fork_mode <- function(x, default = "auto") {
  x <- tolower(Sys.getenv(x, unset = default))
  if (is.na(x) || nchar(x) == 0) return(default)
  if (x %in% c("auto")) return("auto")
  if (x %in% c("1", "true", "t", "yes", "y", "on")) return("on")
  if (x %in% c("0", "false", "f", "no", "n", "off")) return("off")
  stop(sprintf("FASTQR_USE_TIMEOUT_FORK must be auto/on/off or boolean-like, got: %s", x))
}
as_parallel_backend <- function(x, default = "FORK") {
  x <- toupper(Sys.getenv(x, unset = default))
  if (is.na(x) || nchar(x) == 0) return(default)
  if (x %in% c("SEQ", "SEQUENTIAL", "NONE")) return("SEQ")
  if (x %in% c("FORK", "PSOCK")) return(x)
  stop(sprintf("FASTQR_PARALLEL_BACKEND must be FORK, PSOCK, or SEQ, got: %s", x))
}
as_num_vec <- function(x, default) {
  x <- Sys.getenv(x, unset = NA_character_)
  if (is.na(x) || nchar(x) == 0) return(default)
  vals <- strsplit(x, ",", fixed = TRUE)[[1]]
  as.numeric(trimws(vals))
}
as_chr_vec <- function(x, default = character()) {
  x <- Sys.getenv(x, unset = NA_character_)
  if (is.na(x) || nchar(x) == 0) return(default)
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  unique(vals[nzchar(vals)])
}
as_optional_pos_int <- function(x) {
  raw <- trimws(Sys.getenv(x, unset = NA_character_))
  if (is.na(raw) || nchar(raw) == 0 ||
      tolower(raw) %in% c("null", "na", "default")) {
    return(NULL)
  }
  value <- suppressWarnings(as.numeric(raw))
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < 1L || value != floor(value)) {
    stop(sprintf("%s must be empty/default or a positive integer, got: %s", x, raw))
  }
  as.integer(value)
}
parse_rep_id_pool <- function(num_rep) {
  rep_file <- Sys.getenv("FASTQR_REP_ID_FILE", unset = NA_character_)
  rep_list <- Sys.getenv("FASTQR_REP_ID_LIST", unset = NA_character_)

  tokens <- character()
  if (!is.na(rep_file) && nzchar(rep_file)) {
    if (!file.exists(rep_file)) {
      stop(sprintf("FASTQR_REP_ID_FILE does not exist: %s", rep_file))
    }
    file_text <- paste(readLines(rep_file, warn = FALSE), collapse = ",")
    tokens <- c(tokens, strsplit(file_text, "[,[:space:]]+", perl = TRUE)[[1]])
  }
  if (!is.na(rep_list) && nzchar(rep_list)) {
    tokens <- c(tokens, strsplit(rep_list, "[,[:space:]]+", perl = TRUE)[[1]])
  }
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) == 0L) {
    return(NULL)
  }

  vals <- suppressWarnings(as.integer(tokens))
  vals <- sort(unique(vals[!is.na(vals) & vals >= 1L & vals <= num_rep]))
  if (length(vals) == 0L) {
    stop("No valid rep_ids found in FASTQR_REP_ID_FILE / FASTQR_REP_ID_LIST.")
  }
  vals
}
require_pos_int <- function(x, name) {
  if (is.na(x) || !is.finite(x) || x < 1L) {
    stop(sprintf("%s must be a positive integer, got: %s", name, as.character(x)))
  }
}

task_id  <- as_int("SLURM_ARRAY_TASK_ID", 1)
ncores   <- as_int("SLURM_CPUS_PER_TASK", 1)

case     <- as_int("FASTQR_CASE", 1)
tau      <- as_num("FASTQR_TAU", 0.5)
n        <- as_int("FASTQR_N", 500)
num_rep  <- as_int("FASTQR_NUM_REP", 500)

chunk_size <- as_int("FASTQR_CHUNK_SIZE", ncores)
require_pos_int(task_id, "SLURM_ARRAY_TASK_ID")
require_pos_int(ncores, "SLURM_CPUS_PER_TASK")
require_pos_int(n, "FASTQR_N")
require_pos_int(num_rep, "FASTQR_NUM_REP")
require_pos_int(chunk_size, "FASTQR_CHUNK_SIZE")
if (!(case %in% c(1L, 2L))) {
  stop(sprintf("FASTQR_CASE must be 1 or 2, got: %s", as.character(case)))
}
if (is.na(tau) || !is.finite(tau) || tau <= 0 || tau >= 1) {
  stop(sprintf("FASTQR_TAU must be in (0,1), got: %s", as.character(tau)))
}

rep_id_pool <- parse_rep_id_pool(num_rep)
if (is.null(rep_id_pool)) {
  rep_start <- (task_id - 1L) * chunk_size + 1L
  rep_end   <- min(task_id * chunk_size, num_rep)
  rep_ids   <- rep_start:rep_end
  rep_ids   <- rep_ids[rep_ids >= 1 & rep_ids <= num_rep]
  sparse_mode <- FALSE
  total_rep_targets <- num_rep
} else {
  rep_start <- (task_id - 1L) * chunk_size + 1L
  rep_end   <- min(task_id * chunk_size, length(rep_id_pool))
  rep_ids   <- if (rep_start <= rep_end) rep_id_pool[rep_start:rep_end] else integer()
  sparse_mode <- TRUE
  total_rep_targets <- length(rep_id_pool)
}

Mm.factor <- as_num_vec("FASTQR_MM_FACTOR", c(1e-2, 1e-3, 1e-4))
h.factor  <- as_num("FASTQR_H_FACTOR", 1)
tol       <- as_num("FASTQR_TOL", 1e-14)
maxit     <- as_int("FASTQR_MAXIT", 2e6)
bland     <- as_bool("FASTQR_BLAND", FALSE)
track_order <- as_bool("FASTQR_TRACK_ORDER", TRUE)
seed_base <- as_int("FASTQR_SEED_BASE", 2026)
max_attempts_per_rep <- as_int("FASTQR_MAX_ATTEMPTS_PER_REP", 20)
retry_stride <- as_int("FASTQR_RETRY_STRIDE", 1000000)
max_seconds_per_rep <- as_int("FASTQR_MAX_SECONDS_PER_REP", 7200)
timeout_fork_mode <- as_timeout_fork_mode("FASTQR_USE_TIMEOUT_FORK", "auto")
parallel_backend <- as_parallel_backend("FASTQR_PARALLEL_BACKEND", "FORK")
save_h_seq <- as_bool("FASTQR_SAVE_H_SEQ", TRUE)
run_tag <- Sys.getenv("FASTQR_RUN_TAG", unset = "")
llqr_base_dir <- Sys.getenv("FASTQR_LLQR_BASE_DIR", unset = "data/llqr_simu_results")
exact_zero_audit <- as_bool("FASTQR_EXACT_ZERO_AUDIT", FALSE)
exact_zero_methods <- as_chr_vec("FASTQR_EXACT_ZERO_METHODS", character())
exact_zero_block_size <- as_int("FASTQR_EXACT_ZERO_BLOCK_SIZE", 128L)
min_subsample_size <- as_optional_pos_int("FASTQR_MIN_SUBSAMPLE_SIZE")
always_same_h_refit <- as_bool("FASTQR_ALWAYS_SAME_H_REFIT", TRUE)
threshold_lower_bound <- as_bool("FASTQR_THRESHOLD_LOWER_BOUND", TRUE)
threshold_scale_mode <- as_threshold_scale_mode("FASTQR_THRESHOLD_SCALE_MODE", "loglog")
require_pos_int(seed_base, "FASTQR_SEED_BASE")
require_pos_int(max_attempts_per_rep, "FASTQR_MAX_ATTEMPTS_PER_REP")
require_pos_int(retry_stride, "FASTQR_RETRY_STRIDE")
require_pos_int(max_seconds_per_rep, "FASTQR_MAX_SECONDS_PER_REP")
require_pos_int(exact_zero_block_size, "FASTQR_EXACT_ZERO_BLOCK_SIZE")
if (length(Mm.factor) == 0 || any(is.na(Mm.factor))) {
  stop("FASTQR_MM_FACTOR must be a comma-separated numeric list.")
}
if (!threshold_lower_bound && any(!is.finite(Mm.factor) | Mm.factor <= 0)) {
  stop("FASTQR_MM_FACTOR must contain only positive finite values when FASTQR_THRESHOLD_LOWER_BOUND is false.")
}
if (is.na(tol) || !is.finite(tol) || tol <= 0) {
  stop(sprintf("FASTQR_TOL must be > 0, got: %s", as.character(tol)))
}
if (is.na(h.factor) || !is.finite(h.factor) || h.factor <= 0) {
  stop(sprintf("FASTQR_H_FACTOR must be > 0, got: %s", as.character(h.factor)))
}
require_pos_int(maxit, "FASTQR_MAXIT")
if (.Platform$OS.type != "unix" && parallel_backend == "FORK") {
  parallel_backend <- "PSOCK"
}

parallel_chunk_enabled <- ncores > 1 && length(rep_ids) > 1 && parallel_backend != "SEQ"
cluster_type <- if (parallel_chunk_enabled) parallel_backend else "none"
timeout_fork_enabled <- switch(
  timeout_fork_mode,
  on = TRUE,
  off = FALSE,
  auto = !(parallel_chunk_enabled && cluster_type == "FORK")
)

config_base <- list(
  case = case,
  tau = tau,
  n = n,
  num_rep = num_rep,
  h = NULL,
  h.factor = h.factor,
  z = NULL,
  tol = tol,
  maxit = maxit,
  bland = bland,
  track_order = track_order,
  Mm.factor = Mm.factor,
  seed_base = seed_base,
  min_subsample_size = min_subsample_size,
  always_same_h_refit = always_same_h_refit,
  threshold_lower_bound = threshold_lower_bound,
  threshold_scale_mode = threshold_scale_mode,
  run_tag = run_tag,
  llqr_base_dir = llqr_base_dir,
  exact_zero_audit = exact_zero_audit,
  exact_zero_methods = exact_zero_methods,
  exact_zero_block_size = exact_zero_block_size
)

tau_str <- sprintf("tau%02d", as.integer(round(tau * 100)))
partial_dir <- file.path(llqr_base_dir, ".array_tmp",
                         sprintf("case%d_%s_n%d_rep%d", case, tau_str, n, num_rep))
dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== LLQR ARRAY DRIVER ===\n")
cat(sprintf("task_id=%d, ncores=%d, chunk_size=%d\n", task_id, ncores, chunk_size))
cat(sprintf("case=%d, tau=%.2f, n=%d, num_rep=%d\n", case, tau, n, num_rep))
cat(sprintf("project_dir=%s\n", PROJECT_DIR))
cat(sprintf("run_tag=%s\n", if (nzchar(run_tag)) run_tag else "<none>"))
cat(sprintf("llqr_base_dir=%s\n", llqr_base_dir))
if (sparse_mode) {
  cat(sprintf("sparse rep positions: %d-%d of %d (len=%d)\n", rep_start, rep_end, total_rep_targets, length(rep_ids)))
  cat(sprintf("rep ids: %s\n", paste(utils::head(rep_ids, 20L), collapse = ",")))
} else {
  cat(sprintf("rep range: %d-%d (len=%d)\n", rep_start, rep_end, length(rep_ids)))
}
cat(sprintf("partial_dir=%s\n", partial_dir))
cat("Mm.factor:", paste(Mm.factor, collapse = ", "), "\n")
cat("h.factor:", h.factor, "\n")
cat("min_subsample_size:", if (is.null(config_base$min_subsample_size)) "default" else config_base$min_subsample_size, "\n")
cat("always_same_h_refit:", always_same_h_refit, "\n")
cat("threshold_lower_bound:", threshold_lower_bound, "\n")
cat("threshold_scale_mode:", threshold_scale_mode, "\n")
cat("seed_base:", seed_base, "\n\n")
cat("max_attempts_per_rep:", max_attempts_per_rep, "\n")
cat("retry_stride:", retry_stride, "\n\n")
cat("max_seconds_per_rep:", max_seconds_per_rep, "\n\n")
cat("parallel_backend:", parallel_backend, "\n")
cat("parallel_chunk_enabled:", parallel_chunk_enabled, "\n")
cat("timeout_fork_mode:", timeout_fork_mode, "\n")
cat("timeout_fork_enabled:", timeout_fork_enabled, "\n\n")
cat("save_h_seq:", save_h_seq, "\n\n")
cat("exact_zero_audit:", exact_zero_audit, "\n")
cat("exact_zero_methods:", paste(exact_zero_methods, collapse = ","), "\n")
cat("exact_zero_block_size:", exact_zero_block_size, "\n\n")

if (length(rep_ids) == 0) {
  cat("No rep_ids assigned to this task. Exiting.\n")
  quit(save = "no", status = 0)
}

methods <- create_llqr_methods(config_base$Mm.factor)
method_names <- names(methods)
if (exact_zero_audit) {
  missing_audit_methods <- setdiff(exact_zero_methods, method_names)
  if (length(exact_zero_methods) == 0L) {
    stop("FASTQR_EXACT_ZERO_AUDIT=1 requires FASTQR_EXACT_ZERO_METHODS.")
  }
  if (length(missing_audit_methods)) {
    stop("Unknown FASTQR_EXACT_ZERO_METHODS: ", paste(missing_audit_methods, collapse = ", "))
  }
}

run_rep_with_timeout <- function(rep_id, rep_config, methods, timeout_sec, use_timeout_fork) {
  if (.Platform$OS.type != "unix" || !isTRUE(use_timeout_fork)) {
    setTimeLimit(elapsed = timeout_sec, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
    return(run_single_llqr_replication(rep_id, rep_config, methods))
  }

  child <- parallel::mcparallel({
    run_single_llqr_replication(rep_id, rep_config, methods)
  }, silent = TRUE)

  collected <- parallel::mccollect(child, wait = FALSE, timeout = timeout_sec)
  if (is.null(collected)) {
    try(parallel::mckill(child$pid, signal = 9L), silent = TRUE)
    try(parallel::mccollect(child, wait = FALSE, timeout = 0), silent = TRUE)
    stop(sprintf("rep exceeded %d seconds", timeout_sec))
  }

  rr <- collected[[1L]]
  if (inherits(rr, "try-error")) stop(as.character(rr))
  rr
}

run_one <- function(rep_id) {
  rep_config <- config_base
  rep_config$num_rep <- 1
  rep_config$rep_id <- rep_id

  out <- NULL
  for (attempt in seq_len(max_attempts_per_rep)) {
    seed_used <- as.integer(seed_base + rep_id + (attempt - 1L) * retry_stride)
    rep_config$seed_used <- seed_used

    out <- tryCatch({
      rr <- run_rep_with_timeout(rep_id, rep_config, methods, max_seconds_per_rep,
                                 timeout_fork_enabled)

      timing_matrix <- matrix(NA_real_, nrow = 1, ncol = length(method_names),
                              dimnames = list(NULL, method_names))
      for (m in method_names) timing_matrix[1, m] <- as.numeric(rr$timing[[m]])

      estimates_list <- setNames(vector("list", length(method_names)), method_names)
      H_seq_list     <- if (save_h_seq) setNames(vector("list", length(method_names)), method_names) else NULL
      method_metadata <- setNames(vector("list", length(method_names)), method_names)
      exact_zero_results <- setNames(vector("list", length(method_names)), method_names)
      for (m in method_names) {
        estimates_list[[m]] <- list(rr$estimates[[m]])
        if (save_h_seq) {
          H_seq_list[[m]] <- list(rr$H_seq[[m]])
        }
        method_metadata[[m]] <- rr$method_metadata[[m]]
        exact_zero_results[[m]] <- rr$exact_zero_audit[[m]]
      }

      llqr_meta <- rr$method_metadata[["llqr"]]

      partial_results <- list(
        config = rep_config,
        timing_matrix = timing_matrix,
        estimates_list = estimates_list,
        method_metadata = method_metadata,
        exact_zero_audit = exact_zero_results,
        method_names = method_names,
        Mm.factor_mapping = create_llqr_Mm_factor_mapping(method_names, config_base$Mm.factor),
        timestamp = Sys.time(),
        simulation_type = "llqr_partial",
        status = "success",
        error_msg = NULL,
        seed_used = seed_used,
        attempts = attempt,
        max_attempts_per_rep = max_attempts_per_rep,
        max_seconds_per_rep = max_seconds_per_rep,
        save_h_seq = save_h_seq,
        h_used = if (!is.null(llqr_meta$h_used)) llqr_meta$h_used else NA_real_,
        h_retry_factor = if (!is.null(llqr_meta$h_retry_factor)) llqr_meta$h_retry_factor else NA_real_,
        llqr_attempts = if (!is.null(llqr_meta$llqr_attempts)) llqr_meta$llqr_attempts else NA_integer_
      )
      if (save_h_seq) {
        partial_results$H_seq_list <- H_seq_list
      }

      list(ok = TRUE, obj = partial_results)
    }, error = function(e) {
      error_msg <- conditionMessage(e)
      retryable_baseline_error <- grepl("^baseline_error:", error_msg)
      partial_results <- list(
        config = rep_config,
        method_names = method_names,
        timestamp = Sys.time(),
        simulation_type = "llqr_partial",
        status = "error",
        error_msg = error_msg,
        seed_used = seed_used,
        attempts = attempt,
        max_attempts_per_rep = max_attempts_per_rep,
        max_seconds_per_rep = max_seconds_per_rep,
        save_h_seq = save_h_seq,
        retryable_baseline_error = retryable_baseline_error
      )
      list(ok = FALSE, obj = partial_results, retryable = retryable_baseline_error)
    })

    if (isTRUE(out$ok) || !isTRUE(out$retryable)) break
  }
  
  save_path <- file.path(partial_dir, sprintf("rep%04d.RData", rep_id))
  partial_results <- out$obj
  save(partial_results, file = save_path)
  cat(sprintf("[%s] rep=%d attempt=%d seed=%d -> %s (%s)\n",
              format(Sys.time(), "%F %T"), rep_id,
              partial_results$attempts, partial_results$seed_used,
              save_path,
              partial_results$status))
  invisible(TRUE)
}

# Parallelize within the chunk if we have >1 core and >1 rep
if (parallel_chunk_enabled) {
  suppressPackageStartupMessages({
    library(doParallel)
    library(foreach)
  })
  cat(sprintf("parallel cluster type: %s\n", cluster_type))
  cl <- if (cluster_type == "FORK") {
    parallel::makeForkCluster(ncores)
  } else {
    parallel::makeCluster(ncores)
  }
  on.exit(stopCluster(cl), add = TRUE)
  registerDoParallel(cl)
  
  clusterExport(cl, c("PROJECT_DIR", "setup_file"), envir = environment())
  
  clusterEvalQ(cl, {
    setwd(PROJECT_DIR)
    source(setup_file)
    source("scripts/audit_llqr_exact_zero.R")
    suppressPackageStartupMessages(library(microbenchmark))
    NULL
  })
  
  clusterExport(
    cl,
    varlist = c(
      "rep_ids", "run_one", "run_rep_with_timeout", "config_base", "methods",
      "method_names", "max_seconds_per_rep", "max_attempts_per_rep",
      "retry_stride", "seed_base", "partial_dir", "timeout_fork_enabled"
    ),
    envir = environment()
  )
  
  foreach(r = rep_ids) %dopar% {
    run_one(r)
    NULL
  }
} else {
  if (parallel_backend == "SEQ" && ncores > 1 && length(rep_ids) > 1) {
    cat("parallel chunk execution disabled by FASTQR_PARALLEL_BACKEND=SEQ\n")
  }
  for (r in rep_ids) run_one(r)
}


cat("\n=== LLQR ARRAY DRIVER DONE ===\n")
