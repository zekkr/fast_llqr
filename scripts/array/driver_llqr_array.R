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
as_num_vec <- function(x, default) {
  x <- Sys.getenv(x, unset = NA_character_)
  if (is.na(x) || nchar(x) == 0) return(default)
  vals <- strsplit(x, ",", fixed = TRUE)[[1]]
  as.numeric(trimws(vals))
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
tol       <- as_num("FASTQR_TOL", 1e-14)
maxit     <- as_int("FASTQR_MAXIT", 2e6)
bland     <- as_bool("FASTQR_BLAND", FALSE)
track_order <- as_bool("FASTQR_TRACK_ORDER", TRUE)
seed_base <- as_int("FASTQR_SEED_BASE", 2026)
max_attempts_per_rep <- as_int("FASTQR_MAX_ATTEMPTS_PER_REP", 20)
retry_stride <- as_int("FASTQR_RETRY_STRIDE", 1000000)
max_seconds_per_rep <- as_int("FASTQR_MAX_SECONDS_PER_REP", 7200)
require_pos_int(seed_base, "FASTQR_SEED_BASE")
require_pos_int(max_attempts_per_rep, "FASTQR_MAX_ATTEMPTS_PER_REP")
require_pos_int(retry_stride, "FASTQR_RETRY_STRIDE")
require_pos_int(max_seconds_per_rep, "FASTQR_MAX_SECONDS_PER_REP")
if (length(Mm.factor) == 0 || any(is.na(Mm.factor))) {
  stop("FASTQR_MM_FACTOR must be a comma-separated numeric list.")
}
if (is.na(tol) || !is.finite(tol) || tol <= 0) {
  stop(sprintf("FASTQR_TOL must be > 0, got: %s", as.character(tol)))
}
require_pos_int(maxit, "FASTQR_MAXIT")

config_base <- list(
  case = case,
  tau = tau,
  n = n,
  num_rep = num_rep,
  h = NULL,
  z = NULL,
  tol = tol,
  maxit = maxit,
  bland = bland,
  track_order = track_order,
  Mm.factor = Mm.factor,
  seed_base = seed_base
)

tau_str <- sprintf("tau%02d", as.integer(round(tau * 100)))
partial_dir <- file.path("data/llqr_simu_results", ".array_tmp",
                         sprintf("case%d_%s_n%d_rep%d", case, tau_str, n, num_rep))
dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== LLQR ARRAY DRIVER ===\n")
cat(sprintf("task_id=%d, ncores=%d, chunk_size=%d\n", task_id, ncores, chunk_size))
cat(sprintf("case=%d, tau=%.2f, n=%d, num_rep=%d\n", case, tau, n, num_rep))
if (sparse_mode) {
  cat(sprintf("sparse rep positions: %d-%d of %d (len=%d)\n", rep_start, rep_end, total_rep_targets, length(rep_ids)))
  cat(sprintf("rep ids: %s\n", paste(utils::head(rep_ids, 20L), collapse = ",")))
} else {
  cat(sprintf("rep range: %d-%d (len=%d)\n", rep_start, rep_end, length(rep_ids)))
}
cat(sprintf("partial_dir=%s\n", partial_dir))
cat("Mm.factor:", paste(Mm.factor, collapse = ", "), "\n")
cat("seed_base:", seed_base, "\n\n")
cat("max_attempts_per_rep:", max_attempts_per_rep, "\n")
cat("retry_stride:", retry_stride, "\n\n")
cat("max_seconds_per_rep:", max_seconds_per_rep, "\n\n")

if (length(rep_ids) == 0) {
  cat("No rep_ids assigned to this task. Exiting.\n")
  quit(save = "no", status = 0)
}

methods <- create_llqr_methods(config_base$Mm.factor)
method_names <- names(methods)

run_rep_with_timeout <- function(rep_id, rep_config, methods, timeout_sec) {
  if (.Platform$OS.type != "unix") {
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
      rr <- run_rep_with_timeout(rep_id, rep_config, methods, max_seconds_per_rep)

      timing_matrix <- matrix(NA_real_, nrow = 1, ncol = length(method_names),
                              dimnames = list(NULL, method_names))
      for (m in method_names) timing_matrix[1, m] <- as.numeric(rr$timing[[m]])

      estimates_list <- setNames(vector("list", length(method_names)), method_names)
      H_seq_list     <- setNames(vector("list", length(method_names)), method_names)
      method_metadata <- setNames(vector("list", length(method_names)), method_names)
      for (m in method_names) {
        estimates_list[[m]] <- list(rr$estimates[[m]])
        H_seq_list[[m]]     <- list(rr$H_seq[[m]])
        method_metadata[[m]] <- rr$method_metadata[[m]]
      }

      llqr_meta <- rr$method_metadata[["llqr"]]

      partial_results <- list(
        config = rep_config,
        timing_matrix = timing_matrix,
        estimates_list = estimates_list,
        H_seq_list = H_seq_list,
        method_metadata = method_metadata,
        method_names = method_names,
        Mm.factor_mapping = create_llqr_Mm_factor_mapping(method_names, config_base$Mm.factor),
        timestamp = Sys.time(),
        simulation_type = "llqr_partial",
        status = "success",
        error_msg = NULL,
        seed_used = seed_used,
        attempts = attempt,
        max_seconds_per_rep = max_seconds_per_rep,
        h_used = if (!is.null(llqr_meta$h_used)) llqr_meta$h_used else NA_real_,
        h_retry_factor = if (!is.null(llqr_meta$h_retry_factor)) llqr_meta$h_retry_factor else NA_real_,
        llqr_attempts = if (!is.null(llqr_meta$llqr_attempts)) llqr_meta$llqr_attempts else NA_integer_
      )

      list(ok = TRUE, obj = partial_results)
    }, error = function(e) {
      partial_results <- list(
        config = rep_config,
        method_names = method_names,
        timestamp = Sys.time(),
        simulation_type = "llqr_partial",
        status = "error",
        error_msg = conditionMessage(e),
        seed_used = seed_used,
        attempts = attempt,
        max_seconds_per_rep = max_seconds_per_rep
      )
      list(ok = FALSE, obj = partial_results)
    })

    if (isTRUE(out$ok)) break
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
if (ncores > 1 && length(rep_ids) > 1) {
  suppressPackageStartupMessages({
    library(doParallel)
    library(foreach)
  })
  cluster_type <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"
  cat(sprintf("parallel cluster type: %s\n", cluster_type))
  cl <- if (.Platform$OS.type == "unix") {
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
    suppressPackageStartupMessages(library(microbenchmark))
    NULL
  })
  
  clusterExport(cl, varlist = c("rep_ids", "run_one"), envir = environment())
  
  foreach(r = rep_ids) %dopar% {
    run_one(r)
    NULL
  }
} else {
  for (r in rep_ids) run_one(r)
}


cat("\n=== LLQR ARRAY DRIVER DONE ===\n")
