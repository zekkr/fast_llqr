# scripts/array/driver_llqr_array.R
# Run LLQR simulation in SLURM job array (chunked), write partial results per replication.

suppressPackageStartupMessages({
  library(microbenchmark)
})

source("scripts/setup.R")

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

task_id  <- as_int("SLURM_ARRAY_TASK_ID", 1)
ncores   <- as_int("SLURM_CPUS_PER_TASK", 1)

case     <- as_int("FASTQR_CASE", 1)
tau      <- as_num("FASTQR_TAU", 0.5)
n        <- as_int("FASTQR_N", 500)
num_rep  <- as_int("FASTQR_NUM_REP", 500)

chunk_size <- as_int("FASTQR_CHUNK_SIZE", ncores)
if (chunk_size <= 0) chunk_size <- 1

rep_start <- (task_id - 1L) * chunk_size + 1L
rep_end   <- min(task_id * chunk_size, num_rep)
rep_ids   <- rep_start:rep_end
rep_ids   <- rep_ids[rep_ids >= 1 & rep_ids <= num_rep]

Mm.factor <- as_num_vec("FASTQR_MM_FACTOR", c(1e-2, 1e-3, 1e-4))
tol       <- as_num("FASTQR_TOL", 1e-14)
maxit     <- as_int("FASTQR_MAXIT", 2e6)
bland     <- as_bool("FASTQR_BLAND", FALSE)
track_order <- as_bool("FASTQR_TRACK_ORDER", TRUE)
seed_base <- as_int("FASTQR_SEED_BASE", 2026)

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
cat(sprintf("rep range: %d-%d (len=%d)\n", rep_start, rep_end, length(rep_ids)))
cat(sprintf("partial_dir=%s\n", partial_dir))
cat("Mm.factor:", paste(Mm.factor, collapse = ", "), "\n")
cat("seed_base:", seed_base, "\n\n")

if (length(rep_ids) == 0) {
  cat("No rep_ids assigned to this task. Exiting.\n")
  quit(save = "no", status = 0)
}

methods <- create_llqr_methods(config_base$Mm.factor)
method_names <- names(methods)

run_one <- function(rep_id) {
  rep_config <- config_base
  rep_config$num_rep <- 1
  rep_config$rep_id <- rep_id
  
  out <- tryCatch({
    rr <- run_single_llqr_replication(rep_id, rep_config, methods)
    
    timing_matrix <- matrix(NA_real_, nrow = 1, ncol = length(method_names),
                            dimnames = list(NULL, method_names))
    for (m in method_names) timing_matrix[1, m] <- as.numeric(rr$timing[[m]])
    
    estimates_list <- setNames(vector("list", length(method_names)), method_names)
    H_seq_list     <- setNames(vector("list", length(method_names)), method_names)
    for (m in method_names) {
      estimates_list[[m]] <- list(rr$estimates[[m]])
      H_seq_list[[m]]     <- list(rr$H_seq[[m]])
    }
    
    partial_results <- list(
      config = rep_config,
      timing_matrix = timing_matrix,
      estimates_list = estimates_list,
      H_seq_list = H_seq_list,
      method_names = method_names,
      Mm.factor_mapping = create_llqr_Mm_factor_mapping(method_names, config_base$Mm.factor),
      timestamp = Sys.time(),
      simulation_type = "llqr_partial",
      status = "success",
      error_msg = NULL
    )
    
    list(ok = TRUE, obj = partial_results)
  }, error = function(e) {
    partial_results <- list(
      config = rep_config,
      method_names = method_names,
      timestamp = Sys.time(),
      simulation_type = "llqr_partial",
      status = "error",
      error_msg = conditionMessage(e)
    )
    list(ok = FALSE, obj = partial_results)
  })
  
  save_path <- file.path(partial_dir, sprintf("rep%04d.RData", rep_id))
  partial_results <- out$obj
  save(partial_results, file = save_path)
  cat(sprintf("[%s] rep=%d -> %s (%s)\n",
              format(Sys.time(), "%F %T"), rep_id, save_path,
              partial_results$status))
  invisible(TRUE)
}

if (ncores > 1 && length(rep_ids) > 1) {
  suppressPackageStartupMessages({
    library(doParallel)
    library(foreach)
  })
  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl), add = TRUE)
  registerDoParallel(cl)
  
  clusterEvalQ(cl, {
    source("scripts/setup.R")
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
