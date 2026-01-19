# scripts/array/merge_llqr_array.R
# Merge LLQR partial results into a single results object compatible with performance_measurement.R

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

case     <- as_int("FASTQR_CASE", 1)
tau      <- as_num("FASTQR_TAU", 0.5)
n        <- as_int("FASTQR_N", 500)
num_rep  <- as_int("FASTQR_NUM_REP", 500)

Mm.factor <- as_num_vec("FASTQR_MM_FACTOR", c(1e-2, 1e-3, 1e-4))
tol       <- as_num("FASTQR_TOL", 1e-14)
maxit     <- as_int("FASTQR_MAXIT", 2e6)
bland     <- as_bool("FASTQR_BLAND", FALSE)
track_order <- as_bool("FASTQR_TRACK_ORDER", TRUE)
seed_base <- as_int("FASTQR_SEED_BASE", 2026)

config <- list(
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
dir.create("data/llqr_simu_results", recursive = TRUE, showWarnings = FALSE)

cat("=== LLQR MERGE ===\n")
cat(sprintf("partial_dir=%s\n", partial_dir))
cat(sprintf("case=%d, tau=%.2f, n=%d, num_rep=%d\n\n", case, tau, n, num_rep))

methods <- create_llqr_methods(Mm.factor)
method_names <- names(methods)
num_methods <- length(method_names)

timing_matrix <- matrix(NA_real_, nrow = num_rep, ncol = num_methods,
                        dimnames = list(NULL, method_names))

estimates_list <- setNames(vector("list", num_methods), method_names)
H_seq_list     <- setNames(vector("list", num_methods), method_names)
for (m in method_names) {
  estimates_list[[m]] <- vector("list", num_rep)
  H_seq_list[[m]]     <- vector("list", num_rep)
}

missing <- integer()
failed  <- integer()

for (rep_id in 1:num_rep) {
  f <- file.path(partial_dir, sprintf("rep%04d.RData", rep_id))
  if (!file.exists(f)) {
    missing <- c(missing, rep_id)
    next
  }
  env <- new.env()
  load(f, envir = env)
  if (!exists("partial_results", envir = env)) {
    failed <- c(failed, rep_id)
    next
  }
  pr <- get("partial_results", envir = env)
  if (is.null(pr$status) || pr$status != "success") {
    failed <- c(failed, rep_id)
    next
  }
  
  timing_matrix[rep_id, ] <- pr$timing_matrix[1, method_names]
  for (m in method_names) {
    estimates_list[[m]][[rep_id]] <- pr$estimates_list[[m]][[1]]
    H_seq_list[[m]][[rep_id]]     <- pr$H_seq_list[[m]][[1]]
  }
}

results <- list(
  config = config,
  timing_matrix = timing_matrix,
  estimates_list = estimates_list,
  H_seq_list = H_seq_list,
  method_names = method_names,
  Mm.factor_mapping = create_llqr_Mm_factor_mapping(method_names, Mm.factor),
  timestamp = Sys.time(),
  simulation_type = "llqr"
)

final_file <- file.path("data/llqr_simu_results",
                        sprintf("case%d_%s_n%d_rep%d.RData", case, tau_str, n, num_rep))
save(results, file = final_file)

cat(sprintf("\nSaved merged results: %s\n", final_file))
cat(sprintf("Missing reps: %d\n", length(missing)))
cat(sprintf("Failed reps:  %d\n", length(failed)))
if (length(missing) > 0) cat("First missing:", paste(head(missing, 10), collapse = ","), "\n")
if (length(failed)  > 0) cat("First failed: ", paste(head(failed, 10),  collapse = ","), "\n")
cat("=== LLQR MERGE DONE ===\n")
