# scripts/array/merge_tvcqr_array.R
# Merge TVCQR partial results into a single results object compatible with performance_measurement.R

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

case     <- as_int("FASTQR_CASE", 1)
tau      <- as_num("FASTQR_TAU", 0.5)
n        <- as_int("FASTQR_N", 200)
num_rep  <- as_int("FASTQR_NUM_REP", 500)

Mm.factor <- as_num_vec("FASTQR_MM_FACTOR", c(1e-3, 1e-4))
h.factor  <- as_num("FASTQR_H_FACTOR", 1)
tol       <- as_num("FASTQR_TOL", 1e-14)
maxit     <- as_int("FASTQR_MAXIT", 1e6)
bland     <- as_bool("FASTQR_BLAND", FALSE)
eps       <- as_num("FASTQR_EPS", 1e-6)
cpp_helper<- as_bool("FASTQR_CPP_HELPER", FALSE)
seed_base <- as_int("FASTQR_SEED_BASE", 2025)
J         <- as_int("FASTQR_J", 100)
burn_in   <- as_int("FASTQR_BURN_IN", 500)
include_h_seq <- as_bool("FASTQR_MERGE_INCLUDE_H_SEQ", FALSE)
progress_every <- as_int("FASTQR_MERGE_PROGRESS_EVERY", 50)

config <- list(
  case = case,
  tau = tau,
  n = n,
  num_rep = num_rep,
  h = NULL,
  h.factor = h.factor,
  tol = tol,
  maxit = maxit,
  bland = bland,
  Mm.factor = Mm.factor,
  eps = eps,
  cpp_helper = cpp_helper,
  seed_base = seed_base,
  J = J,
  burn_in = burn_in
)
config <- complete_tvcqr_sim_config(config)

tau_str <- sprintf("tau%02d", as.integer(round(tau * 100)))
partial_dir <- file.path("data/tvcqr_simu_results", ".array_tmp",
                         sprintf("case%d_%s_n%d_rep%d", case, tau_str, n, num_rep))
dir.create("data/tvcqr_simu_results", recursive = TRUE, showWarnings = FALSE)

cat("=== TVCQR MERGE ===\n")
cat(sprintf("partial_dir=%s\n", partial_dir))
cat(sprintf("case=%d (%s), tau=%.2f, n=%d, num_rep=%d\n",
            config$case, config$case_label, tau, n, num_rep))
if (config$case == 2L) {
  cat(sprintf("J=%d, burn_in=%d\n", config$J, config$burn_in))
}
cat("\n")

# Determine method_names (create in the same way as simulation)
methods <- create_tvcqr_methods(Mm.factor)
method_names <- names(methods)
num_methods <- length(method_names)

timing_matrix <- matrix(NA_real_, nrow = num_rep, ncol = num_methods,
                        dimnames = list(NULL, method_names))

estimates_list <- setNames(vector("list", num_methods), method_names)
H_seq_list     <- if (include_h_seq) setNames(vector("list", num_methods), method_names) else NULL
method_metadata <- setNames(vector("list", num_methods), method_names)
for (m in method_names) {
  estimates_list[[m]] <- vector("list", num_rep)
  method_metadata[[m]] <- vector("list", num_rep)
  if (include_h_seq) {
    H_seq_list[[m]] <- vector("list", num_rep)
  }
}

missing <- integer()
failed  <- integer()

for (rep_id in 1:num_rep) {
  if (!is.na(progress_every) && progress_every > 0L &&
      (rep_id == 1L || rep_id %% progress_every == 0L || rep_id == num_rep)) {
    cat(sprintf("Merging rep %d / %d\n", rep_id, num_rep))
    flush.console()
  }
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
  
  # Fill timing and estimates
  timing_matrix[rep_id, ] <- pr$timing_matrix[1, method_names]
  for (m in method_names) {
    estimates_list[[m]][[rep_id]] <- pr$estimates_list[[m]][[1]]
    method_metadata[[m]][rep_id] <- list(if (!is.null(pr$method_metadata)) pr$method_metadata[[m]] else NULL)
    if (include_h_seq) {
      H_seq_list[[m]][rep_id] <- list(if (!is.null(pr$H_seq_list)) pr$H_seq_list[[m]][[1]] else NULL)
    }
  }
}

results <- list(
  config = config,
  timing_matrix = timing_matrix,
  estimates_list = estimates_list,
  method_metadata = method_metadata,
  method_metadata_list = method_metadata,
  method_names = method_names,
  Mm.factor_mapping = create_tvcqr_Mm_factor_mapping(method_names, Mm.factor),
  timestamp = Sys.time(),
  case_label = config$case_label,
  case_key = config$case_key,
  simulation_type = "tvcqr"
)
if (include_h_seq) {
  results$H_seq_list <- H_seq_list
}

# Final filename consistent with your existing convention
final_file <- file.path("data/tvcqr_simu_results",
                        sprintf("case%d_%s_n%d_rep%d.RData", case, tau_str, n, num_rep))
save(results, file = final_file)

cat(sprintf("\nSaved merged results: %s\n", final_file))
cat(sprintf("Missing reps: %d\n", length(missing)))
cat(sprintf("Failed reps:  %d\n", length(failed)))
if (length(missing) > 0) cat("First missing:", paste(head(missing, 10), collapse = ","), "\n")
if (length(failed)  > 0) cat("First failed: ", paste(head(failed, 10),  collapse = ","), "\n")
cat("=== TVCQR MERGE DONE ===\n")
