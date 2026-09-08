# Experiment-only orchestration; no solver or public R API changes.
experiment_methods <- c("direct_baseline", "lean_seq", "v41", "v42")

experiment_models <- function() {
  value <- strsplit(Sys.getenv("FASTQR_MODELS", "llqr,tvcqr"), ",", fixed = TRUE)[[1L]]
  value <- unique(trimws(tolower(value)))
  if (!length(value) || any(!value %in% c("llqr", "tvcqr"))) stop("Invalid FASTQR_MODELS")
  value
}

retry_settings <- function() {
  policy <- Sys.getenv("FASTQR_BASELINE_RETRY_POLICY", "none")
  maximum <- suppressWarnings(as.numeric(Sys.getenv("FASTQR_MAX_ATTEMPTS_PER_REP", "20")))
  stride <- suppressWarnings(as.numeric(Sys.getenv("FASTQR_RETRY_STRIDE", "1000000")))
  if (!policy %in% c("none", "paper_baseline_error")) stop("Invalid baseline retry policy")
  for (value in list(maximum, stride)) {
    if (length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value)) {
      stop("Retry limit and stride must be positive integers")
    }
  }
  list(policy = policy, max_attempts = as.integer(maximum), stride = stride)
}

baseline_retryable <- function(rows) {
  b <- rows[rows$method == "direct_baseline", , drop = FALSE]
  nrow(b) == 1L && isTRUE(b$threw_error) && identical(b$error_stage, "fit")
}

run_baseline_retries <- function(rep_id, seed_base, settings, attempt_fun) {
  limit <- if (settings$policy == "paper_baseline_error") settings$max_attempts else 1L
  original_seed <- as.double(seed_base) + rep_id
  if (original_seed + (limit - 1) * settings$stride > .Machine$integer.max) stop("Seed overflow")
  attempts <- list()
  for (attempt in seq_len(limit)) {
    seed <- as.integer(original_seed + (attempt - 1) * settings$stride)
    rows <- attempt_fun(rep_id, seed)
    if (!is.data.frame(rows) || nrow(rows) != 4L || !setequal(rows$method, experiment_methods)) {
      stop("An attempt must return four method records")
    }
    rows$original_seed <- as.integer(original_seed)
    rows$attempt <- attempt
    rows$retry_policy <- settings$policy
    rows$max_attempts <- settings$max_attempts
    rows$retry_stride <- settings$stride
    eligible <- baseline_retryable(rows)
    rows$baseline_retryable <- eligible
    rows$retry_exhausted <- settings$policy == "paper_baseline_error" && eligible && attempt == limit
    attempts[[attempt]] <- rows
    if (settings$policy == "none" || !eligible) break
  }
  all_rows <- do.call(rbind, attempts)
  all_rows$is_final_attempt <- all_rows$attempt == length(attempts)
  rownames(all_rows) <- NULL
  final <- all_rows[all_rows$is_final_attempt, , drop = FALSE]
  rownames(final) <- NULL
  list(schema_version = 2L, final = final, attempts = all_rows)
}

validate_retry_result <- function(value, rep_id, seed_base, identity, settings) {
  # Old one-attempt results remain readable when retry is disabled.
  if (is.data.frame(value)) {
    if (settings$policy != "none") return(FALSE)
    value$threw_error <- FALSE
    value$error_stage <- NA_character_
    value$data_generation_sec <- NA_real_
    value <- run_baseline_retries(rep_id, seed_base, settings, function(...) value)
  }
  required <- c(names(identity), "rep_id", "seed", "method", "method_position",
                "original_seed", "attempt", "retry_policy", "max_attempts", "retry_stride",
                "baseline_retryable", "retry_exhausted", "is_final_attempt", "threw_error", "error_stage")
  if (!is.list(value) || !identical(value$schema_version, 2L) ||
      !is.data.frame(value$final) || !is.data.frame(value$attempts) ||
      !all(required %in% names(value$attempts))) return(FALSE)
  a <- value$attempts
  if (!nrow(a) || anyNA(a[c("rep_id", "seed", "attempt", "method_position", "is_final_attempt")])) return(FALSE)
  for (key in names(identity)) if (!isTRUE(all(a[[key]] == identity[[key]]))) return(FALSE)
  if (!isTRUE(all(a$rep_id == rep_id)) || !isTRUE(all(a$original_seed == seed_base + rep_id)) ||
      !isTRUE(all(a$retry_policy == settings$policy)) ||
      !isTRUE(all(a$max_attempts == settings$max_attempts)) ||
      !isTRUE(all(a$retry_stride == settings$stride))) return(FALSE)
  last <- max(a$attempt)
  limit <- if (settings$policy == "none") 1L else settings$max_attempts
  if (last < 1L || last > limit || !identical(sort(unique(as.integer(a$attempt))), seq_len(last))) return(FALSE)
  for (k in seq_len(last)) {
    x <- a[a$attempt == k, , drop = FALSE]
    expected_seed <- seed_base + rep_id + (k - 1) * settings$stride
    expected_position <- ((match(x$method, experiment_methods) - 1L - (rep_id - 1L) %% 4L) %% 4L) + 1L
    if (nrow(x) != 4L || !setequal(x$method, experiment_methods) ||
        !isTRUE(all(x$seed == expected_seed)) || !isTRUE(all(x$method_position == expected_position)) ||
        !isTRUE(all(x$is_final_attempt == (k == last)))) return(FALSE)
    eligible <- baseline_retryable(x)
    exhausted <- settings$policy == "paper_baseline_error" && eligible && k == limit
    if (!isTRUE(all(x$baseline_retryable == eligible)) || !isTRUE(all(x$retry_exhausted == exhausted))) return(FALSE)
    if (k < last && !eligible) return(FALSE)
    if (k == last && eligible && settings$policy == "paper_baseline_error" && !exhausted) return(FALSE)
  }
  final <- a[a$is_final_attempt, , drop = FALSE]
  rownames(final) <- NULL
  identical(final, value$final)
}
