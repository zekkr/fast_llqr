# Combine the six independently checkpointed local Table S1 configurations.
root <- "results/llqr_multivar_iterations/mst_rep100_local_20260924"
subdirs <- unlist(lapply(1:2, function(ca) sprintf("case%d_n%d", ca, c(300, 500, 800))))
records <- lapply(subdirs, function(id) {
  path <- file.path(root, id)
  manifest <- readRDS(file.path(path, "manifest.rds"))
  result <- readRDS(file.path(path, "results.rds"))
  if (!identical(manifest$config, result$config)) stop("Config mismatch: ", id)
  if (result$config$reps != 100L || length(result$config$cases) != 1L ||
      length(result$config$n_values) != 1L) stop("Unexpected run config: ", id)
  list(id = id, manifest = manifest, result = result)
})
source_hashes <- lapply(records, function(x) x$manifest$source_hashes)
if (!all(vapply(source_hashes, identical, logical(1), source_hashes[[1L]]))) {
  stop("Source hashes differ across configurations.")
}
if (!all(vapply(records, function(x) identical(x$manifest$rng_kind,
                                          records[[1L]]$manifest$rng_kind), logical(1)))) {
  stop("RNG kinds differ across configurations.")
}
raw <- do.call(rbind, lapply(records, function(x) x$result$raw))
diag <- do.call(rbind, lapply(records, function(x) x$result$diagnostics))
if (nrow(raw) != 1200L || nrow(diag) != 600L ||
    !identical(sort(unique(raw$method)), c("cold", "seq_mst"))) {
  stop("Missing or extra method records.")
}
keys <- paste(raw$case, raw$n, raw$replication, raw$method)
dkeys <- paste(diag$case, diag$n, diag$replication)
if (anyDuplicated(keys) || anyDuplicated(dkeys)) stop("Duplicate keys.")
for (ca in 1:2) for (n in c(300, 500, 800)) {
  for (method in c("cold", "seq_mst")) {
    reps <- sort(raw$replication[raw$case == ca & raw$n == n & raw$method == method])
    if (!identical(reps, 1:100)) stop("Incomplete replication set: ", ca, "/", n, "/", method)
  }
  reps <- sort(diag$replication[diag$case == ca & diag$n == n])
  if (!identical(reps, 1:100)) stop("Incomplete diagnostics: ", ca, "/", n)
}
if (any(raw$max_iterations >= 20000L) ||
    any(!is.finite(diag$mst_max_abs_diff)) ||
    any(diag$mst_max_abs_diff > 1e-8)) stop("Fit validation failure.")
expected_seed <- 20260405L + 100000L * diag$case + 1000L * diag$n +
  40L + diag$replication
if (!identical(as.integer(diag$seed), as.integer(expected_seed))) stop("Seed mismatch.")
paired <- merge(raw[raw$method == "cold", c("case", "n", "replication", "total_iterations")],
                raw[raw$method == "seq_mst", c("case", "n", "replication", "total_iterations")],
                by = c("case", "n", "replication"), suffixes = c("_direct", "_mst"))
paired$reduction <- 1 - paired$total_iterations_mst / paired$total_iterations_direct
means <- aggregate(total_iterations ~ case + n + method, raw, mean)
reductions <- aggregate(reduction ~ case + n, paired, mean)
rows <- data.frame()
for (ca in 1:2) for (n in c(300, 500, 800)) {
  direct <- means$total_iterations[means$case == ca & means$n == n & means$method == "cold"]
  mst <- means$total_iterations[means$case == ca & means$n == n & means$method == "seq_mst"]
  reduction <- reductions$reduction[reductions$case == ca & reductions$n == n]
  rows <- rbind(rows, data.frame(case = ca, n = n, direct = direct, mst = mst,
                                 reduction = reduction,
                                 display_direct = format(round(direct), big.mark = ",", trim = TRUE),
                                 display_mst = format(round(mst), big.mark = ",", trim = TRUE),
                                 display_reduction = sprintf("%.1f\\%%", 100 * reduction)))
}
write.csv(raw, file.path(root, "raw_all.csv"), row.names = FALSE)
write.csv(diag, file.path(root, "diagnostics_all.csv"), row.names = FALSE)
write.csv(paired, file.path(root, "paired_all.csv"), row.names = FALSE)
write.csv(rows, file.path(root, "table_s1_all.csv"), row.names = FALSE)
summary <- list(table_s1 = rows, max_fit_difference = max(diag$mst_max_abs_diff),
                mean_parent_child_edge = mean(diag$mst_mean_edge_weight),
                source_hashes = source_hashes[[1L]], raw_rows = nrow(raw),
                diagnostic_rows = nrow(diag))
saveRDS(summary, file.path(root, "aggregate_summary.rds"))
cat("All 600 paired datasets validated.\n")
print(rows, row.names = FALSE)
cat("Maximum fit discrepancy:", format(summary$max_fit_difference, digits = 10), "\n")
cat("Mean parent-child edge:", format(summary$mean_parent_child_edge, digits = 10), "\n")
