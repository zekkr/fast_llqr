#!/usr/bin/env bash
set -euo pipefail

: "${SSQR_PROJECT_ROOT:?}" "${SSQR_EXPERIMENT_DIR:?}" "${SSQR_PREFLIGHT_DIR:?}" "${SSQR_PUSHED_SHA:?}"
cd "$SSQR_PROJECT_ROOT"
[[ "$(git rev-parse HEAD)" == "$SSQR_PUSHED_SHA" ]]
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

Rscript "$SSQR_EXPERIMENT_DIR/test_contract.R"
Rscript "$SSQR_EXPERIMENT_DIR/test_summary.R"
Rscript "$SSQR_EXPERIMENT_DIR/test_retained_diagnostics.R"

out="$SSQR_PREFLIGHT_DIR/hpc_smoke_results"
for model in llqr tvcqr; do
  tag="hpc_rep2_${model}"
  env SSQR_OUTPUT_ROOT="$out" SSQR_ALLOW_SMOKE=1 \
    SSQR_NUM_REP=2 SSQR_SEED_BASE=2025 SSQR_CASE=2 SSQR_TAU=.5 SSQR_N=1000 \
    SSQR_CHUNK_SIZE=2 SLURM_ARRAY_TASK_ID=1 SLURM_CPUS_PER_TASK=2 \
    SSQR_RUN_TAG="$tag" SSQR_MODEL="$model" SSQR_CACHE_FLAGS=27 \
    SSQR_PROVIDER_FLAGS=1 SSQR_BUILD_MODE=optimized \
    Rscript "$SSQR_EXPERIMENT_DIR/driver_array.R"
  env SSQR_OUTPUT_ROOT="$out" SSQR_ALLOW_SMOKE=1 \
    SSQR_NUM_REP=2 SSQR_SEED_BASE=2025 SSQR_CASE=2 SSQR_TAU=.5 SSQR_N=1000 \
    SSQR_RUN_TAG="$tag" SSQR_MODEL="$model" \
    Rscript "$SSQR_EXPERIMENT_DIR/merge_config.R"
done

Rscript -e '
for (model in c("llqr", "tvcqr")) {
  p <- Sys.glob(file.path(Sys.getenv("SSQR_PREFLIGHT_DIR"), "hpc_smoke_results",
                          paste0("hpc_rep2_", model), model, "*", "replication_metrics.rds"))
  stopifnot(length(p) == 1L)
  d <- readRDS(p)
  stopifnot(nrow(d) == 6L, all(d$solver_ok), all(d$accepted_ok),
            all(d$discrepancy_status == "ok"), all(!d$fallback_triggered))
  u <- d[d$method == "unified_u11", ]
  stopifnot(all(u$ierr == 0L), all(u$failed_eval == 0L),
            all(u$h_match_vs_lean), all(u$h_mismatch_eval_count == 0L),
            all(u$retained_decomposition_ok))
  cat(model, "pipeline PASS max discrepancy",
      format(max(d$discrepancy), scientific = TRUE), "\n")
}'
printf '%s\n' "$SSQR_PUSHED_SHA" > "$SSQR_PREFLIGHT_DIR/HPC_PIPELINE_SMOKE_PASS"
