# Local validation — 2026-09-14

Working directory: `/private/tmp/fastllqr-case2-logistic-rep500`.
Base: `61dc4d6b39fab585ca62c590a7024fcc9eb95751`.

## Commands run

```bash
bash experiments/hpc_u11_rep500/build_all.sh
Rscript experiments/hpc_u11_rep500/test_contract.R
Rscript experiments/hpc_u11_rep500/test_interior.R
Rscript experiments/hpc_u11_rep500/test_retained_diagnostics.R
Rscript experiments/hpc_u11_rep500/test_summary.R
bash -n experiments/hpc_u11_rep500/submit_logistic.sh
git diff --check
git diff 61dc4d6 -- R src experiments/hpc_u11_rep500/fortran
```

The compiler completed; macOS SDK deployment-version and duplicate-rpath warnings
were emitted. The frozen core, stable seq and production sources have zero diff.

- Contract tests: 10 checks pass, including baseline-only regeneration,
  candidate/seq failure not regenerating, 20-attempt cap and rotated timing order.
- Interior tests: 14 checks pass, including placeholder, one evaluation point,
  sample IDs greater than n_eval, real single-point and 156-point fitting,
  direct numerical agreement, and transition diagnostic dimensions.
- Retained tests: 2 checks pass, including threshold expansion, basis padding,
  aggregate exclusion and cache modes 25/27/31 preserving outputs.
- Synthetic summary: full 12-config x rep500 x three-method output accepted;
  all requested standalone tables generated. Deliberately deleting one record
  causes status 2, reports 499 present, and leaves its mean runtime undefined.
  This is a synthetic aggregation test, not an actual rep500 experiment.

## Actual six-dataset smoke

```bash
export SSQR_PROJECT_ROOT="$PWD" SSQR_EXPERIMENT_DIR=experiments/hpc_u11_rep500
export SSQR_OUTPUT_ROOT="$PWD/results/local_validation" SSQR_ALLOW_SMOKE=1
export SSQR_NUM_REP=2 SSQR_SEED_BASE=2025 SSQR_CASE=2 SSQR_N=200 SSQR_MODEL=llqr
export SSQR_CHUNK_SIZE=2 SLURM_ARRAY_TASK_ID=1 SLURM_CPUS_PER_TASK=1
export SSQR_CACHE_FLAGS=27 SSQR_PROVIDER_FLAGS=1
for qt in 0.2 0.5 0.8; do
  export SSQR_TAU="$qt" SSQR_RUN_TAG="logistic_smoke_tau${qt}"
  Rscript experiments/hpc_u11_rep500/driver_array.R
  Rscript experiments/hpc_u11_rep500/merge_config.R
done
```

For each tau, seeds 2026 and 2027 used 156 and 163 evaluation points, respectively.
All 18/18 method fits succeeded. Each merge found 2/2 records, no missing,
malformed or extra records. Post-merge checks confirmed:

- Maximum average relative discrepancy over all 18 fits: 4.7585505881147663e-16.
- H row-set mismatches: 0 over all six screen-seq paths.
- Sequential fallback: 0; internal full-active recovery: 1.
- Iteration-limit hits: 0; retained decomposition failures: 0.
- First-pass denominators: 155 and 162, using actual grid sizes.

No objective/KKT proof or full-grid rep500 success is inferred from this smoke.
The new cluster preflight and formal 12-config rep500 run are not run yet.
