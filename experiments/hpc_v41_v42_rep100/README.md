# LLQR/TVCQR v41-v42 rep100 HPC experiment

This isolated harness compares four methods on the current paper DGPs:

- `direct_baseline`: `llqr_local_fit` for LLQR and `tvc_rq` for TVCQR;
- `lean_seq`: the no-history sequential Fortran comparator;
- `v41`: the research-only active-first ppro variant;
- `v42`: the paper/contract-aligned ppro variant.

The approved grid is Cases 1--2, taus 0.2/0.5/0.8, sample sizes
1000/2000/5000/10000, 100 replications, and seed base 2026. Replication `r`
uses seed `2026 + r`, so the realized seeds are 2027--2126. The harness never
substitutes a different seed after failure under its default `none` retry policy.

## LLQR paper-driver regeneration run

The follow-up LLQR-only run selects `FASTQR_MODELS=llqr` and
`FASTQR_BASELINE_RETRY_POLICY=paper_baseline_error`. Only an exception from
the direct fitting call triggers regeneration. The next seed is
`seed_base + rep_id + (attempt - 1) * FASTQR_RETRY_STRIDE`, with defaults
`FASTQR_MAX_ATTEMPTS_PER_REP=20` and `FASTQR_RETRY_STRIDE=1000000`.
Candidate-only failures, returned nonfinite/malformed estimates, and data
generation errors do not trigger regeneration. All four methods run on every
attempt, with the same replication-dependent timing rotation. Their failures
on discarded attempts remain in the audit; they are not described as fixed.

Version-2 partial RDS files contain `final` (four compact method rows) and
`attempts` (all compact attempts), never full estimates or H_seq. The merge
checks identity, actual seeds, four-way pairing, timing positions, and the
baseline-error reason for every regeneration. Final fitting times and numerical
discrepancies use only the last attempt. `retry_summary.csv` reports discarded
fit time and data generation separately, outside the fitting-time comparison.
`seed_map.csv`, `attempt_metrics_all.csv`, and `all_attempt_failures.csv` retain
the complete provenance. A baseline exception on attempt 20 is an exhausted
failure. Final missing/failed calls still invalidate complete-config statistics.

```sh
Rscript experiments/hpc_v41_v42_rep100/test_retry.R
FASTQR_MODELS=llqr FASTQR_BASELINE_RETRY_POLICY=paper_baseline_error \
FASTQR_MAX_ATTEMPTS_PER_REP=20 FASTQR_RETRY_STRIDE=1000000 \
FASTQR_PREVIOUS_RUN_TAG=v41v42_rep100_seed2026_4da4122_20260908_103038 \
FASTQR_RUN_TAG=<new-unique-tag> \
  bash experiments/hpc_v41_v42_rep100/submit_all.sh
```

The user provisionally accepted TVCQR v41 on 2026-09-08; no TVCQR rerun or
production integration is included in this follow-up. Its active-first algorithm
status remains separately documented. The LLQR candidate kernels are unchanged.

## Fixed algorithm settings

- LLQR `Mm.factor=0.1`;
- TVCQR `Mm.factor=1e-5`;
- tolerance `1e-14`;
- `min_subsample_size=1`;
- `threshold_lower_bound=FALSE`;
- `always_same_h_refit=TRUE`;
- log-log threshold scaling;
- `h.factor=1`, Bland off;
- TVCQR Case 2 uses `J=100` and burn-in 500.

The ppro variants are called directly and never call the seq fallback. A v41
active-first failure may still use its internal full-`m` ppro first-point path;
that event is recorded separately and is not a seq fallback.

## Output contract

No full estimates or `H_seq` arrays are persisted. Each replication stores only
elapsed seconds, method/error status, compact ppro diagnostics, and the paper's
relative numerical discrepancy from the direct baseline. For LLQR the metric
averages over `n` fitted conditional quantiles. For TVCQR it averages over the
`n x 4` level-coefficient matrix. A missing, malformed, or non-finite
replication invalidates the corresponding config-method discrepancy summary;
the summary never drops bad replications with `na.rm`.

## Build and submit

```sh
module load compilers/gcc/v12.2.0 soft/R/v4.3.1
bash experiments/hpc_v41_v42_rep100/build_variants.sh
FASTQR_RUN_TAG=<unique-tag> \
  bash experiments/hpc_v41_v42_rep100/submit_all.sh
```

On Linux, the build script reads the LAPACK and BLAS linker flags from the
active R installation. This keeps the experiment DLLs on the same numerical
libraries used by the R fitting session (the cluster R module exposes these as
`Rlapack` and `Rblas`).

The submission script creates one array and one merge job for each of the 48
model/config combinations, then submits one final summary job. It uses 56
workers through `n=5000` and 14 workers at `n=10000`.

Final tables are stored below:

```text
data/v41_v42_rep100/<run-tag>/tables/
```

`candidate_vs_lean_seq.csv` defines its main timing ratio as candidate mean
elapsed time divided by lean-seq mean elapsed time. A value below one means the
candidate is faster.
