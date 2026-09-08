# LLQR/TVCQR v41-v42 rep100 HPC experiment

This isolated harness compares four methods on the current paper DGPs:

- `direct_baseline`: `llqr_local_fit` for LLQR and `tvc_rq` for TVCQR;
- `lean_seq`: the no-history sequential Fortran comparator;
- `v41`: the research-only active-first ppro variant;
- `v42`: the paper/contract-aligned ppro variant.

The approved grid is Cases 1--2, taus 0.2/0.5/0.8, sample sizes
1000/2000/5000/10000, 100 replications, and seed base 2026. Replication `r`
uses seed `2026 + r`, so the realized seeds are 2027--2126. The harness never
substitutes a different seed after failure.

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
bash experiments/hpc_v41_v42_rep100/build_variants.sh
FASTQR_RUN_TAG=<unique-tag> \
  bash experiments/hpc_v41_v42_rep100/submit_all.sh
```

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
