# LLQR formula-interface paired rerun (2026-09-21)

This run changes only direct LLQR fitting to `quantreg::rq(y ~ x, weights=w,
tau=tau, method=method)`, retaining every observation. Default method is `br`;
no stopping tolerance is supplied. Historical `reproduction/frozen` stays immutable.

Paper case 1 is the archived logistic variance-one model on interior observed
points; paper case 2 is the archived normal-error uniform-X model on all sorted
observed points. Both use kernel case 2 and h=n^(-.2). These paper case IDs are
not the old kernel case IDs. All source adapters and Fortran solvers come directly
from `reproduction/frozen`; the live `R/llqr_functions.R` supplies the new baseline.

The three raw method identifiers are direct_baseline, lean_seq, unified_u11.
There are 24 configurations, n=1000/2000/5000/10000, tau=.2/.5/.8, rep=500,
seed_base=2025. Method order rotates by replicate. Each complete call is timed;
generation, pre-call gc, diagnostics, dataset fingerprints and output are excluded.
Only thrown baseline fitting errors regenerate data for all methods (stride
1000000, at most 20 attempts); other failures stop acceptance. Full H paths are
compared in memory and discarded. Compact recovery and fallback records persist.
The discrepancy denominator is pmax(abs(baseline),1e-10), with mean over points
then maximum over replicates. Old archived results are never overwritten.

On the server, after fetching the pushed commit in a separate checkout:

```sh
python3 experiments/llqr_rq_rep500/submit.py preflight --tag UNIQUE_PREFLIGHT --sha PUSHED_SHA
# Only after results/hpc/llqr_rq_runs/UNIQUE_PREFLIGHT/PREFLIGHT_PASS exists:
python3 experiments/llqr_rq_rep500/submit.py formal --tag UNIQUE_REP500 --sha PUSHED_SHA --preflight UNIQUE_PREFLIGHT
```

The preflight builds the frozen solvers, tests the new baseline, runs both small
models, then tests every production configuration once. Formal submission requires
a passing preflight at the exact same commit. Worker records require R4.3.1,
quantreg5.94, Gold6258R, and single-thread numerical libraries. Both stages write
Slurm job IDs and complete effective parameters before launching the arrays.

Validation: `Rscript experiments/llqr_rq_rep500/test_baseline.R`;
`python3 reproduction/check_frozen.py`. New runs do not certify historical TVCQR
environments. The original paper's TVCQR and MST results remain separate sources.
