# LLQR Case 2 logistic/interior rep500 experiment

This branch adapts the paper experiment frozen at
`61dc4d6b39fab585ca62c590a7024fcc9eb95751`. All production R/Fortran sources,
the experiment Fortran kernels, and model adapters remain unchanged.
No manuscript or previous result is overwritten.

## Design and methods

- LLQR Case 2 only: X ~ U(0,1); Y = 1 + 2X^2 + error, where independent
  errors follow logistic(location=0, scale=sqrt(3)/pi), with variance one.
- Fit using ALL n observations. Evaluate at sorted observed X in [0.1,0.9].
  If this set is empty, use the single point 0.5. Save actual n_eval,
  n_interior, n_transitions, and grid_placeholder for each attempt.
- Epanechnikov kernel, h=n^(-1/5), threshold
  0.1*sqrt(log(n)/(n*h))*log(log(n)); n is the full sample size.
- n=1000/2000/5000/10000, tau=.2/.5/.8, 500 replications: 12 configurations.
- Methods: direct_baseline, lean_seq, unified_u11 (paper screen-seq).
- Initial seeds 2026--2525. Only a direct-baseline fit exception triggers
  regeneration, at most 20 attempts, with seed=2025+rep_id+1e6*(attempt-1).
  Failures of seq or screen-seq do not change the dataset. All attempts persist.
- Three method positions rotate by rep_id: 166 or 167 occurrences each.
- Time the complete fitting call, including bandwidth calculation. Exclude
  dataset generation, common grid filtering/sorting, garbage collection,
  post-fit validation, and discarded-attempt calls from the final timing table.

A1's reverse conditional-density bound follows from the logistic log-density
being 1/s-Lipschitz and the range of 1+2x^2 on [0,1] having width 2:
f(X|Y) <= exp(2/s). Take D=[.1,.9] and D+=[.05,.95]. If regeneration occurs,
accepted datasets are selected by the stated numerical acceptance rule.

## Integrity and report definitions

H paths are checked as row sets in memory against lean seq. They have n_eval
rows but contain original observation indices up to n. Full H paths and fit
vectors are not saved. Report maximum over replications of the within-replication
mean |candidate-direct|/max(|direct|,1e-10), using n_eval fitted quantiles.
Derivative estimates are checked for finiteness but excluded from that metric.

First-pass rates use n_eval-1 and exclude the first full fit. A single-point
path has NA first-pass rate and zero repairs and retained maximum; its count
is reported separately. Retained size includes basis padding, excludes aggregate
rows, and is normalized by n*h and n*h*threshold+log(n). Internal full-active
recovery is reported separately from sequential fallback. No speed gate is used.

## Cluster launch

Use a dedicated worktree at the exact pushed SHA. From its root:

```bash
export SSQR_PROJECT_ROOT="$PWD"
export SSQR_PUSHED_SHA="$(git rev-parse HEAD)"
bash experiments/hpc_u11_rep500/submit_logistic.sh
```

The entry point checks the frozen sources and submits a dependency chain:
compiled build -> pipeline smoke -> 12-config cache audit (25/27/31) -> validation
-> dispatcher -> formal rep500 arrays -> merges -> summary/report -> manifest.
Formal arrays are never submitted unless preflight passes for the same SHA.
If preflight fails, the dispatcher remains pending; inspect its preflight logs
and cancel that dispatcher before a corrected new launch. Do not re-run the
entry point merely because formal arrays have not appeared yet.

The launcher prints and saves run/launch/preflight paths and the dispatcher ID.
`FORMAL_SUBMITTED` in the launch directory confirms all formal submissions.
The default result root is `results/hpc/llqr_case2_logistic_runs/` inside the
new checkout. The run directory contains partials, attempts, integrity records,
source hashes, Slurm logs and job IDs. On successful validation,
`paper_staging/REPORT.md` and nine CSV files contain every applicable paper table
plus grid sizes and recovery details. Despite the inherited staging script name,
no files under paper/ are read or written.

Resources remain the paper defaults (cnall/users; GCC12.2/R4.3.1): n<10000 uses
56 workers, 9 array tasks; n=10000 uses 14 workers, 36 tasks, one concurrent
array task per configuration. Numerical libraries use one thread per worker.

## Local validation

See VALIDATION.md for exact commands and outcomes. Cluster preflight has not
been run locally and remains mandatory. Formal rep500 results are not yet available.
