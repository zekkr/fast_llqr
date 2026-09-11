# U11 full-configuration rep500 experiment

This tracked, experiment-only harness compares `direct_baseline`, `lean_seq`,
and `unified_u11` for LLQR and TVCQR. It does not modify the production
solvers, stable sequential solvers, public R API, or paper sources.

## Fixed design

- Models: LLQR and TVCQR.
- Cases: model-specific cases 1 and 2 (paper Cases 1--2 and 3--4,
  respectively).
- Quantiles: 0.2, 0.5, and 0.8.
- Sample sizes: 1,000, 2,000, 5,000, and 10,000.
- Replications: 500 per configuration; `seed_base=2025`, so initial seeds are
  2026--2525.
- Direct-baseline fitting errors trigger the paper driver's regeneration rule:
  at most 20 attempts with stride 1,000,000. Candidate failures do not select a
  replacement data set.
- Method order rotates by replication ID. Each method appears in each position
  166 or 167 times per configuration.

The timed region covers bandwidth selection, method initialization, the full
evaluation-grid fit, and return-object construction. Data generation, garbage
collection, scheduling, and serialization are outside it.

Only compact per-replication metrics and attempt chains are persisted. Full
estimates and H paths exist only in memory for numerical-discrepancy and H
row-set checks.

## Retained-size diagnostics

For each transition after the first evaluation point:

- `initial_threshold_hits` counts individual observations hit by the original
  threshold before any threshold expansion;
- `effective_threshold_hits` counts the threshold-hit observations used by the
  first reduced-solve attempt;
- `basis_forced_rows` counts additional previous-basis observations forced into
  that attempt, including zero-weight padding;
- `first_retained_size` is the number of individual rows in that first attempt;
- `first_aggregate_rows` is the number of aggregate pseudo-rows;
- `first_tableau_rows` is the total first tableau size.

The harness enforces both
`first_retained_size = effective_threshold_hits + basis_forced_rows` and
`first_tableau_rows = first_retained_size + first_aggregate_rows`. The first
evaluation point is recorded as not applicable.

## Execution order

1. Run `build_all.sh` and the local tests.
2. On Explore1000, run `submit_preflight.sh`. It performs the two-model rep2
   pipeline smoke and a one-replication, all-48-configuration comparison of U11
   cache modes 25, 27, and 31.
3. Submit `submit_all.sh` only after the preflight creates `PREFLIGHT_PASS` for
   the exact pushed Git SHA.
4. After all merge jobs finish, `summarize_run.R` evaluates the strict stability
   gate. `build_paper_staging.R` creates paper candidate data and previews only
   when that gate passes.

The formal submitter refuses designs other than rep500/seed_base 2025 and
verifies the checked-out server SHA against the pushed SHA. U11 uses
`cache_flags=27` and `provider_flags=1`; all numerical libraries are restricted
to one thread per worker.
