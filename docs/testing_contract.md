# Testing Contract

This project compares screened ppro solvers against stable direct and sequential baselines for LLQR and TVCQR. The current regression focus is H-sequence agreement with the corresponding stable seq solver.

## Correctness Levels

Level 1: Objective/KKT exactness. The returned solution satisfies the full-sample weighted quantile-regression LP objective and KKT conditions at each evaluation point. This is the strongest correctness target.

Level 2: Estimate match against seq. The returned estimates match the corresponding stable seq solver, within numerical tolerance. For LLQR this means `ll_est` and related coefficients; for TVCQR this means `theta_ll_est`.

Level 3: `H_seq` set match against seq. For each evaluation index, compare the ppro `H_seq` row with the stable seq `H_seq` row as a set. Row order is not significant, but row contents and shape are significant.

Current regression tests require Level 3. Level 3 does not replace Level 1 or Level 2 as algorithmic goals, but it is the active regression gate for ppro basis behavior.

## Expected Regression Grid

Use this grid unless a task states a narrower scope:

- models: LLQR and TVCQR;
- cases: 1 and 2 for both LLQR and TVCQR;
- sample sizes: `n in c(200, 500, 1000)`;
- quantiles: `tau in c(0.2, 0.5, 0.8)`;
- seeds: `2025 + 1` through `2025 + 500`.

Large grids should be run through existing array/integrity scripts or a user-approved compute plan. Do not run large simulations during documentation-only tasks.

## Fallback Accounting

Fallback must be reported separately from ppro success. If a ppro wrapper returns a seq result through fallback, the result must preserve that fact with fields such as `fallback_triggered`, `fallback_reason`, and `returned_backend`.

Unless a test explicitly targets fallback behavior, a fallback result should not be counted as a ppro pass. It can be counted as a recovered run in a separate column.

## Existing Scripts

The repository currently includes these relevant scripts:

- `scripts/check_hseq_set_match.R`: compares saved result files for row-wise set equality of `H_seq`.
- `scripts/array/check_llqr_config_integrity.R`: checks LLQR partial array output for missing or failed replications.
- `scripts/array/check_tvcqr_config_integrity.R`: checks TVCQR partial array output for missing or failed replications.
- `scripts/array/scan_llqr_integrity.R`: scans LLQR array-output grids and writes rerun ID files.
- `scripts/array/scan_tvcqr_integrity.R`: scans TVCQR array-output grids and writes rerun ID files.
- `scripts/array/driver_llqr_array.R`: runs chunked LLQR replications.
- `scripts/array/driver_tvcqr_array.R`: runs chunked TVCQR replications.

## Command Templates

Check saved H-sequence set matches:

```sh
Rscript scripts/check_hseq_set_match.R --models=llqr,tvcqr --cases=1,2 --taus=0.2,0.5,0.8 --ns=200,500,1000 --rep=500
```

Check one LLQR partial-output configuration:

```sh
FASTQR_CASE=1 FASTQR_TAU=0.5 FASTQR_N=500 FASTQR_NUM_REP=500 Rscript scripts/array/check_llqr_config_integrity.R
```

Check one TVCQR partial-output configuration:

```sh
FASTQR_CASE=1 FASTQR_TAU=0.5 FASTQR_N=200 FASTQR_NUM_REP=500 Rscript scripts/array/check_tvcqr_config_integrity.R
```

Scan LLQR partial-output grid:

```sh
FASTQR_CASES=1,2 FASTQR_TAUS=0.2,0.5,0.8 FASTQR_NS=200,500,1000 FASTQR_NUM_REP=500 Rscript scripts/array/scan_llqr_integrity.R
```

Scan TVCQR partial-output grid:

```sh
FASTQR_CASES=1,2 FASTQR_TAUS=0.2,0.5,0.8 FASTQR_NS=200,500,1000 FASTQR_NUM_REP=500 Rscript scripts/array/scan_tvcqr_integrity.R
```

Do not cite or rely on non-existent scripts. If a desired test harness is missing, say so and either use one of the scripts above or propose adding a new test script in a separate implementation task.
