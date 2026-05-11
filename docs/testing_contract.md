# Testing Contract

This project compares screened ppro solvers against stable direct and sequential baselines for LLQR and TVCQR. The current regression focus is H-sequence agreement with the corresponding stable seq solver.

## Scope

This document defines how solver behavior is validated and reported.

It owns:

- correctness levels;
- active regression gates;
- fallback accounting;
- diagnostic output requirements;
- small smoke tests and full regression grids;
- command templates for existing scripts.

It does not define the ppro algorithm itself. Algorithm invariants belong in `docs/algorithm_contract.md`.

## Correctness Levels

Report these separately:

1. Level 1: **Objective/KKT exactness**. The returned solution satisfies the full-sample weighted quantile-regression LP objective and KKT conditions at each evaluation point. This is the strongest correctness target.

2. Level 2: **Estimate match against seq**. The returned estimates match the corresponding stable seq solver, within numerical tolerance. For LLQR this means `ll_est` and related coefficients; for TVCQR this means `theta_ll_est`.

3. Level 3: **`H_seq` set match against seq**. For each evaluation index, compare the ppro `H_seq` row with the stable seq `H_seq` row as a set. Row order is not significant, but row contents and shape are significant.

Current regression tests require Level 3. Level 3 does not replace Level 1 or Level 2 as algorithmic goals, but it is the active regression gate for ppro basis behavior.

## Default Full Regression Grid

Use this grid for final regression validation unless a task states a narrower scope:

- models: LLQR and TVCQR;
- cases: 1 and 2 for both LLQR and TVCQR;
- sample sizes: `n in c(200, 500, 1000)`;
- quantiles: `tau in c(0.2, 0.5, 0.8)`;
- seeds: `2025 + 1` through `2025 + 500`.

Large grids should be run through existing array/integrity scripts or a user-approved compute plan. Do not run large simulations during documentation-only tasks.

For plan-only, diagnosis-only, documentation-only, and small implementation tasks, do not start from the full grid.

Use the escalation ladder below:

1. one smallest reproducer or smoke case;
2. one focused LLQR or TVCQR case;
3. H_seq row-set comparison for the touched model;
4. small multi-case grid if the focused test passes;
5. full regression grid only when explicitly approved or clearly required.

## Fallback Accounting

Fallback must be reported separately from ppro success.

A fallback result may be counted as a recovered run, but it must not be counted as a successful ppro run unless the analysis explicitly targets fallback behavior.

Every table, timing summary, or regression report involving ppro must distinguish at least:

- true ppro success;
- fallback success;
- failure;
- missing or malformed output.

## Existing Scripts

The repository currently includes these relevant scripts:

- `scripts/check_hseq_set_match.R`: compares saved result files for row-wise set equality of `H_seq`.
- `scripts/array/check_llqr_config_integrity.R`: checks LLQR partial array output for missing or failed replications.
- `scripts/array/check_tvcqr_config_integrity.R`: checks TVCQR partial array output for missing or failed replications.
- `scripts/array/scan_llqr_integrity.R`: scans LLQR array-output grids and writes rerun ID files.
- `scripts/array/scan_tvcqr_integrity.R`: scans TVCQR array-output grids and writes rerun ID files.
- `scripts/array/driver_llqr_array.R`: runs chunked LLQR replications.
- `scripts/array/driver_tvcqr_array.R`: runs chunked TVCQR replications.

## Diagnostic Output Schema

Diagnostic scripts or debug-enabled solver calls should print or return at least:

- `model`
- `method`
- `backend`
- `case`
- `n`
- `tau`
- `seed`
- `rep_id`
- `eval_index`
- `h`
- `n_positive_weight`
- `positive_weight_rank`
- `Mm.factor`
- `threshold`
- `n_sl`
- `n_sh`
- `n_S`
- `n_prev_H_forced`
- `n_bad_signs`
- `bad_sign_indices`
- `bad_sign_classes`
- `threshold_doubled`
- `n_resolves`
- `simplex_converged`
- `it_num`
- `ierr` or kernel status code
- `fallback_triggered`
- `fallback_reason`
- `returned_backend`
- `H_seq_ppro`
- `H_seq_seq`

When logs are large, store full vectors in a structured object and print only counts plus the first few indices.

## Command Templates

Rebuild local Fortran shared libraries after editing the corresponding `.f90`
source and before running R-level validation:

```sh
cd src/fortran
gfortran -shared -fPIC -O3 -march=native -funroll-loops -ffast-math -o llqr_seq.so llqr_seq.f90
gfortran -shared -fPIC -O3 -march=native -funroll-loops -ffast-math -o llqr_ppro.so llqr_ppro.f90
gfortran -shared -fPIC -o tvcqr_seq.so tvcqr_seq.f90
gfortran -shared -fPIC -o tvcqr_seq_M_acc.so tvcqr_seq_M_acc.f90 -llapack -lblas
```

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
