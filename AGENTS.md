# AGENTS.md

## Project goal

This repository supports simulations for the paper
“A fast simplex algorithm for local linear quantile regression.”

The main research objects are:
- `llqr_seq_ppro`
- `tvcqr_seq_ppro`
- their Fortran-backed implementations

The `quantreg`/`quantdr`-based solvers (`llqr`, `tvc_rq`) are baseline/oracle solvers. Our sequential solvers (`llqr_seq`, `tvcqr_seq`) and their Fortran-backed implementations are now stable and fast. Their `H_seq` outputs are treated as the current regression-test standard for whether the `seq_ppro` solvers have reached the expected solution path.

Do not replace a ppro implementation by a seq implementation.

## Algorithm contract

The screened sequential ppro solvers must implement:

1. Solve the first evaluation point on the full sample.
2. For each next evaluation point:
   - use previous residuals to form sure-negative `sl`, sure-positive `sh`,
     and uncertain set `S`;
   - force previous H/basis observations into `S`;
   - build weighted aggregates at the current evaluation point:
     - sure-positive residual group uses an `u_H` aggregate and objective `tau`;
     - sure-negative residual group uses a `v_L` aggregate and objective `1 - tau`;
   - build and solve the reduced simplex problem;
   - verify all omitted observations for sign violations;
   - if violations occur, add bad signs back to `S` or enlarge threshold and re-solve;
   - in the worst case, solve the full sample.
3. Never silently accept an unconverged reduced LP.
4. Never silently fall back to seq. If fallback is used, return and report
   `fallback_triggered = TRUE` and `returned_backend = "seq_fallback"`.

## Correctness definitions

Report these separately:
- objective/KKT exactness;
- estimate match against seq;
- H_seq set match against seq.

For current regression tests, H_seq set equality against seq is required.

## Do-not rules

- Do not modify `llqr_seq` or `tvcqr_seq` unless explicitly asked.
- Do not change simulation data-generating mechanisms unless explicitly asked.
- Do not change tests to make a broken implementation pass.
- Do not remove diagnostics that distinguish ppro from seq fallback.
- Do not edit `paper/main_v3.tex` unless the task is explicitly about the paper.

## Files to inspect before editing

- `R/llqr_functions.R`
- `R/tvcqr_functions.R`
- `src/fortran/llqr_ppro.f90`
- `src/fortran/tvcqr_seq_M_acc.f90`
- `scripts/check_hseq_set_match.R`

## Verification expectations

Before proposing a code change, explain:
- which algorithm invariant is being changed or preserved;
- which files and functions are affected;
- how the change will be tested.

After changing code, report:
- exact commands run;
- exact test outcomes;
- whether any fallback occurred;
- whether H_seq set match passed.