# R Wrapper and Reference Instructions

These instructions apply to files under `R/`.

The R code contains reference implementations, method factories, simulation helpers, and wrappers around Fortran kernels. Keep wrapper behavior transparent: callers must be able to tell whether a result came from ppro, seq, or fallback.

## Required Behavior

- Wrappers must expose fallback status. Use fields such as `fallback_triggered`, `fallback_reason`, and `returned_backend` when fallback occurs.
- Wrappers must not hide kernel failures in timing mode. Timing output should preserve failure status instead of silently substituting seq results.
- Certification and fallback should be separate options. A certification failure may trigger a configured fallback, but the returned object must record both the certification outcome and fallback outcome.
- Method factories must not silently swap ppro for seq. A method named `*_ppro*` must run the ppro path unless an explicit, reported fallback branch is taken.
- Diagnostics must preserve provenance. Saved results and diagnostic objects should identify whether estimates, residuals, and `H_seq` came from ppro or from fallback.

## Testing Expectations

For ppro changes, compare against the corresponding stable seq implementation:

- LLQR ppro against `llqr_seq` or `llqr_seq_fortran_wrapper`, as appropriate.
- TVCQR ppro against `tvcqr_seq` or `tvcqr_seq_fortran_wrapper`, as appropriate.

For current regression tests, `H_seq` rows are compared as sets at each evaluation index. Estimate agreement and objective/KKT checks remain stronger correctness targets and should be used when the task scope requires them.
