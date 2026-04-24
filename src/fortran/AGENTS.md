# Fortran Kernel Instructions

These instructions apply to files under `src/fortran/`.

The Fortran kernels implement performance-critical versions of the screened warm-start ppro and stable seq algorithms. Treat them as numerical kernels, not as wrappers. Do not replace a ppro kernel with a call or behavior equivalent to the seq implementation.

## Required Behavior

- No silent seq replacement. A ppro routine must implement the screened ppro algorithm. Any seq fallback must happen through an explicit wrapper or explicit fallback branch that reports status.
- No silent acceptance of nonconvergence. If the simplex loop reaches `maxit`, fails to find a valid pivot, or detects invalid state, return an explicit error/status flag.
- No stale aggregate rows. Rebuild aggregate rows whenever weights, evaluation point, threshold, `sl`, `sh`, `S`, bad signs, or reduced problem dimensions change.
- Recompute counts after changing `sl`, `sh`, or `S`. Do not reuse old reduced dimensions, offsets, or aggregate indices after threshold doubling or bad-sign insertion.
- Return explicit error/status flags. Prefer existing `ierr`-style outputs where present; do not hide kernel failures behind plausible numeric output.
- Preserve R/Fortran interface signatures unless a task explicitly asks for an interface change. If an interface change is required, update the R wrapper and document the compatibility impact.

## Validation

Fortran changes must be validated against both:

- the R reference implementation for the same algorithm path, where available;
- the stable seq solver's `H_seq`, comparing each evaluation row as a set rather than an ordered vector.

Fallback results must be reported separately from ppro kernel successes.
