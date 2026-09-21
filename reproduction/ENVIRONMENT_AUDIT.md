# Environment and direct-fitting audit — 2026-09-21

## Evidence tied to the reported runs

- Original run: `u11_rep500_seed2025_61dc4d6_20260911`, retained paper cases 2–4.
- Logistic run: `llqr_c2_logistic_rep500_seed2025_0863a92_20260914_112526`, paper case 1.
- All 72,000 method records in the retained 48 configurations identify Intel(R)
  Xeon(R) Gold 6258R CPU @ 2.70GHz (54,000 original + 18,000 logistic).
- Both archived build records: Linux 3.10.0-1160.el7.nosig.x86_64, GCC 12.2.0.
  The kernel is recorded; this is not independently an OS distribution version.
- Both submission scripts specify R 4.3.1 and GCC 12.2.0; actual worker
  `sessionInfo()` was not captured. Both use one numerical-library thread.
- Optimized flags: -O3 -march=native -funroll-loops -ffast-math, plus shared-library
  -fPIC -ffree-line-length-none; linkage uses R-configured LAPACK/BLAS.
- Exact September quantreg, KernSmooth, BLAS/LAPACK build and OS distribution
  versions have not been established from the run archives. An August memory
  record mentions quantreg 5.94 and CentOS 7; it is supporting context, not proof
  for September. A read-only SSH attempt on 2026-09-21 was reset before login.

## Frozen direct-fitting calls

LLQR (`frozen/R/llqr_functions.R`, `llqr_local_fit`): fixed design cbind(1,x),
active = w > 0, `quantreg::rq.wfit(..., weights=w[active], tau=tau, method="br")`.
The fitted value is intercept + z*slope. No tolerance argument is supplied.

TVCQR (`frozen/R/tvcqr_functions.R`, `tvc_rq`):
`quantreg::rq(y ~ x2, tau=tau, weights=w)` with the local centered-time design,
formula intercept and full weight vector. No method/tolerance argument and no
explicit zero-weight row removal in the wrapper.

The CRAN quantreg 5.94 source was inspected:
https://cran.r-project.org/src/contrib/Archive/quantreg/quantreg_5.94.tar.gz
Its rq default method is "br"; rq.wfit multiplies design and response by weights,
without removing zero-weight rows; rq.fit.br sets internal tolerance to
`.Machine$double.eps^(2/3)` (approximately 3.666853e-11). These are also the
locally installed 5.97 behaviors. This comparison does not prove that 5.94 or
5.97 was the September installed version. Hence the manuscript reports no
override/default behavior, and does not assert an unverified numerical baseline
tolerance or package version. The U11 tolerance is independently fixed at 1e-14.

## Release validation environment

Local R 4.3.3 / quantreg 5.97 / KernSmooth 2.23.24 / Homebrew GCC 15.1 on macOS.
The generated run environment and compiler command records are retained with the
smoke output. They are new validation metadata, not historical substitutions.
