# Implementation and validation — 2026-09-21

## Changes

- Formal U11 core: `src/fortran/u11/`, byte-identical to both reported rep500
  snapshots. R research API: `R/u11_functions.R` and `R/u11/`.
- fastllqr 0.2.0: shared U11 backend, unchanged default estimate fields,
  optional full audit and diagnostics, always-present provenance attribute.
- Public wrapper safeguard: if the transported interpolation basis contains
  zero-weight rows, check the full positive support rank before returning.
  A deterministic duplicated-x/compact-support test demonstrated why this is
  needed; the archived Fortran kernel is not changed.
- Stable seq, old research entry points, DGPs, historical outputs and user-edited
  AGENTS/contract files were not modified.
- Three advisor reminders replaced by blue text. The original Main package
  sentence and all other red revisions are preserved. Historical metadata
  gaps are explicitly stated, not filled using the local environment.

## Commands actually run

From the research repository root (paths under tmp are retained build evidence):

```sh
python3 reproduction/check_source_sync.py ../fastllqr
python3 reproduction/check_frozen.py
R_MAKEVARS_USER="$PWD/tmp/fastllqr-release/Makevars.local" R CMD INSTALL --preclean --library=tmp/fastllqr-release/library tmp/fastllqr-release/fastllqr
env -u LC_ALL R_LIBS_USER="$PWD/tmp/fastllqr-release/library" Rscript --vanilla -e 'library(fastllqr); testthat::test_dir("tmp/fastllqr-release/fastllqr/tests/testthat",package="fastllqr",reporter="summary")'
env -u LC_ALL R_LIBS_USER="$PWD/tmp/fastllqr-release/library" Rscript --vanilla reproduction/validate_u11.R
env -u LC_ALL R_LIBS_USER="$PWD/tmp/fastllqr-release/library" Rscript --vanilla reproduction/benchmark_profiles.R
Rscript --vanilla reproduction/reproduce.R --mode archived --output tmp/fastllqr-release/archived_complete
Rscript --vanilla reproduction/reproduce.R --mode smoke --output tmp/fastllqr-release/smoke_final
R_MAKEVARS_USER="$PWD/tmp/fastllqr-release/Makevars.local" R CMD build tmp/fastllqr-release/fastllqr
R_MAKEVARS_USER="$PWD/tmp/fastllqr-release/Makevars.local" R_LIBS_USER="$PWD/tmp/fastllqr-release/library" LC_ALL=en_US.UTF-8 R CMD check --no-manual --no-build-vignettes --output=tmp/fastllqr-release fastllqr_0.2.0.tar.gz
python3 -B 'paper/fast_llqr_submission_v2_bundle_20260919_TW 2/paper/scripts/check_bundle_summaries.py' --report reproduction/reports/paper_numeric_check.json
```

The tested tarball was then moved to `output/u11_release_20260921/`.
Both manuscripts were compiled with
`latexmk -pdf -interaction=nonstopmode -halt-on-error <filename>.tex`
from the corrected advisor bundle's paper directory. Main page 3 and Supplement
pages 29–30 were rendered with pdftoppm and visually inspected; pypdf checked URLs.

## Outcomes

- Source synchronization: 9 R/Fortran files identical; both historical core
  fingerprints unchanged. Frozen harness: 32 files; archived inputs: 22 files.
- Package tests: all six test files pass, including unchanged historical H-row
  fixtures and the new compact-support rank-deficiency witness. Tests/examples
  also pass inside R CMD check. Full log: `reports/package_tests.txt`.
- R CMD check: **0 errors, 1 warning**. The warning is the non-portable Fortran
  flag `-ffree-line-length-none`, retained to compile unchanged archived sources.
  This is a GitHub source release candidate, not a clean CRAN submission.
- Cross-implementation gate: **60 datasets, 8,652 evaluation points**, no H-row-set
  mismatch; exact equality of public research/package results and archived/core
  diagnostic paths. Maximum estimate error against direct fit 1.136868e-13;
  maximum normalized objective gap 5.752049e-14; KKT violation 5.662137e-15.
- **1,481 full-active recoveries, zero seq fallbacks** in the small validation
  grid. These recoveries are counted, not hidden or relabeled as fallback.
- Portable-package versus historical-flag core: four n=1000 configurations,
  six alternating timed repetitions per method; estimate/H agreement passes.
  The package includes wrapper overhead and was slower in these small examples.
  See `reports/profile_benchmark.csv`; this is not a replacement for rep500.
- Smoke execution: four configurations × two replications × three methods =
  **24 accepted final fits**, complete merges, no missing/malformed records.
- Archived summaries: 48 configurations / 144 methods / 500 replications each;
  archived U11 totals: **14,060 full-active recoveries, zero seq fallback**.
- Manuscript numerical audit: **486 numeric checks + 720 mapping checks**, zero
  discrepancy. No paper number changed. Main remains 18 pages; Supplement 33.
  Three blue insertions; no pending calls; both GitHub URL annotations present.
  Main has 36 pre-existing overfull-vbox warnings under the Biometrika line-number
  output routine, exactly the same count on a clean compile of the original
  advisor source; zero overfull-hbox and zero undefined references. The edited
  regions render without clipping or overlap.

## Limits and failed intermediate checks

- Initial local compiler configuration referenced an absent /opt/gfortran;
  Homebrew GCC 15.1 with `-g` rejected the internal procedure interface. The
  task-local Makevars workaround is documented; no system config was changed.
- Existing local waldo 0.5.2 failed to load with the local glue version. A newer
  waldo was installed only in the task's isolated R library; tests were rerun.
- An intermediate support-guard implementation dropped matrix dimensions;
  corrected before the final install, full gate and R CMD check above.
- No independent Linux check, full rep500 rerun, MST fitted-data rerun or new HPC
  submission was performed. Full simulation CLI exists and the same driver path
  passed the smoke test; this does not certify a completed 500-replication run.
- September quantreg/BLAS versions and public repository access are unresolved.
  No public release or repository visibility change is asserted.
