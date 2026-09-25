# Reproducing the paper simulations

This directory freezes the U11 implementation underlying the four primary cases
in “A fast simplex algorithm for local linear quantile regression.” The R package
`fastllqr` 0.2.0 provides the corresponding public LLQR/TVCQR interfaces.

## Install the package

The prepared source distribution is `fastllqr_0.2.0.tar.gz`:

```sh
R CMD INSTALL fastllqr_0.2.0.tar.gz
```

Install the fixed GitHub reference (not a moving branch):

```r
install.packages("remotes")
remotes::install_github("zekkr/fastllqr@v0.2.0")
```

Release availability and validation status are recorded in `RELEASE_STATUS.md`.
R, a Fortran compiler, and R packages `quantreg` and `KernSmooth` are needed for
new simulation runs. `testthat` is needed only for package tests. Python 3 is
needed for archived aggregation. No additional R package is installed silently.

## Three reproduction modes

Run commands from the research repository root. Each output path must be new.

```sh
Rscript reproduction/reproduce.R --mode archived --output output/archived-rebuilt
Rscript reproduction/reproduce.R --mode smoke --output output/simulation-smoke
Rscript reproduction/reproduce.R --mode full --reps 500 --seed-base 2025 --output output/simulation-full
```

- `archived` rebuilds 48 configuration / 144 method summaries from included
  frozen aggregate inputs. These are independent of the current machine's
  runtime. It also includes the separately archived MST summaries; it does not
  rerun fitted datasets or re-estimate timing variability.
- `smoke` runs four cases at n=100, tau=0.5, two paired replications. This verifies
  the execution/merge/acceptance pipeline; its timing is not a paper result.
- `full` runs n=1000,2000,5000,10000 and tau=0.2,0.5,0.8, with 500 replications.
  It may be expensive. It is implemented but has not been rerun for this release.
  `--config 1` through `--config 48` selects one configuration for independent
  scheduler jobs, each with a distinct output path. Each job is single-process;
  parallelize configurations, not numerical-library threads.

`--profile portable` (default) uses -O2. `--profile paper-hpc` retains the historical
-O3 -march=native -funroll-loops -ffast-math flags. `--profile checked` enables
runtime checks. Full runs preserve the frozen fitting calls and timing boundary,
including method-order rotation, baseline-error-only regeneration (20 attempts),
seed 2025+r+1000000*(attempt-1), compact status and H-row-set comparisons. New
runs save actual R/package/system information and the compiler commands. They
never overwrite historical results. Timings on a different system are new
measurements, not bitwise reproductions of archived seconds.

## Case mapping and source provenance

| Paper case | Archived source | Raw model/case | Evaluation grid |
|---|---|---|---|
| 1 | 0863a92 logistic run | LLQR / 2 | observed x in [0.1,0.9] |
| 2 | 61dc4d6 original run | LLQR / 2 | all observed x, sorted |
| 3 | 61dc4d6 original run | TVCQR / 1 | i/n |
| 4 | 61dc4d6 original run | TVCQR / 2 | i/n |

The old Gaussian LLQR case is excluded from the paper mapping. Package support
for the Gaussian kernel is still tested. `frozen/` preserves the historical
harness. The `original_` driver and merger are from 61dc4d6; the logistic driver
and merger are from 0863a92. The core and entry are byte-identical in both commits.
`FROZEN_SHA256.json` and `src/fortran/u11/PROVENANCE.json` record fingerprints.
Do not use the historical cluster-specific submission scripts as portable setup
instructions; the top-level reproduction entry handles portable execution.

The canonical core is in `src/fortran/u11/`; the canonical public R adapters are
in `R/u11/`. The package contains independent buildable copies. Check synchronization:

```sh
python3 reproduction/check_source_sync.py ../fastllqr
```

The archived core is unchanged. Public adapters additionally support arbitrary
bandwidths, restore user grid order, expose diagnostics, and reject rank-deficient
positive support hidden by zero-weight interpolation padding. Therefore wrapper
validation/output overhead is included in package timings but not retroactively
in the historical experiment timings.

## Research API and validation

```r
source("R/u11_functions.R")
u11 <- load_u11()
fit <- u11$llqr_seq_ppro(x, y, kernel="epanechnikov", diagnostics=TRUE)
attr(fit, "solver_info")
```

The stable legacy seq and ppro entry points are retained; use the explicit U11
entry above for new research work. A package installation uses
`fastllqr::llqr_seq_ppro` and `fastllqr::tvcqr_seq_ppro` directly.

```sh
Rscript reproduction/validate_u11.R
Rscript reproduction/benchmark_profiles.R
```

Validation compares installed package, research adapter, archived U11 and frozen
stable seq; checks H row sets, estimates, full weighted objective and KKT; and
reports full-active recovery separately from seq fallback. The profile benchmark
is six alternating repetitions of four n=1000 examples, not rep500 evidence.

## macOS compiler configuration

Use a working Fortran toolchain compatible with your R installation. On the
validation machine, R's default `/opt/gfortran` path was absent. Homebrew GCC 15.1
also rejected the internal procedure interface with debug `-g`; without `-g`
the identical sources compiled. `Makevars.macos-homebrew.example` records the
local workaround; use `R_MAKEVARS_USER=/absolute/path/to/that/file` for package
installation/checks if your paths match. It is not a global R configuration
change and is not required by the historical GCC 12.2 profile. Independent
Linux validation remains pending.

## Coverage and unresolved historical metadata

See `ENVIRONMENT_AUDIT.md` for evidence versus missing metadata. The primary
four-case full rerun is implemented. The latest Supplementary Table S1 uses a separate, fully saved 100-replication
run; see `archive/mst_rep100_20260925/`. Its source, checkpoints, seeds and
per-replication statistics are included. Other historical MST summaries remain
archives and should not be described as newly rerun experiments.

## Full-weight formula LLQR baseline (2026-09-21)

Fresh `smoke` and `full` runs now default to the formula interface
`quantreg::rq(y ~ x, weights=w, tau=tau, method="br")`, retaining zero-weight rows.
Use `--llqr-baseline historical` explicitly to replay the earlier filtered
`rq.wfit` baseline. The `archived` mode still reconstructs the immutable original
archive; it must not be described as the new LLQR rerun.

The new paired HPC run and its SHA-gated submission commands are documented in
`experiments/llqr_rq_rep500/README.md`. Its results have paper case IDs 1 (logistic,
interior grid) and 2 (normal errors, full sorted grid), both using kernel case 2.
TVCQR and MST sources remain unchanged. Each run has distinct outputs and actual
worker environment records. For the updated manuscript, the bundle's
`output/llqr_rq_revision/source_manifest.json` selects the new LLQR and retained
TVCQR archives; `verify_manuscript.py --bundle BUNDLE` follows that manifest.

Rebuild the **updated** paper's LLQR summaries from the included fixed compressed
archive, then combine them with the retained TVCQR summaries:

```sh
python3 reproduction/rebuild_current_llqr.py output/current-paper-rebuilt
```

This creates a new directory, reconstructs statistics from all 12,000 paired
replication records, verifies 36,000 fits and compact H/recovery records, and
writes `current_paper_method_summary.csv` (144 rows). It does not rerun fitting.
For a small fresh run use:

```sh
Rscript reproduction/reproduce.R --mode smoke --reps 2 --llqr-baseline formula --output output/current-paper-smoke
```

For a full three-method **LLQR-only** rerun on the matching HPC environment use
`experiments/llqr_rq_rep500/submit.py` as documented there. This runs exactly the
24 LLQR configurations; the general `reproduce.R --mode full` includes all four
cases and is not the command used for this LLQR-only update.

## Supplementary Table S1 (2026-09-25)

The current six-row Table S1 was generated from 100 paired replications per
configuration. The fixed archive, exact DGPs, source snapshots, full checkpoint
set, validation report, and a small fresh-run command are in
[`archive/mst_rep100_20260925/`](archive/mst_rep100_20260925/README.md). The
runner is `scripts/run_llqr_multivar_iteration_simulation.R`; it directly loads
`R/llqr_multivar_functions.R`. The older April 2026 MST summaries elsewhere in
`archive/` are retained only as historical comparisons.
