# Fast local quantile regression research

This repository contains the research implementation and simulation sources for
*A fast simplex algorithm for local quantile regression with exactness and
dimension-reduction guarantees*. The installable R package is
[`fastllqr`](https://github.com/zekkr/fastllqr); the canonical U11 solver sources
are in `src/fortran/u11/` and `R/u11/` here.

## Reproduce the paper

See [reproduction/README.md](reproduction/README.md) for the fixed-source
simulation commands, archived numerical inputs, compiler profiles, and known
limits of historical environment metadata. In particular:

```sh
python3 reproduction/rebuild_current_llqr.py output/current-paper-rebuilt
Rscript reproduction/reproduce.R --mode smoke --reps 2 --llqr-baseline formula --output output/current-paper-smoke
```

Supplementary Table S1 has a separate [100-replication archive and rerun
instructions](reproduction/archive/mst_rep100_20260925/README.md). The archived
results are fixed historical measurements; fresh timings depend on the machine,
compiler, and R environment.

## Source and validation

Use `R/u11_functions.R` for the U11 research interface. The legacy LLQR and TVCQR
functions remain in `R/`. The package's buildable copies of the U11 sources can
be checked with `python3 reproduction/check_source_sync.py ../fastllqr` when the
repositories are checked out beside one another. See
[reproduction/VALIDATION.md](reproduction/VALIDATION.md) for the validation record.
