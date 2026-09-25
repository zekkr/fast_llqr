# Supplementary Table S1: 100 paired replications

This directory publishes the fixed results behind Supplementary Table S1. The
six table rows, 600 paired datasets, and numerical agreement checks are described
in [REPORT.md](REPORT.md). `table_s1_all.csv` contains the displayed summaries;
`raw_all.csv`, `paired_all.csv`, and `diagnostics_all.csv` provide the underlying
per-replication statistics. `complete_run.tar.gz` includes the saved datasets,
checkpoints, source snapshots, manifests, session information, and original
aggregation script. `SHA256.json` records hashes of the files in this directory.

To inspect and independently regenerate the table from the saved checkpoints,
extract the archive in the repository root. It creates
`mst_rep100_local_20260924/`; move that directory to
`results/llqr_multivar_iterations/` before running its `aggregate.R`. The
aggregate script uses that fixed relative path and checks all six configurations,
replication keys, source hashes, seeds, and fit tolerances.

To produce **new** results, run the source in
`scripts/run_llqr_multivar_iteration_simulation.R` from the repository root.
Each configuration needs a separate empty output directory. The exact commands
and seed rule are in [REPORT.md](REPORT.md); a small execution check is:

```sh
Rscript scripts/run_llqr_multivar_iteration_simulation.R --reps 1 --cases 1 --n-values 300 --output-dir results/mst_smoke_case1_n300
```

The experiment compares independent direct fits with full-LP MST warm starts.
It does not use screening or aggregation. Machine-specific iteration counts
should be checked against the archived source and tolerance settings; elapsed
times are not Table S1 outcomes.
