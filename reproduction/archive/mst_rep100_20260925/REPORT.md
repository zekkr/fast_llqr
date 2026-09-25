# MST Table S1 local rep100 result

Date: 2026-09-25 (Asia/Shanghai). Current manuscript was not edited in this run.

## Experiment

- Models: Case 1 mean `1 + 1.2*x1 - 0.9*x2`; Case 2 mean `0.5 + 0.8*x1 + 0.6*x2^2`. The covariate and noise laws are in the saved runner source.
- Six configurations: cases 1/2, n=300/500/800, p=4, tau=0.5, 100 paired replications each.
- Methods: custom R simplex cold start (`direct`) and full-LP parent-tableau warm start (`seq-MST`).
- Seed per dataset: `20260405 + 100000*case + 1000*n + 40 + replication`.

## Table S1 values

| Case | n | direct | seq-MST | MST vs. direct |
|---:|---:|---:|---:|---:|
| 1 | 300 | 64,944 | 21,616 | 66.6% |
| 1 | 500 | 176,403 | 53,762 | 69.5% |
| 1 | 800 | 442,840 | 125,013 | 71.7% |
| 2 | 300 | 48,606 | 15,401 | 68.3% |
| 2 | 500 | 132,601 | 37,570 | 71.6% |
| 2 | 800 | 334,075 | 85,657 | 74.3% |

The first two numeric columns are the means of each method’s total iterations, rounded to integers. The last column is the mean of `100*(1 - MST/direct)` across paired replications, rounded to one decimal place.

## Validation

- 600 saved paired datasets and 1200 method rows; all six configurations have replications 1–100 and no duplicate keys.
- Maximum direct/MST fitted-value difference: `8.302691867e-12` (manuscript display: `8.30 × 10^-12`).
- Mean of the 600 within-dataset mean MST parent-child edge lengths: `0.9110037101820336` (manuscript display: `0.91`).
- Maximum pointwise simplex iterations: 1475, below the configured limit 20,000. No recorded accuracy or max-iteration failures.
- The saved X/Y regenerated bit-for-bit for all 600 seeds from the saved generator source.
- For the overlapping first 10 replications, all 120 direct/MST iteration records matched the April 2026 archive exactly.
- An independent Python calculation from `raw_all.csv` matched all six method means and paired reductions.
- All six configuration manifests have identical source hashes and RNG kinds. R session details and source snapshots are saved in each configuration directory.

## Reproduction

From the repository root, run each configuration into its own empty directory (or add `--resume` for a previously started one):

```sh
for case_id in 1 2; do
  for n in 300 500 800; do
    Rscript scripts/run_llqr_multivar_iteration_simulation.R \
      --reps 100 --cases "$case_id" --n-values "$n" \
      --output-dir "results/llqr_multivar_iterations/mst_rep100_local_20260924/case${case_id}_n${n}"
  done
done
Rscript results/llqr_multivar_iterations/mst_rep100_local_20260924/aggregate.R
```

The aggregate script writes `raw_all.csv`, `diagnostics_all.csv`, `paired_all.csv`, `table_s1_all.csv`, and `aggregate_summary.rds` in this directory. `SHA256SUMS` covers the completed run artifacts.
