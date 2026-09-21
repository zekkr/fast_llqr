# Multivariate Gaussian LLQR MST experiment

This isolated experiment compares independent `quantreg::rq(method="br")`
fits with the experimental compiled `seq-screen-MST` solver. It does not alter
the public package or the archived U11 sources.

Settings: two manuscript data-generating models, `p=4`, `tau=0.5`, all sample
points as evaluation points, `n=500,1000,2000`, 100 paired replications, and
`seed_base=2025`. The Gaussian bandwidth and empirical threshold are documented
in `solver.R`. Timed screened calls include bandwidth selection, MST construction,
weights, all solves, verification, repair/recovery, and output construction.

Local validation:

```sh
bash experiments/multivar_mst/build.sh
MST_BUILD_MODE=checked Rscript experiments/multivar_mst/tests/test_solver.R
MST_BUILD_MODE=optimized Rscript experiments/multivar_mst/tests/test_solver.R
```

On the SHA-matched HPC checkout, run the preflight before the formal run:

```sh
python3 experiments/multivar_mst/submit.py preflight --tag <tag> --sha <sha>
python3 experiments/multivar_mst/submit.py formal --tag <tag> --sha <sha> --preflight <preflight-tag>
python3 experiments/multivar_mst/verify_results.py results/hpc/multivar_mst_runs/<tag>
```
