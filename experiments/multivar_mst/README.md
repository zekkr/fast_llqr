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

On the SHA-matched HPC checkout, run the preflight and the preregistered
threshold pilot before the formal run. The pilot uses seeds 9026--9030, which
do not overlap the formal seeds, and compares threshold factors 0.1, 0.2, 0.5,
and 1.0. The formal submission reads the objectively selected factor from the
pilot archive.

```sh
python3 experiments/multivar_mst/submit.py preflight --tag <tag> --sha <sha>
python3 experiments/multivar_mst/submit.py pilot --tag <tag> --sha <sha> --preflight <preflight-tag>
python3 experiments/multivar_mst/submit.py formal --tag <tag> --sha <sha> --preflight <preflight-tag> --pilot <pilot-tag>
python3 experiments/multivar_mst/verify_results.py results/hpc/multivar_mst_runs/<tag>
```
