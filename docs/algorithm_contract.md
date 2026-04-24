# Algorithm Contract

This repository supports simulations for "A fast simplex algorithm for local linear quantile regression." The active screened warm-start implementations are `llqr_seq_ppro`, `tvcqr_seq_ppro`, and their Fortran-backed kernels. These implementations must remain real screened ppro algorithms. They must not be replaced by calls to the corresponding seq solver.

## Screened Warm-Start Invariant

For every evaluation point, a ppro solver must return the same full-sample quantile-regression solution that the unscreened local LP would return, up to numerical tolerance. Screening is only a dimension-reduction device. It cannot change the accepted solution.

The required flow is:

1. Solve the first evaluation point on the full sample.
2. For each later evaluation point:
   - use residuals from the previous accepted solution to classify observations;
   - force previous H/basis observations into the uncertain set;
   - build current-point weighted aggregates for omitted sure-negative and sure-positive observations;
   - solve the reduced simplex problem;
   - verify every omitted observation against the reduced solution;
   - if verification finds bad signs, expand the reduced problem and solve again;
   - if screening cannot be certified, recover the full-sample problem.
3. Never silently accept an unconverged reduced LP.
4. Never silently fall back to seq.

## Terms

`sl`: The sure-negative set. These observations are predicted, from the previous residuals and the current threshold, to keep negative residual sign at the current evaluation point.

`sh`: The sure-positive set. These observations are predicted to keep positive residual sign at the current evaluation point.

`S`: The uncertain set solved explicitly in the reduced LP. It contains observations near the fitted quantile, observations whose sign is not certified by the screening threshold, any bad signs added after verification, and all previous H/basis observations.

Bad signs: Omitted observations whose residual sign under the reduced solution disagrees with their assigned sure-negative or sure-positive class, or otherwise violates the sign certificate used to aggregate them.

H/basis observations: The observations represented by the active simplex basis/H sequence. Previous H observations must be forced into `S` before solving the next reduced problem, because dropping them can break warm-start validity and cause `H_seq` mismatches.

Weighted aggregates: Current-evaluation pseudo-observations representing the exact objective contribution of omitted `sl` and `sh` observations conditional on their certified signs. Aggregates must be rebuilt whenever the evaluation point, weights, `sl`, `sh`, `S`, threshold, or bad-sign set changes.

Reduced LP: The simplex problem containing explicit rows for `S` plus aggregate rows for omitted certified signs. It is a computational shortcut, not a different estimator.

Verification: The post-solve check over all omitted observations. Verification decides whether the reduced solution is certifiable as the full-sample solution.

## Verification Failure

When verification fails, the solver must not return the failed reduced solution. It must do one of the following:

- add the bad-sign observations back into `S`, rebuild aggregates, and re-solve;
- enlarge or double the screening threshold, recompute `sl`, `sh`, `S`, counts, aggregates, and warm-start state as needed, then re-solve;
- abandon screening for that evaluation point and solve the full-sample problem.

After any threshold or set change, all dependent quantities must be recomputed. Reusing stale aggregate rows or stale counts is a correctness bug.

## Fallback

Production fallback is allowed only as an explicit, reported recovery path. A result produced by fallback must include fields such as:

- `fallback_triggered = TRUE`
- `fallback_reason`
- `returned_backend = "seq_fallback"`

Debugging fallback is different. During diagnosis, fallback may be disabled so kernel failures surface directly. Debugging output should report where certification failed, which rows were added, and whether the kernel converged.

Fallback output must not be counted as ppro success unless the test or analysis explicitly says it is measuring fallback behavior.

## Non-Negotiable Rule

`ppro` cannot be implemented by simply calling `seq`. A ppro method may call seq only through an explicit fallback branch that records and reports the fallback status.
