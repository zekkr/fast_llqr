# Algorithm Contract

This repository supports simulations for "A fast simplex algorithm for local linear quantile regression." The active screened warm-start implementations are `llqr_seq_ppro`, `tvcqr_seq_ppro`, and their Fortran-backed kernels. These implementations must remain real screened ppro algorithms. They must not be replaced by calls to the corresponding seq solver.

## Scope and source of truth

This document is the source of truth for screened ppro algorithm invariants.

It governs:

- `llqr_seq_ppro`;
- `tvcqr_seq_ppro`;
- their Fortran-backed kernels;
- R wrappers or method factories that expose these solvers.

This document does not define task prompt formats, regression grids, or command templates. Those belong in `docs/codex_task_template.md` and `docs/testing_contract.md`.

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

## H_seq is a regression gate, not the estimator definition

`H_seq` row-set agreement with stable seq is the current regression gate for basis-path behavior. It is not the mathematical definition of the estimator.

The estimator-level target remains the full-sample weighted quantile-regression solution. If `H_seq` agrees but objective/KKT or estimate-level checks fail, the implementation is still wrong.

## seq_ppro tableau initialization principle

Both `llqr_seq_ppro` and `tvcqr_seq_ppro` solve the same weighted quantile regression LP as their corresponding `seq` solvers:

`min tau * w^T u + (1 - tau) * w^T v`

subject to `A beta + u - v = y`, `u, v >= 0`, and free `beta`.

The reduced subsample variable order is:

- `1:q` for beta/free variables;
- `q + i` for `u_i`;
- `q + ms + i` for `v_i`.

For LLQR, `q = nvar + 1`, `A = cbind(1, x)`, the evaluation index is `rd`, and the output uses `crossprod(c(1, z[rd]), estimate)` plus `estimate[1 + nvar]`. For TVCQR, `q = 2 * (nvar + 1)`, `A = cbind(1, x, time_index * cbind(1, x))`, the evaluation index is `eva_t`, and the output uses `estimate[1:(nvar + 1)] + (eva_t / m) * estimate[(nvar + 2):(2 * (nvar + 1))]`.

The first evaluation point is solved on the full sample. After convergence, the solver stores the original observation indices of the zero-residual/interpolating rows in `H_seq`: LLQR uses `H <- r1 - 1 - nvar`, and TVCQR uses `H <- r1 - 2 - 2 * nvar`.

For every later evaluation point, including LLQR `rd == 2`, the solver starts from the screened subsample. It uses `r_prev` to build `M`, `sl <- r < -M`, `sh <- r > M`, and `not_jl_or_jh <- !(sl | sh)`. Previous basis rows are forced back into the retained subsample before aggregate rows are built. In code, `H <- match(H_seq[previous, ], idx_not_jl_or_jh)` maps original `H_seq` indices into current subsample row positions.

Low aggregate rows from `sl` are built from `glob.wx = colSums(A[sl, ] * w[sl])` and `glob.wy = sum(y[sl] * w[sl])`. High aggregate rows from `sh` are built from `ghib.wx = colSums(A[sh, ] * w[sh])` and `ghib.wy = sum(y[sh] * w[sh])`. Aggregate rows append `ws <- c(ws, 1)` because the original kernel weights have already been folded into aggregate `A` and `y`.

Inside the reduced problem, `idpos <- which(r[not_jl_or_jh] > 0)` identifies individual positive residual rows and `idneg <- which(r[not_jl_or_jh] < 0)` identifies individual negative residual rows. Current `H` rows are removed from `idpos` and `idneg`. Positive rows use `u` basis variables, negative rows use `v` basis variables, and the nonbasic residual pairs for the `H` rows are `r1 <- H + q` and `r2 <- r1 + ms`.

`Hbar` contains the basis residual rows outside `H`, including any aggregate rows. `IBs` starts with beta variables `1:q`, then the residual basis variables for `Hbar`. `P` records the residual sign of each `Hbar` basis column:

- individual positive residual rows use `P = +1`;
- individual negative residual rows use `P = -1`;
- low aggregate rows from `sl` represent the negative-residual side and use `P = -1`;
- high aggregate rows from `sh` represent the positive-residual side and use `P = +1`.

The objective coefficients for `Hbar` are stored in `lambda`: `tau * ws[idpos]` for positive individual rows, `(1 - tau) * ws[idneg]` for negative individual rows, `1 - tau` for low aggregate rows, and `tau` for high aggregate rows.

The H-based block initialization uses:

- `XH = gammaxs[H, ]`;
- `XHb = gammaxs[Hbar, ]`;
- `yH = bs[H]`;
- `yHb = bs[Hbar]`;
- `xhinv = solve(XH)` or the LLQR two-by-two helper.

Conceptually the basis matrix is `B = [[XH, 0], [XHb, P]]`, so `B^{-1} = [[XH^{-1}, 0], [-P XHb XH^{-1}, P]]`. In TVCQR this is implemented with matrix `P`. In LLQR, `P` is a sign vector and the same operation is implemented by row scaling:

- `Pxhbar <- gammaxs[Hbar, ] * P`;
- `Pxhbarxhinv <- Pxhbar %*% xhinv`;
- `gammaxs <- rbind(xhinv, -Pxhbarxhinv)`;
- `bs <- c(xhinv %*% bs[H], -Pxhbarxhinv %*% bs[H] + bs[Hbar] * P)`.

The final reduced-cost row is `tau * ws[H] + t(lambda) %*% Pxhbarxhinv`. After this initialization, the solver enters the same residual-pair simplex loop used by the sequential method.

For LLQR `rd == 2`, the solver must use the screened subsample and aggregate rows described above, but it must still use explicit H-based block initialization because cached transformed rows do not exist yet. For LLQR `rd > 2`, and analogously for later TVCQR points, the solver reuses cached transformed rows `idx_Hbar_pos`, `idx_Hbar_neg`, `gammaxs.pos`, `gammaxs.neg`, `bs.pos`, and `bs.neg` when rows remain in the same sign class. New individual rows and aggregate rows are recomputed, then the final reduced-cost row is rebuilt.

## Terms

`sl`: The sure-negative set. These observations are predicted, from the previous residuals and the current threshold, to keep negative residual sign at the current evaluation point.

`sh`: The sure-positive set. These observations are predicted to keep positive residual sign at the current evaluation point.

`S`: The uncertain set solved explicitly in the reduced LP. It contains observations near the fitted quantile, observations whose sign is not certified by the screening threshold, any bad signs added after verification, and all previous H/basis observations.

Bad signs: Omitted observations whose residual sign under the reduced solution disagrees with their assigned sure-negative or sure-positive class, or otherwise violates the sign certificate used to aggregate them.

H/basis observations: The observations represented by the active simplex basis/H sequence. Previous H observations must be forced into `S` before solving the next reduced problem, because dropping them can break warm-start validity and cause `H_seq` mismatches.

Weighted aggregates: Current-evaluation pseudo-observations representing the exact objective contribution of omitted `sl` and `sh` observations conditional on their certified signs. Aggregates must be rebuilt whenever the evaluation point, weights, `sl`, `sh`, `S`, threshold, or bad-sign set changes.

Reduced LP: The simplex problem containing explicit rows for `S` plus aggregate rows for omitted certified signs. It is a computational shortcut, not a different estimator.

Verification: The post-solve check over all omitted observations. Verification decides whether the reduced solution is certifiable as the full-sample solution.

## Bounded-kernel local LP degeneracy

LLQR case 2 uses a bounded Epanechnikov kernel and can have an evaluation point with too few positive-weight observations to define the local LP. For one-dimensional LLQR, `q = nvar + 1 = 2`, so an evaluation point with only one positive-weight observation has an effective design that is rank deficient before any ppro screening is applied.

This is a full-sample local LP degeneracy, not a screened-ppro certification failure. Enlarging or doubling the screening threshold can add more explicit rows from the screened problem, but it cannot create positive kernel weights outside the bounded kernel support. Such a case must not be silently accepted as ppro success. A recovery path must either report an explicit failure or use an explicit fallback branch with recorded provenance.

## Verification Failure

When verification fails, the solver must not return the failed reduced solution. It must do one of the following:

- add the bad-sign observations back into `S`, rebuild aggregates, and re-solve;
- enlarge or double the screening threshold, recompute `sl`, `sh`, `S`, counts, aggregates, and warm-start state as needed, then re-solve;
- abandon screening for that evaluation point and solve the full-sample problem.

After any threshold or set change, all dependent quantities must be recomputed. Reusing stale aggregate rows or stale counts is a correctness bug.

When a candidate has no bad signs and the same reported `H` rows are in range,
distinct, and full rank, a ppro implementation may take one numerical rescue
step before rejection: refit the estimate from that same `H`, then run
certification once more. This same-`H` refit is only a stabilization step; it
must not change `H`, relax tolerances, bypass certification, or replace the
existing reject, expand, or fallback path if recertification still fails.

## Certified ppro output

A ppro output is certified only if one of the following holds:

1. the reduced LP passes verification at every evaluation point;
2. after one or more repair steps, the enlarged reduced LP passes verification;
3. screening is abandoned for an evaluation point and the full-sample LP is solved explicitly;
4. an explicit fallback branch returns a seq result and records fallback status.

A failed reduced solve, nonconverged simplex state, stale aggregate solve, or unverified screened solution is not certified output.

## Fallback

Production fallback is allowed only as an explicit, reported recovery path. A result produced by fallback must include fields such as:

- `fallback_triggered = TRUE`
- `fallback_reason`
- `returned_backend = "seq_fallback"`

Debugging fallback is different. During diagnosis, fallback may be disabled so kernel failures surface directly. Debugging output should report where certification failed, which rows were added, and whether the kernel converged.

## Non-Negotiable Rule

`ppro` cannot be implemented by simply calling `seq`. A ppro method may call seq only through an explicit fallback branch that records and reports the fallback status.
