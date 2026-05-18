# Debugging Playbook

Use this workflow when a screened ppro method has an `H_seq` mismatch against the corresponding stable seq method.

## Scope

Use this playbook when a screened ppro method has one of the following symptoms:

- `H_seq` mismatch against the corresponding stable seq method;
- estimate mismatch after ppro screening;
- unexpected fallback;
- missing or malformed fallback status;
- suspected stale aggregate, stale residual, or invalid warm start.

This document is a debugging workflow, not an implementation contract. Algorithm invariants belong in `docs/algorithm_contract.md`; diagnostic output fields belong in `docs/testing_contract.md`.

## Workflow

1. Reproduce the smallest failing case: model, case, `n`, `tau`, seed, method name, `Mm.factor`, evaluation index, and backend.
2. Run the stable seq solver and the ppro solver with fallback disabled if possible, so kernel errors are not hidden.
3. Compare `H_seq` row sets at each evaluation index and identify the first mismatch.
4. Inspect the previous evaluation point, not only the failing one. Most ppro failures come from stale residuals, stale aggregates, or an invalid warm start carried forward.
5. Recompute `sl`, `sh`, and `S` independently for the failing transition and confirm previous H/basis observations are forced into `S`.
6. Rebuild aggregate rows from scratch at the current evaluation point and compare them to the solver's rows.
7. Solve the reduced LP and verify all omitted observations. Record bad signs by original observation index and assigned sign class.
8. If bad signs exist, confirm they are re-added to `S` or that threshold doubling recomputes every dependent quantity before re-solving.
9. Confirm the reduced simplex state converged before accepting it.
10. If the wrapper returns a seq fallback, inspect the original kernel status and fallback reason before treating the output as useful.

## Checks

Stale aggregate rows: Aggregates must use current weights, current `sl/sh/S`, current `tau`, and current local design. Any threshold or set change invalidates prior aggregate rows.

Incorrect `sl/sh/S` updates: Classification must be rebuilt from the intended residual source and threshold. `S` must be the complement of certified omitted signs plus required basis observations.

Previous H not forced into `S`: Previous H/basis observations must stay explicit in the next reduced LP. Missing this step can produce plausible estimates with different basis rows.

Bad signs not re-added correctly: Verification failures must map back to original observation indices. Adding reduced-row positions instead of original indices is a common source of repeated failure.

Threshold doubling not recomputing all dependent quantities: After doubling or otherwise changing the threshold, recompute `sl`, `sh`, `S`, counts, aggregates, reduced design rows, objective coefficients, and warm-start metadata.

Accepting nonconverged simplex states: A reduced solve that reaches `maxit`, has invalid pivots, or returns an error flag must not be accepted as certified output.

Zero-weight observations: Local kernels can produce zero weights. Confirm zero-weight rows do not enter aggregates or basis comparisons in a way that changes the effective LP.

Bounded-kernel positive-support degeneracy: For LLQR case 2, record `sum(w > 0)`, `which(w > 0)`, `qr(A[w > 0, , drop = FALSE])$rank`, `h`, `case`, `seed`, and the traceback line at the first failing `rd`. If the positive-weight count is below `q` or the positive-weight design rank is below `q`, diagnose a full-sample local LP degeneracy before attributing the failure to `Mm.factor`, screening, or `H_seq` mismatch.

Tolerance and tie-breaking differences: Compare tolerances, Bland-rule settings, pivot tie handling, and near-zero residual classification. A row near the threshold should be treated as uncertain rather than certified.

Fortran wrapper fallback hiding a kernel error: Wrappers must expose kernel status and fallback status. If `returned_backend = "seq_fallback"`, debug the kernel error first and do not count the run as ppro success.

## Diagnostic output

For required diagnostic fields, see `docs/testing_contract.md#diagnostic-output-schema`.

## Stop conditions

Stop debugging and report a diagnosis when any of the following is identified:

- the smallest failing configuration is found;
- the first failing evaluation index is found;
- fallback hides the original kernel error;
- a stale aggregate, stale residual, stale threshold, or stale `S/sl/sh` update is confirmed;
- previous H/basis observations are not forced into `S`;
- a reduced LP is accepted despite nonconvergence;
- the failure cannot be reproduced with the available scripts or data.

Do not implement a fix from this playbook unless the current task explicitly allows implementation.
