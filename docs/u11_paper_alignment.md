# U11 paper-alignment staging note

Date: 2026-09-11

This note records manuscript changes that would be required if the shared U11
LLQR/TVCQR implementation becomes the paper implementation. It is staging
documentation only: neither `main_submission_TW_WC_reviewed.tex` nor
`supplement_TW_WC_reviewed.tex` is changed in this experiment branch.

## Shared weighted-quantile-regression algorithm

U11 uses the same screened weighted quantile-regression core for LLQR and
TVCQR. Model adapters supply the design, evaluation grid, kernel weights, and
threshold. The core performs the full first solve, screening, basis padding,
aggregation, reduced solves, deterministic verification or certification,
repair, and explicitly recorded full-active numerical recovery. LLQR/TVCQR
specifics therefore remain outside the generic Algorithm 1 except for their
design and weight definitions.

Aggregate pseudo-rows are protected and do not participate in the simplex
ratio test. This is an accepted implementation convention for the shared core.
The individual observations represented by those rows remain subject to the
same strict-sign certification rule.

## Definition of the implemented first retained set

The current main manuscript defines the theoretical threshold set as
\[
  S_{j,n}^{\mathrm{thr}}
  =\{i\in N_n(s_{j+1}):|\hat e_{i,n}(s_j)|\le\gamma_n\}.
\]
The implemented first retained set should instead be denoted explicitly as
\[
  S_{j,n}^{\mathrm{impl}}
  =S_{j,n}^{\mathrm{eff}}\cup B_{j,n}^{\mathrm{pad}},
\]
where effective refers to the threshold used by the first reduced-solve
attempt after any pre-solve threshold expansion, and
`B_{j,n}^{pad}` contains previous-basis observations forced into that attempt,
including zero-weight padding. Aggregate pseudo-rows are not members of this
set. The first grid point is excluded from `max_j |S_{j,n}^{impl}|`.

If no threshold expansion occurs, `S_eff=S_thr`; moreover `|B_pad|<=q`. Thus
the existing order bound transfers immediately after adding at most `q` rows.
If expansion occurs, the theorem must either control the effective threshold or
condition the displayed first-pass result on no pre-solve expansion. This is an
open theoretical choice and must not be hidden by treating tableau rows as
`|S_{j,n}|`.

The staged screening CSV therefore keeps the current-paper ratio using the
nominal threshold and also reports a separately labelled conservative
sensitivity ratio using the largest first-attempt effective threshold in each
replication. These columns must not be conflated when drafting the theorem or
table note.

The main-paper paragraph currently saying that the program records the first
tableau size should be replaced: the U11 experiment records the actual number
of retained individual observations and separately records basis padding and
zero-to-two aggregate rows. The simulation table should use
`first_retained_size`, not legacy `first_n_sub`/tableau size.

## Exact residual-image certificate (`cache_flags=27`)

Algorithm 1 currently says to recompute the omitted residuals and verify every
strict sign after each reduced solve. U11 may skip that literal scan only under
a deterministic exact-image certificate:

1. the reduced candidate coefficient vector is bitwise identical to the
   preceding accepted coefficient vector whose raw residual image is cached;
2. the current interpolation H row-set equals the preceding H row-set; and
3. the current threshold is finite and nonnegative.

The current screening tags were formed from that cached raw-residual image.
Under these conditions each omitted residual has the same floating-point image
at the candidate, so a tag strictly beyond the nonnegative threshold certifies
its strict sign. This is not probabilistic skipping and does not use an
approximate coefficient-displacement bound. Cache mode 31 audits certificate
hits by recomputing all current active residuals, requiring bitwise agreement
with the cache, and checking the omitted strict signs; mode 25 disables the
certificate and always uses literal verification.

Recommended Algorithm 1 wording: require a deterministic verification
certificate at every acceptance, with literal omitted-row verification as the
default certificate. Put the bitwise exact-image condition and cache audit in
the supplement/implementation section. If the paper continues to say that all
omitted residuals are recomputed on every attempt, U11 mode 27 is not literally
aligned with that statement.

## Numerical solve failure versus mathematical infeasibility

For the LLQR and TVCQR formulations in the paper, the transported preceding
fit supplies a feasible primal point for the first reduced LP. Consequently, a
floating-point solver failure should not be described as mathematical
infeasibility without a separate feasibility certificate.

U11 treats a numerical reduced-solve failure as a finite-precision safeguard:
it may enlarge the threshold and rebuild the reduced problem; after repeated
failure it performs an explicitly recorded full-active solve. These operations
do not change the statistical objective. Algorithm 1 should distinguish:

- certified mathematical infeasibility, which permits immediate full-LP solve;
- numerical nonconvergence, singular initialization, or invalid tableau state,
  which triggers a documented finite-precision retry/recovery policy.

The supplement's generic infeasibility propositions may remain mathematical,
but its implementation prose should state that U11's observed failure flag is
not itself an infeasibility proof. Numerical retry/full-active recovery can be
given in a separate implementation-safeguards paragraph rather than embedded
in the statistical algorithm.

## Strict signs, numerical discrepancy, and timing

- Verification is strict: a sure-positive row with recomputed residual `<=0`
  and a sure-negative row with residual `>=0` fails certification. Exact zero
  is therefore a failure, consistent with the reviewed supplement's current
  implementation paragraph.
- The numerical discrepancy used by this experiment is
  `mean(abs(candidate-direct)/pmax(abs(direct),1e-10))`. Any supplement formula
  or data script must use the same denominator and must not drop near-zero
  direct entries.
- The U11 rep500 timing includes bandwidth selection, initialization, complete
  grid fitting, and returned-object construction. It excludes DGP, explicit
  garbage collection, scheduling, and serialization. This timing boundary must
  replace any older prose that times only an already-initialized kernel call.
- The paper's global Cases 1--2 are LLQR model cases 1--2; global Cases 3--4
  are TVCQR model cases 1--2. Staged output must retain this mapping.

## Manuscript locations to revise after staging approval

Main manuscript:

- Algorithm 1: retained-set union, deterministic certificate interface, and
  separation of mathematical infeasibility from numerical recovery.
- Paragraph immediately after Algorithm 1: replace the `|S|` versus tableau
  approximation with the exact retained-size diagnostic.
- Theorem 2 discussion and Corollary 1: decide how the at-most-`q` basis union
  and any effective-threshold expansion enter the row-count and work bounds.
- Simulation section and retained-size tables: use U11-generated
  `first_retained_size`, and update discrepancy/runtime statements only as one
  complete rep500 data replacement.

Supplement:

- Algorithm/tableau construction: state that aggregate rows are excluded from
  the ratio test and that basis padding is included among retained individual
  rows.
- First-pass theory: distinguish the theoretical threshold set from the
  implemented basis-augmented set, or restate the bound for their union.
- Exactness proof: allow a deterministic certificate whose validity implies
  the literal strict-sign check; specify exact-image certification separately.
- Implementation safeguards: document threshold retry and full-active recovery
  as finite-precision behavior, not evidence of mathematical infeasibility.
- Simulation appendix: replace old data only after all 48 configurations pass
  the declared gate.

## Open questions

1. Restate the dimension-reduction theorem directly for
   `S_thr union B_pad`, or retain the threshold-set theorem and add a corollary
   using `|B_pad|<=q`?
2. If a first attempt uses an expanded threshold, should theory analyze that
   effective random threshold or should first-pass claims be conditional on no
   pre-solve expansion?
3. Should Algorithm 1 expose a generic deterministic-certificate interface, or
   state the exact bitwise residual-image certificate in the main algorithm?
4. Should numerical retry/full-active recovery appear in Algorithm 1 or only
   as a finite-precision implementation safeguard in the supplement?
5. Should the lean-seq ablation report only time ratio and direct-relative
   discrepancy, or also compare basis reuse and recovery mechanisms?
6. After a successful rep500 gate, should U11 replace the current production
   LLQR/TVCQR backends?
