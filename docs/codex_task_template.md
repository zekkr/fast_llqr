# Codex Task Templates

Use these templates to keep future tasks scoped and auditable. Fill in the bracketed fields before starting work.

## Diagnosis-Only Task

```text
Task type: diagnosis only. Do not modify code unless I approve a follow-up implementation task.

Target files to inspect:
- [files]

Files not to edit:
- R/llqr_functions.R
- R/tvcqr_functions.R
- src/fortran/llqr_ppro.f90
- src/fortran/tvcqr_seq_M_acc.f90
- scripts/check_hseq_set_match.R
- paper/main_v3.tex
- tests and simulation scripts

Algorithm invariants:
- ppro must implement screened warm-start logic, not call seq as its implementation.
- previous H/basis observations must be forced into S.
- verification failures must add bad signs or enlarge the threshold and re-solve.
- fallback to seq must be explicit and reported.

Required tests:
- identify the smallest failing model/case/n/tau/seed/evaluation index;
- compare ppro H_seq to stable seq H_seq as row sets;
- report whether any returned result came from fallback.

Final report fields:
- files inspected;
- first failing configuration;
- suspected cause;
- evidence;
- whether fallback was involved;
- recommended next edit scope.
```

## Implementation Task After Approval

```text
Task type: implementation after approval.

Target files:
- [exact files allowed to edit]

Files not to edit:
- [exact files excluded from this task]
- paper/main_v3.tex unless explicitly requested

Algorithm invariants:
- do not replace ppro with seq;
- do not silently accept nonconverged reduced LPs;
- recompute aggregates after any sl/sh/S or threshold change;
- force previous H/basis observations into S;
- report fallback with fallback_triggered, fallback_reason, and returned_backend.

Required tests:
- run the smallest reproducer before and after the change;
- run H_seq row-set comparison against stable seq for the touched model;
- verify fallback counts separately from ppro successes;
- do not run large simulations unless explicitly approved.

Final report fields:
- files changed;
- algorithm behavior changed;
- tests run;
- first failing case before fix;
- result after fix;
- fallback count/status;
- remaining risks.
```

## Test-Only Task

```text
Task type: tests only. Do not change algorithm implementations.

Target files:
- [test or diagnostic files]

Files not to edit:
- R/llqr_functions.R
- R/tvcqr_functions.R
- src/fortran/*.f90 unless explicitly approved
- paper/main_v3.tex

Algorithm invariants to assert:
- ppro H_seq rows match stable seq H_seq rows as sets;
- fallback is not counted as ppro success unless explicitly intended;
- missing H_seq, shape mismatch, and row-set mismatch are separate failure modes.

Required tests:
- use existing scripts when possible, such as scripts/check_hseq_set_match.R;
- include LLQR/TVCQR cases 1/2, n in c(200, 500, 1000), tau in c(0.2, 0.5, 0.8), seeds 2026 through 2525 when compute budget allows;
- document any narrower grid.

Final report fields:
- files changed;
- commands run;
- grid covered;
- failures found;
- fallback handling;
- missing scripts or data.
```

## Fortran-Specific Bugfix Task

```text
Task type: Fortran bugfix.

Target files:
- src/fortran/[file].f90
- R/[wrapper file].R only if interface handling must be updated

Files not to edit:
- unrelated R algorithms;
- tests and simulation scripts unless this task explicitly includes test updates;
- paper/main_v3.tex.

Algorithm invariants:
- no silent seq replacement;
- return explicit status/error flags for kernel failures;
- do not accept nonconvergence;
- do not reuse stale aggregate rows;
- recompute counts after changing sl/sh/S or threshold;
- preserve R/Fortran interface signatures unless this task explicitly asks for an interface change.

Required tests:
- compare Fortran output against R reference for a small reproducer;
- compare Fortran ppro H_seq against stable Fortran seq H_seq as row sets;
- report fallback separately.

Final report fields:
- Fortran routines changed;
- interface changes, if any;
- reproducer;
- tests run;
- kernel status behavior;
- fallback status;
- remaining numerical risks.
```

## R Wrapper/Certificate Task

```text
Task type: R wrapper/certificate.

Target files:
- R/[wrapper/helper file].R

Files not to edit:
- src/fortran/*.f90 unless explicitly approved;
- paper/main_v3.tex;
- unrelated simulation scripts.

Algorithm invariants:
- wrappers must expose fallback status;
- wrappers must not hide kernel failures in timing mode;
- certification and fallback must be separate options;
- method factories must not silently swap ppro for seq;
- diagnostics must preserve whether the returned result came from ppro or fallback.

Required tests:
- force a small successful ppro path and confirm returned_backend identifies ppro;
- force or simulate a fallback path and confirm fallback_triggered and returned_backend are present;
- run H_seq row-set comparison against stable seq for at least the smallest relevant grid.

Final report fields:
- wrapper fields added or changed;
- method factory behavior;
- fallback/certification behavior;
- commands run;
- outputs inspected;
- compatibility risks.
```
