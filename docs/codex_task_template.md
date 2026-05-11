# Codex Task Templates

## Scope

This document defines reusable Codex task prompts.

It does not restate algorithm or testing rules. Instead, every task template should refer to:

- `docs/algorithm_contract.md` for algorithm invariants;
- `docs/testing_contract.md` for validation, fallback accounting, and command templates;
- `docs/debugging_playbook.md` for debugging workflow.

Use these templates to keep future tasks scoped and auditable. Fill in the bracketed fields before starting work.

## Minimal Prompt Skeleton

Use this when the repository contracts already contain the stable constraints.

```text
Task type: [plan only / diagnosis only / implementation after approval / test only].

Before working, read and follow:
- AGENTS.md
- docs/codex_task_template.md
- docs/algorithm_contract.md
- docs/testing_contract.md
- docs/debugging_playbook.md

Use the corresponding task template from docs/codex_task_template.md.

Objective:
[task-specific objective]

Target files to inspect:
- [file 1]
- [file 2]

Files allowed to edit:
- [none for plan-only/diagnosis-only]
- [exact files for implementation tasks]

Files forbidden to edit:
- [explicit exclusions]

Task-specific constraints:
- [constraint 1]
- [constraint 2]

Required deliverable:
[plan-only report / diagnosis report / implementation report / test report]
```

## Common Preamble for All Tasks

Every task should start from this shared preamble:

```text
Before working, read and follow:
- AGENTS.md
- docs/algorithm_contract.md
- docs/testing_contract.md
- docs/debugging_playbook.md when debugging H_seq mismatch, fallback, stale aggregates, or certification failure
- docs/codex_task_template.md

Do not restate repository contracts in the task prompt. Refer to the source-of-truth documents instead.
```

## Plan-Only Code Modification Task

```text
Task type: plan only. Do not modify, create, delete, rename, format, or patch any file. Do not run implementation. Do not produce a diff. Your only deliverable is a concrete, auditable modification plan.

### Objective

I want a code modification plan for:

[Describe the intended change, bug fix, refactor, performance improvement, diagnostic improvement, or simulation/reporting update.]

The goal is not to implement it now. The goal is to decide exactly what should be changed, where, why, and how it will be validated.

### Repository context

Before planning, inspect the relevant repository instructions and contracts:

- AGENTS.md
- R/AGENTS.md if touching R code
- src/fortran/AGENTS.md if touching Fortran code
- docs/algorithm_contract.md
- docs/testing_contract.md
- docs/debugging_playbook.md
- docs/codex_task_template.md

Use the repository’s existing naming conventions, diagnostics, and algorithm contracts. Do not invent a parallel framework unless the existing structure is clearly inadequate.

### Target files to inspect

Inspect only the files needed to produce the plan. Likely relevant files include:

- [file 1]
- [file 2]
- [file 3]

Optional reference files:

- [reference file 1]
- [reference file 2]

If you need to inspect additional files, list them in the final report and explain why they were necessary.

### Edit rule

Files allowed to edit: none.

Plan-only tasks may inspect relevant files but must not modify, create, delete, rename, format, or patch any file.

### Contract references

Do not restate algorithm invariants in the task prompt.

The plan must explicitly say how it preserves the relevant rules in:

- `docs/algorithm_contract.md`;
- `docs/testing_contract.md`.

If the task touches H_seq mismatch, fallback, stale aggregates, or verification failures, also use:

- `docs/debugging_playbook.md`.

### Required plan content

The final plan must include:

1. Problem interpretation:
   - requested change;
   - target behavior;
   - out of scope.

2. Current code path:
   - entry points;
   - affected functions/scripts/backends;
   - relevant existing diagnostics.

3. Proposed edit scope:
   - exact future files/functions to edit;
   - required vs optional vs risky edits.

4. Implementation plan:
   - ordered steps;
   - expected behavior change for each step;
   - no patch and no full replacement function.

5. Validation and risks:
   - smallest smoke test;
   - smallest reproducer if applicable;
   - H_seq/fallback checks when relevant;
   - risks and blocking questions.

### Required final report format

Use exactly this structure:

Plan-only report

1. Problem interpretation
- Requested change:
- Target behavior:
- Out of scope:

2. Files inspected
- ...

3. Current code path
- Entry points:
- Core functions:
- Backend/interface path:
- Existing diagnostics/tests:

4. Proposed edit scope
- Required edits:
- Optional edits:
- Avoided edits:

5. Step-by-step implementation plan
Step 1:
Step 2:
Step 3:
...

6. Validation plan
- Smoke test:
- Reproducer:
- H_seq/set comparison:
- Fallback/status check:
- Regression grid:
- Commands to run:

7. Risks and blocking questions
- Risks:
- Blocking questions, if any:

8. Minimal follow-up implementation prompt
[Copyable implementation task prompt]

### Hard constraints

Do not edit files.
Do not produce a patch.
Do not run long simulations.
Do not hide uncertainty.
Do not recommend broad rewrites unless the current structure makes the requested change unsafe.
Prefer the smallest correct change over a large refactor.
If the task touches ppro correctness, correctness dominates speed.
If the task touches timing summaries, separate true ppro success, fallback success, and failure.
```

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

Contract requirements:
- Follow `docs/algorithm_contract.md`.
- Follow `docs/testing_contract.md`.
- If debugging H_seq mismatch, fallback, or certification failure, follow `docs/debugging_playbook.md`.
- In the final report, state which contract rules were most relevant.

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

Purpose:
Implement the approved change only. Do not broaden the scope.

Target files:
- [exact files allowed to edit]

Files allowed to edit:
- [exact files allowed to edit]

Files forbidden to edit:
- [exact files excluded from this task]
- paper/main_v3.tex unless explicitly requested

Contract checks:
- follow `docs/algorithm_contract.md`;
- follow `docs/testing_contract.md`;
- use `docs/debugging_playbook.md` if the change touches H_seq mismatch, fallback, stale aggregates, or certification failure;
- before editing, state which contract rule is preserved, changed, or relied on.

Required actions:
- run the smallest reproducer before and after the change when applicable;
- run H_seq row-set comparison when ppro/H_seq behavior is touched;
- check fallback accounting when wrapper, certificate, timing, or method-factory behavior is touched;
- do not run large simulations unless explicitly approved.

Final report fields:
- files changed;
- behavior changed;
- commands run;
- exact test outcomes;
- fallback status or not applicable;
- H_seq row-set status or not applicable;
- remaining risks.
```

## Test-Only Task

```text
Task type: tests only. Do not change algorithm implementations.

Purpose:
Add, update, or run tests/diagnostics only.

Target files:
- [test or diagnostic files]

Files allowed to edit:
- [exact test/diagnostic files]

Files forbidden to edit:
- R/llqr_functions.R
- R/tvcqr_functions.R
- src/fortran/*.f90 unless explicitly approved
- paper/main_v3.tex unless explicitly requested

Contract checks:
- follow `docs/testing_contract.md`;
- do not redefine correctness levels, fallback accounting, or regression gates in this task.

Required actions:
- use existing scripts when possible;
- document the grid actually covered;
- separate true ppro success, fallback success, failure, and malformed output;
- do not change tests merely to make a broken implementation pass.

Final report fields:
- files changed;
- commands run;
- grid covered;
- failures found;
- fallback handling;
- missing scripts or data;
- remaining test coverage gaps.
```

## Fortran-Specific Bugfix Task

```text
Task type: Fortran bugfix.

Purpose:
Fix a Fortran-backed solver issue without changing unrelated R algorithms.

Target files:
- src/fortran/[file].f90
- R/[wrapper file].R only if interface handling must be updated

Files allowed to edit:
- [exact files]

Files forbidden to edit:
- unrelated R algorithms
- tests and simulation scripts unless explicitly included
- paper/main_v3.tex unless explicitly requested

Contract checks:
- follow `docs/algorithm_contract.md`;
- follow `docs/testing_contract.md`;
- report any interface change explicitly.

Required actions:
- compare Fortran output against the relevant R/reference path on a small reproducer;
- compare Fortran ppro H_seq against stable seq when ppro/H_seq behavior is touched;
- report kernel status and fallback separately.

Final report fields:
- Fortran routines changed;
- interface changes, if any;
- reproducer;
- commands run;
- exact test outcomes;
- kernel status behavior;
- fallback status;
- remaining numerical risks.
```

## R Wrapper/Certificate Task

```text
Task type: R wrapper/certificate.

Purpose:
Change wrapper, certificate, method-factory, or returned-status behavior without changing the underlying solver unless explicitly approved.

Target files:
- R/[wrapper/helper file].R

Files allowed to edit:
- [exact files]

Files forbidden to edit:
- src/fortran/*.f90 unless explicitly approved
- paper/main_v3.tex unless explicitly requested
- unrelated simulation scripts

Contract checks:
- follow `docs/testing_contract.md` for returned status fields and fallback accounting;
- follow `docs/algorithm_contract.md` if wrapper behavior affects solver certification.

Required actions:
- test one successful ppro path when applicable;
- test or simulate one fallback path when applicable;
- confirm returned_backend, fallback_triggered, and fallback_reason behavior when touched;
- run H_seq row-set comparison when ppro/H_seq behavior is touched.

Final report fields:
- wrapper fields changed;
- method-factory behavior changed;
- certification/fallback behavior changed;
- commands run;
- exact test outcomes;
- compatibility risks.
```

## Documentation-Only Task

```text
Task type: documentation only.

Target documents:
- [doc 1]
- [doc 2]

Files allowed to edit:
- [exact documentation files]

Files forbidden to edit:
- R/*.R
- src/fortran/*.f90
- scripts/*.R unless explicitly requested
- paper/main_v3.tex unless explicitly requested

Contract checks:
- preserve the source-of-truth map in AGENTS.md;
- do not duplicate long contract sections;
- do not introduce algorithm or testing rules in the wrong document.

Required checks:
- confirm which document is the source of truth for each moved rule;
- check for duplicated or conflicting rules after edits;
- report any remaining ambiguity.

Final report fields:
- files changed;
- sections deleted;
- sections moved;
- sections added;
- source-of-truth mapping;
- remaining documentation risks.
