# AGENTS.md

## Project goal

This repository supports simulations for the paper
“A fast simplex algorithm for local linear quantile regression.”

The main research objects are:
- `llqr_seq_ppro`
- `tvcqr_seq_ppro`
- their Fortran-backed implementations

The `quantreg`/`quantdr`-based solvers (`llqr`, `tvc_rq`) are baseline/oracle solvers. Our sequential solvers (`llqr_seq`, `tvcqr_seq`) and their Fortran-backed implementations are now stable and fast. Their `H_seq` outputs are treated as the current regression-test standard for whether the `seq_ppro` solvers have reached the expected solution path.

Do not replace a ppro implementation by a seq implementation.

## Contract documents

Do not restate algorithm, testing, or debugging rules in this file.

The source-of-truth documents are:

- `docs/algorithm_contract.md` for screened ppro algorithm invariants;
- `docs/testing_contract.md` for correctness levels, regression gates, fallback accounting, and command templates;
- `docs/debugging_playbook.md` for H_seq mismatch and fallback debugging workflow;
- `docs/codex_task_template.md` for task-specific prompt templates.

## Task-mode discipline

At the start of every task, identify the task mode:

- plan only;
- diagnosis only;
- implementation after approval;
- test only;
- Fortran-specific bugfix;
- R wrapper/certificate task;
- paper/documentation task.

If the task says **plan only**, do not modify, create, delete, rename, format, or patch any file. Produce only the requested plan.

If the task says **diagnosis only**, do not implement a fix unless a follow-up implementation task explicitly approves it.

If the task says **implementation after approval**, edit only the files explicitly allowed by the prompt.

If the task scope is ambiguous, prefer a narrow plan or diagnosis over broad edits.

## Do-not rules

- Do not modify `llqr_seq` or `tvcqr_seq` unless explicitly asked.
- Do not change simulation data-generating mechanisms unless explicitly asked.
- Do not change tests to make a broken implementation pass.
- Do not remove diagnostics that distinguish ppro from seq fallback.
- Do not edit `paper/main_v3.tex` unless the task is explicitly about the paper.

## File inspection rule

Inspect only files relevant to the current task.

For ppro algorithm changes, the usual core files are:

- `R/llqr_functions.R`
- `R/tvcqr_functions.R`
- `src/fortran/llqr_ppro.f90`
- `src/fortran/tvcqr_seq_M_acc.f90`

For H_seq regression or fallback accounting, the usual diagnostic file is:

- `scripts/check_hseq_set_match.R`

For documentation-only, plan-only, or paper-only tasks, do not inspect implementation files unless needed for the requested plan.

## Code-change reporting floor

For any task that proposes or implements a code change, the response must state:

Before implementation:
- which contract rule is being preserved, changed, or relied on;
- which files and functions are expected to be affected;
- how the change will be validated.

After implementation:
- exact commands run, or `not run` with the reason;
- exact test outcomes, not just "passed";
- whether fallback occurred, or `not applicable`;
- whether `H_seq` row-set match passed when ppro/H_seq behavior is touched, or `not applicable`.

Detailed algorithm rules belong in `docs/algorithm_contract.md`.
Detailed validation rules belong in `docs/testing_contract.md`.
Task-specific report formats belong in `docs/codex_task_template.md`.

## Source-of-truth map

Use this map to avoid duplicated or conflicting instructions:

| Topic | Source of truth |
|---|---|
| ppro must not call seq except explicit fallback | `docs/algorithm_contract.md` |
| `sl`, `sh`, `S`, bad signs, H/basis observations | `docs/algorithm_contract.md` |
| aggregate recomputation and verification failure handling | `docs/algorithm_contract.md` |
| fallback fields and fallback accounting | `docs/testing_contract.md` |
| correctness levels and active regression gate | `docs/testing_contract.md` |
| H_seq mismatch debugging workflow | `docs/debugging_playbook.md` |
| command templates and regression grids | `docs/testing_contract.md` |
| task prompt formats | `docs/codex_task_template.md` |