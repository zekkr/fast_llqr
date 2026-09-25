# Release status — 2026-09-25

The paper code repositories are public:

- R package: [zekkr/fastllqr](https://github.com/zekkr/fastllqr), fixed tag `v0.2.0` (commit `7c97c99`).
- Research and reproduction code: [zekkr/fast_llqr](https://github.com/zekkr/fast_llqr), fixed tag `paper-u11-v0.2.0`.

Anonymous HTTP access to both repository pages, the package tag, the reproduction
directory, and the current Supplementary Table S1 archive returned 200 on
2026-09-25. Existing Git histories were preserved. The release-source backups
and pre-clean manuscript copies are stored locally under the ignored
`fast_llqr/tmp/paper_release_20260925/` directory.

The package source and research U11 files passed the synchronization check.
`R CMD check` returned 0 errors and 1 warning for the unchanged non-portable
Fortran `-ffree-line-length-none` flag. The independent validation covered 60
datasets and 8,652 evaluation points: no H-row-set mismatch, no seq fallback,
and 1,481 full-active recoveries. The current LLQR archive reconstruction
covered 12,000 datasets and 36,000 fits with no H mismatch or seq fallback.
The Table S1 archive contains 600 paired datasets and their source snapshots;
a fresh one-pair execution check completed with maximum fitted-value difference
1.56e-13. See `VALIDATION.md`, `README.md`, and
`archive/mst_rep100_20260925/README.md` for methods and limits.

Historical September quantreg/BLAS metadata and an independent Linux check
remain unavailable. Fresh runtime measurements on other systems need not equal
the archived seconds. No new full four-case 500-replication run was performed
for this release.
