#!/usr/bin/env bash
set -euo pipefail

: "${SSQR_RUN_DIR:?}" "${SSQR_LOG_DIR:?}"
sacct_bin="${SSQR_SACCT:-/rmprog/slurm/v22.05.7/bin/sacct}"
meta="$SSQR_RUN_DIR/_run_meta"
mkdir -p "$meta/slurm_logs"
cp -a "$SSQR_LOG_DIR"/. "$meta/slurm_logs/"

ids="$(
  {
    awk -F'\t' 'NR>1 {print $6; print $7}' "$meta/submission_manifest.tsv"
    awk -F'\t' 'NR>1 && NF==2 {print $2}' "$meta/job_ids.tsv"
  } | awk 'NF' | sort -u | paste -sd, -
)"
if [[ -n "$ids" ]]; then
  "$sacct_bin" -j "$ids" --format=JobIDRaw,JobName%36,Partition,State,ExitCode,Elapsed,AllocCPUS,NodeList%30 \
    -P > "$meta/slurm_sacct.txt" || true
fi

cd "$SSQR_RUN_DIR"
find . -type f ! -name RESULTS_SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
  > RESULTS_SHA256SUMS
printf 'finalized_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$meta/finalized.txt"
# Refresh once so finalized.txt is included in the immutable download manifest.
find . -type f ! -name RESULTS_SHA256SUMS -print0 | sort -z | xargs -0 sha256sum \
  > RESULTS_SHA256SUMS
printf 'finalized %s files\n' "$(wc -l < RESULTS_SHA256SUMS)"
