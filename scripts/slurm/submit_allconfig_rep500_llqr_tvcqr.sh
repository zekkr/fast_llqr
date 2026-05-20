#!/bin/bash
set -euo pipefail

# Submit LLQR / TVCQR full-config rep=500 jobs using the existing per-model
# SLURM wrappers in this repo.
#
# Usage:
#   bash scripts/slurm/submit_allconfig_rep500_llqr_tvcqr.sh
#
# Optional overrides:
#   PROJECT_DIR=/path/to/fast_llqr
#   RUN_LLQR=1 RUN_TVCQR=1
#   LLQR_CASES="1 2"
#   TVCQR_CASES="1 2"
#   TAUS="0.2 0.5 0.8"
#   NS="200 500 1000 2000 5000"
#   FASTQR_NUM_REP=500
#   CPUS_PER_TASK=56
#   PARTITION=cnall
#   ACCOUNT=users
#   WALLTIME_LLQR=08:00:00
#   WALLTIME_TVCQR=08:00:00
#   FASTQR_MM_FACTOR_LLQR="1e-2,1e-3,1e-4"
#   FASTQR_MM_FACTOR_TVCQR="1e-3,1e-4"

PROJECT_DIR="${PROJECT_DIR:-/home/$USER/WORK/Erkang/fast_llqr}"
RUN_LLQR="${RUN_LLQR:-1}"
RUN_TVCQR="${RUN_TVCQR:-1}"

LLQR_CASES="${LLQR_CASES:-1 2}"
TVCQR_CASES="${TVCQR_CASES:-1 2}"
TAUS="${TAUS:-0.2 0.5 0.8}"
NS="${NS:-200 500 1000 2000 5000}"

FASTQR_NUM_REP="${FASTQR_NUM_REP:-500}"
CPUS_PER_TASK="${CPUS_PER_TASK:-56}"
PARTITION="${PARTITION:-cnall}"
ACCOUNT="${ACCOUNT:-users}"

WALLTIME_LLQR="${WALLTIME_LLQR:-08:00:00}"
WALLTIME_TVCQR="${WALLTIME_TVCQR:-08:00:00}"

FASTQR_MM_FACTOR_LLQR="${FASTQR_MM_FACTOR_LLQR:-1e-2,1e-3,1e-4}"
FASTQR_MM_FACTOR_TVCQR="${FASTQR_MM_FACTOR_TVCQR:-1e-3,1e-4}"

FASTQR_MERGE_DEPENDENCY_TYPE="${FASTQR_MERGE_DEPENDENCY_TYPE:-afterany}"

submit_one() {
  local model="$1"
  local case_id="$2"
  local tau="$3"
  local n="$4"

  export PROJECT_DIR
  export FASTQR_CASE="$case_id"
  export FASTQR_TAU="$tau"
  export FASTQR_N="$n"
  export FASTQR_NUM_REP
  export CPUS_PER_TASK
  export PARTITION
  export ACCOUNT
  export FASTQR_MERGE_DEPENDENCY_TYPE

  if [[ "$model" == "llqr" ]]; then
    export FASTQR_MM_FACTOR="$FASTQR_MM_FACTOR_LLQR"
    export WALLTIME="$WALLTIME_LLQR"
    echo "[LLQR] case=${case_id} tau=${tau} n=${n} rep=${FASTQR_NUM_REP}"
    bash "${PROJECT_DIR}/scripts/slurm/submit_llqr_array.sh"
  elif [[ "$model" == "tvcqr" ]]; then
    export FASTQR_MM_FACTOR="$FASTQR_MM_FACTOR_TVCQR"
    export WALLTIME="$WALLTIME_TVCQR"
    echo "[TVCQR] case=${case_id} tau=${tau} n=${n} rep=${FASTQR_NUM_REP}"
    bash "${PROJECT_DIR}/scripts/slurm/submit_tvcqr_array.sh"
  else
    echo "Unknown model: ${model}" >&2
    exit 1
  fi
}

cd "$PROJECT_DIR"

echo "=== submit_allconfig_rep500_llqr_tvcqr ==="
echo "PROJECT_DIR=${PROJECT_DIR}"
echo "LLQR_CASES=${LLQR_CASES}"
echo "TVCQR_CASES=${TVCQR_CASES}"
echo "TAUS=${TAUS}"
echo "NS=${NS}"
echo "FASTQR_NUM_REP=${FASTQR_NUM_REP}"
echo "CPUS_PER_TASK=${CPUS_PER_TASK}"
echo "PARTITION=${PARTITION}"
echo "ACCOUNT=${ACCOUNT}"
echo "RUN_LLQR=${RUN_LLQR} RUN_TVCQR=${RUN_TVCQR}"
echo

if [[ "$RUN_LLQR" == "1" ]]; then
  for case_id in $LLQR_CASES; do
    for tau in $TAUS; do
      for n in $NS; do
        submit_one llqr "$case_id" "$tau" "$n"
      done
    done
  done
fi

if [[ "$RUN_TVCQR" == "1" ]]; then
  for case_id in $TVCQR_CASES; do
    for tau in $TAUS; do
      for n in $NS; do
        submit_one tvcqr "$case_id" "$tau" "$n"
      done
    done
  done
fi

echo
echo "All submissions sent."
