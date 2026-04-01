#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/home/$USER/WORK/Erkang/fast_llqr}"
PARTITION="${PARTITION:-cnall}"
ACCOUNT="${ACCOUNT:-users}"
CHECK_PREFIX="${FASTQR_LLQR_SCAN_PREFIX:-results/table/llqr_rerun}"
SCAN_CORES="${FASTQR_SCAN_CORES:-8}"
export FASTQR_NUM_REP="${FASTQR_NUM_REP:-1000}"

cd "$PROJECT_DIR"
mkdir -p logs results/table

module load compilers/gcc/v12.2.0
export PATH=/apps/soft/R/R-4.3.1/bin:$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64:$LD_LIBRARY_PATH
export R_LIBS_USER="/home/$USER/R/x86_64-pc-linux-gnu-library/4.3"

export FASTQR_LLQR_SCAN_PREFIX="$CHECK_PREFIX"
export FASTQR_SCAN_CORES="$SCAN_CORES"
Rscript scripts/array/scan_llqr_integrity.R

PROBLEM_TSV="${CHECK_PREFIX}_problem_configs.tsv"
if [[ ! -f "$PROBLEM_TSV" ]]; then
  echo "ERROR: problem config file not found: $PROBLEM_TSV" >&2
  exit 1
fi

squeue -h -u "$USER" -o "%i %j %T %R" | awk '$2 ~ /^llqr_/ && $3 == "PENDING" && $4 ~ /DependencyNeverSatisfied/ {print $1}' | xargs -r scancel

tail -n +2 "$PROBLEM_TSV" | while IFS=$'\t' read -r case_id tau n n_missing n_failed missing_rate failed_rate; do
  [[ -z "$case_id" ]] && continue

  export FASTQR_CASE="$case_id"
  export FASTQR_TAU="$tau"
  export FASTQR_N="$n"
  export PARTITION ACCOUNT PROJECT_DIR

  export CPUS_PER_TASK=1
  export FASTQR_CHUNK_SIZE=1
  export FASTQR_MM_FACTOR="1e-3,1e-4"
  export FASTQR_MAX_ATTEMPTS_PER_REP=1
  export FASTQR_RETRY_STRIDE=1000000
  export FASTQR_MAX_SECONDS_PER_REP=7200
  export FASTQR_INCLUDE_LLQR_PPRO=0
  export FASTQR_INCLUDE_LLQR_PPRO_FORTRAN=1
  export FASTQR_LLQR_CHECK_OUT_DIR="results/table"

  if [ "$n" -le 1000 ]; then
    export WALLTIME=04:00:00
  elif [ "$n" -eq 2000 ]; then
    export WALLTIME=06:00:00
  else
    export WALLTIME=08:00:00
  fi

  echo "Resubmitting LLQR case=${case_id} tau=${tau} n=${n} with chunk=1 walltime=${WALLTIME} (missing=${n_missing}, failed=${n_failed})"
  ./scripts/slurm/submit_llqr_rerun_array.sh
  sleep 2
done

echo "Submitted reruns for all problem LLQR configs listed in ${PROBLEM_TSV}."
echo "Monitor with: squeue -u $USER"
