#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/home/$USER/WORK/Erkang/fast_llqr}"

is_pos_int() {
  [[ "$1" =~ ^[0-9]+$ ]] && [ "$1" -gt 0 ]
}

export FASTQR_CASE="${FASTQR_CASE:?FASTQR_CASE is required}"
export FASTQR_TAU="${FASTQR_TAU:?FASTQR_TAU is required}"
export FASTQR_N="${FASTQR_N:?FASTQR_N is required}"
export FASTQR_NUM_REP="${FASTQR_NUM_REP:-1000}"

export FASTQR_MM_FACTOR="${FASTQR_MM_FACTOR:-1e-3,1e-4}"
export FASTQR_H_FACTOR="${FASTQR_H_FACTOR:-1}"
export FASTQR_TOL="${FASTQR_TOL:-1e-14}"
export FASTQR_MAXIT="${FASTQR_MAXIT:-2000000}"
export FASTQR_BLAND="${FASTQR_BLAND:-0}"
export FASTQR_TRACK_ORDER="${FASTQR_TRACK_ORDER:-1}"
export FASTQR_SEED_BASE="${FASTQR_SEED_BASE:-2026}"
export FASTQR_MAX_ATTEMPTS_PER_REP="${FASTQR_MAX_ATTEMPTS_PER_REP:-1}"
export FASTQR_RETRY_STRIDE="${FASTQR_RETRY_STRIDE:-1000000}"
export FASTQR_MAX_SECONDS_PER_REP="${FASTQR_MAX_SECONDS_PER_REP:-7200}"
export FASTQR_INCLUDE_LLQR_PPRO="${FASTQR_INCLUDE_LLQR_PPRO:-0}"
export FASTQR_INCLUDE_LLQR_PPRO_FORTRAN="${FASTQR_INCLUDE_LLQR_PPRO_FORTRAN:-1}"
export FASTQR_REP_ID_FILE="${FASTQR_REP_ID_FILE:-}"
export FASTQR_REP_ID_LIST="${FASTQR_REP_ID_LIST:-}"

CPUS_PER_TASK="${CPUS_PER_TASK:-1}"
WALLTIME="${WALLTIME:-04:00:00}"
PARTITION="${PARTITION:-cnall}"
ACCOUNT="${ACCOUNT:-users}"
CHECK_OUT_DIR="${FASTQR_LLQR_CHECK_OUT_DIR:-results/table}"

export FASTQR_CHUNK_SIZE="${FASTQR_CHUNK_SIZE:-1}"

if [[ "$FASTQR_CASE" != "1" && "$FASTQR_CASE" != "2" ]]; then
  echo "ERROR: FASTQR_CASE must be 1 or 2, got '$FASTQR_CASE'" >&2
  exit 1
fi

for v in FASTQR_N FASTQR_NUM_REP CPUS_PER_TASK FASTQR_CHUNK_SIZE FASTQR_MAX_ATTEMPTS_PER_REP FASTQR_RETRY_STRIDE FASTQR_MAX_SECONDS_PER_REP; do
  val="${!v}"
  if ! is_pos_int "$val"; then
    echo "ERROR: $v must be a positive integer, got '$val'" >&2
    exit 1
  fi
done

if [[ ! -d "$PROJECT_DIR" ]]; then
  echo "ERROR: PROJECT_DIR does not exist: $PROJECT_DIR" >&2
  exit 1
fi

count_sparse_rep_ids() {
  if [[ -n "$FASTQR_REP_ID_FILE" ]]; then
    if [[ ! -f "$FASTQR_REP_ID_FILE" ]]; then
      echo "ERROR: FASTQR_REP_ID_FILE does not exist: $FASTQR_REP_ID_FILE" >&2
      exit 1
    fi
    tr ',[:space:]' '\n' < "$FASTQR_REP_ID_FILE" | awk 'NF && !seen[$1]++ {count++} END{print count+0}'
  elif [[ -n "$FASTQR_REP_ID_LIST" ]]; then
    printf '%s' "$FASTQR_REP_ID_LIST" | tr ',[:space:]' '\n' | awk 'NF && !seen[$1]++ {count++} END{print count+0}'
  else
    printf '%s\n' "$FASTQR_NUM_REP"
  fi
}

RERUN_REP_COUNT="$(count_sparse_rep_ids)"
if ! is_pos_int "$RERUN_REP_COUNT"; then
  echo "ERROR: resolved rerun rep count must be a positive integer, got '$RERUN_REP_COUNT'" >&2
  exit 1
fi

NUM_TASKS=$(( (RERUN_REP_COUNT + FASTQR_CHUNK_SIZE - 1) / FASTQR_CHUNK_SIZE ))
mkdir -p "$PROJECT_DIR/logs" "$PROJECT_DIR/$CHECK_OUT_DIR"
TAU_INT=$(awk -v t="$FASTQR_TAU" 'BEGIN{printf "%02d", int(t*100+0.5)}')
JOB_NAME="llqr_c${FASTQR_CASE}_tau${TAU_INT}_n${FASTQR_N}_rep${FASTQR_NUM_REP}_rerun"

cd "$PROJECT_DIR"

ARRAY_JOB_ID=$(sbatch <<EOF | awk '{print $4}'
#!/bin/bash
#SBATCH -J ${JOB_NAME}
#SBATCH -p ${PARTITION}
#SBATCH -A ${ACCOUNT}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --time=${WALLTIME}
#SBATCH --array=1-${NUM_TASKS}
#SBATCH -o logs/%x.%A_%a.out
#SBATCH -e logs/%x.%A_%a.err

module load compilers/gcc/v12.2.0
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64:\$LD_LIBRARY_PATH
export R_LIBS_USER="/home/$USER/R/x86_64-pc-linux-gnu-library/4.3"

cd ${PROJECT_DIR}
Rscript scripts/array/driver_llqr_array.R
EOF
)

echo "Submitted rerun array job: ${ARRAY_JOB_ID} (tasks: ${NUM_TASKS}, chunk_size: ${FASTQR_CHUNK_SIZE}, rerun_reps: ${RERUN_REP_COUNT})"

MERGE_JOB_ID=$(sbatch --dependency=afterany:${ARRAY_JOB_ID} <<EOF | awk '{print $4}'
#!/bin/bash
#SBATCH -J ${JOB_NAME}_merge
#SBATCH -p ${PARTITION}
#SBATCH -A ${ACCOUNT}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err

module load compilers/gcc/v12.2.0
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64:\$LD_LIBRARY_PATH
export R_LIBS_USER="/home/$USER/R/x86_64-pc-linux-gnu-library/4.3"

cd ${PROJECT_DIR}
Rscript scripts/array/merge_llqr_array.R
EOF
)

echo "Submitted rerun merge job: ${MERGE_JOB_ID} (afterany:${ARRAY_JOB_ID})"

CHECK_JOB_ID=$(sbatch --dependency=afterany:${MERGE_JOB_ID} <<EOF | awk '{print $4}'
#!/bin/bash
#SBATCH -J ${JOB_NAME}_check
#SBATCH -p ${PARTITION}
#SBATCH -A ${ACCOUNT}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err

module load compilers/gcc/v12.2.0
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64:\$LD_LIBRARY_PATH
export R_LIBS_USER="/home/$USER/R/x86_64-pc-linux-gnu-library/4.3"
export FASTQR_LLQR_CHECK_OUT_DIR=${CHECK_OUT_DIR}

cd ${PROJECT_DIR}
Rscript scripts/array/check_llqr_config_integrity.R
EOF
)

echo "Submitted rerun check job: ${CHECK_JOB_ID} (afterany:${MERGE_JOB_ID})"
