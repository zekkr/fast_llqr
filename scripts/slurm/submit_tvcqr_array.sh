#!/bin/bash
set -euo pipefail

# ---- USER CONFIG ----
PROJECT_DIR="${PROJECT_DIR:-/home/$USER/WORK/Erkang/fast_llqr}"

is_pos_int() {
  [[ "$1" =~ ^[0-9]+$ ]] && [ "$1" -gt 0 ]
}

is_nonneg_int() {
  [[ "$1" =~ ^[0-9]+$ ]]
}

export FASTQR_CASE="${FASTQR_CASE:-1}"
export FASTQR_TAU="${FASTQR_TAU:-0.5}"
export FASTQR_N="${FASTQR_N:-200}"
export FASTQR_NUM_REP="${FASTQR_NUM_REP:-500}"

# Mm.factor as comma-separated list
export FASTQR_MM_FACTOR="${FASTQR_MM_FACTOR:-1e-3,1e-4}"

# algorithm knobs
export FASTQR_H_FACTOR="${FASTQR_H_FACTOR:-1}"
export FASTQR_TOL="${FASTQR_TOL:-1e-14}"
export FASTQR_MAXIT="${FASTQR_MAXIT:-1000000}"
export FASTQR_BLAND="${FASTQR_BLAND:-0}"
export FASTQR_EPS="${FASTQR_EPS:-1e-6}"
export FASTQR_CPP_HELPER="${FASTQR_CPP_HELPER:-0}"
export FASTQR_SEED_BASE="${FASTQR_SEED_BASE:-2025}"
export FASTQR_J="${FASTQR_J:-100}"
export FASTQR_BURN_IN="${FASTQR_BURN_IN:-500}"
export FASTQR_MAX_SECONDS_PER_REP="${FASTQR_MAX_SECONDS_PER_REP:-7200}"

# resources
CPUS_PER_TASK="${CPUS_PER_TASK:-56}"
WALLTIME="${WALLTIME:-04:00:00}"
PARTITION="${PARTITION:-cnall}"
ACCOUNT="${ACCOUNT:-users}"
MERGE_DEPENDENCY_TYPE="${FASTQR_MERGE_DEPENDENCY_TYPE:-afterany}"

# chunking: how many reps each array task handles
# recommended: chunk_size = CPUS_PER_TASK (each core runs one rep in parallel)
export FASTQR_CHUNK_SIZE="${FASTQR_CHUNK_SIZE:-$CPUS_PER_TASK}"

if [[ "$FASTQR_CASE" != "1" && "$FASTQR_CASE" != "2" ]]; then
  echo "ERROR: FASTQR_CASE must be 1 or 2, got '$FASTQR_CASE'" >&2
  exit 1
fi

if ! awk -v t="$FASTQR_TAU" 'BEGIN { exit !(t+0==t && t>0 && t<1) }'; then
  echo "ERROR: FASTQR_TAU must be numeric and in (0,1), got '$FASTQR_TAU'" >&2
  exit 1
fi

if [[ "$MERGE_DEPENDENCY_TYPE" != "afterok" && "$MERGE_DEPENDENCY_TYPE" != "afterany" ]]; then
  echo "ERROR: FASTQR_MERGE_DEPENDENCY_TYPE must be 'afterok' or 'afterany', got '$MERGE_DEPENDENCY_TYPE'" >&2
  exit 1
fi

for v in FASTQR_N FASTQR_NUM_REP CPUS_PER_TASK FASTQR_CHUNK_SIZE FASTQR_MAX_SECONDS_PER_REP; do
  val="${!v}"
  if ! is_pos_int "$val"; then
    echo "ERROR: $v must be a positive integer, got '$val'" >&2
    exit 1
  fi
done

for v in FASTQR_J FASTQR_BURN_IN; do
  val="${!v}"
  if ! is_nonneg_int "$val"; then
    echo "ERROR: $v must be a non-negative integer, got '$val'" >&2
    exit 1
  fi
done

# number of array tasks
NUM_TASKS=$(( (FASTQR_NUM_REP + FASTQR_CHUNK_SIZE - 1) / FASTQR_CHUNK_SIZE ))

mkdir -p "$PROJECT_DIR/logs"

TAU_INT=$(awk -v t="$FASTQR_TAU" 'BEGIN{printf "%02d", int(t*100+0.5)}')
JOB_NAME="tvcqr_c${FASTQR_CASE}_tau${TAU_INT}_n${FASTQR_N}_rep${FASTQR_NUM_REP}"

cd "$PROJECT_DIR"

# ---- submit array ----
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
Rscript scripts/array/driver_tvcqr_array.R
EOF
)

echo "Submitted array job: ${ARRAY_JOB_ID}  (tasks: ${NUM_TASKS}, chunk_size: ${FASTQR_CHUNK_SIZE})"

# ---- submit merge with dependency ----
MERGE_JOB_ID=$(sbatch --dependency=${MERGE_DEPENDENCY_TYPE}:${ARRAY_JOB_ID} <<EOF | awk '{print $4}'
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
Rscript scripts/array/merge_tvcqr_array.R
EOF
)

echo "Submitted merge job: ${MERGE_JOB_ID} (${MERGE_DEPENDENCY_TYPE}:${ARRAY_JOB_ID})"
