#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/home/$USER/WORK/Erkang/fast_llqr}"

export FASTQR_CASE="${FASTQR_CASE:-1}"
export FASTQR_TAU="${FASTQR_TAU:-0.5}"
export FASTQR_N="${FASTQR_N:-500}"
export FASTQR_NUM_REP="${FASTQR_NUM_REP:-500}"

export FASTQR_MM_FACTOR="${FASTQR_MM_FACTOR:-1e-2,1e-3,1e-4}"
export FASTQR_TOL="${FASTQR_TOL:-1e-14}"
export FASTQR_MAXIT="${FASTQR_MAXIT:-2000000}"
export FASTQR_BLAND="${FASTQR_BLAND:-0}"
export FASTQR_TRACK_ORDER="${FASTQR_TRACK_ORDER:-1}"
export FASTQR_SEED_BASE="${FASTQR_SEED_BASE:-2026}"

CPUS_PER_TASK="${CPUS_PER_TASK:-56}"
WALLTIME="${WALLTIME:-04:00:00}"
PARTITION="${PARTITION:-cnall}"
ACCOUNT="${ACCOUNT:-users}"

export FASTQR_CHUNK_SIZE="${FASTQR_CHUNK_SIZE:-$CPUS_PER_TASK}"
NUM_TASKS=$(( (FASTQR_NUM_REP + FASTQR_CHUNK_SIZE - 1) / FASTQR_CHUNK_SIZE ))

mkdir -p "$PROJECT_DIR/logs"

TAU_INT=$(awk -v t="$FASTQR_TAU" 'BEGIN{printf "%02d", int(t*100+0.5)}')
JOB_NAME="llqr_c${FASTQR_CASE}_tau${TAU_INT}_n${FASTQR_N}_rep${FASTQR_NUM_REP}"

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

echo "Submitted array job: ${ARRAY_JOB_ID}  (tasks: ${NUM_TASKS}, chunk_size: ${FASTQR_CHUNK_SIZE})"

MERGE_JOB_ID=$(sbatch --dependency=afterok:${ARRAY_JOB_ID} <<EOF | awk '{print $4}'
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

echo "Submitted merge job: ${MERGE_JOB_ID} (afterok:${ARRAY_JOB_ID})"

