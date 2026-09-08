#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/home/wuweic/WORK/Erkang/fast_llqr_codex_test}"
RUN_TAG="${FASTQR_RUN_TAG:-v41v42_rep100_seed2026_$(date +%Y%m%d_%H%M%S)}"
BASE_DIR="${FASTQR_V4142_BASE_DIR:-data/v41_v42_rep100}"
NUM_REP="${FASTQR_NUM_REP:-100}"
SEED_BASE="${FASTQR_SEED_BASE:-2026}"
PARTITION="${PARTITION:-cnall}"
ACCOUNT="${ACCOUNT:-users}"
WALLTIME="${WALLTIME:-08:00:00}"

if [[ "${NUM_REP}" != "100" ]]; then
  printf 'Refusing non-approved replication count: %s\n' "${NUM_REP}" >&2
  exit 2
fi
if [[ "${SEED_BASE}" != "2026" ]]; then
  printf 'Refusing non-approved seed base: %s\n' "${SEED_BASE}" >&2
  exit 2
fi

cd "${PROJECT_DIR}"

for binary in llqr_v41.so llqr_v42.so llqr_lean_seq.so tvcqr_v41.so tvcqr_v42.so tvcqr_lean_seq.so; do
  if [[ ! -s "experiments/hpc_v41_v42_rep100/build/${binary}" ]]; then
    printf 'Missing experiment binary: %s\n' "${binary}" >&2
    exit 2
  fi
done

LOG_DIR="logs/v41_v42_rep100/${RUN_TAG}"
mkdir -p "${LOG_DIR}"
MANIFEST="${LOG_DIR}/submission_manifest.tsv"
printf 'run_tag\tmodel\tcase\ttau\tn\tarray_job_id\tmerge_job_id\tcpus\tchunk_size\n' > "${MANIFEST}"

merge_jobs=()

submit_config() {
  local model="$1"
  local case_id="$2"
  local tau="$3"
  local n="$4"
  local cpus chunk_size array_tasks tau_tag job_name array_job merge_job

  if [[ "${n}" -eq 10000 ]]; then
    cpus=14
    chunk_size=14
  else
    cpus=56
    chunk_size=56
  fi
  array_tasks=$(( (NUM_REP + chunk_size - 1) / chunk_size ))
  tau_tag="$(awk -v value="${tau}" 'BEGIN { printf "%02d", int(100 * value + 0.5) }')"
  job_name="v4142_${model}_c${case_id}_t${tau_tag}_n${n}"

  array_job="$(sbatch --parsable <<EOF
#!/usr/bin/env bash
#SBATCH -J ${job_name}
#SBATCH -p ${PARTITION}
#SBATCH -A ${ACCOUNT}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --time=${WALLTIME}
#SBATCH --array=1-${array_tasks}
#SBATCH -o ${LOG_DIR}/%x.%A_%a.out
#SBATCH -e ${LOG_DIR}/%x.%A_%a.err

set -euo pipefail
module load compilers/gcc/v12.2.0 soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64/R/lib:/apps/soft/R/R-4.3.1/lib64:\${LD_LIBRARY_PATH:-}
export R_LIBS_USER=/home/wuweic/R/x86_64-pc-linux-gnu-library/4.3
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export FASTQR_PROJECT_DIR=${PROJECT_DIR}
export FASTQR_MODEL=${model}
export FASTQR_CASE=${case_id}
export FASTQR_TAU=${tau}
export FASTQR_N=${n}
export FASTQR_NUM_REP=${NUM_REP}
export FASTQR_SEED_BASE=${SEED_BASE}
export FASTQR_CHUNK_SIZE=${chunk_size}
export FASTQR_RUN_TAG=${RUN_TAG}
export FASTQR_V4142_BASE_DIR=${BASE_DIR}
cd ${PROJECT_DIR}
Rscript experiments/hpc_v41_v42_rep100/driver_array.R
EOF
)"

  merge_job="$(sbatch --parsable --dependency="afterany:${array_job}" <<EOF
#!/usr/bin/env bash
#SBATCH -J ${job_name}_merge
#SBATCH -p ${PARTITION}
#SBATCH -A ${ACCOUNT}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH -o ${LOG_DIR}/%x.%j.out
#SBATCH -e ${LOG_DIR}/%x.%j.err

set -euo pipefail
module load soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export FASTQR_PROJECT_DIR=${PROJECT_DIR}
export FASTQR_MODEL=${model}
export FASTQR_CASE=${case_id}
export FASTQR_TAU=${tau}
export FASTQR_N=${n}
export FASTQR_NUM_REP=${NUM_REP}
export FASTQR_SEED_BASE=${SEED_BASE}
export FASTQR_RUN_TAG=${RUN_TAG}
export FASTQR_V4142_BASE_DIR=${BASE_DIR}
cd ${PROJECT_DIR}
Rscript experiments/hpc_v41_v42_rep100/merge_config.R
EOF
)"

  merge_jobs+=("${merge_job}")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${RUN_TAG}" "${model}" "${case_id}" "${tau}" "${n}" \
    "${array_job}" "${merge_job}" "${cpus}" "${chunk_size}" >> "${MANIFEST}"
  printf 'submitted %s case=%s tau=%s n=%s array=%s merge=%s\n' \
    "${model}" "${case_id}" "${tau}" "${n}" "${array_job}" "${merge_job}"
}

for model in llqr tvcqr; do
  for case_id in 1 2; do
    for tau in 0.2 0.5 0.8; do
      for n in 1000 2000 5000 10000; do
        submit_config "${model}" "${case_id}" "${tau}" "${n}"
      done
    done
  done
done

dependency="$(IFS=:; printf '%s' "${merge_jobs[*]}")"
summary_job="$(sbatch --parsable --dependency="afterany:${dependency}" <<EOF
#!/usr/bin/env bash
#SBATCH -J v4142_summary
#SBATCH -p ${PARTITION}
#SBATCH -A ${ACCOUNT}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH -o ${LOG_DIR}/%x.%j.out
#SBATCH -e ${LOG_DIR}/%x.%j.err

set -euo pipefail
module load soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export FASTQR_PROJECT_DIR=${PROJECT_DIR}
export FASTQR_NUM_REP=${NUM_REP}
export FASTQR_SEED_BASE=${SEED_BASE}
export FASTQR_RUN_TAG=${RUN_TAG}
export FASTQR_V4142_BASE_DIR=${BASE_DIR}
cd ${PROJECT_DIR}
Rscript experiments/hpc_v41_v42_rep100/summarize_run.R
EOF
)"

printf 'summary_job_id\t%s\n' "${summary_job}" >> "${MANIFEST}"
printf 'RUN_TAG=%s\nMANIFEST=%s\nSUMMARY_JOB=%s\n' "${RUN_TAG}" "${MANIFEST}" "${summary_job}"
