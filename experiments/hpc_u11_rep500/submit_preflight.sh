#!/usr/bin/env bash
set -euo pipefail

: "${SSQR_PROJECT_ROOT:?}" "${SSQR_PREFLIGHT_ROOT:?}" "${SSQR_PREFLIGHT_TAG:?}" "${SSQR_PUSHED_SHA:?}"
experiment_dir="${SSQR_EXPERIMENT_DIR:-experiments/hpc_u11_rep500}"
partition="${SSQR_PARTITION:-cnall}"
account="${SSQR_ACCOUNT:-users}"
sbatch_bin="${SSQR_SBATCH:-/rmprog/slurm/v22.05.7/bin/sbatch}"
[[ "$SSQR_PREFLIGHT_TAG" =~ ^[A-Za-z0-9_.-]+$ && "$SSQR_PUSHED_SHA" =~ ^[0-9a-f]{40}$ ]] || exit 2
[[ "$(cd "$SSQR_PROJECT_ROOT" && git rev-parse HEAD)" == "$SSQR_PUSHED_SHA" ]] || exit 2

preflight_dir="$SSQR_PREFLIGHT_ROOT/$SSQR_PREFLIGHT_TAG"
[[ ! -e "$preflight_dir" ]] || { printf 'preflight tag already exists\n' >&2; exit 2; }
mkdir -p "$preflight_dir/logs" "$preflight_dir/meta"
printf '%s\n' "$SSQR_PUSHED_SHA" > "$preflight_dir/meta/git_sha.txt"

build_job="$($sbatch_bin --parsable <<EOF
#!/usr/bin/env bash
#SBATCH -J u11_pf_build
#SBATCH -p $partition
#SBATCH -A $account
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH -o $preflight_dir/logs/%x.%j.out
#SBATCH -e $preflight_dir/logs/%x.%j.err
set -euo pipefail
module load compilers/gcc/v12.2.0 soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64/R/lib:/apps/soft/R/R-4.3.1/lib64:\${LD_LIBRARY_PATH:-}
export R_LIBS_USER=/home/wuweic/R/x86_64-pc-linux-gnu-library/4.3
export SSQR_PROJECT_ROOT=$SSQR_PROJECT_ROOT SSQR_EXPERIMENT_DIR=$experiment_dir
cd $SSQR_PROJECT_ROOT
[[ \$(git rev-parse HEAD) == $SSQR_PUSHED_SHA ]]
bash $experiment_dir/build_all.sh
cp $experiment_dir/build/compiler_version.txt $preflight_dir/meta/
cp $experiment_dir/build/compiler_flags.txt $preflight_dir/meta/
cp $experiment_dir/build/build_uname.txt $preflight_dir/meta/
cp $experiment_dir/build/build_lscpu.txt $preflight_dir/meta/ 2>/dev/null || true
cp $experiment_dir/build/source_sha256.txt $preflight_dir/meta/
EOF
)"

smoke_job="$($sbatch_bin --parsable --dependency="afterok:$build_job" <<EOF
#!/usr/bin/env bash
#SBATCH -J u11_pf_smoke
#SBATCH -p $partition
#SBATCH -A $account
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH -o $preflight_dir/logs/%x.%j.out
#SBATCH -e $preflight_dir/logs/%x.%j.err
set -euo pipefail
module load compilers/gcc/v12.2.0 soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64/R/lib:/apps/soft/R/R-4.3.1/lib64:\${LD_LIBRARY_PATH:-}
export R_LIBS_USER=/home/wuweic/R/x86_64-pc-linux-gnu-library/4.3
export SSQR_PROJECT_ROOT=$SSQR_PROJECT_ROOT SSQR_EXPERIMENT_DIR=$experiment_dir
export SSQR_PREFLIGHT_DIR=$preflight_dir SSQR_PUSHED_SHA=$SSQR_PUSHED_SHA
cd $SSQR_PROJECT_ROOT
bash $experiment_dir/hpc_pipeline_smoke.sh
EOF
)"

audit_job="$($sbatch_bin --parsable --dependency="afterok:$smoke_job" <<EOF
#!/usr/bin/env bash
#SBATCH -J u11_pf_audit
#SBATCH -p $partition
#SBATCH -A $account
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-48%8
#SBATCH --time=08:00:00
#SBATCH -o $preflight_dir/logs/%x.%A_%a.out
#SBATCH -e $preflight_dir/logs/%x.%A_%a.err
set -euo pipefail
module load compilers/gcc/v12.2.0 soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64/R/lib:/apps/soft/R/R-4.3.1/lib64:\${LD_LIBRARY_PATH:-}
export R_LIBS_USER=/home/wuweic/R/x86_64-pc-linux-gnu-library/4.3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export SSQR_PROJECT_ROOT=$SSQR_PROJECT_ROOT SSQR_EXPERIMENT_DIR=$experiment_dir
export SSQR_PREFLIGHT_DIR=$preflight_dir
cd $SSQR_PROJECT_ROOT
[[ \$(git rev-parse HEAD) == $SSQR_PUSHED_SHA ]]
Rscript $experiment_dir/audit_cache_modes.R
EOF
)"

validate_job="$($sbatch_bin --parsable --dependency="afterany:$audit_job" <<EOF
#!/usr/bin/env bash
#SBATCH -J u11_pf_check
#SBATCH -p $partition
#SBATCH -A $account
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH -o $preflight_dir/logs/%x.%j.out
#SBATCH -e $preflight_dir/logs/%x.%j.err
set -euo pipefail
module load soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export SSQR_PREFLIGHT_DIR=$preflight_dir SSQR_PUSHED_SHA=$SSQR_PUSHED_SHA
cd $SSQR_PROJECT_ROOT
Rscript $experiment_dir/validate_preflight.R
EOF
)"

printf 'build_job\t%s\nsmoke_job\t%s\naudit_job\t%s\nvalidate_job\t%s\n' \
  "$build_job" "$smoke_job" "$audit_job" "$validate_job" | tee "$preflight_dir/meta/job_ids.tsv"
