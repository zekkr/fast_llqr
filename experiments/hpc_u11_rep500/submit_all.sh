#!/usr/bin/env bash
set -euo pipefail

: "${SSQR_PROJECT_ROOT:?}" "${SSQR_OUTPUT_ROOT:?}" "${SSQR_RUN_TAG:?}" "${SSQR_PUSHED_SHA:?}" "${SSQR_PREFLIGHT_DIR:?}"
experiment_dir="${SSQR_EXPERIMENT_DIR:-experiments/hpc_u11_rep500}"
partition="${SSQR_PARTITION:-cnall}"
account="${SSQR_ACCOUNT:-users}"
walltime="${SSQR_WALLTIME:-08:00:00}"
sbatch_bin="${SSQR_SBATCH:-/rmprog/slurm/v22.05.7/bin/sbatch}"
num_rep="${SSQR_NUM_REP:-500}"
seed_base="${SSQR_SEED_BASE:-2025}"

[[ "$num_rep" == 500 && "$seed_base" == 2025 ]] || {
  printf 'paper run requires rep500 and seed_base=2025\n' >&2
  exit 2
}
[[ "$SSQR_RUN_TAG" =~ ^[A-Za-z0-9_.-]+$ && "$SSQR_PUSHED_SHA" =~ ^[0-9a-f]{40}$ ]] || exit 2
remote_sha="$(cd "$SSQR_PROJECT_ROOT" && git rev-parse HEAD)"
[[ "$remote_sha" == "$SSQR_PUSHED_SHA" ]] || {
  printf 'remote_sha=%s differs from pushed_sha=%s\n' "$remote_sha" "$SSQR_PUSHED_SHA" >&2
  exit 2
}
[[ -f "$SSQR_PREFLIGHT_DIR/PREFLIGHT_PASS" ]] || {
  printf 'missing successful cluster preflight: %s\n' "$SSQR_PREFLIGHT_DIR/PREFLIGHT_PASS" >&2
  exit 2
}
[[ "$(tr -d '[:space:]' < "$SSQR_PREFLIGHT_DIR/PREFLIGHT_PASS")" == "$SSQR_PUSHED_SHA" ]] || {
  printf 'preflight SHA does not match pushed SHA\n' >&2
  exit 2
}

run_dir="$SSQR_OUTPUT_ROOT/$SSQR_RUN_TAG"
log_dir="$SSQR_PROJECT_ROOT/logs/u11_rep500/$SSQR_RUN_TAG"
[[ ! -e "$run_dir" && ! -e "$log_dir" ]] || {
  printf 'run tag already exists\n' >&2
  exit 2
}
mkdir -p "$run_dir/_run_meta" "$log_dir"
printf '%s\n' "$SSQR_PUSHED_SHA" > "$run_dir/_run_meta/git_sha.txt"
printf '%s\n' "$SSQR_PREFLIGHT_DIR" > "$run_dir/_run_meta/preflight_dir.txt"
cp "$SSQR_PREFLIGHT_DIR/cache_mode_audit.csv" "$run_dir/_run_meta/"
cp "$SSQR_PREFLIGHT_DIR/meta/job_ids.tsv" "$run_dir/_run_meta/preflight_job_ids.tsv"
printf 'run_tag=%s\nnum_rep=%s\nseed_base=%s\ncache_flags=27\nprovider_flags=1\n' \
  "$SSQR_RUN_TAG" "$num_rep" "$seed_base" > "$run_dir/_run_meta/run_config.txt"

manifest="$run_dir/_run_meta/submission_manifest.tsv"
printf 'run_tag\tmodel\tcase\ttau\tn\tarray_job\tmerge_job\tcpus\tchunk\tarray_spec\n' > "$manifest"

build_job="$($sbatch_bin --parsable <<EOF
#!/usr/bin/env bash
#SBATCH -J u11_build
#SBATCH -p $partition
#SBATCH -A $account
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH -o $log_dir/%x.%j.out
#SBATCH -e $log_dir/%x.%j.err
set -euo pipefail
module load compilers/gcc/v12.2.0 soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64/R/lib:/apps/soft/R/R-4.3.1/lib64:\${LD_LIBRARY_PATH:-}
export R_LIBS_USER=/home/wuweic/R/x86_64-pc-linux-gnu-library/4.3
export SSQR_PROJECT_ROOT=$SSQR_PROJECT_ROOT SSQR_EXPERIMENT_DIR=$experiment_dir
cd $SSQR_PROJECT_ROOT
[[ \$(git rev-parse HEAD) == $SSQR_PUSHED_SHA ]]
bash $experiment_dir/build_all.sh
cp $experiment_dir/build/compiler_version.txt $run_dir/_run_meta/
cp $experiment_dir/build/compiler_flags.txt $run_dir/_run_meta/
cp $experiment_dir/build/build_uname.txt $run_dir/_run_meta/
cp $experiment_dir/build/build_lscpu.txt $run_dir/_run_meta/ 2>/dev/null || true
cp $experiment_dir/build/source_sha256.txt $run_dir/_run_meta/
EOF
)"

merge_jobs=()
submit_config() {
  local model="$1" case_id="$2" tau="$3" n="$4"
  local cpus chunk tasks spec exclusive tau_tag name array_job merge_job
  if [[ "$n" == 10000 ]]; then
    cpus=14; chunk=14; tasks=36; spec="1-36%1"; exclusive="#SBATCH --exclusive"
  else
    cpus=56; chunk=56; tasks=9; spec="1-9"; exclusive=""
  fi
  tau_tag="$(awk -v x="$tau" 'BEGIN{printf "%02d",int(100*x+.5)}')"
  name="u11_${model}_c${case_id}_t${tau_tag}_n${n}"
  array_job="$($sbatch_bin --parsable --dependency="afterok:$build_job" <<EOF
#!/usr/bin/env bash
#SBATCH -J $name
#SBATCH -p $partition
#SBATCH -A $account
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=$cpus
#SBATCH --time=$walltime
#SBATCH --array=$spec
$exclusive
#SBATCH -o $log_dir/%x.%A_%a.out
#SBATCH -e $log_dir/%x.%A_%a.err
set -euo pipefail
module load compilers/gcc/v12.2.0 soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64/R/lib:/apps/soft/R/R-4.3.1/lib64:\${LD_LIBRARY_PATH:-}
export R_LIBS_USER=/home/wuweic/R/x86_64-pc-linux-gnu-library/4.3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export SSQR_PROJECT_ROOT=$SSQR_PROJECT_ROOT SSQR_EXPERIMENT_DIR=$experiment_dir
export SSQR_OUTPUT_ROOT=$SSQR_OUTPUT_ROOT SSQR_RUN_TAG=$SSQR_RUN_TAG
export SSQR_MODEL=$model SSQR_CASE=$case_id SSQR_TAU=$tau SSQR_N=$n
export SSQR_NUM_REP=$num_rep SSQR_SEED_BASE=$seed_base SSQR_CHUNK_SIZE=$chunk
export SSQR_CACHE_FLAGS=27 SSQR_PROVIDER_FLAGS=1 SSQR_BUILD_MODE=optimized
cd $SSQR_PROJECT_ROOT
[[ \$(git rev-parse HEAD) == $SSQR_PUSHED_SHA ]]
Rscript $experiment_dir/driver_array.R
EOF
)"
  merge_job="$($sbatch_bin --parsable --dependency="afterany:$array_job" <<EOF
#!/usr/bin/env bash
#SBATCH -J ${name}_merge
#SBATCH -p $partition
#SBATCH -A $account
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH -o $log_dir/%x.%j.out
#SBATCH -e $log_dir/%x.%j.err
set -euo pipefail
module load soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export SSQR_PROJECT_ROOT=$SSQR_PROJECT_ROOT SSQR_EXPERIMENT_DIR=$experiment_dir
export SSQR_OUTPUT_ROOT=$SSQR_OUTPUT_ROOT SSQR_RUN_TAG=$SSQR_RUN_TAG
export SSQR_MODEL=$model SSQR_CASE=$case_id SSQR_TAU=$tau SSQR_N=$n
export SSQR_NUM_REP=$num_rep SSQR_SEED_BASE=$seed_base
cd $SSQR_PROJECT_ROOT
Rscript $experiment_dir/merge_config.R
EOF
)"
  merge_jobs+=("$merge_job")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$SSQR_RUN_TAG" "$model" "$case_id" "$tau" "$n" "$array_job" "$merge_job" \
    "$cpus" "$chunk" "$spec" >> "$manifest"
  printf 'submitted %s case=%s tau=%s n=%s array=%s merge=%s\n' \
    "$model" "$case_id" "$tau" "$n" "$array_job" "$merge_job"
}

for model in llqr tvcqr; do
  for case_id in 1 2; do
    for tau in 0.2 0.5 0.8; do
      for n in 1000 2000 5000 10000; do
        submit_config "$model" "$case_id" "$tau" "$n"
      done
    done
  done
done

dependency="$(IFS=:; printf '%s' "${merge_jobs[*]}")"
summary_job="$($sbatch_bin --parsable --dependency="afterany:$dependency" <<EOF
#!/usr/bin/env bash
#SBATCH -J u11_summary
#SBATCH -p $partition
#SBATCH -A $account
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH -o $log_dir/%x.%j.out
#SBATCH -e $log_dir/%x.%j.err
set -euo pipefail
module load soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:\$PATH
export SSQR_PROJECT_ROOT=$SSQR_PROJECT_ROOT SSQR_EXPERIMENT_DIR=$experiment_dir
export SSQR_OUTPUT_ROOT=$SSQR_OUTPUT_ROOT SSQR_RUN_TAG=$SSQR_RUN_TAG
export SSQR_NUM_REP=$num_rep SSQR_SEED_BASE=$seed_base
cd $SSQR_PROJECT_ROOT
Rscript $experiment_dir/summarize_run.R
Rscript $experiment_dir/build_paper_staging.R
EOF
)"
printf 'kind\tjob_id\nbuild_job_id\t%s\nsummary_job_id\t%s\n' "$build_job" "$summary_job" \
  > "$run_dir/_run_meta/job_ids.tsv"

finalize_job="$($sbatch_bin --parsable --dependency="afterany:$summary_job" <<EOF
#!/usr/bin/env bash
#SBATCH -J u11_finalize
#SBATCH -p $partition
#SBATCH -A $account
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH -o $log_dir/%x.%j.out
#SBATCH -e $log_dir/%x.%j.err
set -euo pipefail
export SSQR_RUN_DIR=$run_dir SSQR_LOG_DIR=$log_dir
cd $SSQR_PROJECT_ROOT
bash $experiment_dir/finalize_run.sh
EOF
)"
printf 'finalize_job_id\t%s\n' "$finalize_job" >> "$run_dir/_run_meta/job_ids.tsv"
printf 'RUN_TAG=%s\nBUILD_JOB=%s\nSUMMARY_JOB=%s\nFINALIZE_JOB=%s\nMANIFEST=%s\n' \
  "$SSQR_RUN_TAG" "$build_job" "$summary_job" "$finalize_job" "$manifest"
