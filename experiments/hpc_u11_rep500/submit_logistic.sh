#!/usr/bin/env bash
# Login-node entry point. Submit formal arrays only after successful cluster preflight.
set -euo pipefail
: "${SSQR_PUSHED_SHA:?Set the verified GitHub commit SHA}"
project="${SSQR_PROJECT_ROOT:-$(pwd)}"
cd "$project"
project="$(pwd -P)"
experiment_dir=experiments/hpc_u11_rep500
[[ "$(git rev-parse HEAD)" == "$SSQR_PUSHED_SHA" ]] || { echo 'Wrong checkout SHA' >&2; exit 2; }
git diff --quiet HEAD -- R src "$experiment_dir" || { echo 'Tracked experiment sources have local edits' >&2; exit 2; }
git diff --quiet 61dc4d6b39fab585ca62c590a7024fcc9eb95751 -- R src "$experiment_dir/fortran" || {
  echo 'Frozen baseline or solver sources changed' >&2; exit 2;
}
sbatch_bin="${SSQR_SBATCH:-/rmprog/slurm/v22.05.7/bin/sbatch}"
partition="${SSQR_PARTITION:-cnall}"; account="${SSQR_ACCOUNT:-users}"
run_tag="${SSQR_RUN_TAG:-llqr_c2_logistic_rep500_seed2025_${SSQR_PUSHED_SHA:0:7}_$(date -u +%Y%m%d_%H%M%S)}"
[[ "$run_tag" =~ ^[A-Za-z0-9_.-]+$ ]] || exit 2
output_root="${SSQR_OUTPUT_ROOT:-$project/results/hpc/llqr_case2_logistic_runs}"
launch="$output_root/${run_tag}_launch"
[[ ! -e "$launch" && ! -e "$output_root/$run_tag" ]] || { echo 'Run tag already exists; inspect existing jobs before retrying' >&2; exit 2; }
mkdir -p "$launch"
export SSQR_PROJECT_ROOT="$project" SSQR_EXPERIMENT_DIR="$experiment_dir"
export SSQR_PREFLIGHT_ROOT="$output_root/preflight" SSQR_PREFLIGHT_TAG="${run_tag}_preflight"
export SSQR_OUTPUT_ROOT="$output_root" SSQR_RUN_TAG="$run_tag" SSQR_NUM_REP=500 SSQR_SEED_BASE=2025
export SSQR_PREFLIGHT_DIR="$SSQR_PREFLIGHT_ROOT/$SSQR_PREFLIGHT_TAG"
bash "$experiment_dir/submit_preflight.sh" > "$launch/preflight_submission.log" 2>&1
cat "$launch/preflight_submission.log"
validate_job="$(awk '$1=="validate_job" {print $2}' "$SSQR_PREFLIGHT_DIR/meta/job_ids.tsv")"
[[ "$validate_job" =~ ^[0-9]+$ ]] || { echo 'Missing preflight validation job ID' >&2; exit 2; }
dispatch_job="$($sbatch_bin --parsable --dependency="afterok:$validate_job" <<SBATCH
#!/usr/bin/env bash
#SBATCH -J logistic_dispatch
#SBATCH -p $partition
#SBATCH -A $account
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH -o $launch/dispatch.%j.out
#SBATCH -e $launch/dispatch.%j.err
set -euo pipefail
export SSQR_PROJECT_ROOT=$project SSQR_EXPERIMENT_DIR=$experiment_dir
export SSQR_OUTPUT_ROOT=$output_root SSQR_RUN_TAG=$run_tag
export SSQR_PUSHED_SHA=$SSQR_PUSHED_SHA SSQR_PREFLIGHT_DIR=$SSQR_PREFLIGHT_DIR
export SSQR_NUM_REP=500 SSQR_SEED_BASE=2025
export SSQR_SBATCH=$sbatch_bin SSQR_PARTITION=$partition SSQR_ACCOUNT=$account
cd $project
[[ \$(git rev-parse HEAD) == $SSQR_PUSHED_SHA ]]
bash $experiment_dir/submit_all.sh > $launch/formal_submission.log 2>&1
cat $launch/formal_submission.log
printf 'formal_submission_complete\\n' > $launch/FORMAL_SUBMITTED
SBATCH
)"
[[ "$dispatch_job" =~ ^[0-9]+$ ]] || { echo 'Missing dispatch job ID' >&2; exit 2; }
{
 printf 'RUN_TAG=%s\nRUN_DIR=%s\nPREFLIGHT_DIR=%s\nLAUNCH_DIR=%s\nDISPATCH_JOB=%s\nSHA=%s\n' \
 "$run_tag" "$output_root/$run_tag" "$SSQR_PREFLIGHT_DIR" "$launch" "$dispatch_job" "$SSQR_PUSHED_SHA"
 printf '\nFormal rep500 jobs will only be submitted after preflight passes.\n'
 printf 'If preflight fails, inspect %s/logs and cancel the pending dispatcher %s before a corrected new launch.\n' "$SSQR_PREFLIGHT_DIR" "$dispatch_job"
} | tee "$launch/LAUNCH_INFO.txt"
