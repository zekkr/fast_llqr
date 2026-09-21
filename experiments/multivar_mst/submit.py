#!/usr/bin/env python3
"""Submit SHA-gated multivariate MST preflight or formal arrays."""
import argparse,json,os,pathlib,shlex,subprocess

p=argparse.ArgumentParser()
p.add_argument('stage',choices=['preflight','formal'])
p.add_argument('--tag',required=True);p.add_argument('--sha',required=True)
p.add_argument('--preflight');a=p.parse_args()
root=pathlib.Path.cwd();exp='experiments/multivar_mst'
output=root/'results/hpc/multivar_mst_runs';run=output/a.tag
assert all(c.isalnum() or c in '-_.' for c in a.tag)
assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==a.sha
subprocess.check_call(['git','diff','--quiet','HEAD','--',exp])
if a.stage=='formal':
    assert a.preflight
    marker=output/a.preflight/'PREFLIGHT_PASS'
    assert marker.read_text().strip()==a.sha
    assert (root/exp/'build/build_git_sha.txt').read_text().strip()==a.sha
run.mkdir(parents=True,exist_ok=False);(run/'logs').mkdir();meta=run/'_run_meta';meta.mkdir()
num_rep=3 if a.stage=='preflight' else 100
config=dict(stage=a.stage,git_sha=a.sha,cases=[1,2],ns=[500,1000,2000],tau=.5,
            dimension=4,num_rep=num_rep,seed_base=2025,
            methods=['direct_fit','seq_screen_mst'],include_H_seq=False,
            direct_backend='quantreg::rq(method="br")',threshold_factor=.1,
            preflight=a.preflight)
(meta/'run_config.json').write_text(json.dumps(config,indent=2)+'\n')
print(json.dumps(config),flush=True)

base='''export LC_ALL=C LANG=C
source /etc/profile
set -euo pipefail
module load compilers/gcc/v12.2.0 soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64/R/lib:/apps/soft/R/R-4.3.1/lib64:${LD_LIBRARY_PATH:-}
export R_LIBS_USER=/home/wuweic/R/x86_64-pc-linux-gnu-library/4.3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
'''
env=dict(MST_PROJECT_ROOT=str(root),MST_OUTPUT_ROOT=str(output),MST_RUN_TAG=a.tag,
         MST_NUM_REP=str(num_rep),MST_SEED_BASE='2025',MST_PUSHED_SHA=a.sha,
         MST_REQUIRE_HPC='1',MST_ALLOW_SMOKE='1' if a.stage=='preflight' else '0')
base+='\n'.join('export '+k+'='+shlex.quote(v) for k,v in env.items())
base+='\ncd '+shlex.quote(str(root))+'\n[[ $(git rev-parse HEAD) == '+a.sha+' ]]\n'
manifest=[]
def submit(name,body,deps=None,array=None,time='08:00:00'):
    script=run/'logs'/(name+'.sh');script.write_text('#!/usr/bin/env bash\n'+base+body+'\n')
    cmd=['/rmprog/slurm/v22.05.7/bin/sbatch','--parsable','-J',name,'-p','cnall','-A','users',
         '-N','1','-n','1','--cpus-per-task=1','--time='+time,
         '-o',str(run/'logs/%x.%A_%a.out'),'-e',str(run/'logs/%x.%A_%a.err')]
    if deps:cmd+=['--dependency=afterok:'+':'.join(deps)]
    if array:cmd+=['--array='+array]
    jid=subprocess.check_output(cmd+[str(script)],text=True).strip().split(';')[0]
    assert jid.isdigit();manifest.append(dict(name=name,job_id=jid,array=array))
    (meta/'jobs.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(name,jid,flush=True);return jid

build=[]
if a.stage=='preflight':
    build=[submit('mst_build',f'''bash {exp}/build.sh
printf "%s\\n" "$MST_PUSHED_SHA" > {exp}/build/build_git_sha.txt
MST_BUILD_MODE=checked Rscript {exp}/tests/test_solver.R
MST_BUILD_MODE=optimized Rscript {exp}/tests/test_solver.R''',time='00:30:00')]
else:
    import shutil
    shutil.copytree(root/exp/'build',meta/'build',ignore=shutil.ignore_patterns('*.so'))

merges=[]
for case in [1,2]:
    for n in [500,1000,2000]:
        name=f'mst_c{case}_n{n}'
        setup=f'export MST_CASE={case} MST_N={n} MST_CHUNK_SIZE=1\n'
        spec=f'1-{num_rep}%{min(num_rep,50)}'
        job=submit(name,setup+f'Rscript {exp}/driver_array.R',build,array=spec)
        merges.append(submit(name+'_merge',setup+f'Rscript {exp}/merge_config.R',[job],time='00:20:00'))
summary_body=f'Rscript {exp}/summarize_run.R\n'
if a.stage=='preflight':summary_body+='printf "%s\\n" "$MST_PUSHED_SHA" > "$MST_OUTPUT_ROOT/$MST_RUN_TAG/PREFLIGHT_PASS"\n'
submit('mst_summary',summary_body,merges,time='00:20:00')
print('RUN_DIR='+str(run),flush=True)
