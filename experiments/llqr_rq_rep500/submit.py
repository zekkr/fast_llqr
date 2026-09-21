#!/usr/bin/env python3
"""Submit isolated, SHA-gated preflight or formal arrays; never overwrite a run."""
import argparse, os, subprocess, pathlib, json, shlex
p=argparse.ArgumentParser();p.add_argument('stage',choices=['preflight','formal']);p.add_argument('--tag',required=True);p.add_argument('--sha',required=True);p.add_argument('--preflight');a=p.parse_args()
root=pathlib.Path.cwd(); exp='experiments/llqr_rq_rep500'; output=root/'results/hpc/llqr_rq_runs'; run=output/a.tag
assert all(c.isalnum() or c in '-_.' for c in a.tag)
assert subprocess.check_output(['git','rev-parse','HEAD']).decode().strip()==a.sha
subprocess.check_call(['git','diff','--quiet','HEAD','--','R',exp,'reproduction/frozen'])
if a.stage=='formal':
 assert a.preflight and (output/a.preflight/'PREFLIGHT_PASS').read_text().strip()==a.sha
 assert (root/exp/'build/build_git_sha.txt').read_text().strip()==a.sha
run.mkdir(parents=True,exist_ok=False);(run/'logs').mkdir();meta=run/'_run_meta';meta.mkdir()
num_rep=1 if a.stage=='preflight' else 500
params=dict(stage=a.stage,git_sha=a.sha,paper_cases=[1,2],kernel_case=2,ns=[1000,2000,5000,10000],taus=[.2,.5,.8],num_rep=num_rep,seed_base=2025,methods=['direct_baseline','lean_seq','unified_u11'],cache_flags=27,provider_flags=1,include_H_seq=False,preflight=a.preflight)
(meta/'run_config.json').write_text(json.dumps(params,indent=2)+'\n');print(json.dumps(params),flush=True)
base='''export LC_ALL=C LANG=C
source /etc/profile
set -euo pipefail
module load compilers/gcc/v12.2.0 soft/R/v4.3.1
export PATH=/apps/soft/R/R-4.3.1/bin:$PATH
export LD_LIBRARY_PATH=/apps/soft/R/R-4.3.1/lib64/R/lib:/apps/soft/R/R-4.3.1/lib64:${LD_LIBRARY_PATH:-}
export R_LIBS_USER=/home/wuweic/R/x86_64-pc-linux-gnu-library/4.3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
'''
env=dict(SSQR_PROJECT_ROOT=str(root),SSQR_EXPERIMENT_DIR=exp,SSQR_OUTPUT_ROOT=str(output),SSQR_RUN_TAG=a.tag,SSQR_NUM_REP=str(num_rep),SSQR_SEED_BASE='2025',SSQR_MODEL='llqr',SSQR_CACHE_FLAGS='27',SSQR_PROVIDER_FLAGS='1',SSQR_BUILD_MODE='optimized',SSQR_PUSHED_SHA=a.sha,SSQR_REQUIRE_HPC='1',SSQR_ALLOW_SMOKE='1' if num_rep==1 else '0')
base+='\n'.join('export '+k+'='+shlex.quote(v) for k,v in env.items())+'\ncd '+shlex.quote(str(root))+'\n[[ $(git rev-parse HEAD) == '+a.sha+' ]]\n'
manifest=[]
def submit(name,body,deps=None,cpus=1,array=None,exclusive=False):
 script=run/'logs'/(name+'.sh');script.write_text('#!/usr/bin/env bash\n'+base+body+'\n')
 cmd=['/rmprog/slurm/v22.05.7/bin/sbatch','--parsable','-J',name,'-p','cnall','-A','users','-N','1','-n','1','--cpus-per-task='+str(cpus),'--time=08:00:00','-o',str(run/'logs/%x.%A_%a.out'),'-e',str(run/'logs/%x.%A_%a.err')]
 if deps: cmd+=['--dependency=afterok:'+':'.join(deps)]
 if array: cmd+=['--array='+array]
 if exclusive: cmd+=['--exclusive']
 jid=subprocess.check_output(cmd+[str(script)]).decode().strip().split(';')[0];assert jid.isdigit()
 manifest.append(dict(name=name,job_id=jid,cpus=cpus,array=array));(meta/'jobs.json').write_text(json.dumps(manifest,indent=2)+'\n');print(name,jid,flush=True);return jid
build=[]
if a.stage=='preflight':
 build=[submit('rq_build',f'bash {exp}/build_all.sh\nprintf "%s\\n" "$SSQR_PUSHED_SHA" > {exp}/build/build_git_sha.txt\nRscript {exp}/test_baseline.R\nSSQR_CASE=1 SSQR_TAU=0.5 SSQR_N=200 SSQR_NUM_REP=2 SSQR_CHUNK_SIZE=2 SSQR_RUN_TAG={a.tag}_smoke Rscript {exp}/driver_array.R\nSSQR_CASE=2 SSQR_TAU=0.5 SSQR_N=200 SSQR_NUM_REP=2 SSQR_CHUNK_SIZE=2 SSQR_RUN_TAG={a.tag}_smoke Rscript {exp}/driver_array.R')]
else:
 import shutil
 shutil.copytree(root/exp/'build',meta/'build',ignore=shutil.ignore_patterns('*.so'))
merges=[]
for case in [1,2]:
 for tau in [.2,.5,.8]:
  for n in [1000,2000,5000,10000]:
   cpus=1 if num_rep==1 else (14 if n==10000 else 56);chunk=cpus;tasks=(num_rep+chunk-1)//chunk
   spec='1-'+str(tasks)+('%1' if n==10000 and num_rep>1 else '')
   name='rq_c%d_t%d_n%d'%(case,round(tau*100),n)
   config='export SSQR_CASE=%d SSQR_TAU=%s SSQR_N=%d SSQR_CHUNK_SIZE=%d\n'%(case,tau,n,chunk)
   job=submit(name,config+'Rscript '+exp+'/driver_array.R',build,cpus,spec,n==10000 and num_rep>1)
   merges.append(submit(name+'_merge',config+'Rscript '+exp+'/merge_config.R',[job]))
submit('rq_summary','Rscript '+exp+'/summarize_run.R',merges)
print('RUN_DIR='+str(run),flush=True)
