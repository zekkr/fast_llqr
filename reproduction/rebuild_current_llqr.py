#!/usr/bin/env python3
"""Rebuild current LLQR summary tables from the fixed paired-replication archive."""
from pathlib import Path
import argparse,tarfile,subprocess,os,json,csv,hashlib
p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args()
root=Path(__file__).resolve().parents[1];tag='llqr_rq_rep500_seed2025_4444298_20260921'
archive=root/'reproduction/archive'/(tag+'.tar.gz')
expected=archive.with_suffix(archive.suffix+'.sha256').read_text().split()[0]
assert hashlib.sha256(archive.read_bytes()).hexdigest()==expected,'Archive checksum mismatch'
assert not a.output.exists(),'Output must be a new directory'
a.output.mkdir(parents=True);out=a.output.resolve()
with tarfile.open(archive) as tf:
 for m in tf.getmembers():
  assert m.isfile(),m.name
  dest=out/m.name;dest.resolve().relative_to(out)
  dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(tf.extractfile(m).read())
run=out/tag;sha=(run/'RESULTS_PASS').read_text().strip()
env=dict(os.environ,SSQR_PROJECT_ROOT=str(root),SSQR_OUTPUT_ROOT=str(out),SSQR_RUN_TAG=tag,SSQR_NUM_REP='500',SSQR_SEED_BASE='2025',SSQR_PUSHED_SHA=sha)
subprocess.run(['Rscript','experiments/llqr_rq_rep500/summarize_run.R'],cwd=root,env=env,check=True)
subprocess.run(['python3','experiments/llqr_rq_rep500/verify_results.py',str(run)],cwd=root,check=True)
read=lambda f:list(csv.DictReader(f.open()))
old=read(root/'reproduction/archive/original_summary.csv')
new=read(run/'tables/config_method_summary.csv')
rows=[r for r in old if r['model']=='tvcqr']+[dict(r,paper_case=r['case']) for r in new]
assert len(rows)==144
fields=sorted(set().union(*(r.keys() for r in rows)))
with (out/'current_paper_method_summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
(out/'source_manifest.json').write_text(json.dumps(dict(llqr_tag=tag,llqr_sha=sha,llqr_archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),tvcqr_source='reproduction/archive/original_summary.csv',scope='24 new LLQR + 24 retained TVCQR configurations; no new simulation fitting'),indent=2)+'\n')
print('PASS current-paper archived reconstruction:',out)
