#!/usr/bin/env python3
"""Independent replication-level audit; does not import R aggregation code."""
import argparse,csv,json,math,statistics,itertools,hashlib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('run',type=Path);a=p.parse_args();run=a.run
read=lambda path:list(csv.DictReader(path.open()))
summary=read(run/'tables/config_method_summary.csv');index={(int(r['case']),float(r['tau']),int(r['n']),r['method']):r for r in summary}
methods=['direct_baseline','lean_seq','unified_u11'];stats=[];regen=0;points=0;recovery=0;rawhashes={};checks=0
assert len(index)==72
for case,tau,n in itertools.product((1,2),(.2,.5,.8),(1000,2000,5000,10000)):
 d=run/'llqr'/f'case{case}_tau{round(100*tau):02d}_n{n}'
 raw=read(d/'replication_metrics.csv');attempts=read(d/'attempt_metrics.csv');integ=read(d/'integrity.csv')[0]
 assert integ['complete']=='TRUE' and int(integ['n_complete'])==500 and int(integ['n_missing'])==int(integ['n_malformed'])==int(integ['n_extra'])==0
 assert len(raw)==1500
 assert len({(r['rep_id'],r['method']) for r in raw})==1500
 for rep in range(1,501):
  rows=[r for r in raw if int(r['rep_id'])==rep];aa=[r for r in attempts if int(r['rep_id'])==rep]
  assert len(rows)==3 and {r['method'] for r in rows}==set(methods)
  assert len({(r['seed'],r['dataset_md5'],r['n_eval']) for r in rows})==1
  assert all(len(r['dataset_md5'])==32 for r in rows)
  last=max(int(r['attempt']) for r in aa);assert len(aa)==3*last and 1<=last<=20
  assert all(int(r['seed'])==2025+rep+(int(r['attempt'])-1)*1000000 for r in aa)
  for at in range(1,last):
   b=[r for r in aa if int(r['attempt'])==at and r['method']=='direct_baseline'];assert len(b)==1 and b[0]['threw_error']=='TRUE' and b[0]['error_stage']=='fit'
  assert all(int(r['attempt'])==last for r in rows)
  regen+=last>1;points+=int(rows[0]['n_eval'])
 for m in methods:
  x=[r for r in raw if r['method']==m];s=index[case,tau,n,m]
  assert all(r['solver_ok']==r['accepted_ok']==r['finite_estimate']=='TRUE' and r['discrepancy_status']=='ok' and r['fallback_triggered']=='FALSE' and r['iteration_limit_hit']=='FALSE' for r in x)
  assert {int(r['rep_id']) for r in x}==set(range(1,501))
  if m!='direct_baseline':assert all(r['h_check_status']=='match' and r['h_mismatch_eval_count']=='0' for r in x)
  times=[float(r['elapsed_sec']) for r in x];dis=[float(r['discrepancy']) for r in x]
  assert all(math.isfinite(v) and v>=0 for v in times+dis)
  values={'time_mean_sec':statistics.mean(times),'time_min_sec':min(times),'time_max_sec':max(times),'max_average_relative_discrepancy':max(dis)}
  for k,v in values.items():assert math.isclose(v,float(s[k]),rel_tol=1e-12,abs_tol=1e-15),(case,tau,n,m,k);checks+=1
  if m=='unified_u11':
   rr=sum(int(r['internal_recovery_total']) for r in x);recovery+=rr
   assert rr==int(s['internal_recovery_total'])
  stats.append(dict(case=case,tau=tau,n=n,method=m,**values))
 for path in (d/'replication_metrics.csv',d/'attempt_metrics.csv'):
  rawhashes[str(path.relative_to(run))]=hashlib.sha256(path.read_bytes()).hexdigest()
result=dict(configurations=24,datasets=12000,fits=36000,compared_H_points=points,H_mismatches=0,seq_fallbacks=0,full_active_recoveries=recovery,regenerated_replications=regen,independent_numeric_checks=checks,source_sha=(run/'RESULTS_PASS').read_text().strip(),raw_sha256=rawhashes)
(run/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='raw_sha256'},indent=2))
