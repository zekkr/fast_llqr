#!/usr/bin/env python3
import argparse,csv,json,math,pathlib
p=argparse.ArgumentParser();p.add_argument('run',type=pathlib.Path);a=p.parse_args()
summary=list(csv.DictReader((a.run/'tables/runtime_summary.csv').open()))
issues=[];raw_count=0
for s in summary:
    case,n=int(s['case']),int(s['n'])
    rows=list(csv.DictReader((a.run/f'case{case}_n{n}'/'replication_metrics.csv').open()))
    raw_count+=len(rows)
    if len(rows)!=100 or [int(r['rep_id']) for r in rows]!=list(range(1,101)):
        issues.append(f'bad replication IDs case={case} n={n}')
    for r in rows:
        if r['direct_ok']!='TRUE' or r['screen_ok']!='TRUE' or r['audit_ok']!='TRUE' or r['h_set_match']!='TRUE':
            issues.append(f'failed row case={case} n={n} rep={r["rep_id"]}')
    checks={
      'direct_mean':sum(float(r['direct_time']) for r in rows)/100,
      'direct_min':min(float(r['direct_time']) for r in rows),
      'direct_max':max(float(r['direct_time']) for r in rows),
      'screen_mean':sum(float(r['screen_time']) for r in rows)/100,
      'screen_min':min(float(r['screen_time']) for r in rows),
      'screen_max':max(float(r['screen_time']) for r in rows),
      'max_objective_gap':max(float(r['objective_gap']) for r in rows),
      'max_kkt_residual':max(float(r['kkt_residual']) for r in rows),
      'max_fit_discrepancy':max(float(r['fit_discrepancy']) for r in rows)}
    for key,value in checks.items():
        if not math.isclose(float(s[key]),value,rel_tol=1e-12,abs_tol=1e-14):issues.append(f'{key} mismatch case={case} n={n}')
result={'consistent':not issues,'raw_replications':raw_count,'summary_rows':len(summary),'issues':issues,
        'scope':'Independent aggregation check of saved replication metrics; does not rerun fits.'}
print(json.dumps(result,indent=2))
if issues:raise SystemExit(1)
