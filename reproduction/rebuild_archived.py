#!/usr/bin/env python3
"""Rebuild selected paper summaries from immutable configuration-level archives.
This is an aggregation reproduction, not a rerun of fitted datasets.
"""
import csv,sys,json,statistics
from pathlib import Path
root=Path(__file__).resolve().parent
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
def read(name):
    with (root/'archive'/name).open() as f:return list(csv.DictReader(f))
a=[r for r in read('original_summary.csv') if int(r['paper_case']) in (2,3,4)]
a += [dict(r,paper_case='1') for r in read('logistic_summary.csv')]
a.sort(key=lambda r:(int(r['paper_case']),float(r['tau']),int(r['n']),r['method']))
assert len(a)==144 and len({(r['paper_case'],r['tau'],r['n'],r['method']) for r in a})==144
assert all(int(r['n_present'])==500 for r in a)
def write(name,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (out/name).open('w') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)
write('mapped_summary.csv',a)
write('runtime.csv',[{k:r[k] for k in ('paper_case','tau','n','method','time_mean_sec','time_min_sec','time_max_sec')} for r in a])
write('ablation.csv',[dict(paper_case=c,**{m:statistics.mean(float(r['time_mean_sec']) for r in a if int(r['paper_case'])==c and r['method']==m) for m in ('lean_seq','unified_u11')}) for c in range(1,5)])
# Main screening tables use logistic Case 1 and original TVCQR Case 3.
names=dict(gamma='gamma_mean',first_mean='first_pass_prop_mean',first_min='first_pass_prop_min',first_max='first_pass_prop_max',F_mean='Fr_mean',F_median='Fr_median',F_q90='Fr_q90',F_max='Fr_max',R_mean='total_repair_mean',R_min='total_repair_min',R_max='total_repair_max',all_grid_no_repair='uniform',mean_max_S='mean_max_Sj',over_nb='max_Sj_over_nb',over_nb_gamma_log='max_Sj_over_nb_gamma_logn')
screen=[dict(r,paper_case='1') for r in read('logistic_screening.csv') if float(r['tau'])==.5]
screen += [dict(paper_case='3',n=r['n'],tau='.5',**{k:r[v] for k,v in names.items()}) for r in read('original_screening.csv') if r['paper_case']=='3']
write('screening_tau05.csv',screen)
u=[r for r in a if r['method']=='unified_u11']
write('screened_diagnostics.csv',u)
# MST is a separate experiment: preserve its archived summaries, not U11 timings.
for name in ('llqr_multivar_iteration_summary_20260406_033218.csv','llqr_multivar_iteration_reduction_20260406_033218.csv'):
    write(name,read(name))
status=dict(configurations=48,method_summaries=144,replications_per_configuration=500,
            full_active_recovery=sum(int(r['internal_recovery_total']) for r in u),
            seq_fallback=sum(int(r['fallback_count']) for r in u),scope='four primary cases; MST archived summaries only')
(out/'verification.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status))
