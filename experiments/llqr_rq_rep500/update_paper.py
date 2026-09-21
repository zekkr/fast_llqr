#!/usr/bin/env python3
"""Update manuscript numeric blocks from accepted full-weight rq LLQR and archived TVCQR.

By default, update only the supplement to protect the independently edited main.
Use --check for read-only verification; --target main or both explicitly opts in
to processing the main manuscript.
Original run case IDs and all input files remain unchanged.
"""
import argparse, csv, hashlib, json, math, re, statistics
from pathlib import Path
from decimal import Decimal, ROUND_HALF_UP

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bundle',type=Path,required=True)
parser.add_argument('--run',type=Path,required=True)
parser.add_argument('--check',action='store_true')
parser.add_argument('--target',choices=('both','main','supplement'),default='both')
args=parser.parse_args()
ROOT=args.bundle.resolve()
PAPER=ROOT/'paper'
OUT=ROOT/'output/llqr_rq_revision'
OLD_TAG='u11_rep500_seed2025_61dc4d6_20260911'
OLD=ROOT/'results/hpc/u11_rep500_runs'/OLD_TAG
NEW=args.run.resolve();NEW_TAG=NEW.name
assert (NEW/'RESULTS_PASS').exists(), 'Formal rerun has not passed'
NS = (1000, 2000, 5000, 10000)
TAUS = (.2, .5, .8)
SOURCES = []

def read(path):
    SOURCES.append(path)
    with path.open() as f:
        return list(csv.DictReader(f))

def num(r, k):
    return float(r[k])

def sig(x):
    return format(x, f'.{max(0,3-math.floor(math.log10(abs(x))))}f') if x else '0.000'

def half_up(value, places=2):
    """Decimal half-up rounding; keep trailing zeros at the requested precision."""
    return format(Decimal(str(value)).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP), 'f')

def runtime_time(value):
    """Two decimals at >=0.1; two significant digits for positive times <0.1."""
    x = Decimal(str(value))
    assert x >= 0
    if x == 0 or x >= Decimal('0.1'):
        return half_up(x)
    rounded = x.quantize(Decimal(1).scaleb(x.adjusted()-1), rounding=ROUND_HALF_UP)
    # A carry may change the decimal position needed for two significant digits.
    return half_up(x, 1-rounded.adjusted())

def old_time(x):
    return f'{x:.4f}' if abs(x) < 1 else sig(x)

def sci(x):
    a, b = f'{x:.2e}'.split('e')
    return a + r'\times 10^{' + str(int(b)) + '}'

def replace_block(text, name, fn):
    pattern = r'(% BEGIN '+name+r'\n)(.*?)(% END '+name+r')'
    result, count = re.subn(pattern, lambda m: m[1]+fn(m[2])+m[3], text, flags=re.S)
    assert count == 1, name
    return result

rows = {}
for base, tag in ((OLD, OLD_TAG), (NEW, NEW_TAG)):
    for r in read(base/'tables/config_method_summary.csv'):
        cid = int(r['case']) + (2 if r['model']=='tvcqr' else 0)
        if tag == OLD_TAG and cid in (1,2):
            continue
        if tag == NEW_TAG:
            assert r['model']=='llqr' and int(r['case']) in (1,2)
            cid = int(r['case'])
        key = (cid, float(r['tau']), int(r['n']), r['method'])
        assert key not in rows
        assert int(r['n_present']) == int(r['n_solver_ok']) == int(r['n_accepted']) == 500
        assert all(num(r,k)==0 for k in ('n_nonfinite','ierr_nonzero','failed_eval_nonzero','h_mismatch_count','fallback_count'))
        assert all(r[k]=='TRUE' for k in ('solver_stable','path_accepted','timing_order_balanced','discrepancy_valid'))
        rows[key] = dict(r, source_run=tag, source_case=r['case'], paper_case=str(cid))
assert len(rows)==144
assert set(rows)=={(c,t,n,m) for c in range(1,5) for t in TAUS for n in NS for m in ('direct_baseline','lean_seq','unified_u11')}
# Old TVCQR keeps its own verified sources. New LLQR must pass its independent gate.
for r in read(OLD/'tables/integrity.csv'):
    if r['model'] != 'tvcqr': continue
    assert r['complete']=='TRUE' and int(r['n_complete'])==500
for r in read(NEW/'tables/integrity.csv'):
    assert r['complete']=='TRUE' and int(r['n_complete'])==500 and int(r['n_missing'])==int(r['n_malformed'])==0
regen=read(NEW/'tables/baseline_regeneration_attempts.csv')
regenerated_reps=len({(r['case'],r['tau'],r['n'],r['rep_id']) for r in regen})
names=dict(gamma='gamma_mean', first_mean='first_pass_prop_mean',first_min='first_pass_prop_min',first_max='first_pass_prop_max',F_mean='Fr_mean',F_median='Fr_median',F_q90='Fr_q90',F_max='Fr_max',R_mean='total_repair_mean',R_min='total_repair_min',R_max='total_repair_max',all_grid_no_repair='uniform',mean_max_S='mean_max_Sj',over_nb='max_Sj_over_nb',over_nb_gamma_log='max_Sj_over_nb_gamma_logn')
new_screen={int(r['n']):{k:r[v] for k,v in names.items()} for r in read(NEW/'tables/screening_all_taus.csv') if int(r['case'])==1 and float(r['tau'])==.5}
old_screen={int(r['n']):r for r in read(PAPER/'data/u11_rep500/rep500_seed2025_screening_tau05.csv') if int(r['paper_case'])==3}
assert all(int(r['threshold_expansion_total'])==0 for r in rows.values() if r['method']=='unified_u11'), 'Review threshold-expansion prose before publication'

def row(c,t,n,m): return rows[c,t,n,m]
def tm(c,t,n,m): return num(row(c,t,n,m),'time_mean_sec')

speeds = [tm(c,.5,n,'direct_baseline')/tm(c,.5,n,'unified_u11') for c in range(1,5) for n in NS]
dis = [num(r,'max_average_relative_discrepancy') for k,r in rows.items() if k[3]=='unified_u11']
seq_dis = [num(r,'max_average_relative_discrepancy') for k,r in rows.items() if k[3]=='lean_seq']
macros = {'AblationScreenDiscrepancyMedian':sci(statistics.median(dis)), 'AblationSeqDiscrepancyMedian':sci(statistics.median(seq_dis)), 'DiscrepancyOverallMin':sci(min(dis)), 'DiscrepancyOverallMax':sci(max(dis)),
          'RuntimeReductionMin':f'{100*(1-1/min(speeds)):.1f}', 'RuntimeReductionMax':f'{100*(1-1/max(speeds)):.1f}',
          'RuntimeSpeedupMin':f'{min(speeds):.1f}', 'RuntimeSpeedupMax':f'{max(speeds):.1f}'}

runtime = [r'\begin{table}[htbp]',r'\centering',r'\small',r'\setlength{\tabcolsep}{4pt}',r'\renewcommand{\arraystretch}{1.05}',r'\begin{tabular}{lrrr}',r'\toprule',r'Case & $n$ & direct local fit & \texttt{screen-seq}\\',r'\midrule']
for c in range(1,5):
    for n in NS:
        cells=[]
        for m in ('direct_baseline','unified_u11'):
            z=row(c,.5,n,m)
            cells.append(f"{runtime_time(z['time_mean_sec'])} [{runtime_time(z['time_min_sec'])}, {runtime_time(z['time_max_sec'])}]")
        runtime.append(f'Case {c} & {n} & '+ ' & '.join(cells)+r'\\')
    if c<4: runtime.append(r'\midrule')
runtime += [r'\bottomrule',r'\end{tabular}','\\caption[Computation times at the median quantile.]{Computation times at $\\tau=0.5$ over 500 replications. Entries give the mean [minimum, maximum] in seconds; brackets show the observed range, not a confidence interval.}',r'\label{tab:runtime_tau05_main}',r'\end{table}','']


def screening_body(body):
    lines=[]
    for line in body.splitlines():
        if not re.match(r'Case [123] &',line):
            lines.append(line);continue
        c=int(line[5]); n=int(line.split('&')[1])
        if c in (1,2):
            c=1;z=new_screen[n]
            rate=[num(z,k) for k in ('first_mean','first_min','first_max')]
            fs=[num(z,k) for k in ('F_mean','F_median','F_q90','F_max')]
            rs=[num(z,k) for k in ('R_mean','R_min','R_max')]
            gamma=num(z,'gamma'); uniform=num(z,'all_grid_no_repair')
        else:
            z=old_screen[n]
            rate=[num(z,k) for k in ('first_pass_prop_mean','first_pass_prop_min','first_pass_prop_max')]
            fs=[num(z,k) for k in ('Fr_mean','Fr_median','Fr_q90','Fr_max')]
            rs=[num(z,k) for k in ('total_repair_mean','total_repair_min','total_repair_max')]
            gamma=num(z,'gamma_mean');uniform=num(z,'uniform')
        ratecell='/'.join(half_up(x) for x in rate)
        # Preserve the old Case 3 count rounding; Case 1 needs more precision near zero.
        meanfmt=lambda x:f'{x:.3f}' if c==1 else f'{x:.0f}'
        fcell=meanfmt(fs[0])+'/'+ '/'.join((f'{x:g}' if c==1 else f'{x:.0f}') for x in fs[1:])
        rcell=meanfmt(rs[0])+'/'+ '/'.join(f'{x:g}' for x in rs[1:])
        lines.append(f'Case {c} & {n} & {half_up(gamma)} & {ratecell} & {fcell} & {rcell} & {half_up(uniform)}'+r'\\')
    return '\n'.join(lines)+'\n'

def retained_body(body):
    def change(m):
        n=int(m[1]);r=new_screen[n]
        return f"Case 1 & {n} & {num(r,'gamma'):.3f} & {num(r,'mean_max_S'):.0f} & {num(r,'over_nb'):.3f} & {num(r,'over_nb_gamma_log'):.2f}"+r'\\'
    return re.sub(r'Case [12] & (\d+) & [^\n]+',change,body)


def supp_runtime(body):
    start=body.index('Case 1 &');end=body.index('Case 3 &')
    lines=[]
    for c in (1,2):
        for t in TAUS:
            for m in ('direct_baseline','unified_u11'):
                cells=[]
                for n in NS:
                    r=row(c,t,n,m)
                    cells.append(f"{old_time(num(r,'time_mean_sec'))} ({old_time(num(r,'time_max_sec')-num(r,'time_min_sec'))})")
                prefix=(f'Case {c}' if t==.2 else '')+f' & {t:.1f} & direct local fit' if m=='direct_baseline' else r' &  & \texttt{screen-seq}'
                lines.append(prefix+' & '+' & '.join(cells)+(r'\\*' if m=='direct_baseline' else r'\\'))
            lines.append(r'\midrule')
    return body[:start]+'\n'.join(lines)+'\n'+body[end:]

case_means=[]
for c in range(1,5):
    seq=statistics.mean(tm(c,t,n,'lean_seq') for t in TAUS for n in NS)
    screen=statistics.mean(tm(c,t,n,'unified_u11') for t in TAUS for n in NS)
    case_means.append(dict(case=c,seq=seq,screen=screen,change=100*(screen/seq-1)))


def unblue(text):
    # Remove only the independently added finalRev macro, respecting nested braces.
    needle=r'\finalRev{'
    while needle in text:
        i=text.index(needle);j=i+len(needle);level=1;k=j
        while level:
            if text[k]=='{':level+=1
            elif text[k]=='}':level-=1
            k+=1
        text=text[:i]+text[j:k-1]+text[k:]
    return text

def blue_cells(old,new):
    oldrows=[l for l in old.splitlines() if '&' in l and r'\\' in l]
    nr=iter(oldrows);lines=[]
    for line in new.splitlines():
        if '&' in line and r'\\' in line:
            previous=next(nr)
            cells=line.split('&');oldcells=previous.split('&');assert len(cells)==len(oldcells)
            for k,(a,b) in enumerate(zip(oldcells,cells)):
                if unblue(a).strip()==b.strip():cells[k]=a
                elif k>0:
                    ending=r'\\*' if b.rstrip().endswith(r'\\*') else r'\\' if b.rstrip().endswith(r'\\') else ''
                    value=b.strip()[:-len(ending)].strip() if ending else b.strip()
                    cells[k]=' '+r'\finalRev{'+value+'}'+ending
            line='&'.join(cells)
        lines.append(line)
    return '\n'.join(lines)+'\n'

def update_block(text,name,fn):
    return replace_block(text,name,lambda old:blue_cells(old,fn(unblue(old))))

changed=[]
targets=('main','supplement') if args.target=='both' else (args.target,)
for target in targets:
    filename=target+'_submission_new_v2_TW.tex'
    path=PAPER/filename;before=path.read_text();text=before
    for name,value in macros.items():
        text=re.sub(r'\\newcommand\{\\'+name+r'\}\{[^\n]+\}',lambda m: m[0] if unblue(m[0])=='\\newcommand{\\'+name+'}{'+value+'}' else '\\newcommand{\\'+name+'}{\\finalRev{'+value+'}}',text)
    text=text.replace('% Auto-generated by paper/scripts/build_u11_rep500_assets.R.','% Generated by paper/scripts/build_v2_logistic_case1_assets.py.')
    if filename.startswith('main'):
        text=update_block(text,'V2_RUNTIME_TAU05',lambda body: re.sub(r'Case [1-4] & [^\n]+',lambda m: next(l for l in runtime if l.startswith(' & '.join(m[0].split(' & ')[:2])+' & ')),body))
        text=update_block(text,'U11_REP500_MAIN_SCREENING_VERIFICATION',screening_body)
        text=update_block(text,'U11_REP500_MAIN_SCREENING_SIZE',retained_body)
    else:
        text=update_block(text,'U11_REP500_SUPPLEMENT_RUNTIME_TABLE',supp_runtime)
        def ablation(body):
            for c in (1,2):
                r=case_means[c-1]
                body=re.sub(r'Case '+str(c)+r' & [^\n]+',lambda _:f"Case {c} & {r['seq']:.3f} & {r['screen']:.3f} & ${r['change']:.1f}\\%$"+r'\\',body)
            return body
        text=update_block(text,'U11_REP500_SUPPLEMENT_LEAN_ABLATION',ablation)
        for t in TAUS:
            def discrepancy(body,t=t):
                for c in (1,2):
                    cells=['$'+sci(num(row(c,t,n,'unified_u11'),'max_average_relative_discrepancy'))+'$' for n in NS]
                    body=re.sub(r'Case '+str(c)+r' & [^\n]+',lambda _: f'Case {c} & '+' & '.join(cells)+r'\\',body)
                return body
            text=update_block(text,'U11_REP500_SUPPLEMENT_DISCREPANCY_TAU'+str(round(t*10)).zfill(2),discrepancy)
    if target=='supplement':
        text=text.replace(r'R 4.3.1 and GNU Fortran 12.2.0.}',r'R 4.3.1 and GNU Fortran 12.2.0. The LLQR simulations used \texttt{quantreg} 5.94.}')
        old_direct=r'\finalRev{Direct LLQR fitting uses \texttt{quantreg::rq.wfit} with \texttt{method="br"}, after removing zero-weight observations. Direct TVCQR fitting uses \texttt{quantreg::rq} with its default method, without explicitly removing zero-weight observations. Both calls use the default stopping tolerances.}'
        text=text.replace(old_direct,r'\finalRev{Direct LLQR and TVCQR fitting use the formula interface of \texttt{quantreg::rq}, with all observations and their local kernel weights, the Barrodale--Roberts method, and default stopping tolerances.}')
        if regenerated_reps:
            text=text.replace('No direct-fitting errors or regeneration occurred in the 48 reported configurations; all final fits had finite output.',r'\finalRev{Direct-fitting errors required dataset regeneration in '+str(regenerated_reps)+' LLQR replications; all final fits had finite output.}')
    if target=='main':
        for old,new in ((r'$0.092$ to $0.041$ in Case~1', f"${float(new_screen[1000]['over_nb']):.3f}$ to ${float(new_screen[10000]['over_nb']):.3f}$ in Case~1"),('from 23 to 66 over the reported sample sizes',f"from {float(new_screen[1000]['mean_max_S']):.0f} to {float(new_screen[10000]['mean_max_S']):.0f} over the reported sample sizes")):
            if old!=new:text=text.replace(old,r'\finalRev{'+new+'}')
        # The two percentages outside tables are also generated from the new Case 1 data.
        start=100*float(new_screen[1000]['all_grid_no_repair']);end=100*float(new_screen[10000]['all_grid_no_repair'])
        if (start,end)!=(8.4,99.6):
            text=text.replace(r'$8.4\%$ to $99.6\%$',r'\finalRev{$'+f'{start:.1f}'+r'\%$ to $'+f'{end:.1f}'+r'\%$}')

    if text!=before:
        changed.append(filename)
        if not args.check:path.write_text(text)
if args.check:
    assert not changed, 'Stale generated manuscript blocks: '+str(changed)
else:
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'source_manifest.json').write_text(json.dumps(dict(llqr=str(NEW.relative_to(ROOT)),tvcqr=str(OLD.relative_to(ROOT)),llqr_cases=[1,2],tvcqr_cases=[3,4],llqr_sha=(NEW/'RESULTS_PASS').read_text().strip()),indent=2)+'\n')
    fields=['paper_case','source_run','source_case','model','tau','n','method','n_present','time_mean_sec','time_min_sec','time_max_sec','max_average_relative_discrepancy']
    with (OUT/'mapped_method_summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows[k] for k in sorted(rows))
    (OUT/'data_validation.json').write_text(json.dumps(dict(configurations=48,method_rows=144,replications_per_configuration=500,regenerated_replications=regenerated_reps,macros=macros,case_ablation=case_means,screen_discrepancy_median=statistics.median(dis),seq_discrepancy_median=statistics.median(seq_dis),sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCES}),indent=2)+'\n')
print('PASS: 48 configurations / 144 method summaries; regenerated LLQR replications:',regenerated_reps)
print('PASS: generated v2 numeric blocks '+('match sources.' if args.check else 'updated.'))
print('Discrepancy medians:',statistics.median(dis),statistics.median(seq_dis))
