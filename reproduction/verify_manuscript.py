"""Audit embedded manuscript numbers from CSV inputs without importing the generator.

This checks aggregation, mappings and displayed rounding, not raw simulation fits.
"""
from pathlib import Path
import csv
import json
import math
import re
import statistics
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--report', type=Path, help='Optionally save the full JSON result at this path.')
parser.add_argument('--llqr-run',type=Path,help='Accepted new LLQR paired run; omit to verify historical manuscript.')
parser.add_argument('--bundle', type=Path, required=True)
args = parser.parse_args()
BUNDLE = args.bundle
manifest=BUNDLE/'output/llqr_rq_revision/source_manifest.json'
if args.llqr_run is None and manifest.exists():
    args.llqr_run=BUNDLE/json.loads(manifest.read_text())['llqr']
def read_csv(p):
    with (BUNDLE / p).open() as f:
        return list(csv.DictReader(f))

old = read_csv('results/hpc/u11_rep500_runs/u11_rep500_seed2025_61dc4d6_20260911/tables/config_method_summary.csv')
new = read_csv('results/hpc/llqr_case2_logistic_runs/llqr_c2_logistic_rep500_seed2025_0863a92_20260914_112526/tables/config_method_summary.csv')
# The old LLQR Case 1 is excluded. New logistic LLQR raw Case 2 maps to paper Case 1.
source = [dict(r) for r in old if int(r['paper_case']) in (2, 3, 4)]
source.extend(dict(r, paper_case='1') for r in new)
if args.llqr_run:
    assert (args.llqr_run/'RESULTS_PASS').exists()
    with (args.llqr_run/'tables/config_method_summary.csv').open() as f:
        source=[dict(r) for r in old if int(r['paper_case']) in (3,4)]+[dict(r,paper_case=r['case']) for r in csv.DictReader(f)]
index = {(int(r['paper_case']), float(r['tau']), int(r['n']), r['method']): r for r in source}
main = (BUNDLE / 'paper/main_submission_new_v2_TW.tex').read_text()
supp = (BUNDLE / 'paper/supplement_submission_new_v2_TW.tex').read_text()
def strip_blue(text):
    needle=r'\finalRev{'
    while needle in text:
        i=text.index(needle);j=i+len(needle);depth=1;k=j
        while depth:
            if text[k]=='{':depth+=1
            elif text[k]=='}':depth-=1
            k+=1
        text=text[:i]+text[j:k-1]+text[k:]
    return text
main=strip_blue(main);supp=strip_blue(supp)
checks = []
issues = []

def block(text, name):
    return text.split('% BEGIN ' + name, 1)[1].split('% END ' + name, 1)[0]

def case_rows(text):
    return [line.split('\\\\', 1)[0].split('&') for line in text.splitlines() if re.match(r'^Case [1-4] &', line)]

def check_num(label, printed, actual):
    printed = printed.strip().replace('$', '').replace(r'\%', '')
    sci = re.fullmatch(r'(-?[\d.]+)\\times 10\^\{(-?\d+)\}', printed)
    if sci:
        mant, exp = sci.groups()
        decimals = len(mant.split('.')[1]) if '.' in mant else 0
        value = float(mant) * 10 ** int(exp)
        tol = .5 * 10 ** (int(exp) - decimals)
    else:
        decimals = len(printed.split('.')[1]) if '.' in printed else 0
        value = float(printed)
        tol = .5 * 10 ** -decimals
    passed = abs(value - actual) <= tol + 1e-12 * max(tol, abs(actual))
    record = dict(label=label, printed=printed, actual=actual, tolerance=tol, consistent=passed)
    checks.append(record)
    if not passed:
        issues.append(record)

# Main Table 1, 96 mean/minimum/maximum numbers.
for row in case_rows(block(main, 'V2_RUNTIME_TAU05')):
    case, n = int(row[0].strip().split()[1]), int(row[1])
    for method, cell in zip(('direct_baseline', 'unified_u11'), row[2:]):
        vals = re.fullmatch(r'\s*([\d.]+) \[([\d.]+), ([\d.]+)\]\s*', cell).groups()
        for key, val in zip(('time_mean_sec', 'time_min_sec', 'time_max_sec'), vals):
            check_num(f'main T1 case={case} n={n} {method} {key}', val, float(index[case, .5, n, method][key]))

# Main Tables 2 and 3 draw Case 1 from the new screening summary, Case 3 from the old one.
new_screen = read_csv('output/llqr_case2_logistic_rep500/reproduced/screening_retained.csv')
old_screen = read_csv('paper/data/u11_rep500/rep500_seed2025_screening_tau05.csv')
screen = {(1, int(r['n'])): r for r in new_screen if float(r['tau']) == .5}
names = dict(gamma='gamma_mean', first_mean='first_pass_prop_mean', first_min='first_pass_prop_min', first_max='first_pass_prop_max', F_mean='Fr_mean', F_median='Fr_median', F_q90='Fr_q90', F_max='Fr_max', R_mean='total_repair_mean', R_min='total_repair_min', R_max='total_repair_max', all_grid_no_repair='uniform', mean_max_S='mean_max_Sj', over_nb='max_Sj_over_nb', over_nb_gamma_log='max_Sj_over_nb_gamma_logn')
if args.llqr_run:
    with (args.llqr_run/'tables/screening_all_taus.csv').open() as f:
        screen={(1,int(r['n'])):{k:r[v] for k,v in names.items()} for r in csv.DictReader(f) if r['case']=='1' and float(r['tau'])==.5}
for r in old_screen:
    if r['paper_case'] == '3':
        screen[3, int(r['n'])] = {k: r[v] for k, v in names.items()}
for row in case_rows(block(main, 'U11_REP500_MAIN_SCREENING_VERIFICATION')):
    case, n = int(row[0].strip().split()[1]), int(row[1])
    cells = [row[2]] + row[3].split('/') + row[4].split('/') + row[5].split('/') + [row[6]]
    keys = ('gamma', 'first_mean', 'first_min', 'first_max', 'F_mean', 'F_median', 'F_q90', 'F_max', 'R_mean', 'R_min', 'R_max', 'all_grid_no_repair')
    for val, key in zip(cells, keys):
        check_num(f'main T2 case={case} n={n} {key}', val, float(screen[case, n][key]))
for row in case_rows(block(main, 'U11_REP500_MAIN_SCREENING_SIZE')):
    case, n = int(row[0].strip().split()[1]), int(row[1])
    for val, key in zip(row[2:], ('gamma', 'mean_max_S', 'over_nb', 'over_nb_gamma_log')):
        check_num(f'main T3 case={case} n={n} {key}', val, float(screen[case, n][key]))

# Supplement Table S2, including independent max-minus-min calculation.
case = tau = None
for line in block(supp, 'U11_REP500_SUPPLEMENT_RUNTIME_TABLE').splitlines():
    if not ('direct local fit &' in line or r'\texttt{screen-seq} &' in line):
        continue
    row = line.split('\\\\', 1)[0].split('&')
    if row[0].strip():
        case = int(row[0].strip().split()[1])
    if row[1].strip():
        tau = float(row[1])
    method = 'direct_baseline' if row[2].strip() == 'direct local fit' else 'unified_u11'
    for n, cell in zip((1000, 2000, 5000, 10000), row[3:]):
        mean, range_ = re.fullmatch(r'\s*([\d.]+) \(([\d.]+)\)\s*', cell).groups()
        r = index[case, tau, n, method]
        check_num(f'supp S2 case={case} tau={tau} n={n} {method} mean', mean, float(r['time_mean_sec']))
        check_num(f'supp S2 case={case} tau={tau} n={n} {method} range', range_, float(r['time_max_sec']) - float(r['time_min_sec']))

# Supplement Table S3, case averages use 12 equally weighted configuration means.
ablation = []
for row in case_rows(block(supp, 'U11_REP500_SUPPLEMENT_LEAN_ABLATION')):
    case = int(row[0].strip().split()[1])
    means = {m: statistics.mean(float(r['time_mean_sec']) for r in source if int(r['paper_case']) == case and r['method'] == m) for m in ('lean_seq', 'unified_u11')}
    change = 100 * (means['unified_u11'] / means['lean_seq'] - 1)
    for val, key, actual in zip(row[1:], ('seq_mean', 'screen_mean', 'change_percent'), (means['lean_seq'], means['unified_u11'], change)):
        check_num(f'supp S3 case={case} {key}', val, actual)
    ablation.append(dict(case=case, **means, change_percent=change))

# Supplement Tables S4-S6.
for tau, suffix in ((.2, '02'), (.5, '05'), (.8, '08')):
    for row in case_rows(block(supp, 'U11_REP500_SUPPLEMENT_DISCREPANCY_TAU' + suffix)):
        case = int(row[0].strip().split()[1])
        for n, val in zip((1000, 2000, 5000, 10000), row[1:]):
            check_num(f'supp discrepancy case={case} tau={tau} n={n}', val, float(index[case, tau, n, 'unified_u11']['max_average_relative_discrepancy']))

speedups = [float(index[c, .5, n, 'direct_baseline']['time_mean_sec']) / float(index[c, .5, n, 'unified_u11']['time_mean_sec']) for c in range(1, 5) for n in (1000, 2000, 5000, 10000)]
derived = dict(speedup_min=min(speedups), speedup_max=max(speedups), reduction_min=100 * (1 - 1/min(speedups)), reduction_max=100 * (1 - 1/max(speedups)))
for macro, key in (('RuntimeSpeedupMin','speedup_min'), ('RuntimeSpeedupMax','speedup_max'), ('RuntimeReductionMin','reduction_min'), ('RuntimeReductionMax','reduction_max')):
    val = re.search(r'\\newcommand\{\\' + macro + r'\}\{([^}]+)\}', main).group(1)
    check_num(macro, val, derived[key])
for text_name, text in (('main',main), ('supp',supp)):
    for macro, method, agg in (('DiscrepancyOverallMax', 'unified_u11', max), ('DiscrepancyOverallMin', 'unified_u11', min), ('AblationScreenDiscrepancyMedian', 'unified_u11', statistics.median), ('AblationSeqDiscrepancyMedian', 'lean_seq', statistics.median)):
        match = re.search(r'\\newcommand\{\\' + macro + r'\}\{([^\n]+)\}', text)
        if match:
            actual = agg(float(r['max_average_relative_discrepancy']) for r in source if r['method'] == method)
            check_num(text_name + ' ' + macro, match.group(1), actual)
            derived[macro] = actual

# Check mapped values against source, never use mapped values as primary evidence.
mapped = read_csv('output/'+('llqr_rq_revision' if args.llqr_run else 'llqr_case1_paper_revision')+'/mapped_method_summary.csv')
map_issues = []
for r in mapped:
    s = index[int(r['paper_case']), float(r['tau']), int(r['n']), r['method']]
    for key in ('n_present', 'time_mean_sec', 'time_min_sec', 'time_max_sec', 'max_average_relative_discrepancy'):
        if not math.isclose(float(r[key]), float(s[key]), rel_tol=1e-14, abs_tol=0):
            map_issues.append((r, key))

recovery = [{k:r[k] for k in ('paper_case','tau','n','internal_recovery_total','repaired_points_total','repair_total')} for r in source if r['method']=='unified_u11' and int(r['internal_recovery_total']) > 0]
slower_screen = []
for c in range(1,5):
    for tau in (.2,.5,.8):
        for n in (1000,2000,5000,10000):
            seq = float(index[c,tau,n,'lean_seq']['time_mean_sec'])
            scr = float(index[c,tau,n,'unified_u11']['time_mean_sec'])
            if scr >= seq:
                slower_screen.append(dict(case=c,tau=tau,n=n,seq_mean=seq,screen_mean=scr,change_percent=100*(scr/seq-1)))
result = dict(
    scope='Summary-level arithmetic and displayed-rounding check; no simulation fits or raw replication reconstruction',
    source_rows=len(source), unique_source_keys=len(index), all_current_source_rows_have_500_present=all(int(r['n_present'])==500 for r in source),
    printed_numeric_checks=len(checks), printed_numeric_mismatches=issues,
    mapped_value_checks=len(mapped)*5, mapped_value_mismatches=map_issues,
    checks=checks, derived=derived, ablation=ablation,
    screen_not_faster_than_seq_configurations=slower_screen,
    nonzero_full_active_recoveries=recovery,
    total_full_active_recoveries=sum(int(r['internal_recovery_total']) for r in source if r['method']=='unified_u11'),
    total_seq_fallbacks=sum(int(r['fallback_count']) for r in source if r['method']=='unified_u11'),
    initial_threshold_expansions=sum(int(r['threshold_expansion_total']) for r in source if r['method']=='unified_u11'),
)
if args.report:
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2) + '\n')
failed = bool(issues or map_issues or len(source) != 144 or len(index) != 144
              or len(checks) != 486 or len(mapped) * 5 != 720
              or not result['all_current_source_rows_have_500_present'])
print(json.dumps(dict(consistent=not failed, printed_numeric_checks=len(checks),
                      printed_numeric_mismatches=len(issues), mapped_value_checks=len(mapped)*5,
                      mapped_value_mismatches=len(map_issues), scope=result['scope']), indent=2))
raise SystemExit(1 if failed else 0)
