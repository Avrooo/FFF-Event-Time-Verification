#!/usr/bin/env python3
"""Verify archived numerical summaries against manuscript Tables 2 and 3.
No thermal simulation is executed, and no external Python package is required.
"""
from pathlib import Path
import csv, json, math
ROOT=Path(__file__).resolve().parents[1]
def rows(path):
 with (ROOT/path).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def near(a,b,precision=4):
 # Decimal rendering in manuscript uses the given number of decimal places.
 if round(float(a),precision)!=round(float(b),precision):
  raise AssertionError(f'{a} != manuscript {b} to {precision} decimal places')
n0=rows('N0/N0_results/N0_full_print_summary.csv')
assert len(n0)==8, f'N0 entries: {len(n0)}'
expect={'A0':{0.25:2.3664,0.1:0.9008},'A2':{0.25:55.8542,0.1:57.0967},'A6':{0.25:48.1436,0.1:52.5634},'R2':{0.25:8.3669,0.1:3.3386}}
for r in n0:
 case=r['case'];thr=float(r['storage_threshold_K']);assert case in expect and thr in expect[case]
 near(r['abs_residual_pct_heat_input'],expect[case][thr]);assert int(r['n_elements'])==51204
 near(r['physical_deposition_s'],929.96,2)
assert len({(r['case'],r['storage_threshold_K']) for r in n0})==8
print('PASS: N0 8 cases, Table 2 values, 51,204 cells and 929.96 s command clock')
a=rows('N1/N1_results/N1_energy_audit_summary.csv')
p=rows('N1/N1_results/paired_full_solver_results.csv')
assert len(a)==36 and len(p)==18, (len(a),len(p))
by={}
for r in a:
 key=(r['toolpath'],float(r['speed_mm_s']),float(r['dt_s']))
 by.setdefault(key,{})[r['timing_mode']]=float(r['energy_residual_pct'])
assert len(by)==18
for k,v in by.items():assert set(v)=={'corrected','legacy_timing'} and v['corrected']<v['legacy_timing'],k
# Table 3 displays the 40 mm/s, 0.050 s reference values, not the 9-setting medians.
ref={('b01_six_layer_segment','corrected'):0.6002,('b01_six_layer_segment','legacy_timing'):54.3119,('independent_serpentine_six_layers','corrected'):0.0447,('independent_serpentine_six_layers','legacy_timing'):42.4282}
for (path,mode),expected in ref.items():
 value=next(float(r['energy_residual_pct']) for r in a if r['toolpath']==path and r['timing_mode']==mode and float(r['speed_mm_s'])==40 and float(r['dt_s'])==0.05)
 near(value,expected,4)
# Verify the manuscript's Table 3 median/ranges, allowing its display rounding.
from statistics import median
for path,correct_med,legacy_med in [('b01_six_layer_segment',0.600,54.31),('independent_serpentine_six_layers',0.045,42.43)]:
 c=[float(r['energy_residual_pct']) for r in a if r['toolpath']==path and r['timing_mode']=='corrected']
 l=[float(r['energy_residual_pct']) for r in a if r['toolpath']==path and r['timing_mode']=='legacy_timing']
 near(median(c),correct_med,3);near(median(l),legacy_med,2)
print('PASS: N1 36 cases, 18 pairs, four Section 4.4 values, Table 3 medians')
# Manuscript conversion A0 0.10 K, check exact algebra rather than only rounded values.
with (ROOT/'N0/N0_results/A0_md0p1_audit.json').open() as f:q=json.load(f)['audit']
h=q['heat_in_J']; delta=q['stored_change_J']; bed=q['bed_net_J']; conv=q['natural_conv_J']; rad=q['radiation_J']
new_residual=h-bed-conv-rad-(h+delta)
assert math.isclose(new_residual,-q['balance_error_J'],abs_tol=1e-7)
assert math.isclose(100*abs(new_residual)/h,0.9007993788850952,abs_tol=1e-7)
print('PASS: manuscript-vs-archived bed/residual sign conversion')
print('Archived-result checks PASSED. No numerical simulations were rerun.')
