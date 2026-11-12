# NAMRC N1 FINAL: 36 full Rust simulations + independent energy audits + ZIP.
# Run from an extracted NAMRC_N1_FINAL folder. Keeps raw results for audits/resume.
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import os,sys,subprocess,shutil,json,csv,hashlib,traceback,time

ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)
sys.path.insert(0,str(ROOT))
os.environ['PATH']=str(Path.home()/'.cargo'/'bin')+os.pathsep+os.environ.get('PATH','')
import run_full_solver_N1 as N1

RES=ROOT/'N1_results'; RUN=ROOT/'runs'
RES.mkdir(exist_ok=True); RUN.mkdir(exist_ok=True)
CASES=list(csv.DictReader((ROOT/'experiment_matrix.csv').open(newline='')))
assert len(CASES)==36, f'Expected 36 frozen cases, got {len(CASES)}'
assert len({(c['toolpath'],c['speed_mm_s'],c['dt_s'],c['timing_mode']) for c in CASES})==36
assert (ROOT/'solver'/'Cargo.toml').exists()
assert (ROOT/'toolpaths'/'b01_six_layer_segment.gcode').exists()
assert (ROOT/'toolpaths'/'independent_serpentine_six_layers.gcode').exists()
print('Frozen matrix: 2 paths × 3 speeds × 3 time steps × 2 modes = 36 cases',flush=True)
print('Cargo executable:',shutil.which('cargo'),flush=True)
if not shutil.which('cargo'):raise RuntimeError('Rust/Cargo unavailable in this Jupyter kernel')

# Compile once. Fail clearly rather than generating guessed results.
EXE=N1.build()
print('Rust solver:',EXE,flush=True)

failures=[]; raw_missing=[]
for i,c in enumerate(CASES,1):
    name=f"{c['toolpath']}__v{c['speed_mm_s']}__dt{c['dt_s']}__{c['timing_mode']}"
    rp=RES/(name+'.json')
    rd=RUN/name
    # Reuse a completed case ONLY if its raw histories remain available for the audit.
    if rp.exists() and (rd/'outputfiles'/'elem_temps.csv').exists():
        try:
            cached=json.loads(rp.read_text())
            if cached.get('status')=='complete':
                print(f'[{i}/36] REUSE {name}',flush=True)
                continue
        except Exception:pass
    if rp.exists() and not (rd/'outputfiles'/'elem_temps.csv').exists():
        rp.unlink() # Rerun if older package deleted raw output.
    print(f'[{i}/36] RUN {name}',flush=True)
    try:
        r=N1.execute_case(c,EXE)
        print(f"     {r['n_elements']} elements; mean final T={r['mean_final_C']:.2f} C; {r['execution_wall_s']:.1f} s",flush=True)
    except Exception as exc:
        print('     FAILED:',repr(exc),flush=True)
        failures.append(dict(case=name,stage='simulation',error=str(exc)))
        continue

N1.paired_report()

# Energy audits are intentionally independent of solver integration. Same laws for
# both timing modes. No fan; fixed 25 C, M2 bed, 0.1 K history deadband.
try:
    import numpy as np
    from N1_energy_audit import audit
except Exception as e:
    raise RuntimeError('Install numpy in this Jupyter environment before auditing') from e

P=dict(h=3.0,Tdep=488.15,Ta=298.15,eps=0.95,Tbed=333.15,Rc=0.005)
energy=[]
for i,c in enumerate(CASES,1):
    name=f"{c['toolpath']}__v{c['speed_mm_s']}__dt{c['dt_s']}__{c['timing_mode']}"
    rd=RUN/name; out=rd/'outputfiles'/'elem_temps.csv'
    audit_path=RES/('energy_audit_'+name+'.json')
    if not out.exists():
        failures.append(dict(case=name,stage='audit',error='Missing raw histories'))
        continue
    try:
        if audit_path.exists():
            obj=json.loads(audit_path.read_text()); a=obj['audit']
            print(f'[{i}/36] AUDIT REUSE {name} = {100*a["rel_err_vs_heat_in"]:.3f}%',flush=True)
        else:
            print(f'[{i}/36] AUDIT {name}',flush=True)
            a=audit(str(rd),P)
            obj=dict(case=name,audit=a,storage_threshold_K=0.1,assumptions=P,
                     legacy_energy_audited_against_true_clock=True)
            audit_path.write_text(json.dumps(obj,indent=2))
            print(f'     residual {100*a["rel_err_vs_heat_in"]:.3f}%',flush=True)
        energy.append(dict(toolpath=c['toolpath'],speed_mm_s=int(c['speed_mm_s']),
                           dt_s=float(c['dt_s']),timing_mode=c['timing_mode'],
                           case=name,energy_residual_pct=100*a['rel_err_vs_heat_in'],
                           signed_balance_error_J=a['balance_error_J'],
                           deposited_heat_J=a['heat_in_J'],
                           natural_convection_J=a['natural_conv_J'],
                           radiation_J=a['radiation_J'],
                           bed_net_J=a['bed_net_J'],
                           stored_change_J=a['stored_change_J'],
                           audit_t_end_s=a['t_end_s']))
    except Exception as exc:
        print('     AUDIT FAILED:',repr(exc),flush=True)
        failures.append(dict(case=name,stage='audit',error=str(exc)))

if energy:
    with (RES/'N1_energy_audit_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(energy[0])); w.writeheader();w.writerows(energy)

# Bring numerical and energy results into one manuscript-focused paired summary.
paired=[]
for c in CASES:
    if c['timing_mode']!='corrected':continue
    other={**c,'timing_mode':'legacy_timing'}
    def nm(q):return f"{q['toolpath']}__v{q['speed_mm_s']}__dt{q['dt_s']}__{q['timing_mode']}"
    cp=RES/(nm(c)+'.json'); lp=RES/(nm(other)+'.json')
    ca=RES/('energy_audit_'+nm(c)+'.json');la=RES/('energy_audit_'+nm(other)+'.json')
    if not all(p.exists() for p in (cp,lp,ca,la)):continue
    cr=json.loads(cp.read_text());lr=json.loads(lp.read_text())
    ac=json.loads(ca.read_text())['audit']; al=json.loads(la.read_text())['audit']
    paired.append(dict(toolpath=c['toolpath'],speed_mm_s=int(c['speed_mm_s']),dt_s=float(c['dt_s']),
                       corrected_elements=cr['n_elements'],legacy_elements=lr['n_elements'],
                       corrected_final_mean_C=cr['mean_final_C'],legacy_final_mean_C=lr['mean_final_C'],
                       legacy_minus_corrected_mean_C=lr['mean_final_C']-cr['mean_final_C'],
                       corrected_energy_residual_pct=100*ac['rel_err_vs_heat_in'],
                       legacy_energy_residual_pct=100*al['rel_err_vs_heat_in'],
                       corrected_runtime_s=cr['execution_wall_s'],legacy_runtime_s=lr['execution_wall_s']))
if paired:
    with (RES/'N1_final_paired_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(paired[0]));w.writeheader();w.writerows(paired)

(RES/'failures_final.json').write_text(json.dumps(failures,indent=2))
status=dict(expected_cases=36,completed_simulations=sum((RES/(f"{c['toolpath']}__v{c['speed_mm_s']}__dt{c['dt_s']}__{c['timing_mode']}"+'.json')).exists() for c in CASES),
            completed_energy_audits=len(energy),complete_pairs=len(paired),failures=len(failures),
            source_commit=(ROOT/'R8_4_source_commit.txt').read_text().strip(),
            note='All energy audits use fixed baseline physical laws and 0.1 K output reconstruction; not physical FLIR validation.')
(RES/'N1_run_status_final.json').write_text(json.dumps(status,indent=2))

# Manuscript review ZIP: avoid enormous raw histories but keep logs and case settings.
zip_path=ROOT/'N1_results_for_review.zip'
with ZipFile(zip_path,'w',ZIP_DEFLATED,compresslevel=6) as z:
    for p in [ROOT/'experiment_matrix.csv',ROOT/'toolpath_audit.json',ROOT/'R8_4_source_commit.txt',
              ROOT/'path_level_clock_diagnostics.csv',ROOT/'run_full_solver_N1.py',
              ROOT/'N1_energy_audit.py']:
        if p.exists():z.write(p,p.relative_to(ROOT))
    for p in RES.rglob('*'):
        if p.is_file() and p.suffix in ('.csv','.json','.txt'):
            z.write(p,p.relative_to(ROOT))
    for d in RUN.iterdir():
        if not d.is_dir():continue
        for rel in ['inputfiles/Input_file.txt','solver_stdout.log']:
            p=d/rel
            if p.exists():z.write(p,p.relative_to(ROOT))

print('\n',json.dumps(status,indent=2),flush=True)
print('Results ready:',zip_path,flush=True)
if failures:print('Some cases failed: inspect N1_results/failures_final.json; do NOT claim 36 complete.',flush=True)
