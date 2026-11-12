"""N0 full-B01 energy audit for R8.4, portable Windows/Python 3.10+.
Run 'python n0_runner.py --help'.  No data invented; never deletes a finished case.
Requires Rust/Cargo to compile actual solver and numpy for audit.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np
from n0_energy_audit import audit

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'solver'
GCODE=ROOT/'b01_solver.gcode'
RUNS=ROOT/'N0_runs'
RESULTS=ROOT/'N0_results'
EXPERIMENTS={
    'A0':{},
    'A2':{'legacy_timing':1},
    'A6':{'specific_heat_capacity_file_temp_step':5,'legacy_cp_step':1,'legacy_timing':1,'legacy_bed_flag':1,'legacy_mesh':1,'bed_contact_model':0,'convection_model':1},
    'R2':{'specific_heat_capacity_file_temp_step':5,'legacy_cp_step':1,'legacy_bed_flag':1,'legacy_mesh':1,'bed_contact_model':0,'convection_model':1},
}
BASE={
'Number_of_threads':1,'GCode_file':'inputfiles/b01_solver.gcode',
'divisions_per_bead_width':1,'divisions_per_bead_height':1,
'beadwidth':0.55e-3,'beadheight':0.2e-3,'model_type':2,
'specific_heat_capacity_file':'inputfiles/sp_heat_cap_data.csv',
'specific_heat_capacity_file_temp_step':1,
'kx':0.35,'ky':0.35,'kz':0.35,'density':1271.1,
'convective_heat_transfer_film_coefficient':3.0,
'bed_temperature':333.15,'conductivity_at_bed_surface_interface':0.17,
'ambient_temperature':298.15,'extrusion_temperature':488.15,
'emissivity':0.95,'time_step':0.05,'cooldown_period':60.0,
'nodes_file_output':'outputfiles/nodefile.csv',
'element_file_output':'outputfiles/elementfile.csv',
'activation_times_output':'outputfiles/activation_times.csv',
'abaqus_input_mesh_file':'outputfiles/abaqus_input.inp',
'element_temperature_file_output':'outputfiles/elem_temps.csv',
'element_min_temp_diff_store':0.25,
'nodal_temperature_file_output':'outputfiles/node_temps.csv',
'turn_off_layers_at':1000,
'bed_contact_model':2,'bed_contact_resistance':0.005,'convection_model':0,
}

def shasum(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def check():
 for f in [GCODE,SOURCE/'Cargo.toml',SOURCE/'src/main.rs',SOURCE/'inputfiles/sp_heat_cap_data.csv',SOURCE/'inputfiles/conductivity.csv',SOURCE/'inputfiles/Input_file.txt']:
  if not f.is_file(): raise FileNotFoundError(f)
 sample=(SOURCE/'src/legacy.rs').read_text()
 if 'timing' not in sample: raise RuntimeError('Solver source lacks R8.4 legacy timing switch')
 print('G-code SHA256',shasum(GCODE))
 print('Solver source available, original G-code adapter available')

def binary():
 name='DES_thermal_simulation.exe' if sys.platform.startswith('win') else 'DES_thermal_simulation'
 return SOURCE/'target/release'/name

def compile_solver():
 import shutil as sh
 if not sh.which('cargo'):
  raise RuntimeError('Rust Cargo not found. Install Rust from https://rustup.rs/, reopen Anaconda Prompt, verify `cargo --version`.')
 print('Compiling R8.4 solver...')
 subprocess.run(['cargo','build','--release','--locked'],cwd=SOURCE,check=True)
 if not binary().exists():raise RuntimeError(f'Compiled executable missing: {binary()}')
 print('Compiled:',binary())

def prepare(case,md):
 name=f'{case}_md{str(md).replace(".","p")}'
 d=RUNS/name; d.mkdir(parents=True,exist_ok=True)
 for sub in ['inputfiles','outputfiles']: (d/sub).mkdir(exist_ok=True)
 cfg=BASE|EXPERIMENTS[case]|{'element_min_temp_diff_store':md}
 for f in ['sp_heat_cap_data.csv','conductivity.csv']:shutil.copy2(SOURCE/'inputfiles'/f,d/'inputfiles'/f)
 shutil.copy2(GCODE,d/'inputfiles/b01_solver.gcode')
 # Preserve documented input-key order, add only keys missing from template (legacy switches).
 template=(SOURCE/'inputfiles/Input_file.txt').read_text().splitlines()
 result=[]; seen=set()
 for line in template:
  if ':' in line:
   key=line.split(':',1)[0].strip()
   if key in cfg: result.append(f'{key}: {cfg[key]} //');seen.add(key);continue
  result.append(line)
 for key,val in cfg.items():
  if key not in seen:result.append(f'{key}: {val} //')
 (d/'inputfiles/Input_file.txt').write_text('\n'.join(result)+'\n')
 meta={'case':case,'storage_threshold_K':md,'parameters':cfg,'gcode_sha256':shasum(GCODE),'solver_source':str(SOURCE)}
 (d/'case_manifest.json').write_text(json.dumps(meta,indent=2))
 return d

def run_case(d):
 record=d/'run_status.json'
 if record.exists():
  state=json.loads(record.read_text())
  if state.get('returncode')==0 and (d/'outputfiles/elem_temps.csv').exists():
   print('SKIP existing success:',d.name);return state
 if not binary().exists():raise FileNotFoundError('Build the Rust solver first (N0 notebook: compile cell)')
 print('Running',d.name,flush=True)
 tic=time.perf_counter()
 with open(d/'solver_stdout.log','w',encoding='utf-8') as log:
  proc=subprocess.run([str(binary().resolve())],cwd=d,stdout=log,stderr=subprocess.STDOUT)
 state={'returncode':proc.returncode,'wall_s':round(time.perf_counter()-tic,3),'executable':str(binary())}
 record.write_text(json.dumps(state,indent=2))
 if proc.returncode:raise RuntimeError(f'{d.name} failed. See {d/"solver_stdout.log"}')
 return state

def clock_summary(d,case):
 # Independent source-based reconstruction of integrated physics intervals.
 # The G-code physical activation schedule is common to all four configurations.
 a=np.atleast_2d(np.loadtxt(d/'outputfiles/activation_times.csv',delimiter=','))
 ta=a[:,0];diff=np.diff(ta);dt=BASE['time_step'];cool=BASE['cooldown_period']
 if not np.all(diff>=-1e-9):raise ValueError('Activation times not monotonic')
 # Per source: first interval integrated by run_model(time, dt), subsequent positive intervals
 # integrated as floor(time/dt)*dt plus either remainder (corrected) or a FULL dt (legacy).
 good=diff[1:]>1e-4; rem=diff[1:][good]
 time_correct=float(diff[0]+rem.sum())
 step_count=np.floor(rem/dt).astype('int64')
 time_legacy=float(diff[0]+np.sum((step_count+1)*dt))
 # Original D2 cooldown uses stale FIRST interval, not dt (from source model_orthotropic_td_shc_variable_h.rs).
 cooldown_steps=int(cool/dt);cool_legacy=cooldown_steps*float(diff[0]);cool_correct=cooldown_steps*dt
 islegacy=bool(EXPERIMENTS[case].get('legacy_timing',0))
 return {'activation_count':len(ta),'first_activation_s':float(ta[0]),'last_activation_s':float(ta[-1]),
         'physical_deposition_elapsed_s':float(ta[-1]-ta[0]),
         'reconstructed_deposition_integrated_corrected_s':time_correct,
         'reconstructed_deposition_integrated_legacy_s':time_legacy,
         'reconstructed_exposure_ratio_legacy_vs_corrected':time_legacy/time_correct if time_correct else None,
         'requested_cooldown_s':cool,'legacy_cooldown_integrated_s':cool_legacy,
         'corrected_cooldown_integrated_s':cool_correct,
         'case_uses_legacy_timing':islegacy,
         'warning':'Reconstructed from activation schedule and reviewed Rust source, not an instrumented internal step counter. The first interval and <=1e-4s skipped intervals are treated as coded.'}

def audit_case(d,case,md):
 out=RESULTS/f'{d.name}_audit.json'
 if out.exists():
  print('SKIP existing audit',d.name);return json.loads(out.read_text())
 if not (d/'outputfiles/elem_temps.csv').exists(): raise FileNotFoundError(f'Missing raw output: {d}')
 p=dict(h=BASE['convective_heat_transfer_film_coefficient'],Tdep=BASE['extrusion_temperature'],Ta=BASE['ambient_temperature'],eps=BASE['emissivity'],Tbed=BASE['bed_temperature'],Rc=BASE['bed_contact_resistance'],kbed=BASE['conductivity_at_bed_surface_interface'],size=(0.55e-3,0.55e-3,0.2e-3))
 if EXPERIMENTS[case].get('bed_contact_model',2)==0:p['bed_law']='M0'
 if EXPERIMENTS[case].get('convection_model',0)==1:p['hz']=True
 tic=time.perf_counter()
 a=audit(str(d),p)
 extra=clock_summary(d,case)
 info={'case':case,'storage_threshold_K':md,'audit_wall_s':time.perf_counter()-tic,'audit':a,'clock_reconstruction':extra,'parameters':p,'gcode_sha256':shasum(GCODE),
       'important':'4 configurations share fixed material settings and G-code. M0/M2 and legacy cp differ by definition; common accounting method is used with model-specific documented boundary law. Underlying energy audit uses CORRECT physical clock and cp table for all models. No fan.'}
 out.write_text(json.dumps(info,indent=2)); print('AUDITED',d.name,'relative error',round(100*a['rel_err_vs_heat_in'],3),'%',flush=True)
 return info

def summary():
 RESULTS.mkdir(exist_ok=True)
 rows=[]
 for p in sorted(RESULTS.glob('*_audit.json')):
  j=json.loads(p.read_text());a=j['audit'];c=j['clock_reconstruction'];r={'case':j['case'],'storage_threshold_K':j['storage_threshold_K'],'n_elements':a['n'],
 'physical_deposition_s':c['physical_deposition_elapsed_s'],
 'integrated_deposition_s':c['reconstructed_deposition_integrated_legacy_s'] if c['case_uses_legacy_timing'] else c['reconstructed_deposition_integrated_corrected_s'],
 'deposition_clock_ratio':c['reconstructed_exposure_ratio_legacy_vs_corrected'] if c['case_uses_legacy_timing'] else 1.0,
 'nominal_cooldown_s':c['requested_cooldown_s'],
 'integrated_cooldown_s':c['legacy_cooldown_integrated_s'] if c['case_uses_legacy_timing'] else c['corrected_cooldown_integrated_s'],
 'deposited_thermal_excess_J':a['heat_in_J'],'stored_change_since_deposition_J':a['stored_change_J'],
 'natural_convection_J':a['natural_conv_J'],'forced_convection_J':a['forced_conv_J'],
 'radiation_J':a['radiation_J'],'bed_net_J':a['bed_net_J'],
 'total_loss_J':a['loss_total_J'],'signed_energy_balance_residual_J':a['balance_error_J'],
 'abs_residual_pct_heat_input':100*a['rel_err_vs_heat_in'],
 'audit_ends_at_s':a['t_end_s'],'min_temp_C':a['Tmin_C'],'max_temp_C':a['Tmax_C']};rows.append(r)
 if not rows: print('No audits yet');return
 with open(RESULTS/'N0_full_print_summary.csv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 # Differences across thresholds are *numerical storage sensitivity*, not sampling uncertainty.
 deltas=[]
 bycase={c:[r for r in rows if r['case']==c] for c in EXPERIMENTS}
 for c,vals in bycase.items():
  vals=sorted(vals,key=lambda x:x['storage_threshold_K'])
  if len(vals)>1:
   best=vals[0]
   for v in vals[1:]:
    deltas.append({'case':c,'reference_threshold_K':best['storage_threshold_K'],'test_threshold_K':v['storage_threshold_K'],
 'delta_signed_residual_J':v['signed_energy_balance_residual_J']-best['signed_energy_balance_residual_J'],
 'delta_abs_residual_pct_points':v['abs_residual_pct_heat_input']-best['abs_residual_pct_heat_input']})
 if deltas:
  with open(RESULTS/'N0_storage_sensitivity.csv','w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=deltas[0].keys());w.writeheader();w.writerows(deltas)
 print('Wrote',RESULTS/'N0_full_print_summary.csv')
 for r in rows:print(r['case'],f"threshold={r['storage_threshold_K']} K",f"energy residual={r['abs_residual_pct_heat_input']:.3f}%",f"clock ratio={r['deposition_clock_ratio']:.3f}")

if __name__=='__main__':
 parser=argparse.ArgumentParser()
 parser.add_argument('--cases',nargs='+',choices=EXPERIMENTS,default=list(EXPERIMENTS))
 parser.add_argument('--thresholds',nargs='+',type=float,default=[0.25,0.1])
 parser.add_argument('--prepare-only',action='store_true')
 parser.add_argument('--audit-only',action='store_true')
 parser.add_argument('--no-build',action='store_true')
 args=parser.parse_args();check()
 RESULTS.mkdir(exist_ok=True)
 if not(args.prepare_only or args.audit_only or args.no_build):compile_solver()
 for md in args.thresholds:
  if not (0<=md<=1):raise ValueError('Storage threshold must be between 0 and 1 K')
  for case in args.cases:
   d=RUNS/f'{case}_md{str(md).replace(".","p")}'
   if args.audit_only:
    if not d.exists():raise FileNotFoundError(d)
   else:
    d=prepare(case,md)
   if not (args.prepare_only or args.audit_only):run_case(d)
   if not args.prepare_only:audit_case(d,case,md)
 summary()
