#!/usr/bin/env python3
"""NAMRC N1: 36 actual Rust thermal-solver simulations.
Run from this folder: python run_full_solver_N1.py
Requires Python 3 and Rust Cargo. No packages beyond Python standard library required.
Uses frozen R8.4 code, independent six-layer B01 and serpentine toolpaths.
DO NOT use path_level_clock_diagnostics.csv as thermal-solver outcomes.
"""
import argparse,csv,json,math,os,platform,shutil,statistics,subprocess,sys,time,traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parent; SRC=ROOT/'solver'; RUN=ROOT/'runs'; RES=ROOT/'N1_results'

# Required physical settings: R8.4 frozen no-fan baseline. T=25 C to isolate numerical clock effects.
SETTINGS={
 'Number_of_threads':'1','GCode_file':'inputfiles/experiment.gcode',
 'divisions_per_bead_width':'1','divisions_per_bead_height':'1',
 'beadwidth':'0.00055','beadheight':'0.0002','model_type':'2',
 'specific_heat_capacity_file_temp_step':'1',
 'kx':'0.35','ky':'0.35','kz':'0.35','density':'1271.1',
 'convective_heat_transfer_film_coefficient':'3.0',
 'bed_temperature':'333.15','conductivity_at_bed_surface_interface':'0.17',
 'ambient_temperature':'298.15','extrusion_temperature':'488.15',
 'emissivity':'0.95','cooldown_period':'60.0',
 'element_min_temp_diff_store':'0.1','turn_off_layers_at':'1000',
 'nodal_temperature_file_output':'/dev/null' if os.name!='nt' else 'NUL',
 'bed_contact_model':'2','bed_contact_resistance':'0.005','convection_model':'0',
 'fan_model':'0','legacy_cp_step':'0','legacy_bed_flag':'0','legacy_mesh':'0'
}

def write_settings(input_path,speed,dt,mode):
    vals=dict(SETTINGS);vals['time_step']=str(dt);vals['legacy_timing']='1' if mode=='legacy_timing' else '0'
    template=(SRC/'inputfiles'/'Input_file.txt').read_text().splitlines()
    output=[]; seen=set()
    for line in template:
        key=line.split(':',1)[0].strip()
        if key in vals:
            output.append(f'{key}: {vals[key]} // N1 fixed experiment');seen.add(key)
        else: output.append(line)
    for k,v in vals.items():
        if k not in seen:output.append(f'{k}: {v} // N1 fixed experiment')
    input_path.write_text('\n'.join(output)+'\n')

def set_speed(src,dst,speed):
    # Change only extrusion feed. All non-extruding travel feeds are retained.
    # G1 reader treats E present as extrusion; each extrusion line lists X/Y before E.
    n=0;lines=[]
    for line in src.read_text().splitlines():
        if line.startswith('G1 ') and any(p.startswith('E') for p in line.split()[1:]):
            if 'F' not in line: raise ValueError('Extrusion missing F: '+line)
            parts=line.split();parts=[('F'+str(speed*60)) if p.startswith('F') else p for p in parts]
            line=' '.join(parts);n+=1
        lines.append(line)
    if n<50:raise ValueError('Too few extrusion paths')
    dst.write_text('\n'.join(lines)+'\n')
    return n

def read_final_temps(path):
    # Last stored temperature of every active element, Celsius.
    out=[];tmax=0.0;tmin=float('inf');temperature_max=-float('inf')
    with path.open() as f:
        for ln in f:
            a=ln.strip().strip(',').split(',')
            if len(a)<2: continue
            try:
                t=float(a[-2]);v=float(a[-1]);out.append(v-273.15)
                tmax=max(tmax,t);tmin=min(tmin,v-273.15);temperature_max=max(temperature_max,v-273.15)
            except ValueError:continue
    if not out: raise RuntimeError('No element-temperature histories')
    return out,tmax

def build():
    if not shutil.which('cargo'):
        raise RuntimeError('cargo not found. Install Rust from https://rustup.rs/ and open a new terminal.')
    subprocess.run(['cargo','build','--release'],cwd=SRC,check=True)
    exe=SRC/'target'/'release'/('DES_thermal_simulation.exe' if os.name=='nt' else 'DES_thermal_simulation')
    if not exe.exists():raise RuntimeError(f'Binary missing: {exe}')
    return exe

def execute_case(case,exe,limit=None):
    name=f"{case['toolpath']}__v{case['speed_mm_s']}__dt{case['dt_s']}__{case['timing_mode']}"
    rd=RUN/name; rp=RES/(name+'.json')
    if rp.exists():
        old=json.loads(rp.read_text());
        if old.get('status')=='complete':return old
    (rd/'inputfiles').mkdir(parents=True,exist_ok=True)
    (rd/'outputfiles').mkdir(exist_ok=True)
    src=ROOT/'toolpaths'/(case['toolpath']+'.gcode')
    n=set_speed(src,rd/'inputfiles'/'experiment.gcode',int(case['speed_mm_s']))
    for f in ['sp_heat_cap_data.csv','conductivity.csv']:
        shutil.copy2(SRC/'inputfiles'/f,rd/'inputfiles'/f)
    write_settings(rd/'inputfiles'/'Input_file.txt',int(case['speed_mm_s']),float(case['dt_s']),case['timing_mode'])
    t=time.perf_counter()
    with (rd/'solver_stdout.log').open('w') as log:
        p=subprocess.run([str(exe.resolve())],cwd=rd,stdout=log,stderr=subprocess.STDOUT,timeout=limit)
    elapsed=time.perf_counter()-t
    if p.returncode!=0:
        raise RuntimeError(f'{name}: solver exit={p.returncode}; read {rd / "solver_stdout.log"}')
    temps,tmax=read_final_temps(rd/'outputfiles'/'elem_temps.csv')
    result={**case,'status':'complete','case':name,'n_extrusion_moves':n,
         'n_elements':len(temps),'time_last_stored_s':tmax,
         'mean_final_C':statistics.fmean(temps),'min_final_C':min(temps),'max_final_C':max(temps),
         'execution_wall_s':elapsed,
         'full_rust_solver':True,
         'numerical_energy_audit_completed':False,
         'executable_commit':(ROOT/'R8_4_source_commit.txt').read_text().strip() if (ROOT/'R8_4_source_commit.txt').exists() else 'unknown'}
    rp.write_text(json.dumps(result,indent=2));return result

def paired_report():
    pairs={}
    for p in RES.glob('*.json'):
        r=json.loads(p.read_text());
        if r.get('status')!='complete':continue
        k=(r['toolpath'],int(r['speed_mm_s']),float(r['dt_s']))
        pairs.setdefault(k,{})[r['timing_mode']]=r
    out=[]
    for (tp,v,dt),p in sorted(pairs.items()):
        if not {'legacy_timing','corrected'}<=set(p):continue
        c=p['corrected']; l=p['legacy_timing']
        if c['n_elements'] != l['n_elements']:
            raise RuntimeError(f'Element count mismatch at {(tp,v,dt)}: {c["n_elements"]} != {l["n_elements"]}')
        out.append(dict(toolpath=tp,speed_mm_s=v,dt_s=dt,n_elements=c['n_elements'],
                        corrected_final_mean_C=c['mean_final_C'],legacy_final_mean_C=l['mean_final_C'],
                        legacy_minus_corrected_final_C=l['mean_final_C']-c['mean_final_C'],
                        runtime_corrected_s=c['execution_wall_s'],runtime_legacy_s=l['execution_wall_s']))
    RES.mkdir(exist_ok=True)
    with (RES/'paired_full_solver_results.csv').open('w',newline='') as f:
        if out:
            w=csv.DictWriter(f,out[0].keys());w.writeheader();w.writerows(out)
    print(f'Completed pairs {len(out)}/18; cases {sum(len(v) for v in pairs.values())}/36')

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--max-cases',type=int,default=None,help='Run only N cases for a smoke test')
    ap.add_argument('--no-build',action='store_true')
    ap.add_argument('--only',default='',help='Filter case names with substring e.g. serpentine')
    ap.add_argument('--timeout',type=float,default=None)
    ap.add_argument('--keep-outputs',action='store_true',help='Keep raw thermal output directories (large)')
    args=ap.parse_args()
    RUN.mkdir(exist_ok=True); RES.mkdir(exist_ok=True)
    if not (SRC/'Cargo.toml').exists():sys.exit('Missing solver source directory')
    exe=(SRC/'target'/'release'/('DES_thermal_simulation.exe' if os.name=='nt' else 'DES_thermal_simulation')) if args.no_build else build()
    with (ROOT/'experiment_matrix.csv').open() as f:rows=list(csv.DictReader(f))
    if args.only:rows=[r for r in rows if args.only in r['toolpath'] or args.only in str(r)]
    if args.max_cases is not None:rows=rows[:args.max_cases]
    errors=[]
    for i,r in enumerate(rows,1):
        try:
            ans=execute_case(r,exe,args.timeout)
            print(f'[{i}/{len(rows)}] {ans["case"]}: {ans["n_elements"]} elements; final mean {ans["mean_final_C"]:.2f} C; {ans["execution_wall_s"]:.1f} s',flush=True)
            if not args.keep_outputs:
                shutil.rmtree(RUN/ans['case'],ignore_errors=True)
        except Exception as e:
            print('FAILED:',r,e,flush=True);errors.append(dict(case=r,error=str(e)))
    paired_report()
    (RES/'failures.json').write_text(json.dumps(errors,indent=2))
    if errors:sys.exit(f'{len(errors)} cases failed; see N1_results/failures.json')
