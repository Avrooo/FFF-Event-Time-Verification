"""N1 representative raw-history verification for Windows/Jupyter.

Default: INSPECT and re-audit available raw histories; never re-run missing solver cases.
Explicit --rerun-missing: run ONLY the four representative Rust cases with missing raw histories,
back up/restore their existing official results JSON, then re-audit using the frozen N1 auditor.

This is an independent *re-execution* of the published postprocessing, plus algebraic
cross-checks, NOT a new independently coded energy quadrature method.
"""
from pathlib import Path
import argparse, csv, hashlib, json, os, shutil, sys, traceback, zipfile

SPECS=[
    ('b01_six_layer_segment', 40, '0.05', 'corrected'),
    ('b01_six_layer_segment', 40, '0.05', 'legacy_timing'),
    ('independent_serpentine_six_layers', 40, '0.05', 'corrected'),
    ('independent_serpentine_six_layers', 40, '0.05', 'legacy_timing'),
]

def find_root():
    locs=[Path.cwd(),Path(__file__).resolve().parent,Path.home()/'Downloads']
    for loc in locs:
        if (loc/'run_full_solver_N1.py').exists() and (loc/'N1_energy_audit.py').exists():return loc
    for loc in locs:
        if loc.exists():
            for f in loc.rglob('run_N1_complete.py'):
                if (f.parent/'N1_energy_audit.py').exists():return f.parent
    raise FileNotFoundError('Cannot locate extracted NAMRC_N1_FINAL; place this script in the experiment folder.')

def hashfile(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def name(spec):
    tp,v,dt,mode=spec
    return f'{tp}__v{v}__dt{dt}__{mode}'

def case_dict(spec):
    tp,v,dt,mode=spec
    return dict(toolpath=tp,speed_mm_s=str(v),dt_s=dt,timing_mode=mode)

def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--rerun-missing',action='store_true',help='Explicitly re-run only missing raw-history cases')
    parser.add_argument('--include-raw',action='store_true',help='Add four raw elem_temps.csv to ZIP (may be very large)')
    args=parser.parse_args()
    ROOT=find_root();os.chdir(ROOT);sys.path.insert(0,str(ROOT))
    os.environ['PATH']=str(Path.home()/'.cargo'/'bin')+os.pathsep+os.environ.get('PATH','')
    import run_full_solver_N1 as runner
    from N1_energy_audit import audit
    RES=ROOT/'N1_results'; RUN=ROOT/'runs'; DEST=ROOT/'N1_representative_verification';DEST.mkdir(exist_ok=True)
    P=dict(h=3.0,Tdep=488.15,Ta=298.15,eps=0.95,Tbed=333.15,Rc=0.005)
    errors=[]; rows=[];manifest=[]
    print('N1 root:',ROOT,flush=True)
    missing=[s for s in SPECS if not (RUN/name(s)/'outputfiles'/'elem_temps.csv').exists()]
    print('Raw histories available:',len(SPECS)-len(missing),'/4',flush=True)
    for s in SPECS:
        p=RUN/name(s)/'outputfiles'/'elem_temps.csv'
        print(('PRESENT ' if p.exists() else 'MISSING ')+name(s),flush=True)
    if missing and not args.rerun_missing:
        print('\nMissing raw history for',len(missing),'case(s). No simulation will be started.',flush=True)
        print('To regenerate ONLY missing cases, run this script with --rerun-missing.',flush=True)
    if missing and args.rerun_missing:
        exe=runner.SRC/'target'/'release'/('DES_thermal_simulation.exe' if os.name=='nt' else 'DES_thermal_simulation')
        if not exe.exists():exe=runner.build()
        for s in missing:
            nm=name(s);rp=RES/(nm+'.json');backup=DEST/(nm+'.original_report_backup.json')
            if rp.exists():
                shutil.copy2(rp,backup)
                rp.unlink()  # disable cache; restore below in finally
            print('REGENERATING',nm,flush=True)
            try:
                outcome=runner.execute_case(case_dict(s),exe)
                # preserve new result separately; restore the original metadata below
                (DEST/(nm+'.regenerated_run.json')).write_text(json.dumps(outcome,indent=2))
                if not (RUN/nm/'outputfiles'/'elem_temps.csv').exists():raise RuntimeError('Solver did not produce raw history')
            except Exception as e:
                errors.append(dict(case=nm,stage='regenerate',error=str(e)))
                print('FAILED',nm,e,flush=True)
            finally:
                if backup.exists():shutil.copy2(backup,rp)
    for s in SPECS:
        nm=name(s);rd=RUN/nm;raw=rd/'outputfiles'/'elem_temps.csv'
        old=RES/('energy_audit_'+nm+'.json')
        item=dict(case=nm,raw_present=raw.exists(),raw_bytes=raw.stat().st_size if raw.exists() else 0,
                  raw_sha256=hashfile(raw) if raw.exists() else '')
        if old.exists():item['archived_audit_sha256']=hashfile(old)
        manifest.append(item)
        if not raw.exists():
            errors.append(dict(case=nm,stage='raw_presence',error='Raw history unavailable; not verified'))
            continue
        print('AUDITING',nm,flush=True)
        try:
            a=audit(str(rd),P)
            # Separate arithmetic recheck of the aggregate energy identity.
            identity=a['stored_change_J']+a['natural_conv_J']+a['forced_conv_J']+a['radiation_J']+a['bed_net_J']
            if abs(identity-a['balance_error_J'])>1e-6:raise AssertionError('Aggregate energy identity failure')
            if abs(abs(a['balance_error_J'])/a['heat_in_J']-a['rel_err_vs_heat_in'])>1e-9:
                raise AssertionError('Relative energy residual identity failure')
            o=json.loads(old.read_text())['audit'] if old.exists() else None
            row=dict(case=nm,toolpath=s[0],speed_mm_s=s[1],dt_s=float(s[2]),timing_mode=s[3],
                     n_elements=a['n'],deposited_heat_J=a['heat_in_J'],
                     stored_change_J=a['stored_change_J'],
                     natural_convection_J=a['natural_conv_J'],forced_convection_J=a['forced_conv_J'],
                     radiation_J=a['radiation_J'],bed_net_J=a['bed_net_J'],
                     signed_residual_J=a['balance_error_J'],recomputed_residual_pct=100*a['rel_err_vs_heat_in'],
                     original_residual_pct=(100*o['rel_err_vs_heat_in'] if o else ''),
                     residual_difference_percentage_points=(100*(a['rel_err_vs_heat_in']-o['rel_err_vs_heat_in']) if o else ''),
                     raw_sha256=item['raw_sha256'])
            rows.append(row)
            (DEST/(nm+'.verified_audit.json')).write_text(json.dumps(dict(audit=a,comparison=row,assumptions=P),indent=2))
            if o is not None:
                for key in ('heat_in_J','stored_change_J','natural_conv_J','radiation_J','bed_net_J','balance_error_J'):
                    # Archive values are generated with the same numerical postprocessor, so expect close agreement.
                    tolerance=max(0.05,0.002*max(abs(o[key]),1))
                    if abs(a[key]-o[key])>tolerance:
                        errors.append(dict(case=nm,stage='compare',error=f'{key}: recalculated {a[key]}, archive {o[key]}, tolerance {tolerance}'))
            print(f'    residual {row["recomputed_residual_pct"]:.5f}% | archived {row["original_residual_pct"]}',flush=True)
        except Exception as e:
            errors.append(dict(case=nm,stage='audit',error=str(e)))
            print('AUDIT FAILED',nm,e,flush=True)
    if rows:
        with (DEST/'N1_four_case_reaudited.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (DEST/'N1_raw_history_manifest.json').write_text(json.dumps(manifest,indent=2))
    status=dict(cases_requested=4,raw_available=sum(x['raw_present'] for x in manifest),cases_reaudited=len(rows),
                failures=errors,all_four_reaudited=(len(rows)==4 and not errors),
                qualifier='Re-execution of same independent-of-solver energy postprocessor; not independent numerical implementation of energy quadrature.')
    (DEST/'N1_verification_status.json').write_text(json.dumps(status,indent=2))
    # Include input settings and original records for transparency, but not raw files unless opted in.
    zpath=ROOT/'N1_representative_evidence_for_review.zip'
    with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=3) as z:
        for p in DEST.glob('*'):
            if p.is_file():z.write(p,p.relative_to(ROOT))
        for s in SPECS:
            nm=name(s)
            for p in [RUN/nm/'inputfiles'/'Input_file.txt',RUN/nm/'inputfiles'/'experiment.gcode',RES/(nm+'.json'),RES/('energy_audit_'+nm+'.json')]:
                if p.exists():z.write(p,p.relative_to(ROOT))
            if args.include_raw:
                p=RUN/nm/'outputfiles'/'elem_temps.csv'
                if p.exists():z.write(p,p.relative_to(ROOT))
    print('\nSTATUS:',json.dumps(status,indent=2),flush=True)
    print('Upload this file:',zpath,flush=True)
    if len(rows)!=4 or errors:
        print('WARNING: Verification incomplete; do not claim all 4 verified.',flush=True)
    return status

if __name__=='__main__':run()
