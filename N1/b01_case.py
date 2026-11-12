"""Build a B01 case directory for the R6 solver. All physical values are the documented baseline (see R8_1_report.md); nothing is fitted."""
import os, shutil, re
FORK='/home/claude/work/r4/corrected'
BASE=dict(  # SI units as the solver expects
  beadwidth=0.55e-3, beadheight=0.2e-3, divisions_per_bead_width=1, divisions_per_bead_height=1, model_type=2,
  kx=0.35, ky=0.35, kz=0.35, density=1271.1, convective_heat_transfer_film_coefficient=3.0,
  bed_temperature=333.15, conductivity_at_bed_surface_interface=0.17, ambient_temperature=298.15,
  extrusion_temperature=488.15, emissivity=0.95, time_step=0.05, cooldown_period=60.0,
  element_min_temp_diff_store=1.0, turn_off_layers_at=1000, Number_of_threads=1,
  specific_heat_capacity_file_temp_step=1, nodal_temperature_file_output='/dev/null', bed_contact_model=2, bed_contact_resistance=5e-3, convection_model=0)
def make(name, root='/home/claude/work/r8/runs', max_layer=None, **over):
    kv=dict(BASE); kv.update(over); d=f'{root}/{name}'; shutil.rmtree(d,ignore_errors=True)
    os.makedirs(d+'/inputfiles'); os.makedirs(d+'/outputfiles')
    for f in ('sp_heat_cap_data.csv','conductivity.csv'): shutil.copy(f'{FORK}/inputfiles/{f}',d+'/inputfiles/')
    if max_layer is None: shutil.copy('/home/claude/work/r8/b01/inputfiles/b01_solver.gcode',d+'/inputfiles/b01_solver.gcode')
    else:
        import sys; sys.path.insert(0,'/home/claude/work/r8'); from gcode_audit import parse,write_solver_gcode
        write_solver_gcode(parse()[0],d+'/inputfiles/b01_solver.gcode',max_layer=max_layer)
    kv['GCode_file']='inputfiles/b01_solver.gcode'
    out=[]; seen=set()
    for l in open(f'{FORK}/inputfiles/Input_file.txt').read().splitlines():
        k=l.split(':')[0]
        if k in kv: out.append(f'{k}: {kv[k]} //'); seen.add(k)
        else: out.append(l)
    out=[l for l in out if not l.startswith(('bed_contact_model','bed_contact_resistance','convection_model'))]
    for k in ('bed_contact_model','bed_contact_resistance','convection_model'): out.append(f'{k}: {kv[k]} //')
    open(d+'/inputfiles/Input_file.txt','w').write('\n'.join(out)+'\n'); return d
if __name__=='__main__':
    import sys; print(make(sys.argv[1] if len(sys.argv)>1 else 'baseline_test'))
