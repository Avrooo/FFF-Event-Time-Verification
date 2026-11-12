"""R8.4 controlled ablation. A0 = corrected R8.2 reference (no fan). Each A* restores ONE original behaviour; A6 restores all (must equal the original solver);
R* = A6 with ONE behaviour corrected again (interaction analysis). Same B01 G-code, inputs, dt=0.05, storage 0.25 K, hold extraction."""
import sys, json, shutil, numpy as np
sys.path.insert(0,'/home/claude/work/r8')
from b01_case import make
from run_case import run_case
from b01_extract import extract
CP=dict(specific_heat_capacity_file_temp_step=5)   # the published input declares 5 K for the 1 K table
LEG=dict(legacy_cp_step=1,legacy_timing=1,legacy_bed_flag=1,legacy_mesh=1)
ALL=dict(CP,bed_contact_model=0,convection_model=1,**LEG)
def minus(*keys):
    d=dict(ALL)
    for k in keys:
        if k=='cp': d.pop('legacy_cp_step')
        if k=='timing': d.pop('legacy_timing')
        if k=='bedflag': d.pop('legacy_bed_flag')
        if k=='bedlaw': d['bed_contact_model']=2
        if k=='conv': d['convection_model']=0
        if k=='mesh': d.pop('legacy_mesh')
    return d
EXP={'A0':{}, 'A1':dict(CP,legacy_cp_step=1), 'A2':dict(legacy_timing=1), 'A3a':dict(bed_contact_model=0), 'A3b':dict(legacy_bed_flag=1),
     'A3':dict(bed_contact_model=0,legacy_bed_flag=1), 'A4':dict(convection_model=1), 'A5':dict(legacy_mesh=1), 'A6':ALL,
     'R1':minus('cp'),'R2':minus('timing'),'R3':minus('bedlaw','bedflag'),'R4':minus('conv'),'R5':minus('mesh')}
times=np.arange(0,991,1.0)
def build(name,max_layer=None,md=0.25,dt=0.05):
    kv=dict(EXP[name]); legacy={k:v for k,v in kv.items() if k.startswith('legacy_')}; over={k:v for k,v in kv.items() if not k.startswith('legacy_')}
    d=make('abl_'+name,max_layer=max_layer,element_min_temp_diff_store=md,time_step=dt,**over)
    p=d+'/inputfiles/Input_file.txt'; L=open(p).read().splitlines()+[f'{k}: {v} //' for k,v in legacy.items()]; open(p,'w').write('\n'.join(L)+'\n'); return d
if __name__=='__main__':
    for name in sys.argv[1:]:
        d=build(name); info=run_case(d)
        if info['returncode']!=0: json.dump(dict(name=name,failed=True,**info),open(f'results/abl_{name}.json','w')); print(name,'FAILED',flush=True); continue
        F=extract(d,times,faces=('ymin','ymax')); np.savez_compressed(f'results/abl_{name}_field.npz',times=times,layers=F['ymin']['layers'],**{f+'_T':F[f]['T'] for f in F},**{f+'_coord_mm':F[f]['coord'] for f in F})
        # last-stored temperatures of every element (for exact-equality checks) and global extrema
        import numpy as _n
        Tfin=[]; Tmax=0; Tmin=1e9
        for line in open(d+'/outputfiles/elem_temps.csv'):
            v=_n.fromstring(line.strip().strip(','),sep=','); T=v[1::2]; Tfin.append(T[-1]); Tmax=max(Tmax,T.max()); Tmin=min(Tmin,T.min())
        _n.save(f'results/abl_{name}_final_T.npy',_n.array(Tfin))
        json.dump(dict(name=name,keys={k:str(v) for k,v in EXP[name].items()},Tmin_C=Tmin-273.15,Tmax_C=Tmax-273.15,**info),open(f'results/abl_{name}.json','w')); print(name,round(info['wall_s'],1),flush=True)
        shutil.rmtree(d,ignore_errors=True)
