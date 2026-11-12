from pathlib import Path
import csv, math, json, re, shutil
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'toolpaths';OUT.mkdir(exist_ok=True)
SRC=ROOT/'b01/inputfiles/b01_solver.gcode'
lines=SRC.read_text().splitlines(); out=[]; current_layer=-1; layer_z=[]
for ln in lines:
    if ln.startswith(';'): continue
    match=re.search(r'\bZ([\d.]+)',ln)
    if match:
        z=float(match.group(1))
        if not layer_z or abs(z-layer_z[-1])>1e-4:
            if z>=1.4-1e-4: break  # first six complete deposition layers, 0.2-1.2 mm
            layer_z.append(z)
    out.append(ln)
assert len(layer_z)==6,layer_z
b01=[';Printing starts here']+out+['; layer end']
(OUT/'b01_six_layer_segment.gcode').write_text('\n'.join(b01)+'\n')
# Independent geometry: 18x10 mm short-edge alternating zigzag with boundary per layer
# The 0.55 mm road spacing yields different event-length distribution than B01 brim and solid-floor paths.
g=[';Printing starts here','G1 Z0.2000 F5400','G1 X60.0000 Y60.0000 F5400']
for li in range(6):
    z=.2*(li+1)
    g += [f'G1 Z{z:.4f} F5400',f'G1 X60.0000 Y60.0000 F5400']
    # outer 18x10 perimeter, followed by serpentine hatch
    for x,y in [(78,60),(78,70),(60,70),(60,60)]:g.append(f'G1 X{x:.4f} Y{y:.4f} E1.00000 F2400')
    for j in range(1,18):
        y=60+j*.55
        xstart=60 if j%2 else 78
        xend=78 if j%2 else 60
        g.append(f'G1 X{xstart:.4f} Y{y:.4f} F5400')
        g.append(f'G1 X{xend:.4f} Y{y:.4f} E1.00000 F2400')
g+=['; layer end']
(OUT/'independent_serpentine_six_layers.gcode').write_text('\n'.join(g)+'\n')

def audit(path):
    events=[]; x=y=z=0.; curF=0.; prev=None; n_extr=0; layerz=set(); dists=[]
    for line in path.read_text().splitlines():
        if not line.startswith('G1 '): continue
        tokens={p[0]:float(p[1:]) for p in line.split()[1:]}
        nx,ny,nz=(tokens.get(k,v) for k,v in [('X',x),('Y',y),('Z',z)])
        f=tokens.get('F',curF); length=math.dist((nx,ny,nz),(x,y,z))
        if 'E' in tokens and length>1e-8:
            n_extr+=1; layerz.add(round(nz,4));dists.append(length)
        if f<=0:raise ValueError('zero feed rate')
        events.append(dict(x0=x,y0=y,z0=z,x1=nx,y1=ny,z1=nz,length_mm=length,extruding='E' in tokens,base_feed_mm_min=f))
        x,y,z,curF=nx,ny,nz,f
    if len(layerz)!=6 or n_extr<50: raise ValueError((len(layerz),n_extr))
    return events,dict(name=path.name,layers=len(layerz),extruding_moves=n_extr,median_extrusion_length_mm=sorted(dists)[len(dists)//2],total_extrusion_length_mm=sum(dists),all_motion_moves=len(events),z_levels=sorted(layerz))
summary={}
for path in sorted(OUT.glob('*.gcode')):
    e,st=audit(path);summary[path.stem]=st
(ROOT/'toolpath_audit.json').write_text(json.dumps(summary,indent=2))
# 36-run matrix, 3 speeds x 3 dt x 2 modes x 2 toolpaths
rows=[]
for name in summary:
    for speed in [20,40,80]:
        for dt in [.025,.05,.10]:
            for mode in ['corrected','legacy_timing']:
                rows.append(dict(toolpath=name,speed_mm_s=speed,dt_s=dt,timing_mode=mode))
with (ROOT/'experiment_matrix.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
# independent, deterministic movement-level clock diagnostic (not full solver or element activation)
clock=[]
for row in rows:
    events,_=audit(OUT/(row['toolpath']+'.gcode'))
    dt=row['dt_s']; nominal=0.; legacy=0.; n_part=0
    for ev in events:
        f=60*row['speed_mm_s'] if ev['extruding'] else ev['base_feed_mm_min']
        d=ev['length_mm']/(f/60.)
        if d<=0:continue
        nominal+=d
        if d<dt:
            n_part+=1
        # Illustrative clock-only diagnostic: full-step legacy integration on sub-dt intervals;
        # differs from Rust's element-level activation and is NOT a full-solver measurement.
        legacy+= max(dt,d) if row['timing_mode']=='legacy_timing' else d
    clock.append(dict(**row,path_level_nominal_time_s=round(nominal,6),path_level_legacy_proxy_s=round(legacy,6),path_level_proxy_inflation=round(legacy/nominal,6),n_sub_dt_moves=n_part))
with (ROOT/'path_level_clock_diagnostics.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=clock[0]);w.writeheader();w.writerows(clock)
print(json.dumps(summary,indent=2));print('Experiment matrix:',len(rows),'cases; path-level diagnostics prepared, NOT full solver runs.')
