"""B01 G-code audit + adapter for the R6 solver reader, and an independent nominal deposition timeline.
Adapter rules (documented in the report):
  * the R6/Bhandari reader only parses lines between ';Printing starts here' and '; layer end' and only G1 lines whose
    first word is X/Y/Z -> the adapter wraps the executed moves in those markers;
  * the reader IGNORES G0, but B01 uses G0 for every travel and layer change: ignoring them would turn each
    travel into a fictitious extruding bead -> every G0 is rewritten as a non-extruding G1 (same X/Y/Z/F);
  * extrusion is inferred from the presence of an E word (absolute E, only present on extruding moves in B01);
  * non-motion lines (START macro, M109/M82/G90/G92, M106/M107, JOB_COMPLETE, comments) are removed and audited separately."""
import re, sys, json, collections, numpy as np
SRC='/home/claude/work/r8/b01/inputfiles/B01_ADHESION_RETRY_PREVIEW_ONLY.gcode'
def parse(src=SRC):
    ev=[]; cnt=collections.Counter(); layer=-1; x=y=z=0.0; f=0.0; e=0.0; fan=0; other=[]; fan_events=[]
    for ln,raw in enumerate(open(src),1):
        s=raw.split(';')[0].strip()
        if raw.startswith(';LAYER:'): layer=int(raw.split(':')[1])
        if not s: continue
        w=s.split(); c=w[0]; cnt[c]+=1
        if c in ('G0','G1'):
            p={k[0]:float(k[1:]) for k in w[1:]}
            nx,ny,nz=p.get('X',x),p.get('Y',y),p.get('Z',z); f=p.get('F',f); ext=('E' in p)
            de=p['E']-e if ext else 0.0
            if ext: e=p['E']
            ev.append(dict(line=ln,cmd=c,layer=layer,x0=x,y0=y,z0=z,x1=nx,y1=ny,z1=nz,f=f,extr=ext,de=de))
            x,y,z=nx,ny,nz
        elif c=='M106': fan_events.append((ln,layer,s))
        elif c in ('M107','M82','G90','M109','G92','START','JOB_COMPLETE'): pass
        else: other.append((ln,s))
    return ev,cnt,fan_events,other
def timeline(ev):
    t=0.0; out=[]
    for m in ev:
        L=np.sqrt((m['x1']-m['x0'])**2+(m['y1']-m['y0'])**2+(m['z1']-m['z0'])**2); dt=L/(m['f']/60.0) if L>0 else 0.0
        out.append((t,t+dt,L)); t+=dt
    return out,t
def write_solver_gcode(ev,path,max_layer=None):
    with open(path,'w') as g:
        g.write(';Printing starts here\n')
        for m in ev:
            if max_layer is not None and m['layer']>max_layer: break
            w=['G1']
            if abs(m['x1']-m['x0'])>0 or m is ev[0] or True: pass
            w+= [f"X{m['x1']:.4f}",f"Y{m['y1']:.4f}"] if (m['x1']!=m['x0'] or m['y1']!=m['y0'] or m['z1']==m['z0']) else []
            if m['z1']!=m['z0']: w.append(f"Z{m['z1']:.4f}")
            if len(w)==1: continue
            if m['extr']: w.append(f"E{m['de']:.5f}")
            w.append(f"F{m['f']:.0f}")
            g.write(' '.join(w)+'\n')
        g.write('; layer end\n')
if __name__=='__main__':
    ev,cnt,fan,other=parse(); tl,T=timeline(ev)
    print('commands',dict(cnt)); print('fan',fan); print('unrecognised',other[:5])
    ex=[m for m in ev if m['extr']]; print('moves',len(ev),'extruding',len(ex),'path mm',round(sum(np.hypot(m['x1']-m['x0'],m['y1']-m['y0']) for m in ex),1))
    print('nominal time s',round(T,1))
    write_solver_gcode(ev,'/home/claude/work/r8/b01/inputfiles/b01_solver.gcode')
