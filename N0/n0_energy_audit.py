"""R8.2 independent energy audit, STREAMING per element (no global time grid, constant memory), usable on the full 100-layer print.
Uses only solver OUTPUT files + documented laws. Each stored history is sample-and-hold (value valid until the next stored sample), which
is within `element_min_temp_diff_store` of the solver state at all times. Per element, loss is integrated EXACTLY over the piecewise-constant
segments delimited by: its own stored times, its neighbours' activation times (faces become interior), and the fan schedule steps.
  H = rho V int c dT (table integral);  exposed face: [h_eff(t)(T-Ta) + eps sigma (T^4-Ta^4)] A ;  bed face: G (T-T_bed), G(M2)=A/(dz/(2kz)+Rc)
  h_eff = (h_nat^3 + h_forced(t)^3)^(1/3), h_forced = h_full (S/255)^n (documented fan model) evaluated at segment mid-time.
Reports balance error, split of losses (natural conv., forced conv. (fan), radiation, bed) and a negative control with the fan ignored."""
import sys, json, csv, numpy as np
RHO=1271.1; SIG=5.67e-8; KZ=0.35
tab=np.loadtxt(str(__import__('pathlib').Path(__file__).resolve().parent/'solver/inputfiles/sp_heat_cap_data.csv'),delimiter=','); cT,cY=tab[:,1],tab[:,0]
cum=np.concatenate([[0],np.cumsum(0.5*(cY[1:]+cY[:-1])*np.diff(cT))])
def Hs(T):
    T=np.asarray(T,float); i=np.clip(((T-cT[0])/1.0).astype(int),0,len(cT)-2); f=T-cT[i]
    return cum[i]+cY[i]*f+0.5*(cY[i+1]-cY[i])*f*f
def audit(d,P,fan=None):
    """fan = (list of (t,S), h_full, exponent) or None"""
    nodes=np.loadtxt(d+'/outputfiles/nodefile.csv',delimiter=',')[:,1:]
    el=np.atleast_2d(np.loadtxt(d+'/outputfiles/elementfile.csv',delimiter=',',dtype=int))[:,1:]
    bb=nodes[el-1]; cen=bb.mean(axis=1); sz=bb.max(axis=1)-bb.min(axis=1); dx,dy,dz=sz.mean(axis=0)
    if 'size' in P: dx,dy,dz=P['size']      # solver's element size (differs from the node spacing under the legacy D8 mesh)
    ta=np.loadtxt(d+'/outputfiles/activation_times.csv',delimiter=',')[:,0]; n=len(ta)
    lat=np.rint((cen-cen.min(axis=0))/np.array(sz.mean(axis=0))).astype(int); idx={tuple(l):i for i,l in enumerate(lat)}
    offs=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
    nb=np.array([[idx.get((l[0]+a,l[1]+b,l[2]+c),-1) for a,b,c in offs] for l in lat])
    bot=np.abs(bb.min(axis=1)[:,2]-nodes[:,2].min())<1e-9
    A=np.array([dy*dz,dy*dz,dx*dz,dx*dz,dx*dy,dx*dy]); V=dx*dy*dz; G=(dx*dy)/(dz/(2*KZ)+P['Rc']) if P.get('bed_law','M2')=='M2' else P['kbed']*dx*dy/dz
    hz=(3.0+12.0*cen[:,2]**2) if P.get('hz') else None
    ft=np.array([a for a,_ in fan[0]]) if fan else np.array([]); fS=np.array([b for _,b in fan[0]]) if fan else np.array([])
    def hforced(tm):
        if not fan: return np.zeros_like(tm)
        k=np.searchsorted(ft,tm,side='right')-1; sp=np.where(k>=0,fS[np.clip(k,0,None)],0.0)
        return np.where(sp>0,fan[1]*(sp/255.0)**fan[2],0.0)
    tot=dict(nat=0.0,forced=0.0,rad=0.0,bed=0.0,forced_ignored=0.0); Tend=np.empty(n); tmax_all=0.0; Tmin=1e9; Tmax=0.0
    with open(d+'/outputfiles/elem_temps.csv') as fh:
        for i,line in enumerate(fh):
            v=np.fromstring(line.strip().strip(','),sep=','); t=v[0::2]; T=v[1::2]; Tend[i]=T[-1]; tmax_all=max(tmax_all,t[-1]); Tmin=min(Tmin,T.min()); Tmax=max(Tmax,T.max())
            nbt=[ta[j] for j in nb[i] if j>=0 and ta[j]>ta[i]]
            B=np.unique(np.concatenate([[ta[i]],t[t>ta[i]],np.array(nbt),ft[(ft>ta[i])&(ft<t[-1])]])); B=B[B<=t[-1]]
            if len(B)<2: continue
            s0=B[:-1]; dt=np.diff(B); mid=s0+0.5*dt
            k=np.searchsorted(t,s0,side='right')-1; Ts=np.where(k>=0,T[np.clip(k,0,None)],P['Tdep'])
            hf=hforced(mid); hnat=P['h'] if hz is None else hz[i]; heff=np.where(hf>0,np.cbrt(hnat**3+hf**3),hnat)
            rad=P['eps']*SIG*(Ts**4-P['Ta']**4); dT=Ts-P['Ta']
            for f in range(6):
                if f==5 and bot[i]:
                    tot['bed']+=float((G*(Ts-P['Tbed'])*dt).sum()); continue
                j=nb[i,f]; exp=np.ones(len(s0),bool) if j<0 else (ta[j]>s0+1e-9)
                w=A[f]*dt*exp
                tot['nat']+=float((hnat*dT*w).sum()); tot['forced']+=float(((heff-hnat)*dT*w).sum()); tot['rad']+=float((rad*w).sum())
                tot['forced_ignored']+=0.0
    loss=tot['nat']+tot['forced']+tot['rad']+tot['bed']
    e_end=RHO*V*Hs(Tend).sum(); e_dep=RHO*V*Hs(np.full(n,P['Tdep'])).sum(); heat_in=e_dep-RHO*V*Hs(np.full(n,P['Ta'])).sum()
    err=e_end-e_dep+loss; neg=err-tot['forced']          # negative control: pretend the fan did not remove heat
    return dict(n=n,heat_in_J=heat_in,stored_change_J=e_end-e_dep,loss_total_J=loss,natural_conv_J=tot['nat'],forced_conv_J=tot['forced'],radiation_J=tot['rad'],bed_net_J=tot['bed'],
                balance_error_J=err,rel_err_vs_heat_in=abs(err)/heat_in,rel_err_vs_loss=abs(err)/abs(loss),
                negative_control_fan_ignored=dict(balance_error_J=neg,rel_err_vs_heat_in=abs(neg)/heat_in),Tmin_C=float(Tmin-273.15),Tmax_C=float(Tmax-273.15),t_end_s=tmax_all)

