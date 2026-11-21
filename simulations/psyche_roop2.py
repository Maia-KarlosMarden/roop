"""
Integrated ROOP loop repeated on Psyche's next transfers (B1 -> D descent, D -> C raise).
The same relational deviation drives both fronts:
  DeltaF_geom = G'r - G p_geom = h^2 (1/M_hat - 1/M)  -> mass, inside an adjustable admissible band
  DeltaF      = -(v/2)(<G'r>/(G p_ref) - 1)            -> thrust command with the current M_hat
At low altitude the oblateness dominates, so the loop also estimates J2 on line: each quarter-orbit window
gives mu_eff = mu (1 + J2 g), g = 1.5 <(R/r)^2 (1 - 3 z^2/r^2)>, and a running linear regression over the
altitudes visited during the transfer returns mu and J2 (prior: shape-model C20 with 20 % error).
Variants:
  GROUND    : waits T_ground days at the departure orbit for a ground gravity solution (GM 0.1 %, J2 1 %), then transfers
  ROOP_J2   : leaves immediately with M_hat carried from the previous phase (+0.15 %, band +/-0.2 %), mass and J2 on line
  ROOP_FIX  : as ROOP_J2 but J2 frozen at the shape-model prior (ablation)
"""
import numpy as np, json, sys
sys.path.insert(0,'.')
import sim_bodies as S
B=S.BODIES['psyche']; GM,R=B['GM'],B['R']; ws=2*np.pi/(B['Pspin']*3600); G=S.G; J2T=-B['C20']
E_DAMP=0.04
AT=1.0e-4; DT=60.0   # SPT-140-class thrust authority (~0.27 N / 2700 kg)

def rk4(f,t,x,h):
    k1=f(t,x); k2=f(t+h/2,x+h/2*k1); k3=f(t+h/2,x+h/2*k2); k4=f(t+h,x+h*k3)
    return x+h/6*(k1+2*k2+2*k3+k4)

INC_DEG=150.0   # retrograde w.r.t. Psyche spin: stable crossing of the 3:2 and 1:1 spin-orbit resonances
def run(variant, r0, r1, days, seed=0, T_ground=7.0):
    rng=np.random.default_rng(100+seed)
    inc=np.radians(INC_DEG+rng.uniform(-3,3)); node=rng.uniform(0,2*np.pi)
    Rx=np.array([[1,0,0],[0,np.cos(inc),-np.sin(inc)],[0,np.sin(inc),np.cos(inc)]])
    Rz=np.array([[np.cos(node),-np.sin(node),0],[np.sin(node),np.cos(node),0],[0,0,1]]); Q=Rz@Rx
    fd=rng.normal(size=3); f=B['srp']*fd/np.linalg.norm(fd)
    v=np.sqrt(GM/r0); x=np.r_[Q@[r0,0,0],Q@[0,v,0]]
    base=S.accel_fn(GM,R,B['C20'],B['C22'],ws,f); thr=np.zeros(3)
    def rhs(t,x): d=base(t,x); d[3:]+=thr; return d
    if variant=='GROUND': gm_hat=GM*1.0015; j2_hat=1.2*J2T; phase='wait'
    else:                 gm_hat=GM*1.0015; j2_hat=1.2*J2T; phase='transfer'
    lo,hi=gm_hat*0.998,gm_hat*1.002; recent=[]; pairs=[]
    p_ref=r1; up=r1>r0
    buf_r=[];buf_v=[];buf_a=[]; dv=0.0; t=0.0; pulse=0; firing_now=False; basis=None; u_star=0.0; armed_until=0.0; pulse_sign=0
    truth=[]; log=[]; orbit_r=[]
    N=int(days*86400/DT)
    for k in range(N):
        rm=x[:3]+rng.normal(0,B['sr'],3); vm=x[3:]+rng.normal(0,B['sv'],3)
        buf_r.append(rm); buf_v.append(vm); buf_a.append(thr.copy()+f)
        rnow=np.linalg.norm(rm); P=2*np.pi*np.sqrt(rnow**3/gm_hat); n_win=max(20,int(P/4/DT))
        if len(buf_r)>=n_win:
            Rw=np.array(buf_r); Vw=np.array(buf_v); Aw=np.array(buf_a); n=len(Rw); m=n//2
            dv_=np.zeros((n,3)); dr_=np.zeros((n,3))
            for j in range(m+1,n): dv_[j]=dv_[j-1]+0.5*(Aw[j]+Aw[j-1])*DT; dr_[j]=dr_[j-1]+0.5*(dv_[j]+dv_[j-1])*DT
            for j in range(m-1,-1,-1): dv_[j]=dv_[j+1]-0.5*(Aw[j]+Aw[j+1])*DT; dr_[j]=dr_[j+1]-0.5*(dv_[j]+dv_[j+1])*DT
            rc=Rw-dr_; vc=Vw-dv_; rn=np.linalg.norm(rc,axis=1)
            mu_eff=S.rel(rc,vc); g=1.5*np.mean((R/rn)**2*(1-3*(rc[:,2]/rn)**2))
            h2=np.mean(np.linalg.norm(np.cross(rc,vc),axis=1))**2
            if variant!='GROUND':
                pairs.append((mu_eff,g)); pairs=pairs[-200:]
                if variant=='ROOP_J2' and len(pairs)>=8:
                    P_=np.array(pairs); gs=P_[:,1]
                    if gs.max()-gs.min()>0.2*abs(gs.mean()):
                        A_=np.vstack([np.ones(len(P_)),gs]).T; c,*_=np.linalg.lstsq(A_,P_[:,0],rcond=None)
                        j2_hat=float(np.clip(c[1]/c[0],0.5*J2T*1.2,2.0*J2T*1.2))
                mu_win=mu_eff/(1+j2_hat*g)
                # DeltaF_geom mass correction inside the adjustable band (band may re-open to 3 sigma if data disagree)
                recent.append(mu_win); recent=recent[-12:]
                centre=float(np.median(recent))
                sig=1.4826*np.median(np.abs(np.array(recent)-centre))/np.sqrt(len(recent)) if len(recent)>=3 else 0.002*centre
                half=max(3*sig,0.002*centre); lo,hi=centre-half,centre+half; gm_hat=centre
            elif phase=='wait' and t>=T_ground*86400:
                gm_hat=GM*1.001; j2_hat=J2T*1.01; phase='transfer'
            # DeltaF thrust command with current gm_hat (J2-aware reference: mean orbit at r1)
            p_est=h2/gm_hat/(1+j2_hat*g); dlt=p_est/p_ref-1; vbar=np.sqrt(gm_hat/p_ref)
            orbit_r.append(rc); orbit_r=orbit_r[-4:]; Ro=np.vstack(orbit_r)        # last ~1 orbit of corrected positions
            k_=np.cross(rc,vc).mean(0); k_/=np.linalg.norm(k_); i1=Ro[0]-k_*(Ro[0]@k_); i1/=np.linalg.norm(i1); j1=np.cross(k_,i1)
            th=np.arctan2(Ro@j1,Ro@i1); A_=np.vstack([np.ones_like(th),np.cos(th),np.sin(th)]).T
            cf,*_=np.linalg.lstsq(A_,1/np.linalg.norm(Ro,axis=1),rcond=None); ex,ey=cf[1]/cf[0],cf[2]/cf[0]; e_mag=np.hypot(ex,ey)
            cap=min(0.15,AT*P/8)                                  # pulse limited to 1/8 orbit of arc
            if phase in('transfer','damp'):
                if e_mag>E_DAMP: phase='damp'
                elif phase=='damp' and e_mag<E_DAMP/2: phase='transfer'
            if phase=='transfer' and ((not up and dlt<0.002) or (up and dlt>-0.002)): phase='hold'
            if phase=='damp' and pulse<=0:
                pulse=int(min(vbar*e_mag/4,cap)/AT/DT)+1; pulse_sign=1.0
                u_star=np.arctan2(-ey,-ex); basis=(i1,j1); armed_until=t+P
            if phase=='hold':
                dv_req=-(vbar/2)*dlt
                if abs(dlt)>0.003 or e_mag>0.03:
                    dvp=abs(dv_req) if abs(dlt)>0.003 else vbar*e_mag/4
                    pulse=int(min(dvp,cap)/AT/DT)+1
                    pulse_sign=np.sign(dv_req) if abs(dlt)>0.003 else 1.0
                    u_star=np.arctan2(-ey,-ex) if pulse_sign>0 else np.arctan2(ey,ex)
                    basis=(i1,j1); armed_until=t+P
                else: pulse=0; pulse_sign=0
            log.append((t/86400,gm_hat,lo,hi,j2_hat,phase,float(dlt),float(np.linalg.norm(x[:3])),float(g)))
            buf_r=[];buf_v=[];buf_a=[]
        vv=x[3:]; vhat=vv/np.linalg.norm(vv)
        if phase=='transfer': thr=(AT if up else -AT)*vhat
        elif phase in('hold','damp') and pulse>0:
            ang=np.arctan2(x[:3]@basis[1],x[:3]@basis[0]); arc=pulse*DT*2*np.pi/P
            if abs((ang-(u_star-arc/2)+np.pi)%(2*np.pi)-np.pi)<np.radians(4) or firing_now or t>armed_until:
                firing_now=True; thr=pulse_sign*AT*vhat; pulse-=1
                if pulse<=0: firing_now=False
            else: thr=np.zeros(3)
        else: thr=np.zeros(3)
        dv+=np.linalg.norm(thr)*DT
        x=rk4(rhs,t,x,DT); t+=DT
        if k%15==0: truth.append((t/86400,np.linalg.norm(x[:3]),dv))
    tr=np.array(truth)
    Pb=2*np.pi*np.sqrt(r1**3/GM)/86400; w=max(1,int(round(Pb/(15*DT/86400))))
    ma=np.convolve(tr[:,1],np.ones(w)/w,mode='valid'); tol=0.01*r1
    ok=np.where(np.abs(ma-r1)<tol)[0]; ia=ok[0]+w-1 if len(ok) else None
    t_arr=None if ia is None else float(tr[ia,0]); dv_arr=None if ia is None else float(tr[ia,2])
    last=ma[-int(2/(15*DT/86400)):]
    return dict(variant=variant,seed=seed,r0_km=r0/1e3,r1_km=r1/1e3,t_arrive_d=t_arr,dv_at_arrival=dv_arr,
                r_err_km=float(last.mean()/1e3-r1/1e3),gm_err_pct=float(100*(gm_hat/GM-1)),j2_err_pct=float(100*(j2_hat/J2T-1))), log, tr[::4].tolist()

if __name__=='__main__':
    leg=sys.argv[1]
    r0,r1,days={'B1_D':(413.2e3,188.2e3,25.0),'D_C':(188.2e3,303.2e3,20.0)}[leg]
    out={}
    for var in ['GROUND','ROOP_FIX','ROOP_J2']:
        res=[]
        for s in range(3):
            r,log,tr=run(var,r0,r1,days,seed=s); res.append(r)
            if s==0: out[var+'_log']=log; out[var+'_truth']=tr
            print(json.dumps(r),flush=True)
        out[var]=res
    json.dump(out,open(f'psyche_roop2_{leg}.json','w'))
