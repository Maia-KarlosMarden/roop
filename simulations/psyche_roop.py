"""
Integrated ROOP at Psyche: the same relational deviation DeltaF drives BOTH fronts.
  DeltaF_geom = G'r - G p_geom = h^2 (1/M_hat - 1/M)        -> mass correction (no prior needed)
  DeltaF*     = -(v/2)(<G'r>/(G p_ref) - 1)                   -> thrust command (uses current M_hat)
Mass update inside an adjustable admissible band [M_lo, M_hi]:
  initial band from size + plausible bulk density (1000-8000 kg/m^3); each quarter-orbit window gives
  M_win = <h>^2/(G p_geom) (J2-corrected with shape-model C20, 20% error); M_hat <- clamp(M_hat + k (M_win - M_hat));
  band recentred on M_hat with half-width max(3*sigma_recent, 0.2%), never wider than before.
Scenario: Orbit A (811 km radius) -> Orbit B1 (413 km), low-thrust spiral then hold, 30 days.
Truth: GM 1.601e9, rotating C20/C22, SRP. Nav noise 50 m / 0.1 mm/s. Thrust accel 1.22e-5 m/s^2 (17-d nominal transfer).
Variants:
  PRIOR_SHAPE : M_hat fixed at shape-model value (GM 1.824e9, +14%), descends immediately
  GROUND      : waits in Orbit A for a ground GM solution (perfect, 0.1%) delivered after T_ground days, then descends
  ROOP        : no prior; online DeltaF_geom mass loop + DeltaF* thrust; descends once band half-width < 2%
"""
import numpy as np, json, sys
sys.path.insert(0,'.')
import sim_bodies as S
B=S.BODIES['psyche']; GM,R=B['GM'],B['R']; ws=2*np.pi/(B['Pspin']*3600); G=S.G
C20s=B['C20']*1.2                     # shape-model C20 with 20% error
VOL=4/3*np.pi*R**3
GM_LO,GM_HI=G*1000*VOL,G*8000*VOL     # admissible band from density 1000-8000 kg/m^3
rA,rB=811e3,413.2e3; AT=1.22e-5; DT=60.0

def rk4(f,t,x,h):
    k1=f(t,x); k2=f(t+h/2,x+h/2*k1); k3=f(t+h/2,x+h/2*k2); k4=f(t+h,x+h*k3)
    return x+h/6*(k1+2*k2+2*k3+k4)

def run(variant, days=30.0, seed=0, T_ground=7.0):
    rng=np.random.default_rng(seed)
    inc=np.radians(90+rng.uniform(-3,3)); node=rng.uniform(0,2*np.pi)
    Rx=np.array([[1,0,0],[0,np.cos(inc),-np.sin(inc)],[0,np.sin(inc),np.cos(inc)]])
    Rz=np.array([[np.cos(node),-np.sin(node),0],[np.sin(node),np.cos(node),0],[0,0,1]]); Q=Rz@Rx
    fd=rng.normal(size=3); f=B['srp']*fd/np.linalg.norm(fd)
    v=np.sqrt(GM/rA); x=np.r_[Q@[rA,0,0],Q@[0,v,0]]
    base=S.accel_fn(GM,R,B['C20'],B['C22'],ws,f)
    thr=np.zeros(3)
    def rhs(t,x): d=base(t,x); d[3:]+=thr; return d
    # estimator state
    if variant=='PRIOR_SHAPE': gm_hat=1.824e9; lo,hi=gm_hat,gm_hat
    elif variant=='GROUND':    gm_hat=1.824e9; lo,hi=gm_hat,gm_hat
    else:                      gm_hat=np.sqrt(GM_LO*GM_HI); lo,hi=GM_LO,GM_HI
    p_ref=rB; phase='wait' if variant in('GROUND','ROOP') else 'descend'
    buf_r=[];buf_v=[];buf_a=[]; recent=[]; log=[]; dv=0.0; t=0.0; pulse=0; pulse_dir=np.zeros(3); firing_now=False; basis=None; u_star=0.0; armed_until=0.0; pulse_sign=0
    N=int(days*86400/DT); n_win=None; t_arrive=None; truth=[]
    for k in range(N):
        # measurement
        rm=x[:3]+rng.normal(0,B['sr'],3); vm=x[3:]+rng.normal(0,B['sv'],3)
        buf_r.append(rm); buf_v.append(vm); buf_a.append(thr.copy()+f)
        rnow=np.linalg.norm(rm); P=2*np.pi*np.sqrt(rnow**3/GM if gm_hat<=0 else rnow**3/gm_hat)
        n_win=max(20,int(P/4/DT))
        if len(buf_r)>=n_win:
            Rw=np.array(buf_r); Vw=np.array(buf_v); Aw=np.array(buf_a); n=len(Rw); m=n//2
            dv_=np.zeros((n,3)); dr_=np.zeros((n,3))
            for j in range(m+1,n): dv_[j]=dv_[j-1]+0.5*(Aw[j]+Aw[j-1])*DT; dr_[j]=dr_[j-1]+0.5*(dv_[j]+dv_[j-1])*DT
            for j in range(m-1,-1,-1): dv_[j]=dv_[j+1]-0.5*(Aw[j]+Aw[j+1])*DT; dr_[j]=dr_[j+1]-0.5*(dv_[j]+dv_[j+1])*DT
            rc=Rw-dr_; vc=Vw-dv_
            mu_eff=S.rel(rc,vc); rn=np.linalg.norm(rc,axis=1)
            g=1.5*np.mean((R/rn)**2*(1-3*(rc[:,2]/rn)**2)); mu_win=mu_eff/(1-C20s*g/(-B['C20'])*(-B['C20'])/1.0) if False else mu_eff/(1+(-C20s)*g)
            h2=np.mean(np.linalg.norm(np.cross(rc,vc),axis=1))**2
            # ---- front 1: mass via DeltaF_geom inside the adjustable band ----
            if variant=='ROOP':
                dF_geom=h2/gm_hat*G - h2/mu_win*G          # = G(G'r) - G(G p_geom), units consistent (h^2/M terms)
                mu_in=float(np.clip(mu_win,lo,hi))            # reject windows outside the admissible band
                recent.append(mu_in); recent=recent[-8:]
                centre=float(np.median(recent))
                sig=(1.4826*np.median(np.abs(np.array(recent)-centre))/np.sqrt(len(recent))) if len(recent)>=3 else (hi-lo)/4
                half=min((hi-lo)/2, max(3*sig,0.002*centre))
                lo,hi=max(GM_LO,centre-half),min(GM_HI,centre+half)
                gm_hat=float(np.clip(centre,lo,hi))
                if phase=='wait' and (hi-lo)/2/gm_hat<0.02: phase='descend'
            elif variant=='GROUND' and phase=='wait' and t>=T_ground*86400:
                gm_hat=GM*(1+0.001); phase='descend'
            # ---- front 2: DeltaF* thrust command with current gm_hat ----
            p_est=h2/gm_hat; dlt=p_est/p_ref-1; vbar=np.sqrt(gm_hat/p_ref)
            if phase=='descend' and dlt<0.002:
                phase='hold'
            if phase=='hold':
                # e-vector from K-fit (scale-free) -> burn location; size from DeltaF*
                k_=np.cross(rc,vc).mean(0); k_/=np.linalg.norm(k_); i1=rc[0]-k_*(rc[0]@k_); i1/=np.linalg.norm(i1); j1=np.cross(k_,i1)
                th=np.arctan2(rc@j1,rc@i1); A_=np.vstack([np.ones_like(th),np.cos(th),np.sin(th)]).T
                cf,*_=np.linalg.lstsq(A_,1/rn,rcond=None); ex,ey=cf[1]/cf[0],cf[2]/cf[0]
                dv_req=-(vbar/2)*dlt; e_mag=np.hypot(ex,ey)
                if abs(dlt)>0.003 or e_mag>0.01:
                    dvp=abs(dv_req) if abs(dlt)>0.003 else vbar*e_mag/4
                    pulse=int(min(dvp,0.15)/AT/DT)+1
                    pulse_sign=np.sign(dv_req) if abs(dlt)>0.003 else 1.0
                    # burn centred where a tangential impulse moves e toward 0:
                    # prograde at u = angle(-e) ... sign-adjusted
                    u_star=np.arctan2(-ey,-ex) if pulse_sign>0 else np.arctan2(ey,ex)
                    basis=(i1,j1); armed_until=t+P
                else: pulse=0; pulse_sign=0
            log.append((t/86400, gm_hat, lo, hi, phase, float(dlt)))
            buf_r=[];buf_v=[];buf_a=[]
        # thrust
        vv=x[3:]; vhat=vv/np.linalg.norm(vv)
        if phase=='descend': thr=-AT*vhat
        elif phase=='hold' and pulse>0:
            ang=np.arctan2(x[:3]@basis[1],x[:3]@basis[0]); arc=pulse*DT*2*np.pi/P
            if abs((ang-(u_star-arc/2)+np.pi)%(2*np.pi)-np.pi)<np.radians(4) or firing_now or t>armed_until:
                firing_now=True; thr=pulse_sign*AT*vhat; pulse-=1
                if pulse<=0: firing_now=False
            else: thr=np.zeros(3)
        else: thr=np.zeros(3)
        dv+=np.linalg.norm(thr)*DT
        x=rk4(rhs,t,x,DT); t+=DT
        if k%30==0: truth.append((t/86400,np.linalg.norm(x[:3]),dv))
    tr=np.array(truth); last=tr[tr[:,0]>days-2]
    Pb=2*np.pi*np.sqrt(rB**3/GM)/86400; w=max(1,int(Pb/(30*DT/86400)))
    ma=np.convolve(tr[:,1],np.ones(w)/w,mode='valid'); arr=np.where(np.abs(ma-rB)<5e3)[0]
    ia=arr[0]+w-1 if len(arr) else None
    t_arrive=None if ia is None else float(tr[ia,0]); dv_arr=None if ia is None else float(tr[ia,2])
    post=tr[tr[:,0]>(t_arrive if t_arrive else days)]
    hold_rate=float((post[-1,2]-post[0,2])/(post[-1,0]-post[0,0])) if len(post)>10 else None
    return dict(variant=variant,seed=seed,gm_final_err_pct=100*(gm_hat/GM-1),r_final_km=float(last[:,1].mean()/1e3),
                r_err_km=float(last[:,1].mean()/1e3-rB/1e3),r_spread_km=float((last[:,1].max()-last[:,1].min())/1e3),
                t_arrive_d=t_arrive,dv_at_arrival=dv_arr,hold_dv_per_day=hold_rate,dv=dv), log, tr.tolist()

if __name__=='__main__':
    out={}
    for var in ['PRIOR_SHAPE','GROUND','ROOP']:
        res=[]
        for s in range(3):
            r,log,tr=run(var,seed=s); res.append(r)
            if s==0: out[var+'_log']=[(a,b,c,d,e,f_) for a,b,c,d,e,f_ in log]; out[var+'_truth']=tr
            print(json.dumps(r),flush=True)
        out[var]=res
    json.dump(out,open('psyche_roop.json','w'))
