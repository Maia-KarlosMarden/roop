"""
Onboard relational estimation of GM and J2 during Psyche's low-thrust spiral Orbit A -> Orbit B1.
Truth: point mass + rotating degree-2 field (C20, C22), constant SRP, continuous anti-velocity thrust
sized so that 811 km -> 413 km radius takes ~17 days (Transfer to B1 = 17 d in the mission plan).
Onboard: noisy position/velocity (50 m, 0.1 mm/s); commanded thrust known but with 1% magnitude error.
Per quarter-orbit window: remove known accelerations (thrust + SRP) by double integration, REL -> mu_eff(r).
Running linear regression mu_eff = mu + (mu J2) g(r, i)  ->  GM and J2 as a function of transfer day.
"""
import numpy as np, json, sys
sys.path.insert(0,'.')
import sim_bodies as S
B=S.BODIES['psyche']; GM,R=B['GM'],B['R']; ws=2*np.pi/(B['Pspin']*3600); J2=-B['C20']
def run(seed, thrust_err=0.01, days=17.0, dt=120.0):
    rng=np.random.default_rng(seed)
    rA=811e3; rB=413e3; dv=np.sqrt(GM/rB)-np.sqrt(GM/rA); aT=dv/(days*86400)
    inc=np.radians(90+rng.uniform(-3,3)); node=rng.uniform(0,2*np.pi)
    Rx=np.array([[1,0,0],[0,np.cos(inc),-np.sin(inc)],[0,np.sin(inc),np.cos(inc)]])
    Rz=np.array([[np.cos(node),-np.sin(node),0],[np.sin(node),np.cos(node),0],[0,0,1]]); Q=Rz@Rx
    fd=rng.normal(size=3); f=B['srp']*fd/np.linalg.norm(fd)
    v=np.sqrt(GM/rA); x0=np.r_[Q@[rA,0,0],Q@[0,v,0]]
    base=S.accel_fn(GM,R,B['C20'],B['C22'],ws,f)
    def rhs(t,x):
        d=base(t,x); vv=x[3:]; d[3:]+=-aT*vv/np.linalg.norm(vv); return d
    from scipy.integrate import solve_ivp
    t=np.arange(0,days*86400,dt)
    sol=solve_ivp(rhs,(0,t[-1]),x0,t_eval=t,rtol=1e-10,atol=1e-6,method='DOP853'); X=sol.y.T
    rm=X[:,:3]+rng.normal(0,B['sr'],(len(t),3)); vm=X[:,3:]+rng.normal(0,B['sv'],(len(t),3))
    # onboard known accel: commanded thrust (with magnitude error) along measured velocity + SRP model
    an=-(aT*(1+thrust_err))*vm/np.linalg.norm(vm,axis=1)[:,None]+f
    wins=[]; k=0
    while k<len(t)-10:
        r_now=np.linalg.norm(rm[k]); P=2*np.pi*np.sqrt(r_now**3/GM); n=max(20,int(P/4/dt))
        if k+n>=len(t): break
        sl=slice(k,k+n); tt=t[sl]; a=an[sl]; m=n//2
        # integrate known accel from window midpoint
        dv_=np.zeros((n,3)); dr_=np.zeros((n,3))
        for j in range(m+1,n): dv_[j]=dv_[j-1]+0.5*(a[j]+a[j-1])*dt; dr_[j]=dr_[j-1]+0.5*(dv_[j]+dv_[j-1])*dt
        for j in range(m-1,-1,-1): dv_[j]=dv_[j+1]-0.5*(a[j]+a[j+1])*dt; dr_[j]=dr_[j+1]-0.5*(dv_[j]+dv_[j+1])*dt
        rc=rm[sl]-dr_; vc=vm[sl]-dv_
        mu_eff=S.rel(rc,vc); hv=np.cross(rc,vc).mean(0); ii=np.arccos(abs(hv[2])/np.linalg.norm(hv))
        rn=np.linalg.norm(rc,axis=1); rbar=rn.mean(); g=1.5*np.mean((R/rn)**2*(1-3*(rc[:,2]/rn)**2))
        wins.append((tt.mean()/86400, rbar, mu_eff, g)); k+=n
    W=np.array(wins)
    # running regression (weighted equally), report at day marks
    out=[]
    for day in [1,2,4,8,12,17]:
        sel=W[:,0]<=day
        if sel.sum()<3: continue
        A=np.vstack([np.ones(sel.sum()),W[sel,3]]).T; c,*_=np.linalg.lstsq(A,W[sel,2],rcond=None)
        out.append((day, int(sel.sum()), W[sel,1][-1]/1e3-R/1e3, 100*(c[0]/GM-1), 100*(c[1]/c[0]/J2-1)))
    return out
if __name__=='__main__':
    allr={}
    for terr in [0.0,0.01,0.03]:
        rows=[run(s,terr) for s in range(8)]
        print(f'--- thrust magnitude knowledge error {terr*100:.0f}% (8 runs) ---')
        res=[]
        for i in range(len(rows[0])):
            day=rows[0][i][0]; alt=np.mean([r[i][2] for r in rows]); nw=rows[0][i][1]
            gm=np.array([r[i][3] for r in rows]); j2=np.array([r[i][4] for r in rows])
            print(f'day {day:4.0f}  windows {nw:3d}  alt {alt:6.0f} km  GM {gm.mean():+.4f}±{gm.std():.4f}%  J2 {j2.mean():+.2f}±{j2.std():.2f}%')
            res.append(dict(day=day,alt_km=alt,gm_bias=gm.mean(),gm_sd=gm.std(),j2_bias=j2.mean(),j2_sd=j2.std()))
        allr[terr]=res
    json.dump(allr,open('psyche_spiral.json','w'))
