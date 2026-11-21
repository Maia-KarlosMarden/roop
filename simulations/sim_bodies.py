"""
Relational mass estimator on real-parameter Bennu and Psyche.
Gravity: point mass + rotating degree-2 field (C20, C22, unnormalized, ref radius R), spin about z.
SRP: constant inertial acceleration (assumed magnitude per body).
Estimators:
  REL    : M = <h>^2 / (G p_geom), SRP removed to first order (r - f dt^2/2, v - f dt)
  OD_PM  : batch LSQ, mu + 6 states, point-mass dynamics (+SRP known)
  OD_J2  : batch LSQ, mu + C20 + C22 + 6 states, full degree-2 rotating field (+SRP known, spin known)
Usage: python sim_bodies.py <body> <orbit_label> <arcs_in_orbits comma> <trials>
"""
import numpy as np, json, sys, time
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
G=6.6743e-11
BODIES={
 # GM [m^3/s^2], R [m], C20, C22, spin period [h], SRP [m/s^2], pos noise [m], vel noise [m/s], prior GM
 'bennu': dict(GM=4.892, R=244.0, C20=-0.0334, C22=0.00276, Pspin=4.2978, srp=1e-7, sr=2.0, sv=1e-4, prior=5.2,
               orbits={'PS_flyby':('flyby',7250.0), 'OrbA':('circ',1750.0), 'OrbB':('circ',1000.0)}),
 'psyche': dict(GM=1.601e9, R=113.2e3, C20=-0.1027, C22=0.0190, Pspin=4.195948, srp=2.6e-8, sr=50.0, sv=1e-4, prior=1.824e9,
               orbits={'OrbA':('circ',113.2e3+700e3),'OrbB':('circ',113.2e3+300e3),'OrbC':('circ',113.2e3+190e3),'OrbD':('circ',113.2e3+75e3)}),
}
def accel_fn(mu,R,C20,C22,wspin,f):
    def rhs(t,x):
        r=x[:3]; lam=wspin*t; c,s=np.cos(lam),np.sin(lam)
        xb= c*r[0]+s*r[1]; yb=-s*r[0]+c*r[1]; zb=r[2]
        r2=xb*xb+yb*yb+zb*zb; rr=np.sqrt(r2); r5=r2*r2*rr; r7=r5*r2
        q=3*zb*zb-r2; d=xb*xb-yb*yb; k=mu*R*R
        ax=-mu*xb/(r2*rr)+k*(C20*(-xb/r5-2.5*xb*q/r7)+3*C22*(2*xb/r5-5*xb*d/r7))
        ay=-mu*yb/(r2*rr)+k*(C20*(-yb/r5-2.5*yb*q/r7)+3*C22*(-2*yb/r5-5*yb*d/r7))
        az=-mu*zb/(r2*rr)+k*(C20*(2*zb/r5-2.5*zb*q/r7)+3*C22*(-5*zb*d/r7))
        aix=c*ax-s*ay; aiy=s*ax+c*ay
        return np.array([x[3],x[4],x[5],aix+f[0],aiy+f[1],az+f[2]])
    return rhs
def prop(x0,t,mu,R,C20,C22,wspin,f):
    s=solve_ivp(accel_fn(mu,R,C20,C22,wspin,f),(t[0],t[-1]),x0,t_eval=t,rtol=1e-10,atol=1e-9*max(1,np.linalg.norm(x0[:3])/1e3),method='DOP853')
    return s.y.T
def rel(rm,vm):
    Rr=np.linalg.norm(rm,axis=1); hv=np.cross(rm,vm); h=np.linalg.norm(hv,axis=1)
    k=hv.mean(0); k/=np.linalg.norm(k); i1=rm[0]-k*(rm[0]@k); i1/=np.linalg.norm(i1); j1=np.cross(k,i1)
    th=np.arctan2(rm@j1,rm@i1); A=np.vstack([np.ones_like(th),np.cos(th),np.sin(th)]).T
    c,*_=np.linalg.lstsq(A,1/Rr,rcond=None)
    return h.mean()**2*c[0]   # = GM estimate
def case(body,orb,arc_orbits,trials,seed=0,N=60):
    B=BODIES[body]; kind,rad=B['orbits'][orb]; GM,R=B['GM'],B['R']; ws=2*np.pi/(B['Pspin']*3600)
    rng=np.random.default_rng(seed); out={'REL':[],'OD_PM':[],'OD_J2':[]}; tim={'REL':[],'OD_PM':[],'OD_J2':[]}
    for tr in range(trials):
        inc=rng.uniform(np.radians(70),np.radians(110)); node=rng.uniform(0,2*np.pi)
        Rx=np.array([[1,0,0],[0,np.cos(inc),-np.sin(inc)],[0,np.sin(inc),np.cos(inc)]])
        Rz=np.array([[np.cos(node),-np.sin(node),0],[np.sin(node),np.cos(node),0],[0,0,1]]); Q=Rz@Rx
        fdir=rng.normal(size=3); fdir/=np.linalg.norm(fdir); f=B['srp']*fdir
        if kind=='circ':
            T=2*np.pi*np.sqrt(rad**3/GM); u0=rng.uniform(0,2*np.pi); v=np.sqrt(GM/rad)
            x0=np.r_[Q@[rad*np.cos(u0),rad*np.sin(u0),0], Q@[-v*np.sin(u0),v*np.cos(u0),0]]
            t=np.linspace(0,arc_orbits*T,N)
        else:  # flyby, close approach rad, v_inf 0.1 m/s; arc_orbits = window in hours
            vp=np.sqrt(0.01+2*GM/rad); xca=np.r_[Q@[rad,0,0],Q@[0,vp,0]]; half=arc_orbits*3600/2
            x0=prop(xca,np.array([0,-half]),GM,R,B['C20'],B['C22'],ws,np.zeros(3))[-1]
            t=np.linspace(0,2*half,N)
        tr_=prop(x0,t,GM,R,B['C20'],B['C22'],ws,f)
        rm=tr_[:,:3]+rng.normal(0,B['sr'],(N,3)); vm=tr_[:,3:]+rng.normal(0,B['sv'],(N,3)); dt=(t-t.mean())[:,None]
        a=time.perf_counter(); g=rel(rm-0.5*f*dt**2, vm-f*dt); tim['REL'].append(time.perf_counter()-a); out['REL'].append(g/GM-1)
        sc=np.r_[np.full(3,B['sr']*5),np.full(3,B['sv']*5)]
        def res_pm(x): return ((prop(x[:6],t,x[6]*GM,R,0,0,ws,f)[:,:3]-rm)/B['sr']).ravel()
        a=time.perf_counter()
        s1=least_squares(res_pm,np.r_[rm[0],vm[0],B['prior']/GM],x_scale=np.r_[sc,0.05],max_nfev=40)
        tim['OD_PM'].append(time.perf_counter()-a); out['OD_PM'].append(s1.x[6]-1)
        def res_j2(x): return ((prop(x[:6],t,x[6]*GM,R,x[7],x[8],ws,f)[:,:3]-rm)/B['sr']).ravel()
        a=time.perf_counter()
        s2=least_squares(res_j2,np.r_[s1.x,0.0,0.0],x_scale=np.r_[sc,0.05,0.05,0.01],max_nfev=40)
        tim['OD_J2'].append(time.perf_counter()-a); out['OD_J2'].append(s2.x[6]-1)
    st={k:dict(bias=float(100*np.mean(v)),sd=float(100*np.std(v)),rms=float(100*np.sqrt(np.mean(np.square(v))))) for k,v in out.items()}
    for k in tim: st[k]['ms']=float(1e3*np.median(tim[k]))
    return st
if __name__=='__main__':
    body,orb=sys.argv[1],sys.argv[2]; arcs=[float(a) for a in sys.argv[3].split(',')]; trials=int(sys.argv[4])
    B=BODIES[body]; kind,rad=B['orbits'][orb]
    T_h=(2*np.pi*np.sqrt(rad**3/B['GM'])/3600) if kind=='circ' else None
    R={}
    for a in arcs:
        st=case(body,orb,a,trials); R[a]=st
        lab=f'{a} orb ({a*T_h:.1f} h)' if T_h else f'{a} h'
        print(json.dumps(dict(body=body,orbit=orb,arc=lab,T_h=T_h,**{k:{kk:round(vv,4) for kk,vv in v.items()} for k,v in st.items()})),flush=True)
    json.dump(dict(T_h=T_h,res=R),open(f'bd_{body}_{orb}.json','w'))
