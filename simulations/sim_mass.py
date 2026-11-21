"""
Mass estimation around a small body (Psyche-like, point mass) from noisy tracking.
Estimators:
  STD : time-average variant  M = mean(w^2 r^3)/G
  ENH : rms-radius variant  M = mean(w^2) * r_rms^3 / (Ke G), e from (rmax-rmin)/(rmax+rmin)
  REL : relational invariant G'r = Gp with p from a purely geometric conic fit (positions only)
        and h from kinematics:  M = <h>^2 / (G p_geom)
  OD  : classical batch orbit determination, 7 parameters (6 elements + mu), positions only
"""
import numpy as np, json
from scipy.optimize import least_squares
G=6.6743e-11; Mtrue=2.3e19; MU=G*Mtrue
SR=10.0; SV=1e-3

def kep2rv(a,e,inc,O,w,M0,mu,t):
    n=np.sqrt(mu/a**3); M=M0+n*t; E=M.copy()
    for _ in range(30): E=E-(E-e*np.sin(E)-M)/(1-e*np.cos(E))
    cE,sE=np.cos(E),np.sin(E); b=np.sqrt(1-e*e)
    xp=a*(cE-e); yp=a*b*sE; r=a*(1-e*cE)
    vxp=-np.sqrt(mu*a)/r*sE; vyp=np.sqrt(mu*a)/r*b*cE
    cO,sO,cw,sw,ci,si=np.cos(O),np.sin(O),np.cos(w),np.sin(w),np.cos(inc),np.sin(inc)
    R=np.array([[cO*cw-sO*sw*ci,-cO*sw-sO*cw*ci],[sO*cw+cO*sw*ci,-sO*sw+cO*cw*ci],[sw*si,cw*si]])
    return (R@np.vstack([xp,yp])).T, (R@np.vstack([vxp,vyp])).T

def estimators(rm,vm,t,guess):
    R=np.linalg.norm(rm,axis=1); hv=np.cross(rm,vm); h=np.linalg.norm(hv,axis=1); w=h/R**2
    std=np.mean(w**2*R**3)/G
    e_est=(R.max()-R.min())/(R.max()+R.min()); Ke=(1-e_est**2)**0.75
    enh=np.mean(w**2)*np.sqrt(np.mean(R**2))**3/(Ke*G)
    # geometric conic fit in the orbit plane
    k=hv.mean(0); k/=np.linalg.norm(k); i1=rm[0]-k*(rm[0]@k); i1/=np.linalg.norm(i1); j1=np.cross(k,i1)
    th=np.arctan2(rm@j1, rm@i1)
    A=np.vstack([np.ones_like(th),np.cos(th),np.sin(th)]).T
    c,*_=np.linalg.lstsq(A,1/R,rcond=None); p=1/c[0]
    rel=h.mean()**2/(G*p)
    # OD
    def res(x):
        r,_=kep2rv(x[0],x[1],x[2],x[3],x[4],x[5],x[6]*1e9,t); return (r-rm).ravel()/SR
    x0=list(guess); x0[6]*=1.2   # 20% a-priori mass error
    sol=least_squares(res,x0,x_scale=[1e4,0.01,0.01,0.01,0.01,0.01,0.1],max_nfev=200)
    od=sol.x[6]*1e9/G
    return dict(STD=std,ENH=enh,REL=rel,OD=od)

def mc(e,frac,N=60,trials=150,seed=0):
    rng=np.random.default_rng(seed)
    a=3.0e5; inc,O,w=0.3,0.2,0.5; T=2*np.pi*np.sqrt(a**3/MU)
    out={k:[] for k in['STD','ENH','REL','OD']}
    for tr in range(trials):
        M0=rng.uniform(0,2*np.pi)
        t=np.linspace(0,frac*T,N)
        r,v=kep2rv(a,e,inc,O,w,M0,MU,t)
        rm=r+rng.normal(0,SR,r.shape); vm=v+rng.normal(0,SV,v.shape)
        est=estimators(rm,vm,t,[a,e,inc,O,w,M0,MU/1e9])
        for k in out: out[k].append(est[k]/Mtrue-1)
    return {k:(100*np.mean(x),100*np.std(x)) for k,x in out.items()}

if __name__=='__main__':
    R={}
    for e in [0.02,0.1,0.2]:
        for frac in [0.25,0.5,1.0]:
            r=mc(e,frac); R[f"{e}_{frac}"]=r
            print(f"e={e} arc={frac} orb | "+" | ".join(f"{k}: bias {b:+.3f}% sd {s:.3f}%" for k,(b,s) in r.items()),flush=True)
    json.dump(R,open('mass_results.json','w'))
