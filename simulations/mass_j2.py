import numpy as np, sys
sys.path.insert(0,'.')
import sim_mass as S
from scipy.integrate import solve_ivp
G=S.G; MU=S.MU; Rb=1.13e5
def f(t,x,J2):
    r=x[:3]; R=np.linalg.norm(r); z2=(r[2]/R)**2
    a=-MU*r/R**3; k=1.5*J2*MU*Rb**2/R**5
    a=a+k*r*np.array([5*z2-1,5*z2-1,5*z2-3]); return np.r_[x[3:],a]
rng=np.random.default_rng(1)
for J2 in [0.0,0.01,0.05,0.1]:
  for e in [0.02,0.1]:
    out={k:[] for k in['STD','ENH','REL','OD']}
    for tr in range(30):
        a=3e5; M0=rng.uniform(0,2*np.pi); T=2*np.pi*np.sqrt(a**3/MU)
        r0,v0=S.kep2rv(a,e,0.3,0.2,0.5,np.array([M0]),MU,np.array([0.0]))
        t=np.linspace(0,T,60)
        sol=solve_ivp(f,(0,T),np.r_[r0[0],v0[0]],t_eval=t,args=(J2,),rtol=1e-11,atol=1e-6)
        r=sol.y[:3].T; v=sol.y[3:].T
        rm=r+rng.normal(0,S.SR,r.shape); vm=v+rng.normal(0,S.SV,v.shape)
        est=S.estimators(rm,vm,t,[a,e,0.3,0.2,0.5,M0,MU/1e9])
        for k in out: out[k].append(est[k]/S.Mtrue-1)
    print(f"J2={J2} e={e} | "+" | ".join(f"{k}: {100*np.mean(x):+.3f}±{100*np.std(x):.3f}%" for k,x in out.items()),flush=True)
