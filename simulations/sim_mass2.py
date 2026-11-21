"""
Objective 2 - relational mass estimator on small bodies (Bennu-like), orbit and hyperbolic flyby.
GM = 4.89 m^3/s^2 (Bennu, OSIRIS-REx). Optional constant SRP acceleration on the spacecraft.
Noise: position 2 m (optical/landmark nav), velocity 0.1 mm/s, per axis, 1 sigma.
Estimators:
  REL      : M = <h>^2/(G p_geom), conic fit 1/r = A + B cos th + C sin th (no mu, no time, no propagation)
  REL_SRP  : REL after first-order removal of a known SRP acceleration f (r - f dt^2/2, v - f dt)
  OD       : batch least squares, mu + 6 initial-state parameters, numerical propagation with SRP known
"""
import numpy as np, json, sys
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
G=6.6743e-11; GM=4.89; Mtrue=GM/G
SR=2.0; SV=1e-4

def prop(x0,t,mu,f):
    def rhs(_,x):
        r=x[:3]; R=np.linalg.norm(r); return np.r_[x[3:], -mu*r/R**3+f]
    s=solve_ivp(rhs,(t[0],t[-1]),x0,t_eval=t,rtol=1e-10,atol=1e-9,method='DOP853')
    return s.y.T

def rel(rm,vm):
    R=np.linalg.norm(rm,axis=1); hv=np.cross(rm,vm); h=np.linalg.norm(hv,axis=1)
    k=hv.mean(0); k/=np.linalg.norm(k); i1=rm[0]-k*(rm[0]@k); i1/=np.linalg.norm(i1); j1=np.cross(k,i1)
    th=np.arctan2(rm@j1,rm@i1)
    A=np.vstack([np.ones_like(th),np.cos(th),np.sin(th)]).T
    c,*_=np.linalg.lstsq(A,1/R,rcond=None)
    return h.mean()**2/(G/c[0]) if c[0]>0 else np.nan   # = <h>^2 * A / G

def run_case(kind, window_h, f_srp, trials=40, seed=0, N=60):
    rng=np.random.default_rng(seed); out={'REL':[],'REL_SRP':[],'OD':[]}
    f=np.array(f_srp)
    for tr in range(trials):
        if kind=='orbit':
            a=1750.0; e=0.05; th0=rng.uniform(0,2*np.pi)
            p=a*(1-e*e); r0=p/(1+e*np.cos(th0))
            vr=np.sqrt(GM/p)*e*np.sin(th0); vt=np.sqrt(GM/p)*(1+e*np.cos(th0))
            x0=np.array([r0*np.cos(th0),r0*np.sin(th0),0, vr*np.cos(th0)-vt*np.sin(th0), vr*np.sin(th0)+vt*np.cos(th0),0])
            inc=rng.uniform(0.3,1.2); Rx=np.array([[1,0,0],[0,np.cos(inc),-np.sin(inc)],[0,np.sin(inc),np.cos(inc)]])
            x0=np.r_[Rx@x0[:3],Rx@x0[3:]]
            t=np.linspace(0,window_h*3600,N)
        else:  # hyperbolic flyby, close approach 7.25 km at t=0, v_inf=0.1 m/s; window centred on CA
            rp=7250.0; vinf=0.1; vp=np.sqrt(vinf**2+2*GM/rp)
            x_ca=np.array([rp,0,0,0,vp,0]); half=window_h*3600/2
            back=prop(x_ca,np.array([0,-half]),GM,np.zeros(3))[-1]
            x0=back; t=np.linspace(0,2*half,N)
        tr_=prop(x0,t,GM,f)
        rm=tr_[:,:3]+rng.normal(0,SR,(N,3)); vm=tr_[:,3:]+rng.normal(0,SV,(N,3))
        out['REL'].append(rel(rm,vm)/Mtrue-1)
        dt=(t-t.mean())[:,None]
        out['REL_SRP'].append(rel(rm-0.5*f*dt**2, vm-f*dt)/Mtrue-1)
        def res(x):
            return ((prop(x[:6],t,x[6],f)[:,:3]-rm)/SR).ravel()
        g=np.r_[rm[0],vm[0],GM*1.2]
        sol=least_squares(res,g,x_scale=np.r_[np.full(3,10.),np.full(3,1e-3),0.5],max_nfev=60)
        out['OD'].append(sol.x[6]/GM-1)
    return {k:(float(100*np.nanmean(v)),float(100*np.nanstd(v))) for k,v in out.items()}

if __name__=='__main__':
    kind=sys.argv[1]; srp=float(sys.argv[2]); wins=[float(w) for w in sys.argv[3].split(',')]
    R={}
    for w in wins:
        r=run_case(kind,w,[srp,0.3*srp,0.0]); R[w]=r
        print(json.dumps(dict(kind=kind,srp=srp,window_h=w,**{k:[round(x,3) for x in v] for k,v in r.items()})),flush=True)
    json.dump(R,open(f'm2_{kind}_{srp:g}.json','w'))
