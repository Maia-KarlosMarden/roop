"""REL with closed-form J2 correction on Psyche; two-altitude joint (mu, J2) solve."""
import sys, numpy as np, json; sys.path.insert(0,'.')
import sim_bodies as S
B=S.BODIES['psyche']; GM,R=B['GM'],B['R']; J2=-B['C20']; ws=2*np.pi/(B['Pspin']*3600)
def run(orb,arc,trials=25,seed=7,N=60):
    rng=np.random.default_rng(seed); rad=B['orbits'][orb][1]; T=2*np.pi*np.sqrt(rad**3/GM)
    out=[]
    for tr in range(trials):
        inc=rng.uniform(np.radians(70),np.radians(110)); node=rng.uniform(0,2*np.pi)
        Rx=np.array([[1,0,0],[0,np.cos(inc),-np.sin(inc)],[0,np.sin(inc),np.cos(inc)]])
        Rz=np.array([[np.cos(node),-np.sin(node),0],[np.sin(node),np.cos(node),0],[0,0,1]]); Q=Rz@Rx
        fd=rng.normal(size=3); f=B['srp']*fd/np.linalg.norm(fd)
        u0=rng.uniform(0,2*np.pi); v=np.sqrt(GM/rad)
        x0=np.r_[Q@[rad*np.cos(u0),rad*np.sin(u0),0],Q@[-v*np.sin(u0),v*np.cos(u0),0]]
        t=np.linspace(0,arc*T,N); tr_=S.prop(x0,t,GM,R,B['C20'],B['C22'],ws,f)
        rm=tr_[:,:3]+rng.normal(0,B['sr'],(N,3)); vm=tr_[:,3:]+rng.normal(0,B['sv'],(N,3)); dt=(t-t.mean())[:,None]
        rc=rm-0.5*f*dt**2; vc=vm-f*dt
        mu_eff=S.rel(rc,vc)
        hv=np.cross(rc,vc).mean(0); i_est=np.arccos(abs(hv[2])/np.linalg.norm(hv))   # inclination wrt spin axis
        rbar=np.linalg.norm(rc,axis=1).mean()
        g=1.5*(R/rbar)**2*(1-1.5*np.sin(i_est)**2)
        out.append((mu_eff,g))
    return out
res={}
for orb in ['OrbA','OrbB','OrbC','OrbD']:
    o=run(orb,1.0); mu=np.array([x[0] for x in o]); g=np.array([x[1] for x in o])
    raw=100*(mu/GM-1); corr=100*(mu/(1+J2*g)/GM-1); corr_shape=100*(mu/(1+1.2*J2*g)/GM-1)
    res[orb]=(mu,g)
    print(f"{orb}: REL raw {raw.mean():+.3f}±{raw.std():.3f}% | J2 exact {corr.mean():+.3f}±{corr.std():.3f}% | J2 +20% err {corr_shape.mean():+.3f}±{corr_shape.std():.3f}%")
# two-altitude joint solve, closed form: mu_eff = mu + (mu J2) g, pair trial-by-trial
for a,b in [('OrbA','OrbC'),('OrbB','OrbD'),('OrbA','OrbD')]:
    ma,ga=res[a]; mb,gb=res[b]; est_mu=[];est_j2=[]
    for k in range(len(ma)):
        A=np.array([[1,ga[k]],[1,gb[k]]]); x=np.linalg.solve(A,[ma[k],mb[k]]); est_mu.append(x[0]); est_j2.append(x[1]/x[0])
    em=100*(np.array(est_mu)/GM-1); ej=100*(np.array(est_j2)/J2-1)
    print(f"2-alt {a}+{b}: mu {em.mean():+.4f}±{em.std():.4f}%  J2 {ej.mean():+.2f}±{ej.std():.2f}%")
