"""
Relational mass estimator applied to the REAL reconstructed OSIRIS-REx trajectory around Bennu
(JPL Horizons, spacecraft -64 relative to Bennu centre, ecliptic J2000, Orbital A, 2019-01-03..07).
Reference value: GM = 4.892 +/- 0.006 m^3/s^2 (Scheeres et al. 2019, Nature Astronomy).

Steps (all closed-form):
 1. Orbit plane by least squares; out-of-plane offset d (SRP pushes a terminator orbit off the centre).
 2. In-plane relational estimate: mu_fit = <h_in>^2 * A, from 1/rho = A + B cos th + C sin th.
 3. Offset correction: in-plane radial gravity is mu*rho/r^3, so mu = mu_fit * <(r/rho)^3>.
 4. Oblateness correction (Bennu C20 from best-fit ellipsoid, spin axis ~ ecliptic pole):
    mu = mu / (1 + J2 * g), g = 1.5 <(R/r)^2 (1 - 3 z^2/r^2)>.
"""
import numpy as np, json
GM_REF, GM_SIG = 4.892, 0.006
J2, RREF = 0.0334, 244.0
d=np.loadtxt('orex_bennu_vectors.csv',delimiter=',',comments='#',skiprows=2)
t=(d[:,0]-d[0,0])*86400; r=d[:,1:4]*1e3; v=d[:,4:7]*1e3
sun=np.array([-6.518e7,-1.168e8,-1.210e7]); s_hat=sun/np.linalg.norm(sun)
# contiguous segments (2 h cadence)
seg_id=np.r_[0,np.cumsum(np.diff(t)>7300)]

def plane(idx):
    R=r[idx]; V=v[idx]; c=R.mean(0); U,S,Wt=np.linalg.svd(R-c); n=Wt[2]
    if np.dot(np.cross(R,V).mean(0),n)<0: n=-n
    return n, np.dot(R,n).mean()
NG,DG=plane(np.arange(len(t)))          # plane and offset from the whole arc
def estimate(idx, corr_offset=True, corr_j2=True, global_plane=True):
    R=r[idx]; V=v[idx]
    n,doff=(NG,DG) if global_plane else plane(idx)                       # out-of-plane offset of the orbit
    Rin=R-np.outer(np.dot(R,n),n); Vin=V-np.outer(np.dot(V,n),n)
    rho=np.linalg.norm(Rin,axis=1); rr=np.linalg.norm(R,axis=1)
    i1=Rin[0]/rho[0]; j1=np.cross(n,i1); th=np.arctan2(Rin@j1,Rin@i1)
    A_=np.vstack([np.ones_like(th),np.cos(th),np.sin(th)]).T
    cf,*_=np.linalg.lstsq(A_,1/rho,rcond=None)
    h=np.linalg.norm(np.cross(Rin,Vin),axis=1)
    mu=np.mean(h)**2*cf[0]
    out={'mu_raw':mu}
    if corr_offset: mu=mu*np.mean((rr/rho)**3)
    out['mu_off']=mu
    g=1.5*np.mean((RREF/rr)**2*(1-3*(R[:,2]/rr)**2))
    if corr_j2: mu=mu/(1+J2*g)
    out.update(mu=mu,doff=doff,g=g,n=n,ecc=np.hypot(cf[1],cf[2])/cf[0],
               sun_angle=np.degrees(np.arccos(abs(np.dot(n,s_hat)))),hours=(t[idx][-1]-t[idx][0])/3600,
               a_srp_normal=mu*abs(doff)/np.mean(rr)**3, sun_dot=float(np.dot(n,s_hat)))
    return out

rows=[]
segs=[np.where(seg_id==k)[0] for k in np.unique(seg_id)]
print('segments (h):',[round((t[s][-1]-t[s][0])/3600,1) for s in segs])
P_h=60.0
for frac,label in [(0.25,'1/4 de órbita'),(0.5,'1/2 órbita'),(1.0,'1 órbita')]:
    L=int(round(frac*P_h/2))+1   # points per window at 2 h cadence
    ests=[]
    for s in segs:
        for k0 in range(0,len(s)-L+1,max(1,L//2)):
            idx=s[k0:k0+L]; ests.append(estimate(idx))
    for key in ['mu_raw','mu_off','mu']:
        a=np.array([e[key] for e in ests])
        rows.append(dict(arco=label,etapa=key,n=len(a),media=float(a.mean()),dp=float(a.std()),
                         erro_pct=float(100*(a.mean()/GM_REF-1)),dp_pct=float(100*a.std()/GM_REF)))
# full contiguous longest segment + all data jointly
full=estimate(np.concatenate(segs[:1]) if len(segs)==1 else segs[0])
allidx=np.arange(len(t)); tot=estimate(allidx)
print(json.dumps(rows,indent=0,ensure_ascii=False))
for name,e in [('segmento contínuo 0',full),('todos os 49 pontos',tot)]:
    print(name,{k:(round(float(x),5) if np.ndim(x)==0 else None) for k,x in e.items() if np.ndim(x)==0})
json.dump(dict(rows=rows,all=dict((k,float(x)) for k,x in tot.items() if np.ndim(x)==0)),open('orex_real.json','w'),ensure_ascii=False,indent=1)
