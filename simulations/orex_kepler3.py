import numpy as np, json
d=np.loadtxt('orex_bennu_vectors.csv',delimiter=',',comments='#',skiprows=2)
t=(d[:,0]-d[0,0])*86400; r=d[:,1:4]*1e3; v=d[:,4:7]*1e3
seg=np.where(t<=72*3600+1)[0]          # 72 h contiguous segment
R=r[seg]; V=v[seg]; T=t[seg]
n=np.cross(R,V).mean(0); n/=np.linalg.norm(n); i1=R[0]-n*(R[0]@n); i1/=np.linalg.norm(i1); j1=np.cross(n,i1)
th=np.unwrap(np.arctan2(R@j1,R@i1))
# period: time for true-anomaly-like angle to advance 2*pi (interpolated)
Tp=np.interp(th[0]+2*np.pi,th,T)-T[0]
rho=np.linalg.norm(R-np.outer(R@n,n),axis=1)
A=np.vstack([np.ones_like(th),np.cos(th),np.sin(th)]).T; c,*_=np.linalg.lstsq(A,1/rho,rcond=None)
p=1/c[0]; e=np.hypot(c[1],c[2])/c[0]; a=p/(1-e*e)
mu=4*np.pi**2*a**3/Tp**2
print(f'Kepler III: period {Tp/3600:.2f} h, a {a:.1f} m, e {e:.3f}, GM {mu:.4f} ({100*(mu/4.892-1):+.2f}%), needs >= {Tp/3600:.0f} h of data')
json.dump(dict(period_h=Tp/3600,a=a,e=e,gm=mu),open('orex_kepler3.json','w'))
