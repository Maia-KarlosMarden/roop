import sys,time,numpy as np; sys.path.insert(0,'.')
import sim_mass2 as S
from scipy.optimize import least_squares
rng=np.random.default_rng(3); GM=S.GM; f=np.array([1e-7,3e-8,0])
rp=7250.0; vinf=0.1; vp=np.sqrt(vinf**2+2*GM/rp); half=12*3600
x0=S.prop(np.array([rp,0,0,0,vp,0]),np.array([0,-half]),GM,np.zeros(3))[-1]; t=np.linspace(0,2*half,60)
tr=S.prop(x0,t,GM,f)
tR=[];tO=[];fails={0.2:0,0.5:0,1.0:0,2.0:0}; n=15
for k in range(n):
    rm=tr[:,:3]+rng.normal(0,S.SR,(60,3)); vm=tr[:,3:]+rng.normal(0,S.SV,(60,3)); dt=(t-t.mean())[:,None]
    a=time.perf_counter(); m=S.rel(rm-0.5*f*dt**2,vm-f*dt); tR.append(time.perf_counter()-a)
    for g in fails:
        def res(x): return ((S.prop(x[:6],t,x[6],f)[:,:3]-rm)/S.SR).ravel()
        a=time.perf_counter()
        sol=least_squares(res,np.r_[rm[0],vm[0],GM*(1+g)],x_scale=np.r_[np.full(3,10.),np.full(3,1e-3),0.5],max_nfev=60)
        if g==0.2: tO.append(time.perf_counter()-a)
        if abs(sol.x[6]/GM-1)>0.05: fails[g]+=1
print('REL time per estimate (ms):',round(1e3*np.median(tR),3))
print('OD  time per estimate (s):',round(np.median(tO),2),' ratio ~',int(np.median(tO)/np.median(tR)))
print('OD failures (>5% error) by a-priori mass error:',{f'{int(g*100)}%':f'{v}/{n}' for g,v in fails.items()})
