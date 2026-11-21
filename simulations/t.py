import numpy as np
G=6.6743e-11; M=5.972e24; mu=G*M
# 1) dF with osculating elements from a perturbed state (random state)
rng=np.random.default_rng(0)
for _ in range(3):
    r=np.array([7.0e6,1e5,2e5])*(1+0.01*rng.standard_normal(3)); v=np.array([0,7.6e3,1e3])*(1+0.05*rng.standard_normal(3))
    R=np.linalg.norm(r); h=np.linalg.norm(np.cross(r,v)); w=h/R**2
    E=v@v/2-mu/R; a=-mu/(2*E); evec=np.cross(v,np.cross(r,v))/mu-r/R; e=np.linalg.norm(evec)
    Gp=G*a*(1-e**2); Gpr=w**2*R**3/M*R
    print("dF/Gp =",(Gpr-Gp)/Gp)
# 2) mass estimators
def est(e,mode,N=20000):
    a=1e5; n=np.sqrt(mu/a**3)
    if mode=='theta': th=np.linspace(0,2*np.pi,N,endpoint=False)
    else:
        Mn=np.linspace(0,2*np.pi,N,endpoint=False); Ea=Mn.copy()
        for _ in range(50): Ea=Mn+e*np.sin(Ea)
        th=2*np.arctan2(np.sqrt(1+e)*np.sin(Ea/2),np.sqrt(1-e)*np.cos(Ea/2))
    p=a*(1-e**2); r=p/(1+e*np.cos(th)); h=np.sqrt(mu*p); w=h/r**2
    std=np.mean(w**2*r**3/G)/M
    Ke=(1-e**2)**0.75; enh=(np.mean(w**2)*np.mean(r**2)**1.5/(Ke*G))/M
    return std,enh
for e in [0,0.05,0.1,0.2,0.4]:
    print(e,"theta-unif std/enh:",["%.4f"%x for x in est(e,'theta')]," time-unif:",["%.4f"%x for x in est(e,'time')])
