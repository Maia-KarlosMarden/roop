import sys,math,time,random; sys.path.insert(0,'.')
import sim_market as S
rnd=random.Random(0); a0=S.RE+S.H0; x=S.init_state(a0,1e-3,1.0)
samples=[]
for _ in range(20000):
    r=[x[i]+rnd.gauss(0,5) for i in range(3)]; v=[x[3+i]+rnd.gauss(0,.05) for i in range(3)]; samples.append((r,v))
G,M=S.G,S.M; pref=a0
def conv(r,v):           # conventional onboard sensing: osculating a, e-vector in node frame
    a,(ex,ey),h,R,u=S.osc_ae(r,v); return a,ex,ey
def rel(r,v):            # relational sensing: G', Q, K and u
    hv=S.cross(r,v); h=S.norm(hv); R=S.norm(r)
    hh=[c/h for c in hv]; n=S.cross((0,0,1),hh); nn=S.norm(n); n=[c/nn for c in n]; m=S.cross(hh,n)
    u=math.atan2(S.dot(r,m),S.dot(r,n))
    w=h/R**2; rho=3*M/(4*math.pi*R**3); Gp=3*w*w/(4*math.pi*rho)
    return Gp*R/(G*pref), Gp/G, math.cos(u), math.sin(u)
for f in (conv,rel,conv,rel):
    t=time.perf_counter()
    for _ in range(5):
        for r,v in samples: f(r,v)
    print(f.__name__, round((time.perf_counter()-t)/(5*len(samples))*1e6,3),'us/sample')
