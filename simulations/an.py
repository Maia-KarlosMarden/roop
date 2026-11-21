import json,math,random,sys
sys.path.insert(0,'.')
import sim_leo as S
d=json.load(open('out_NONE.json'))['hist']
# out_every=60 steps=600 s; first 2 orbits
dp=[h[3] for h in d[:20]]
print("p osc range over ~2 orbits (m):",min(dp),max(dp))
# noise on instantaneous p from GNSS
a0=S.RE+S.H0; v0=math.sqrt(S.MU/a0); rnd=random.Random(3); vals=[]
for _ in range(20000):
    r=[a0+rnd.gauss(0,S.SIG_R),rnd.gauss(0,S.SIG_R),rnd.gauss(0,S.SIG_R)]
    v=[rnd.gauss(0,S.SIG_V),v0+rnd.gauss(0,S.SIG_V),rnd.gauss(0,S.SIG_V)]
    h=S.norm(S.cross(r,v)); vals.append(h*h/S.MU-a0)
m=sum(vals)/len(vals); sd=math.sqrt(sum((x-m)**2 for x in vals)/len(vals))
print("GNSS noise on instantaneous p (1sigma, m):",sd)
print("fixed dF_lim 1e-10 -> equiv p (m):",1e-10/S.G)
