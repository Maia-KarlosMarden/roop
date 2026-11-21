"""
ROOP* : relational deviation recomputed to be actionable for thrusting.

Per sample, two relational quantities (from G' = 3w^2/4*pi*rho and rho = 3M/4*pi*r^3):
   Q_i = G'_i r_i / (G p_ref)    (= p_i/p_ref : size)
   K_i = G'_i / G                (= p/r = 1 + e cos(u - w) : shape)
Over a sliding window of one orbit:
   dF_p   = <Q> - 1                         (orbit-averaged relational deviation, J2/noise filtered)
   e_vec  = (B/A, C/A) from LSQ  K = A + B cos u + C sin u
   dF*    = Dv_req = -(v/2) dF_p            [m/s]  -> IF and HOW MUCH
   u*     = atan2(de_needed)                 -> WHERE (argument of latitude of the burn)
Trigger: Dv_req > Dv_min ; burn transverse for Dv_req/a_max seconds at u ~ u* ; then 1-orbit lockout.
"""
import math, random, json, sys
sys.path.insert(0,'.')
from sim_leo import G,M,MU,RE,H0,AMAX,SIG_R,SIG_V,DT,rk4,elements,cross,norm

def arg_lat(r,hv):
    hn=norm(hv); hh=[hv[i]/hn for i in range(3)]
    n=cross((0,0,1),hh); nn=norm(n); n=[n[i]/nn for i in range(3)]
    m=cross(hh,n)
    return math.atan2(sum(r[i]*m[i] for i in range(3)), sum(r[i]*n[i] for i in range(3)))

def solve3(Mx,b):
    def det(m): return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    D=det(Mx); out=[]
    for j in range(3):
        mj=[[b[i] if c==j else Mx[i][c] for c in range(3)] for i in range(3)]
        out.append(det(mj)/D)
    return out

def run(ctrl, days=15, seed=1, e0=0.001, dv_min=0.03, e_ref=(0.0,0.0)):
    rnd=random.Random(seed)
    a0=RE+H0; inc=math.radians(51.6)
    rp=a0*(1-e0); vp=math.sqrt(MU*(1+e0)/rp)
    x=[rp,0,0, 0,vp*math.cos(inc),vp*math.sin(inc)]
    N=int(days*86400/DT); P=2*math.pi*math.sqrt(a0**3/MU); W=int(P/DT)
    xs=list(x); sa=0; sp=0
    for _ in range(W):
        aa,_,hh,_,_=elements(xs[:3],xs[3:]); sa+=aa; sp+=hh*hh/MU
        xs,_=rk4(xs,(0,0,0))
    a_ref=sa/W; p_ref=sp/W; vbar=math.sqrt(MU/p_ref)
    # ring buffers
    keys=['c','s','cc','ss','cs','K','Kc','Ks','Q','raw','a']
    buf={k:[0.0]*W for k in keys}; S={k:0.0 for k in keys}; idx=0; filled=0
    lock=0; pulse=0; firing=0; dv=0; dv_drag=0; nfire=0; wait=0; armed=False; dvreq=0; ustar=0
    lim_W=G*100.0; hist=[]; truth=[]; log=[]
    for k in range(N):
        r=[x[i]+rnd.gauss(0,SIG_R) for i in range(3)]
        v=[x[3+i]+rnd.gauss(0,SIG_V) for i in range(3)]
        a,e,h,hv,R=elements(r,v)
        w=h/R**2; rho=3*M/(4*math.pi*R**3); Gp=3*w*w/(4*math.pi*rho)
        u=arg_lat(r,hv); c,s=math.cos(u),math.sin(u)
        K=Gp/G; Q=Gp*R/(G*p_ref)
        vals=dict(c=c,s=s,cc=c*c,ss=s*s,cs=c*s,K=K,Kc=K*c,Ks=K*s,Q=Q,raw=Gp*R-G*p_ref,a=G*(a-a_ref))
        for kk in keys:
            S[kk]+=vals[kk]-buf[kk][idx]; buf[kk][idx]=vals[kk]
        idx=(idx+1)%W; filled=min(filled+1,W)
        if lock>0: lock-=1
        if firing:
            pulse-=1
            if pulse<=0: firing=0; lock=W
        elif lock==0 and filled==W:
            if ctrl in('ROOP_star','ROOP_star_noloc'):
                if not armed:
                    dFp=S['Q']/W-1
                    dvreq=-(vbar/2)*dFp
                    if dvreq>dv_min:
                        A,B,C=solve3([[W,S['c'],S['s']],[S['c'],S['cc'],S['cs']],[S['s'],S['cs'],S['ss']]],[S['K'],S['Kc'],S['Ks']])
                        ex,ey=B/A,C/A; dex,dey=e_ref[0]-ex,e_ref[1]-ey
                        de=math.hypot(dex,dey); dep=2*dvreq/vbar
                        armed=True; wait=0
                        ustar=math.atan2(dey,dex) if (de>0.3*dep and ctrl=='ROOP_star') else None
                        log.append((k*DT/86400,dvreq,ex,ey,ustar))
                if armed:
                    wait+=1
                    ok = ustar is None or abs((u-ustar+math.pi)%(2*math.pi)-math.pi)<math.radians(10) or wait>W
                    if ok:
                        firing=1; pulse=max(1,round(dvreq/AMAX/DT)); nfire+=1; armed=False
            elif ctrl in('ROOP_avg','BASE_a'):
                sig=S['raw' if ctrl=='ROOP_avg' else 'a']/W
                if sig<-lim_W: firing=1; pulse=int(60/DT); nfire+=1
        rt=x[:3]; vt=x[3:]
        if firing:
            if ctrl=='BASE_a':
                vn=norm(vt); d=[vt[i]/vn for i in range(3)]
            else:
                hvt=cross(rt,vt); hn=norm(hvt); rn=norm(rt)
                d=cross([hvt[i]/hn for i in range(3)],[rt[i]/rn for i in range(3)])
            th=tuple(AMAX*d[i] for i in range(3)); dv+=AMAX*DT
        else: th=(0,0,0)
        x,dr=rk4(x,th); dv_drag+=dr*DT
        if k%6==0:
            at,_,ht,hvt,Rt=elements(x[:3],x[3:])
            truth.append((k*DT/86400, at-a_ref, Rt, arg_lat(x[:3],hvt)))
    # per-orbit geometric fit of true trajectory: mean a error and e
    n=int(W/6); orb=[]
    for j in range(0,len(truth)-n,n):
        seg=truth[j:j+n]
        Sx=[[0]*3 for _ in range(3)]; Sb=[0]*3
        for t_,da_,Rt,ut in seg:
            f=(1,math.cos(ut),math.sin(ut))
            for i in range(3):
                Sb[i]+=f[i]/Rt
                for l in range(3): Sx[i][l]+=f[i]*f[l]
        A,B,C=solve3(Sx,Sb)
        orb.append((seg[0][0], sum(z[1] for z in seg)/n, math.hypot(B,C)/A))
    ma=[o[1] for o in orb]; ee=[o[2] for o in orb]
    res=dict(ctrl=ctrl,dv=round(dv,3),dv_drag=round(dv_drag,3),nfire=nfire,
             a_err_final=round(ma[-1],1),a_err_rms=round(math.sqrt(sum(m*m for m in ma)/len(ma)),1),
             e_initial=round(ee[0],6),e_final=round(ee[-1],6),e_mean_last5d=round(sum(ee[-len(ee)//3:])/(len(ee)//3),6))
    return res, orb, log

if __name__=='__main__':
    ctrl=sys.argv[1]; e0=float(sys.argv[2]) if len(sys.argv)>2 else 0.001
    seed=int(sys.argv[3]) if len(sys.argv)>3 else 1
    res,orb,log=run(ctrl,e0=e0,seed=seed); res['seed']=seed
    res['e0']=e0
    print(json.dumps(res))
    json.dump(dict(res=res,orb=orb,log=log),open(f'star_{ctrl}_{e0}_{seed}.json','w'))
