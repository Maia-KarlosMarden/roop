"""
Objective 1 - fair test of the relational trigger dF* against market practice.
LEO 350 km, i=51.6 deg. Forces: two-body + J2 + J3 + drag with stochastic solar activity
(log-normal Ornstein-Uhlenbeck density factor, sigma=0.35, tau=2 d). GNSS noise 5 m / 5 cm/s.
Thruster 0.1 N, 500 kg, transverse. Target: frozen eccentricity vector (J2/J3) and reference mean a.

Controllers
  NONE         : no control
  GROUND_NAIVE : weekly ground cycle (plan with perfect OD + persistence drag forecast, 1-day uplink latency);
                 drag burn at first opportunity; separate e-correction pair when |de|>box
  GROUND_SMART : same cadence; drag burn placed at u* to fold e-correction (best ground practice)
  AUTO_CONV    : onboard, orbit-averaged osculating a and e-vector; same decision logic as ROOP_STAR
  ROOP_STAR    : onboard, relational dF* (Q = G'r/Gp_ref, K = G'/G fit); same decision logic
"""
import math, random, json, sys
G=6.6743e-11; M=5.972e24; MU=G*M; RE=6378137.0; J2=1.08263e-3; J3=-2.53266e-6
CDAm=2.2/500.0; RHO0=6e-12; H0=350e3; HS=50e3
AMAX=0.1/500.0; SIG_R=5.0; SIG_V=0.05; DT=10.0
INC=math.radians(51.6)

def accel(x, th, fden):
    rx,ry,rz,vx,vy,vz=x
    r2=rx*rx+ry*ry+rz*rz; r=math.sqrt(r2); k=-MU/(r2*r)
    ax,ay,az=k*rx,k*ry,k*rz
    z2=rz*rz/r2
    f=1.5*J2*MU*RE*RE/(r2*r2*r)
    ax+=f*rx*(5*z2-1); ay+=f*ry*(5*z2-1); az+=f*rz*(5*z2-3)
    # J3
    g=-2.5*J3*MU*RE**3/(r2**3*r)   # -(5/2) J3 mu R^3 / r^7
    zr=rz/r
    ax+=g*rx*(3*rz-7*rz**3/r2); ay+=g*ry*(3*rz-7*rz**3/r2)
    az+=g*(6*rz*rz-7*rz**4/r2-0.6*r2)
    h=r-RE; rho=RHO0*fden*math.exp(-(h-H0)/HS); v=math.sqrt(vx*vx+vy*vy+vz*vz)
    d=-0.5*rho*CDAm*v
    ax+=d*vx+th[0]; ay+=d*vy+th[1]; az+=d*vz+th[2]
    return (vx,vy,vz,ax,ay,az), abs(d)*v

def rk4(x,th,fden):
    k1,dr=accel(x,th,fden)
    x2=[x[i]+0.5*DT*k1[i] for i in range(6)]; k2,_=accel(x2,th,fden)
    x3=[x[i]+0.5*DT*k2[i] for i in range(6)]; k3,_=accel(x3,th,fden)
    x4=[x[i]+DT*k3[i] for i in range(6)];     k4,_=accel(x4,th,fden)
    return [x[i]+DT/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]) for i in range(6)], dr

def cross(a,b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def dot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def norm(a): return math.sqrt(dot(a,a))

def frame(r,v):
    hv=cross(r,v); h=norm(hv); hh=[c/h for c in hv]
    n=cross((0,0,1),hh); nn=norm(n); n=[c/nn for c in n]; m=cross(hh,n)
    return h,hv,hh,n,m

def osc_ae(r,v):
    R=norm(r); V2=dot(v,v); a=1/(2/R-V2/MU)
    h,hv,hh,n,m=frame(r,v)
    ev=cross(v,hv); ev=[ev[i]/MU-r[i]/R for i in range(3)]
    return a,(dot(ev,n),dot(ev,m)),h,R,math.atan2(dot(r,m),dot(r,n))

def solve3(Mx,b):
    def det(m): return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    D=det(Mx); return [det([[b[i] if c==j else Mx[i][c] for c in range(3)] for i in range(3)])/D for j in range(3)]

def init_state(a,e,w,nu=0.0):
    p=a*(1-e*e); r=p/(1+e*math.cos(nu)); u=w+nu
    rp=(r*math.cos(u), r*math.sin(u)*math.cos(INC), r*math.sin(u)*math.sin(INC))
    vr=math.sqrt(MU/p)*e*math.sin(nu); vt=math.sqrt(MU/p)*(1+e*math.cos(nu))
    rh=[c/r for c in rp]; th=(-math.sin(u), math.cos(u)*math.cos(INC), math.cos(u)*math.sin(INC))
    return [rp[0],rp[1],rp[2], vr*rh[0]+vt*th[0], vr*rh[1]+vt*th[1], vr*rh[2]+vt*th[2]]

def density_series(days,seed,sigma=0.35,tau=2.0):
    rnd=random.Random(1000+seed); n=int(days*24)+2; x=0.0; out=[]
    a=math.exp(-1/(24*tau)); s=sigma*math.sqrt(1-a*a)
    for _ in range(n): out.append(math.exp(x-0.5*sigma*sigma)); x=a*x+s*rnd.gauss(0,1)
    return out

E_FROZEN=-J3*RE*math.sin(INC)/(2*J2*(RE+H0))   # frozen e (w=90 deg), node frame (0, e_f)

def run(ctrl, days=60, seed=1, de0=3e-4, dv_min=0.03, e_box=1e-4, plan_days=7.0, latency_days=1.0):
    rnd=random.Random(seed); den=density_series(days+1,seed)
    a0=RE+H0
    # start near frozen orbit with an acquisition error de0 along e_x
    ex0,ey0=de0,E_FROZEN; e0=math.hypot(ex0,ey0); w0=math.atan2(ey0,ex0)
    x=init_state(a0,e0,w0)
    P=2*math.pi*math.sqrt(a0**3/MU); W=int(P/DT); N=int(days*86400/DT)
    # reference mean a from one uncontrolled orbit (truth)
    xs=list(x); sa=0
    for _ in range(W):
        aa,_,_,_,_=osc_ae(xs[:3],xs[3:]); sa+=aa; xs,_=rk4(xs,(0,0,0),1.0)
    a_ref=sa/W; p_ref=a_ref*(1-E_FROZEN**2); vbar=math.sqrt(MU/p_ref)
    e_ref=(0.0,E_FROZEN)
    keys=['c','s','cc','ss','cs','K','Kc','Ks','Q','a','ex','ey']
    buf={k:[0.0]*W for k in keys}; S={k:0.0 for k in keys}; idx=0; filled=0
    # burn queue: list of (u_target or None, dv_signed, earliest_step)
    queue=[]; firing=0; pulse=0; lock=0; dv=0.0; dv_drag=0.0; nburn=0; ground_events=0
    # ground bookkeeping
    next_plan=int(plan_days*86400/DT); decay_hist=[]; last_plan_a=None; last_plan_k=0
    truth=[]; t_mean_a=[]; sumA=0; cntA=0
    for k in range(N):
        fden=den[int(k*DT/3600)]
        r=[x[i]+rnd.gauss(0,SIG_R) for i in range(3)]
        v=[x[3+i]+rnd.gauss(0,SIG_V) for i in range(3)]
        a,(exo,eyo),h,R,u=osc_ae(r,v)
        if ctrl in('AUTO_CONV','ROOP_STAR'):
            w_=h/R**2; rho=3*M/(4*math.pi*R**3); Gp=3*w_*w_/(4*math.pi*rho)
            c,s=math.cos(u),math.sin(u)
            K=Gp/G; Q=Gp*R/(G*p_ref)
            vals=dict(c=c,s=s,cc=c*c,ss=s*s,cs=c*s,K=K,Kc=K*c,Ks=K*s,Q=Q,a=a,ex=exo,ey=eyo)
            for kk in keys: S[kk]+=vals[kk]-buf[kk][idx]; buf[kk][idx]=vals[kk]
            idx=(idx+1)%W; filled=min(filled+1,W)
        if lock>0: lock-=1
        # ---------- decision ----------
        if not firing and not queue and lock==0:
            if ctrl in('AUTO_CONV','ROOP_STAR') and filled==W:
                if ctrl=='ROOP_STAR':
                    dvreq=-(vbar/2)*(S['Q']/W-1)
                    A,B,C=solve3([[W,S['c'],S['s']],[S['c'],S['cc'],S['cs']],[S['s'],S['cs'],S['ss']]],[S['K'],S['Kc'],S['Ks']])
                    ex,ey=B/A,C/A
                else:
                    am=S['a']/W; dvreq=(vbar/2)*(a_ref-am)/a_ref
                    ex,ey=S['ex']/W,S['ey']/W
                dex,dey=e_ref[0]-ex,e_ref[1]-ey; de=math.hypot(dex,dey); us=math.atan2(dey,dex)
                if dvreq>dv_min:
                    add_drag(queue, us if de>0.3*2*dvreq/vbar else None, dvreq, k, W)
                    res=de-2*dvreq/vbar
                    if res>e_box: add_pair(queue,us,vbar*res/4,k,W)
                elif de>e_box:
                    add_pair(queue,us,vbar*de/4,k,W)
            elif ctrl.startswith('GROUND') and k>=next_plan:
                # ground: perfect OD of last orbit (truth mean a, e) at planning time
                ground_events+=1
                seg=truth[-int(W/6):]
                am=sum(z[1] for z in seg)/len(seg)+a_ref
                ex,ey=fit_e(seg)
                # persistence forecast of decay rate from last cycle (truth mean a history)
                if last_plan_a is not None:
                    rate=(am-last_plan_a)/((k-last_plan_k)*DT)
                else:
                    rate=-2.4e-3   # first cycle: nominal model value [m/s]
                last_plan_a=am; last_plan_k=k
                lat=int(latency_days*86400/DT); cyc=plan_days*86400
                a_exec=am+rate*latency_days*86400
                dA=(a_ref-a_exec)+(-rate)*cyc/2
                dvreq=max(0.0,(vbar/2)*dA/a_ref)
                dex,dey=e_ref[0]-ex,e_ref[1]-ey; de=math.hypot(dex,dey); us=math.atan2(dey,dex)
                if ctrl=='GROUND_SMART':
                    if dvreq>0: add_drag(queue,us,dvreq,k+lat,W)
                    res=de-2*dvreq/vbar
                else:
                    if dvreq>0: add_drag(queue,None,dvreq,k+lat,W)
                    res=de
                if res>e_box: add_pair(queue,us,vbar*res/4,k+lat,W)
                next_plan+=int(plan_days*86400/DT)
        # ---------- execution ----------
        if not firing and queue and k>=queue[0][2]:
            ut,dvs,_,deadline=queue[0]
            np_=max(1,round(abs(dvs)/AMAX/DT)); arc=np_*DT*2*math.pi/P
            ok= ut is None or abs((u-(ut-arc/2)+math.pi)%(2*math.pi)-math.pi)<math.radians(3) or k>deadline
            if ok:
                firing=1 if dvs>0 else -1; pulse=np_; nburn+=1; queue.pop(0)
        if firing:
            pulse-=1
            hv=cross(x[:3],x[3:]); hn=norm(hv); rn=norm(x[:3])
            d=cross([c/hn for c in hv],[c/rn for c in x[:3]])
            th=tuple(firing*AMAX*d[i] for i in range(3)); dv+=AMAX*DT
            if pulse<=0:
                firing=0
                if not queue: lock=W
        else: th=(0,0,0)
        x,dr=rk4(x,th,fden); dv_drag+=dr*DT
        if k%6==0:
            at,_,_,Rt,ut=osc_ae(x[:3],x[3:])
            truth.append((k*DT/86400, at-a_ref, Rt, ut, dv, dv_drag))
    # per-orbit metrics
    n=int(W/6); orb=[]
    for j in range(0,len(truth)-n,n):
        seg=truth[j:j+n]; ex,ey=fit_e(seg)
        orb.append((seg[0][0], sum(z[1] for z in seg)/n, ex, ey, math.hypot(ex-e_ref[0],ey-e_ref[1])))
    ma=[o[1] for o in orb]; de_=[o[4] for o in orb]
    burn_in=len(orb)//6   # ignore first 10 days for steady-state e metric
    res=dict(ctrl=ctrl,seed=seed,days=days,dv=round(dv,3),dv_drag=round(dv_drag,3),
             de_initial=round(de_[0],6),dv_e_lower=round(vbar*de_[0]/2,3),
             nburn=nburn,ground_events=ground_events,
             a_rms=round(math.sqrt(sum(m*m for m in ma)/len(ma)),1),a_maxabs=round(max(abs(m) for m in ma),1),
             de_rms_steady=round(math.sqrt(sum(d*d for d in de_[burn_in:])/len(de_[burn_in:])),7),
             de_final=round(de_[-1],7))
    def at(day):
        for z in truth:
            if z[0]>=day: return z
        return truth[-1]
    z1,z2=at(20.0),truth[-1]; span=z2[0]-z1[0]
    res['dv_steady_per_yr']=round((z2[4]-z1[4])/span*365.25,2)
    res['drag_steady_per_yr']=round((z2[5]-z1[5])/span*365.25,2)
    res['dv_first20']=round(z1[4],3)
    return res,orb

DV_CAP=0.1   # m/s per pulse (~500 s arc at 2e-4 m/s^2, ~33 deg of arc)
def add_drag(queue,ut,dv,start,W):
    m=max(1,math.ceil(dv/DV_CAP))
    for j in range(m): queue.append([ut,dv/m,start,start+(j+1)*W])
def add_pair(queue,ut,d1,start,W):
    m=max(1,math.ceil(d1/DV_CAP))
    for j in range(m):
        queue.append([ut,d1/m,start,start+(2*j+1)*W]); queue.append([ut+math.pi,-d1/m,start,start+(2*j+2)*W])

def fit_e(seg):
    Sx=[[0]*3 for _ in range(3)]; Sb=[0]*3
    for z in seg:
        Rt,ut=z[2],z[3]; f=(1,math.cos(ut),math.sin(ut))
        for i in range(3):
            Sb[i]+=f[i]/Rt
            for l in range(3): Sx[i][l]+=f[i]*f[l]
    A,B,C=solve3(Sx,Sb); return B/A,C/A

if __name__=='__main__':
    ctrl=sys.argv[1]; seed=int(sys.argv[2]) if len(sys.argv)>2 else 1
    days=float(sys.argv[3]) if len(sys.argv)>3 else 60
    pdays=float(sys.argv[4]) if len(sys.argv)>4 else 7.0
    ebox=float(sys.argv[5]) if len(sys.argv)>5 else 1e-4
    res,orb=run(ctrl,days=days,seed=seed,plan_days=pdays,e_box=ebox)
    res['plan_days']=pdays; res['e_box']=ebox
    print(json.dumps(res),flush=True)
    tag='' if (pdays==7.0 and ebox==1e-4) else f'_p{pdays:g}_b{ebox:g}'
    json.dump(dict(res=res,orb=orb),open(f'mk_{ctrl}_{seed}{tag}.json','w'))
