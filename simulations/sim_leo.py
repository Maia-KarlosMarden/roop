"""
ROOP test bench - LEO station-keeping (altitude maintenance) under J2 + drag + GNSS noise.
Controllers:
  ROOP_osc   : literal reading. dF = G'r - G p_osc, p_osc from osculating elements of same state.
  ROOP_ref   : dF = G'r - G p_ref (p_ref = target), instantaneous, fixed threshold in G'r units.
  ROOP_ref_W : same, wide threshold (100 m equivalent in p).
  ROOP_avg   : orbit-averaged dF (one-period running mean), deadband 100 m equivalent in p.
  BASE_a     : conventional: orbit-averaged semi-major-axis error, tangential thrust, same deadband.
Thrust: on/off micro-thruster, 0.1 N, 500 kg sat. ROOP fires transverse (theta-hat, d|h|/dt = r*a_theta);
BASE fires along velocity.
"""
import math, random, json, sys
G=6.6743e-11; M=5.972e24; MU=G*M; RE=6378137.0; J2=1.08263e-3
CDAm=2.2*1.0/500.0           # Cd*A/m  [m^2/kg]
RHO0=6e-12; H0=350e3; HS=50e3
AMAX=0.1/500.0               # 2e-4 m/s^2
SIG_R=5.0; SIG_V=0.05        # GNSS 1-sigma per axis
DT=10.0

def accel(x, th):
    rx,ry,rz,vx,vy,vz=x
    r2=rx*rx+ry*ry+rz*rz; r=math.sqrt(r2)
    k=-MU/(r2*r)
    ax,ay,az=k*rx,k*ry,k*rz
    # J2
    f=1.5*J2*MU*RE*RE/(r2*r2*r); z2=rz*rz/r2
    ax+=f*rx*(5*z2-1); ay+=f*ry*(5*z2-1); az+=f*rz*(5*z2-3)
    # drag (non-rotating atmosphere)
    h=r-RE; rho=RHO0*math.exp(-(h-H0)/HS); v=math.sqrt(vx*vx+vy*vy+vz*vz)
    d=-0.5*rho*CDAm*v
    ad=(d*vx,d*vy,d*vz)
    ax+=ad[0]+th[0]; ay+=ad[1]+th[1]; az+=ad[2]+th[2]
    return (vx,vy,vz,ax,ay,az), math.sqrt(ad[0]**2+ad[1]**2+ad[2]**2)

def rk4(x,th):
    k1,dr=accel(x,th)
    x2=[x[i]+0.5*DT*k1[i] for i in range(6)]; k2,_=accel(x2,th)
    x3=[x[i]+0.5*DT*k2[i] for i in range(6)]; k3,_=accel(x3,th)
    x4=[x[i]+DT*k3[i] for i in range(6)];     k4,_=accel(x4,th)
    return [x[i]+DT/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]) for i in range(6)], dr

def cross(a,b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def norm(a): return math.sqrt(a[0]**2+a[1]**2+a[2]**2)

def elements(r,v):
    R=norm(r); V2=v[0]**2+v[1]**2+v[2]**2
    hv=cross(r,v); h=norm(hv)
    a=1/(2/R-V2/MU)
    ev=cross(v,hv); ev=[ev[i]/MU-r[i]/R for i in range(3)]; e=norm(ev)
    return a,e,h,hv,R

def run(ctrl, days=15, seed=1, out_every=60):
    rnd=random.Random(seed)
    a0=RE+H0; e0=0.001; inc=math.radians(51.6)
    rp=a0*(1-e0); vp=math.sqrt(MU*(1+e0)/rp)
    x=[rp,0,0, 0,vp*math.cos(inc),vp*math.sin(inc)]
    N=int(days*86400/DT); P=2*math.pi*math.sqrt(a0**3/MU); W=int(P/DT)
    # references = orbit-mean of one uncontrolled orbit (removes J2 short-period offset)
    xs=list(x); sa=0; sp=0
    for _ in range(W):
        aa,_,hh,_,_=elements(xs[:3],xs[3:]); sa+=aa; sp+=hh*hh/MU
        xs,_=rk4(xs,(0,0,0))
    a_ref=sa/W; p_ref=sp/W
    Gp_ref=G*p_ref
    buf=[0.0]*W; s=0.0; idx=0; filled=0
    lock=0; pulse=0
    firing=0; dv=0.0; dv_drag=0.0; nfire=0; th=(0,0,0)
    # thresholds
    db_p=100.0                      # m  (equivalent)
    lim_fixed=1e-10                # fixed dF_lim in G'r units (negative test)
    lim_W=G*db_p                    # dF equivalent of 100 m in p
    hist=[]
    for k in range(N):
        # measurement
        r=[x[i]+rnd.gauss(0,SIG_R) for i in range(3)]
        v=[x[3+i]+rnd.gauss(0,SIG_V) for i in range(3)]
        a,e,h,hv,R=elements(r,v)
        w=h/R**2                                 # instantaneous angular velocity
        rho_orb=3*M/(4*math.pi*R**3)             # 'orbital density'
        Gprime=3*w*w/(4*math.pi*rho_orb)         # = w^2 R^3 / M
        GpR=Gprime*R                             # = h^2/M
        if ctrl=='ROOP_osc':
            sig=GpR-G*a*(1-e*e); lim=lim_fixed
        elif ctrl in('ROOP_ref','ROOP_ref_W'):
            sig=GpR-Gp_ref; lim=lim_fixed if ctrl=='ROOP_ref' else lim_W
        elif ctrl=='ROOP_avg':
            raw=GpR-Gp_ref
            s+=raw-buf[idx]; buf[idx]=raw; idx=(idx+1)%W; filled=min(filled+1,W)
            sig=s/filled; lim=lim_W
        elif ctrl=='BASE_a':
            raw=G*(a-a_ref)                       # scale to same units for uniform logic
            s+=raw-buf[idx]; buf[idx]=raw; idx=(idx+1)%W; filled=min(filled+1,W)
            sig=s/filled; lim=lim_W
        elif ctrl=='NONE':
            sig=0; lim=1
        # control logic. One-sided (drag makeup: thrust only raises the orbit).
        if ctrl in('ROOP_avg','BASE_a'):
            # pulse-and-wait: 60 s pulse, then one-orbit lockout so the filter refreshes
            if lock>0: lock-=1
            if firing!=0:
                pulse-=1
                if pulse<=0: firing=0; lock=W
            elif lock==0 and filled==W and sig<-lim:
                firing=1; pulse=int(60/DT); nfire+=1
        else:
            # literal reading: fire when |dF|>lim, stop when dF -> 0 (sign crossing)
            if firing==0 and sig<-lim:
                firing=1; nfire+=1
            elif firing!=0 and sig>=0:
                firing=0
        # thrust direction from TRUE state (attitude assumed perfect)
        rt=x[:3]; vt=x[3:]
        if firing!=0:
            if ctrl=='BASE_a':
                vn=norm(vt); d=[vt[i]/vn for i in range(3)]
            else:
                hvt=cross(rt,vt); hn=norm(hvt); rn=norm(rt)
                hh=[hvt[i]/hn for i in range(3)]; rr=[rt[i]/rn for i in range(3)]
                d=cross(hh,rr)
            th=tuple(firing*AMAX*d[i] for i in range(3)); dv+=AMAX*DT
        else:
            th=(0,0,0)
        x,dr=rk4(x,th); dv_drag+=dr*DT
        if k%out_every==0:
            at,et,ht,_,Rt=elements(x[:3],x[3:])
            hist.append((k*DT/86400, at-a_ref, et, ht*ht/MU-p_ref, Rt-RE, sig/G, firing))
    # summary on true state, orbit-mean a over last part
    arr=hist
    da=[hh[1] for hh in arr]; dp=[hh[3] for hh in arr]
    # orbit-mean (smooth) of a error to remove J2 short-period
    wnd=max(1,int(P/(DT*out_every)))
    ma=[sum(da[i:i+wnd])/wnd for i in range(0,len(da)-wnd)]
    res=dict(ctrl=ctrl, dv=dv, dv_drag=dv_drag, nfire=nfire,
             mean_a_err_final=ma[-1], mean_a_err_rms=math.sqrt(sum(m*m for m in ma)/len(ma)),
             e_final=arr[-1][2], duty=dv/(AMAX*N*DT))
    return res, hist

if __name__=='__main__':
    ctrl=sys.argv[1]; days=float(sys.argv[2]) if len(sys.argv)>2 else 15
    res,hist=run(ctrl,days)
    print(json.dumps(res))
    json.dump(dict(res=res,hist=hist),open(f'out_{ctrl}.json','w'))
