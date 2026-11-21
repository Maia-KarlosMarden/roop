import sys,math,random,statistics as st; sys.path.insert(0,'.')
import sim_market as S
rnd=random.Random(5); a0=S.RE+S.H0; x=S.init_state(a0,1e-3,1.0)
A=[];P=[];EX=[];KX=[]
a_t,(ex_t,ey_t),h_t,R_t,u_t=S.osc_ae(x[:3],x[3:])
for _ in range(20000):
    r=[x[i]+rnd.gauss(0,5) for i in range(3)]; v=[x[3+i]+rnd.gauss(0,.05) for i in range(3)]
    a,(ex,ey),h,R,u=S.osc_ae(r,v); A.append(a); P.append(h*h/S.MU)
print('sigma a (vis-viva) m:',round(st.pstdev(A),1),' sigma p=h^2/mu (relational Q*p_ref) m:',round(st.pstdev(P),1))
