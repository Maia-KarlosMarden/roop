import json, math, statistics as st, glob
g0=9.80665
# ---------- Objective 1 ----------
def load(pattern):
    R=[json.load(open(f))['res'] for f in sorted(glob.glob(pattern))]
    return R
S={
 'Solo semanal, ingênuo':load('mk_GROUND_NAIVE_[123].json'),
 'Solo semanal, inteligente':load('mk_GROUND_SMART_[123].json'),
 'Solo a cada 3 dias, inteligente':load('mk_GROUND_SMART_[123]_p3_b0.0001.json'),
 'Autônomo convencional':load('mk_AUTO_CONV_[123].json'),
 'Autônomo ΔF* (relacional)':load('mk_ROOP_STAR_[123].json'),
}
out={}
for k,R in S.items():
    dv=[r['dv_steady_per_yr'] for r in R]; dr=[r['drag_steady_per_yr'] for r in R]
    eff=[d/x for d,x in zip(dr,dv)]
    out[k]=dict(n=len(R),dv_yr=st.mean(dv),dv_sd=st.pstdev(dv),drag_yr=st.mean(dr),waste=100*(st.mean(dv)/st.mean(dr)-1),
                a_rms=st.mean(r['a_rms'] for r in R),de=st.mean(r['de_rms_steady'] for r in R),
                ev_yr=st.mean(r['ground_events'] for r in R)/60*365.25, burns_yr=st.mean(r['nburn'] for r in R)/60*365.25)
auto=out['Autônomo ΔF* (relacional)']['dv_yr']
for k,v in out.items():
    v['saving_vs_this_%']=100*(1-auto/v['dv_yr'])
    for isp,lab in((220,'hidrazina'),(1500,'elétrica')):
        v[f'prop_kg_yr_{lab}']=500*(1-math.exp(-v['dv_yr']/(isp*g0)))
    print(k,{kk:round(vv,3) if isinstance(vv,float) else vv for kk,vv in v.items()})
json.dump(out,open('econ_obj1.json','w'),ensure_ascii=False,indent=1)
# ---------- Objective 2 ----------
water_opt=2.3e-4*150*86400*0.05; pt_opt=0.35*150*86400*1e-5
print('Hein optimistic: water kg/day',round(water_opt,1),' Pt kg/day',round(pt_opt,1))
prices={'água (US$ 20 mil/kg, cislunar)':(water_opt,20000),'platina (US$ 43,7 mil/kg, spot 19/11/2026)':(pt_opt,43666)}
rows=[]
for name,(kgd,p) in prices.items():
    for thr,lab in ((1.0,'otimista'),(0.1,'conservador (10%)')):
        vday=kgd*thr*p
        for D in (3,9,15):
            rows.append((name,lab,round(kgd*thr,1),round(vday/1e6,2),D,round(D*vday/1e6,1),round(D,1)))
for r in rows: print(r)
json.dump(rows,open('econ_obj2.json','w'),ensure_ascii=False)
