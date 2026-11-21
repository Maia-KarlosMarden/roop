import json, math, datetime as dt
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9.5,'axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False,'axes.edgecolor':'#c3c2b7','figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb'})
# Official baseline (Psyche Mission Description and Design Rationale): durations in days
base=[('Orbit A',56),('Transf. B1',17),('Orbit B1',92),('Transf. D',98),('Orbit D',163),('Transf. C',90),('Orbit C',100),('Transf. B2',23),('Orbit B2',100)]
A_MIN=40; D_MIN=100; ORB_PER_DAY_D=684/100
sc={'Conservador':dict(A=0,marg=0.5),'Central':dict(A=8,marg=1.0),'Otimista':dict(A=16,marg=2.0)}
rows=[]
for k,v in sc.items():
    save=v['A']+4*v['marg']; Dn=163+save
    rows.append(dict(cen=k,A_cut=v['A'],marg=4*v['marg'],dias=save,D=Dn,pct_planned=100*save/163,pct_min=100*save/D_MIN,
                     orbitas=save*ORB_PER_DAY_D,precisao=100*(1-math.sqrt(163/Dn))))
for r in rows: print({k:(round(x,1) if isinstance(x,float) else x) for k,x in r.items()})
json.dump(rows,open('psyche_time.json','w'),ensure_ascii=False)
col={'Orbit':'#c3c2b7','Transf':'#e1e0d9','Orbit D':'#1baf7a','ganho':'#0d7a52'}
fig,ax=plt.subplots(figsize=(10,2.9))
lines=[('Plano oficial',0,0)]+[(f"{r['cen']} (+{r['dias']:.0f} d na Orbit D)",sc[r['cen']]['A'],4*sc[r['cen']]['marg']) for r in rows]
for yi,(lab,acut,marg) in enumerate(lines[::-1]):
    y=yi; x=0
    for name,d in base:
        dd=d
        if name=='Orbit A': dd=d-acut
        if name.startswith('Transf') and marg>0 and name!='Transf. C' and name!='Transf. B2' or (name in('Transf. C','Transf. B2') and marg>0):
            dd=d-marg/4
        c=col['Orbit D'] if name=='Orbit D' else (col['Transf'] if name.startswith('Transf') else col['Orbit'])
        if name=='Orbit D':
            extra=acut+marg
            ax.barh(y,163,left=x,color=c,height=0.6,edgecolor='#fcfcfb',lw=1.5)
            if extra>0: ax.barh(y,extra,left=x+163,color=col['ganho'],height=0.6,edgecolor='#fcfcfb',lw=1.5)
            dd=163+extra
        else:
            ax.barh(y,dd,left=x,color=c,height=0.6,edgecolor='#fcfcfb',lw=1.5)
        if yi==len(lines)-1 and name in('Orbit A','Orbit B1','Orbit D','Orbit C','Orbit B2'):
            ax.text(x+dd/2,y+0.42,name,ha='center',fontsize=8,color='#1f1f1e')
        x+=dd
ax.set_yticks(range(len(lines))); ax.set_yticklabels([l[0] for l in lines[::-1]])
ax.set_xlabel('dias desde o início da Orbit A (ago/2029)'); ax.set_xlim(0,745)
ax.grid(axis='x',color='#e1e0d9',lw=0.6); ax.set_axisbelow(True); ax.tick_params(axis='y',length=0)
ax.set_title('Psyche: dias de Orbit A e de margens de transferência realocados para a Orbit D (mapeamento elementar, 75 km)',fontsize=9.5,loc='left',pad=16)
fig.tight_layout(); fig.savefig('fig8_psyche_time.png',dpi=200)
