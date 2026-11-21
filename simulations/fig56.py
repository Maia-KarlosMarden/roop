import json, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9.5,'axes.spines.top':False,'axes.spines.right':False,
  'axes.edgecolor':'#c3c2b7','figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb'})
E=json.load(open('econ_obj1.json'))
order=['Solo semanal, ingênuo','Solo semanal, inteligente','Solo a cada 3 dias, inteligente','Autônomo convencional','Autônomo ΔF* (relacional)']
LBL=['Solo semanal, ingênuo','Solo semanal, inteligente','Solo a cada 3 dias, inteligente','Autônomo convencional','Autônomo ΔF (relacional)']
col=['#c3c2b7','#c3c2b7','#c3c2b7','#2a78d6','#1baf7a']
fig,ax=plt.subplots(1,3,figsize=(10,3.6),sharey=True)
y=range(len(order))[::-1]
vals=[E[k]['dv_yr'] for k in order]; sd=[E[k]['dv_sd'] for k in order]
ax[0].barh(list(y),vals,xerr=sd,color=col,height=0.6,ecolor='#5f5e5a',capsize=3)
ax[0].axvline(E['Autônomo convencional']['drag_yr'],color='#e34948',lw=1.2,ls='--')
ax[0].text(E['Autônomo convencional']['drag_yr']+0.4,4.42,'mínimo físico (arrasto)',fontsize=8,color='#e34948',va='bottom'); ax[0].set_ylim(-0.5,4.85)
ax[0].set_xlabel('Δv em regime (m/s por ano)'); ax[0].set_yticks(list(y)); ax[0].set_yticklabels(LBL)
for i,v in zip(y,vals): ax[0].text(1,i,f'{v:.1f}',va='center',fontsize=8,color='#1f1f1e')
a=[E[k]['a_rms'] for k in order]
ax[1].barh(list(y),a,color=col,height=0.6); ax[1].set_xlabel('erro de altitude média, RMS (m)')
for i,v in zip(y,a): ax[1].text(v+8,i,f'{v:.0f}',va='center',fontsize=8)
ev=[E[k]['ev_yr'] for k in order]
ax[2].barh(list(y),ev,color=col,height=0.6); ax[2].set_xlabel('ciclos de solo por ano')
for i,v in zip(y,ev): ax[2].text(v+2,i,f'{v:.0f}',va='center',fontsize=8)
for a_ in ax: a_.grid(axis='x',color='#e1e0d9',lw=0.6); a_.set_axisbelow(True)
fig.suptitle('Objetivo 1 — LEO 350 km, 60 dias × 3 sementes: autonomia atinge o mínimo físico com altitude 6–11× mais justa e zero ciclos de solo',fontsize=9.5,x=0.01,ha='left')
fig.tight_layout(); fig.savefig('fig5_market.png',dpi=200)

fig,ax=plt.subplots(1,2,figsize=(10,3.5))
for j,(kind,xs,title) in enumerate([('orbit',[2,4,8,16],'Órbita a 1,75 km (período 58 h)'),('flyby',[6,12,24,48],'Sobrevoo hiperbólico a 7,25 km (tipo OSIRIS-REx)')]):
    A=json.load(open(f'm2_{kind}_1e-07.json')); B=json.load(open(f'm2_{kind}_0.json'))
    w=[str(float(x)) for x in xs]
    ax[j].plot(xs,[A[k]['REL_SRP'][1] for k in w],'-o',color='#1baf7a',lw=2,ms=5,label='Relacional M = ⟨h⟩²/(G·p_geom), SRP compensada')
    ax[j].plot(xs,[A[k]['OD'][1] for k in w],'--s',color='#2a78d6',lw=2,ms=5,label='OD em lote (μ + 6 estados, SRP modelada)')
    if kind=='orbit':
        ax[j].plot(xs,[A[k]['REL'][1] for k in w],':^',color='#eb6834',lw=2,ms=5,label='Relacional sem compensar SRP')
    ax[j].axhline(1,color='#e34948',lw=1,ls='--'); ax[j].text(xs[0],1.08,'1 % (requisito OSIRIS-REx p/ inserção)',fontsize=7.5,color='#e34948')
    ax[j].set_xscale('log',base=2); ax[j].set_yscale('log'); ax[j].set_xticks(xs); ax[j].set_xticklabels([f'{x} h' for x in xs])
    ax[j].set_xlabel('duração do arco de dados'); ax[j].set_title(title,fontsize=9,loc='left')
    ax[j].grid(color='#e1e0d9',lw=0.6)
ax[0].set_ylabel('erro da massa, 1σ (%)')
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,frameon=False,fontsize=8,loc='lower center',ncol=3,bbox_to_anchor=(0.5,0))
fig.suptitle('Objetivo 2 — corpo tipo Bennu (GM 4,89 m³/s²), SRP 1×10⁻⁷ m/s², 40 Monte Carlo por ponto',fontsize=9.5,x=0.01,ha='left')
fig.tight_layout(rect=(0,0.08,1,1)); fig.savefig('fig6_mass.png',dpi=200)
