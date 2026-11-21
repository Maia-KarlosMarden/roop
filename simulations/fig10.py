import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9.5,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#c3c2b7','figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb'})
D=json.load(open('psyche_roop.json')); GM=1.601e9
fig,ax=plt.subplots(1,2,figsize=(10.5,3.9))
L=np.array([(l[0],l[1],l[2],l[3]) for l in D['ROOP_log']]); m=L[:,0]<=4
ax[0].fill_between(L[m,0],L[m,2]/1e9,L[m,3]/1e9,color='#1baf7a',alpha=0.25,lw=0,label='faixa admissível (ajustável)')
ax[0].plot(L[m,0],L[m,1]/1e9,'-o',color='#1baf7a',ms=3.5,lw=2,label='GM estimado (ΔF_geom → 0)')
ax[0].axhline(GM/1e9,color='#5f5e5a',lw=1,ls='--'); ax[0].text(2.2,GM/1e9*0.91,'GM verdadeiro 1,601',fontsize=8,color='#5f5e5a')
ax[0].axhline(1.824,color='#eb6834',lw=1,ls=':'); ax[0].text(2.2,1.824*1.04,'modelo de forma 1,824 (+14 %)',fontsize=8,color='#eb6834')
ax[0].set_yscale('log'); ax[0].set_ylim(0.35,3.4); from matplotlib.ticker import FixedLocator,FixedFormatter,NullFormatter; ax[0].yaxis.set_major_locator(FixedLocator([0.4,0.6,1.0,1.6,2.0,3.0])); ax[0].yaxis.set_major_formatter(FixedFormatter(['0,4','0,6','1,0','1,6','2,0','3,0'])); ax[0].yaxis.set_minor_formatter(NullFormatter()); ax[0].set_xlim(0,4)
ax[0].set_xlabel('dias desde o início da Orbit A'); ax[0].set_ylabel('GM (10⁹ m³/s²)')
ax[0].set_title('Massa sem informação prévia: faixa de densidade 1–8 g/cm³\nse fecha em ±0,2 % no 1º dia',fontsize=9,loc='left')
ax[0].legend(frameon=False,fontsize=8,loc='lower right'); ax[0].grid(color='#e1e0d9',lw=0.6)
cols={'PRIOR_SHAPE':('#eb6834','massa prévia fixa (+14 %)'),'GROUND':('#2a78d6','espera massa do solo (7 dias)'),'ROOP':('#1baf7a','ROOP integrado (sem prévia)')}
for k,(c,lab) in cols.items():
    tr=np.array(D[k+'_truth']); w=int((2*np.pi*np.sqrt(413.2e3**3/GM)/86400)/(30*60/86400))
    ma=np.convolve(tr[:,1],np.ones(w)/w,mode='valid'); tt=tr[w-1:,0]
    ax[1].plot(tt,ma/1e3,color=c,lw=2,label=lab)
ax[1].axhline(413.2,color='#5f5e5a',lw=1,ls='--'); ax[1].text(0.5,418,'alvo Orbit B1 (413 km de raio)',fontsize=8,color='#5f5e5a')
ax[1].set_xlabel('dias desde o início da Orbit A'); ax[1].set_ylabel('raio médio de uma órbita (km)')
ax[1].set_title('Descida A→B1: ROOP chega ~6 dias antes; massa errada\nestaciona 60 km acima do alvo',fontsize=9,loc='left')
ax[1].legend(frameon=False,fontsize=8,loc='upper right'); ax[1].grid(color='#e1e0d9',lw=0.6)
fig.tight_layout(); fig.savefig('fig10_psyche_roop.png',dpi=200)
