import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9.5,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#c3c2b7','figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb'})
GM=1.601e9
cols={'GROUND':('#2a78d6','wait for ground gravity (7 days)'),'ROOP_FIX':('#eda100','ROOP, J2 fixed from shape model'),'ROOP_J2':('#1baf7a','integrated ROOP (mass + J2 in flight)')}
fig,ax=plt.subplots(1,2,figsize=(10.5,3.8))
for a,(leg,r1,title) in zip(ax,[('B1_D',188.2e3,'B1 → D (413 → 188 km), retrograde descent'),('D_C',303.2e3,'D → C (188 → 303 km)')]):
    D=json.load(open(f'psyche_roop2_{leg}.json'))
    for k,(c,lab) in cols.items():
        tr=np.array(D[k+'_truth']); dt=tr[1,0]-tr[0,0]
        w=max(1,int(round((2*np.pi*np.sqrt(r1**3/GM)/86400)/dt)))
        ma=np.convolve(tr[:,1],np.ones(w)/w,mode='valid'); a.plot(tr[w-1:,0],ma/1e3,color=c,lw=2,label=lab,ls='--' if k=='ROOP_FIX' else '-')
    a.axhline(r1/1e3,color='#5f5e5a',lw=1,ls=':'); a.set_xlabel('days since start of transfer'); a.set_ylabel('one-orbit mean radius (km)')
    a.set_title(title,fontsize=9.5,loc='left'); a.grid(color='#e1e0d9',lw=0.6)
ax[0].text(20,196,'D target',fontsize=8,color='#5f5e5a'); ax[1].text(14,295,'C target',fontsize=8,color='#5f5e5a')
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,frameon=False,fontsize=8,loc='lower center',ncol=3,bbox_to_anchor=(0.5,0))
fig.tight_layout(rect=(0,0.07,1,1)); fig.savefig('fig11_psyche_legs_en.png',dpi=200)
