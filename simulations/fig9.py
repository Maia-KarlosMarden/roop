import json, numpy as np, sys
sys.path.insert(0,'.')
import importlib.util
spec=importlib.util.spec_from_file_location('orx','orex_real.py')
src=open('orex_real.py').read().split("rows=[]")[0]       # reuse definitions only
ns={}; exec(src,ns)
est=ns['estimate']; t=ns['t']; seg_id=ns['seg_id']
segs=[np.where(seg_id==k)[0] for k in np.unique(seg_id)]
pts={}
for frac,lab in [(0.25,'1/4 de órbita\n(16 h)'),(0.5,'1/2 órbita\n(30 h)'),(1.0,'1 órbita\n(60 h)')]:
    L=int(round(frac*60/2))+1; raw=[];fin=[]
    for s in segs:
        for k0 in range(0,len(s)-L+1,max(1,L//2)):
            e=est(s[k0:k0+L]); raw.append(e['mu_raw']); fin.append(e['mu'])
    pts[lab]=(raw,fin)
k3=json.load(open('orex_kepler3.json'))['gm']
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9.5,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#c3c2b7','figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb'})
fig,ax=plt.subplots(figsize=(8.6,4.3))
ax.axhspan(4.892-0.006,4.892+0.006,color='#c3c2b7',alpha=0.45,lw=0)
ax.axhline(4.892,color='#5f5e5a',lw=1)
ax.text(2.55,4.892+0.0068,'referência OSIRIS-REx: 4,892 ± 0,006 m³/s²',fontsize=8,color='#5f5e5a',ha='center',va='bottom')
rng=np.random.default_rng(1)
for i,(lab,(raw,fin)) in enumerate(pts.items()):
    j=rng.uniform(-0.12,0.12,len(raw))
    ax.scatter(i-0.18+j*0.6,raw,s=28,facecolors='none',edgecolors='#2a78d6',lw=1.4,label='relacional bruto' if i==0 else None,zorder=3)
    ax.scatter(i+0.18+j*0.6,fin,s=30,color='#1baf7a',label='relacional + correções fechadas (offset de SRP, J2)' if i==0 else None,zorder=3)
ax.scatter([3],[k3],marker='D',s=40,color='#eb6834',label='3ª lei de Kepler (exige 1 órbita)',zorder=3)
ax.set_xticks([0,1,2,3]); ax.set_xticklabels(list(pts.keys())+['Kepler III\n(61 h)'])
ax.set_ylabel('GM de Bennu (m³/s²)'); ax.set_ylim(4.855,4.935)
sec=ax.secondary_yaxis('right',functions=(lambda y:100*(y/4.892-1),lambda p:4.892*(1+p/100))); sec.set_ylabel('erro (%)')
ax.text(-0.3,4.857,'requisito da missão antes da Orbital A: 2 % (fora da escala, ±0,10 m³/s²)',fontsize=7.5,color='#e34948')
ax.grid(axis='y',color='#e1e0d9',lw=0.6); fig.legend(*ax.get_legend_handles_labels(),frameon=False,fontsize=8,loc='lower center',ncol=3,bbox_to_anchor=(0.5,0))
ax.set_title('Efemérides reais da OSIRIS-REx (JPL Horizons, Orbital A, 3–7 jan 2019): GM de Bennu pelo estimador relacional',fontsize=9.5,loc='left')
fig.tight_layout(rect=(0,0.07,1,1)); fig.savefig('fig9_orex_real.png',dpi=200)
print({k:(np.round(v[0],4).tolist(),np.round(v[1],4).tolist()) for k,v in pts.items()})
