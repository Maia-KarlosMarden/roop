import json, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
  'axes.edgecolor':'#c3c2b7','axes.grid':True,'grid.color':'#e1e0d9','grid.linewidth':0.6,
  'figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb'})
C=['#eda100','#2a78d6','#1baf7a']
S=[('BASE_a','Convencional (a médio, pulso fixo)',C[0],'-'),
   ('ROOP_star_noloc','ΔF* só com dimensionamento',C[1],'-'),
   ('ROOP_star','ΔF* com dimensionamento + local u*',C[2],'-')]
fig,(a1,a2)=plt.subplots(2,1,figsize=(7.6,5.9),sharex=True)
for k,lab,c,ls in S:
    o=json.load(open(f'star_{k}_0.001_1.json'))['orb']
    t=[x[0] for x in o]
    a1.plot(t,[x[1] for x in o],color=c,lw=2,label=lab)
    a2.plot(t,[x[2]*1e4 for x in o],color=c,lw=2,label=lab)
a1.axhline(0,color='#c3c2b7',lw=1); a1.set_ylabel('erro de a médio (m)')
a1.set_title('Altitude: ΔF* dispara 2,4× menos e mantém a órbita 3× mais justa',fontsize=9.5,loc='left')
a2.set_ylabel('excentricidade média (×10⁻⁴)'); a2.set_xlabel('tempo (dias)')
a2.set_title('Excentricidade: escolher onde disparar circulariza a órbita sem Δv extra',fontsize=9.5,loc='left')
fig.legend(*a1.get_legend_handles_labels(),frameon=False,fontsize=8.5,loc='lower center',ncol=3,bbox_to_anchor=(0.5,0))
fig.tight_layout(rect=(0,0.05,1,1)); fig.savefig('fig4_star.png',dpi=200)
