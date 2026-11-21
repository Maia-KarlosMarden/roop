import json, math, sys
sys.path.insert(0,'.')
import sim_leo as S
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
  'axes.edgecolor':'#c3c2b7','axes.grid':True,'grid.color':'#e1e0d9','grid.linewidth':0.6,
  'figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb','text.color':'#1f1f1e','axes.labelcolor':'#1f1f1e'})
C=['#2a78d6','#eb6834','#1baf7a','#eda100']
def ma(xs,w):
    return [sum(xs[i:i+w])/w for i in range(len(xs)-w)]
W=9  # ~1 orbit at 600 s sampling
fig,ax=plt.subplots(figsize=(7,3.8))
series=[('NONE','Sem controle = ROOP literal (p osculador)',C[0],'-'),
        ('ROOP_ref','ROOP instantâneo (p de referência)',C[1],'-'),
        ('ROOP_avg','ROOP com média orbital',C[2],'-'),
        ('BASE_a','Convencional (a médio)',C[3],'--')]
for k,lab,c,ls in series:
    h=json.load(open(f'out_{k}.json'))['hist']
    t=[x[0] for x in h]; da=[x[1]/1000 for x in h]
    ax.plot(t[:len(t)-W],ma(da,W),color=c,lw=2,ls=ls,label=lab)
ax.axhline(0,color='#c3c2b7',lw=1)
ax.set_xlabel('tempo (dias)'); ax.set_ylabel('erro do semieixo maior médio (km)')
ax.set_title('Manutenção de altitude em LEO 350 km (J2 + arrasto + ruído GNSS)',fontsize=10,loc='left')
ax.legend(frameon=False,fontsize=8.5,loc='center left')
ax.annotate('verde e amarelo coincidem',(13,-0.1),xytext=(10.5,1.0),fontsize=8,color='#5f5e5a',arrowprops=dict(arrowstyle='-',color='#5f5e5a',lw=0.8))
fig.tight_layout(); fig.savefig('fig1_leo.png',dpi=200)

# fig2: instantaneous p over 2 orbits, uncontrolled
res,h=S.run('NONE',days=0.13,out_every=1)
t=[x[0]*24 for x in h]; dp=[x[3]/1000 for x in h]
fig,ax=plt.subplots(figsize=(7,3.6))
ax.plot(t,dp,color=C[0],lw=2,label='p instantâneo − p de referência (J2)')
ax.axhspan(-0.088,0.088,color=C[1],alpha=0.35,lw=0,label='ruído GNSS 1σ em p (±88 m)')
ax.axhline(0,color='#1f1f1e',lw=0.8)
ax.set_xlabel('tempo (horas)',labelpad=2); ax.set_ylabel('desvio de p (km)')
ax.set_title('Sinal ΔF/G instantâneo sem controle: o limiar do pedido (±1,5 m) é invisível nesta escala',fontsize=9.5,loc='left')
ax.legend(frameon=False,fontsize=8.5,loc='upper center',bbox_to_anchor=(0.5,-0.22),ncol=2)
fig.tight_layout(); fig.savefig('fig2_signal.png',dpi=200)

# fig3: mass estimator bias, e=0.1, one orbit (from MC table) 
R=json.load(open('mass_results.json'))
fig,ax=plt.subplots(figsize=(7,3.2))
keys=['STD','ENH','REL','OD']; labs=['Standard\n(pedido)','Enhanced Kₑ\n(pedido)','Relacional\nh²/(G·p geom.)','OD clássica\n(μ livre)']
for j,(arc,c) in enumerate([('0.25',C[0]),('1.0',C[1])]):
    b=[R[f'0.1_{arc}'][k][0] for k in keys]; s=[R[f'0.1_{arc}'][k][1] for k in keys]
    xs=[i+(j-0.5)*0.36 for i in range(4)]
    ax.bar(xs,b,width=0.34,color=c,yerr=s,capsize=3,ecolor='#5f5e5a',
           label=f'arco de {"1/4 de órbita" if arc=="0.25" else "1 órbita"}')
ax.set_xticks(range(4)); ax.set_xticklabels(labs,fontsize=8.5)
ax.axhline(0,color='#1f1f1e',lw=0.8)
ax.set_ylabel('erro da massa (%)  — média ± 1σ')
ax.set_title('Estimativa de massa, corpo tipo Psyche, e = 0,1, 60 medidas, 150 Monte Carlo',fontsize=9.5,loc='left')
for i in (2,3):
    ax.annotate('< 0,01 %\n(ambos os arcos)',(i,0),xytext=(i,1.6),ha='center',fontsize=8,color='#5f5e5a')
ax.legend(frameon=False,fontsize=8.5)
fig.tight_layout(); fig.savefig('fig3_mass.png',dpi=200)
