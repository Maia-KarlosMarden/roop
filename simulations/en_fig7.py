import json, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9.5,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#c3c2b7','figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb'})
rms=lambda d:np.sqrt(d['bias']**2+d['sd']**2)
fig,ax=plt.subplots(1,2,figsize=(10,3.8))
# Bennu Orbit A, arcs
d=json.load(open('bd_bennu_OrbA.json')); T=d['T_h']; arcs=['0.125','0.25','0.5','1.0']; x=[float(a)*T for a in arcs]
for k,c,ls,lab in [('REL','#1baf7a','-o','Relational (SRP compensated)'),('OD_PM','#2a78d6','--s','Point-mass OD'),('OD_J2','#6250d6',':^','OD with free C20, C22')]:
    ax[0].plot(x,[rms(d['res'][a][k]) for a in arcs],ls,color=c,lw=2,ms=5,label=lab)
ax[0].axhline(1,color='#e34948',lw=1,ls='--'); ax[0].text(x[0],1.1,'1 %',color='#e34948',fontsize=8)
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xticks(x); ax[0].set_xticklabels([f'{v:.0f} h' for v in x])
ax[0].set_xlabel('data arc in Orbit A (r = 1.75 km)'); ax[0].set_ylabel('RMS mass error (%)')
ax[0].set_title('Bennu: relational ≈ OD up to half an orbit; long arcs need windows',fontsize=9,loc='left')
ax[0].legend(frameon=False,fontsize=7.5,loc='lower left'); from matplotlib.ticker import NullFormatter; ax[0].xaxis.set_minor_formatter(NullFormatter()); ax[0].grid(color='#e1e0d9',lw=0.6)
# Psyche vs altitude, full orbit
orbs=['OrbA','OrbB','OrbC','OrbD']; alt=[700,300,190,75]
od=[json.load(open(f'bd_psyche_{o}.json'))['res']['1.0'] for o in orbs]
rel_raw=[0.129,0.492,0.908,2.44]; rel_corr=[0.0032,0.0139,0.2008,0.827]
ax[1].plot(alt,rel_raw,'-o',color='#eb6834',lw=2,ms=5,label='Relational, raw')
ax[1].plot(alt,rel_corr,'-o',color='#1baf7a',lw=2,ms=5,label='Relational + closed-form J2 correction')
ax[1].plot(alt,[rms(o['OD_PM']) for o in od],'--s',color='#2a78d6',lw=2,ms=5,label='Point-mass OD')
ax[1].plot(alt,[rms(o['OD_J2']) for o in od],':^',color='#6250d6',lw=2,ms=5,label='OD with free C20, C22')
ax[1].set_yscale('log'); ax[1].invert_xaxis(); ax[1].set_xticks(alt); ax[1].set_xticklabels(['A\n700 km','B\n300 km','C\n190 km','D\n75 km'])
ax[1].set_xlabel('Psyche mission orbit (altitude)'); ax[1].set_title('Psyche: oblateness bias grows with descent; the J2 correction removes it',fontsize=9,loc='left')
ax[1].legend(frameon=False,fontsize=7.5,loc='upper left'); ax[1].grid(color='#e1e0d9',lw=0.6)
fig.tight_layout(); fig.savefig('fig7_bodies_en.png',dpi=200)
