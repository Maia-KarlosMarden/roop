import sys,json; sys.path.insert(0,'.')
import sim_market as S
for de0 in [0.0,3e-4]:
    res,orb=S.run('NONE',days=40,seed=1,de0=de0)
    print(de0,res)
    for o in orb[::60]: print("  t=%.1f ex=%.2e ey=%.2e |de|=%.2e"%(o[0],o[2],o[3],o[4]))
