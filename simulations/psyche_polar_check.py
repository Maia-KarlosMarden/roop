"""Polar vs retrograde B1->D descent check (Section 5.8): min radius, arrival, GM error.
Usage: python3 psyche_polar_check.py INC_DEG SEED [E_DAMP]   (E_DAMP=9 disables eccentricity damping)"""
import sys,json,numpy as np; sys.path.insert(0,'.')
import psyche_roop2 as P
inc=float(sys.argv[1]); seed=int(sys.argv[2]); P.INC_DEG=inc
if len(sys.argv)>3: P.E_DAMP=float(sys.argv[3])
r,log,tr=P.run('ROOP_J2',413.2e3,188.2e3,25.0,seed=seed)
tr=np.array(tr); i=np.argmin(tr[:,1])
print(json.dumps(dict(inc=inc,seed=seed,min_r_km=round(tr[i,1]/1e3,1),t_min_d=round(tr[i,0],2),arrive=r['t_arrive_d'],gm_err=round(r['gm_err_pct'],2))))
