"""Private generator entry. Truth, source and modal data never passed to RX."""
import csv,time,subprocess
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from hla_rx1_common import ROOT,OUT,LOCAL,read,dump,sha,table,state,platform,guard

def parse_mod(path,depths):
 b=Path(path).read_bytes();rec=4*int(np.frombuffer(b[:4],'<i4')[0]);h=np.frombuffer(b[84:108],'<i4');nt,nm=int(h[2]),int(h[3])
 z=np.frombuffer(b[4*rec:5*rec],'<f4')[:nt].astype(float);n=int(np.frombuffer(b[5*rec:5*rec+4],'<i4')[0])
 ph=np.stack([np.frombuffer(b[(7+j)*rec:(8+j)*rec],'<c8')[:nm] for j in range(n)],axis=1).astype(complex)
 k=np.frombuffer(b[(7+n)*rec:(7+n)*rec+n*8],'<c8').astype(complex)
 if nt!=nm or n<1 or not np.array_equal(z,depths) or not np.isfinite(k).all() or not np.isfinite(ph).all():raise ValueError('invalid exact-depth provider')
 return z,ph,k

def solve(config,start):
 directory=LOCAL/'modal';directory.mkdir(exist_ok=True);rows=[];mods={}
 for c in config['propagation_plan']:
  guard(start);stem=c['stem'];dest=directory/stem
  source=OUT/'environment'/f'{stem}.env';dest.with_suffix('.env').write_bytes(source.read_bytes());began=time.monotonic()
  try:
   p=subprocess.run([config['kraken_executable'],stem],cwd=directory,capture_output=True,timeout=90)
   dest.with_suffix('.stdout.txt').write_bytes(p.stdout+p.stderr)
   if p.returncode:raise RuntimeError('solver return '+str(p.returncode))
   mod=parse_mod(dest.with_suffix('.mod'),config['modal_depths_m']);mods[(c['line'],c['mesh'],c['offset_hz'])]=mod
   rows.append(dict(**c,execution_valid=True,seconds=time.monotonic()-began,sha256=sha(dest.with_suffix('.mod')),error=''))
  except Exception as e:
   rows.append(dict(**c,execution_valid=False,seconds=time.monotonic()-began,sha256='',error=repr(e)));table(OUT/'SOLVER_CALLS.csv',rows);raise
  table(OUT/'SOLVER_CALLS.csv',rows)
  print('PROPAGATION',len(rows),'/',len(config['propagation_plan']),flush=True)
 return mods

def group_provider(center,minus,plus,step):
 z,ph,k=center;indices=[];quality=[]
 for other in [minus,plus]:
  if other[1].shape!=ph.shape:raise ValueError('DRIFT_PROVIDER_NOT_READY: mode counts differ')
  a=abs(ph.conj().T@other[1])/(np.linalg.norm(ph,axis=0)[:,None]*np.linalg.norm(other[1],axis=0)[None,:]+1e-300)
  cost=1-a+.01*abs(k[:,None]-other[2][None,:])/np.maximum(np.median(np.diff(np.sort(k.real))),1e-8)
  rr,cc=linear_sum_assignment(cost);per=np.empty(len(k),int);per[rr]=cc
  quality.append(float(np.min(a[np.arange(len(k)),per])));indices.append(per)
 dk=(plus[2][indices[1]].real-minus[2][indices[0]].real)/(2*step)
 if min(quality)<.95 or not np.all(dk>0):raise ValueError('DRIFT_PROVIDER_NOT_READY: group mode identity not supported')
 return dk,quality

def clean_field(row,line,source_family,config,center,dk):
 # Adiabatic narrowband moving-source modal kernel. Solve each mode emission
 # u=t-r(u,t)/vg. Carrier correction prevents phase/delay double counting.
 z,ph,k=center;p0,v=state(row);zs=float(row['z_label']);iz=np.flatnonzero(z==zs)[0];ir=np.flatnonzero(z==200)[0]
 f0=config['private_source_frequency_hz'][line];lo=config['receiver_metadata']['digital_reference_hz'][line]
 slope=config['private_drift_hz_per_s'][line] if source_family=='DRIFT' else 0.
 t=np.arange(0,1200,1/config['receiver_metadata']['sample_rate_hz']);rec=platform(t)[:,None,:]+np.stack([np.array(config['receiver_metadata']['offsets_m']),np.zeros(8)],axis=-1)[None,:,:]
 weight=ph[iz]*ph[ir];vg=2*np.pi/dk;out=np.empty((len(t),8),complex)
 for a in range(0,len(t),128):
  tt=t[a:a+128,None,None];rr=rec[a:a+128,:,None,:]
  r=np.linalg.norm(p0+v*tt[...,None]-rr,axis=-1)
  u=tt-r/vg
  for _ in range(5):r=np.linalg.norm(p0+v*u[...,None]-rr,axis=-1);u=tt-r/vg
  phase=2*np.pi*(f0-lo)*tt-k[None,None,:]*r+np.pi*slope*u**2
  out[a:a+128]=np.sum(weight[None,None,:]*np.sqrt(2*np.pi/(k[None,None,:]*r))*np.exp(1j*phase-1j*np.pi/4),axis=-1)
 return out

def main():
 cfg=read(OUT/'DESIGN_FREEZE.json');start=time.monotonic();mods=solve(cfg,start)
 dump(LOCAL/'PROVIDER_INDEX.json',{'paths':[str(p) for p in (LOCAL/'modal').glob('*.mod')]})
if __name__=='__main__':main()
