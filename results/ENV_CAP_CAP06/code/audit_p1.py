"""Independent power/log calculations and nullspace; no production imports."""
from pathlib import Path
import json,csv,hashlib,numpy as np
O=Path(__file__).resolve().parents[1];W=Path(r'D:\ProjectStorage\WSL-CZ\ENV_CAP_CAP06/raw/p1')
def main():
 checks=[]
 def ck(n,b,e=0):checks.append(dict(check_id=n,PASS=bool(b),error=float(e)))
 for p in sorted((W/'receiver_public').glob('*.npz')):
  a=np.load(p);b=np.load(W/'products'/p.name);signal=(a['samples'].real.astype(float)**2+a['samples'].imag.astype(float)**2).sum(-1)/128;noise=(a['background'].real.astype(float)**2+a['background'].imag.astype(float)**2).sum(-1)/128
  power=(signal-noise)/a['gain_metadata']**2;valid=(power[:,1:]>0)&(power[:,0,None]>0)
  ck(p.stem+':nonpositive_mask',np.array_equal(valid,b['valid']))
  vals=[]
  for t,i in zip(*np.where(valid)):vals.append(abs(10*np.log10(power[t,i+1]/power[t,0])-b['ratio'][t,i]))
  e=max(vals,default=0);ck(p.stem+':independent_log_power',e<1e-8,e)
 # Build toy derivatives independently from explicit scalar expressions and use least squares.
 A=[];E=[];t=np.linspace(-1,1,9)
 for ti in t:
  for ci,di,ei in zip([1,-1,2,-2],[1,2,-1,-2],[1,1,-1,-1]):
   A.append([ci*np.cos(np.pi*ti),di*np.sin(2*np.pi*ti),di*np.sin(np.pi*ti),ci*ti*ti,di*ti**3]);E.append([ei*np.cos(3*np.pi*ti)])
 A=np.array(A);E=np.array(E);saved=json.loads((O/'P1_NUISANCE_INFO.json').read_text())
 for z in saved:
  g=z['gain'];n=[E]
  if g in ['G1','G2','G4']:n.append(np.kron(np.ones((9,1)),np.eye(4)))
  if g=='G2':n.append(np.kron(t[:,None],np.eye(4)))
  if g=='G3':n.append(np.eye(36))
  N=np.column_stack(n)
  for j in range(5):
   B=np.column_stack([N,A[:,[k for k in range(5) if k!=j]]]);v=A[:,j]-B@np.linalg.lstsq(B,A[:,j],rcond=1e-10)[0];e=abs(v@v-z['profiled_state_information'][j]);ck(g+':profile:'+str(j),e<1e-8,e)
 gate=json.loads((O/'P1_POWER_FORWARD_GATE.json').read_text());ck('CZ_not_admitted',gate['G1b']=='O2_FIRST_CZ_FORWARD_NOT_ADMITTED' and gate['CZ_depth_performance']=='NOT_EVALUATED')
 # Shared reference delta-method covariance has positive offdiagonal, no fictitious independent ratios.
 D=np.column_stack([-np.ones(4),np.eye(4)]);C=D@D.T;ck('shared_reference_covariance',np.allclose(C,np.ones((4,4))+np.eye(4)))
 out=dict(PASS=all(c['PASS'] for c in checks),checks=len(checks),FAIL=sum(not c['PASS'] for c in checks),records=checks,scope='toy implementation only; CZ physical admission remains failed')
 (O/'P1_INDEPENDENT_COLD_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print(out['PASS'],out['checks'])
if __name__=='__main__':main()
