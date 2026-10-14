#!/usr/bin/env python3
"""Cold verifier: independent PRT/MOD parsing and Prüfer IVP roots for fixed 201-Hz modes.
Does not import the main analysis or call KRAKEN / FIELD / SCOOTER.
"""
from pathlib import Path
import argparse,csv,hashlib,json,math,struct,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def readcsv(p):
 with Path(p).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for x in iter(lambda:f.read(1048576),b''):h.update(x)
 return h.hexdigest()
def do_prt(p):
 s=Path(p).read_text(encoding='utf-8');seen=False;out=[]
 for line in s.splitlines():
  if 'k (1/m)' in line and 'Phase Speed' in line:seen=True;continue
  if not seen:continue
  f=line.strip().split()
  if len(f)==5 and f[0].isdigit():
   try:out.append((int(f[0]),float(f[1]),float(f[2])))
   except ValueError:pass
 return out
def do_mod(p):
 buf=Path(p).read_bytes();step=struct.unpack_from('<i',buf)[0]*4
 m=struct.unpack_from('<i',buf,step*5)[0]
 start=(7+m)*step
 if start+8*m>len(buf):raise ValueError('unexpected mod size')
 return np.frombuffer(buf,dtype='<c8',count=m,offset=start).copy()
def env_speed(p):
 lines=Path(p).read_text(encoding='utf-8-sig').splitlines();idx=next(i for i,l in enumerate(lines) if l.startswith("'R'"));
 a=np.asarray([list(map(float,s.split()[:6])) for s in lines[5:idx]],float)
 return a[:,0],a[:,1]
def phase_at_end(f,z,c,k,step,rtol):
 omega=2*np.pi*f;scale=omega/1500
 def fun(depth,theta):
  t=theta[0];s=math.sin(t);v=math.cos(t);cz=float(np.interp(depth,z,c))
  return [scale*v*v+((omega/cz)**2-k*k)/scale*s*s]
 ans=solve_ivp(fun,(0.,float(z[-1])),[0.],method='DOP853',max_step=step,rtol=rtol,atol=rtol/100)
 if not ans.success or not np.isfinite(ans.y[0,-1]):raise ValueError('PRUFER_IVP_FAILED')
 return float(ans.y[0,-1])
def prufer_root(f,z,c,mode,initial,step,rtol,wall,started):
 target=(mode-.5)*np.pi
 def objective(k):
  if time.monotonic()-started>wall:raise RuntimeError('COLD_ROOT_WALL_BUDGET_STOP')
  return phase_at_end(f,z,c,k,step,rtol)-target
 width=3e-7
 for n in range(16):
  a=initial-width;b=initial+width
  if a<=0:raise RuntimeError('root bracket crossed zero')
  va=objective(a);vb=objective(b)
  if np.isfinite(va) and np.isfinite(vb) and va*vb<0:
   return float(brentq(objective,a,b,xtol=1e-14,rtol=5e-15,maxiter=80)),n+1
  width*=2
 raise RuntimeError(f'PRUFER_ROOT_NOT_BRACKETED mode={mode}, initial={initial}')
def audit(repo,config,out):
 began=time.monotonic();checks=[]
 def ch(name,error,tol,detail=''):
  checks.append({'check':name,'error':float(error),'tolerance':float(tol),'PASS':bool(np.isfinite(error) and error<=tol),'detail':detail})
 # Separate binary and text parsers; byte hashes pinned to immutable N2 archive.
 for entry in config['archive_sets']:
  sources={x:repo/entry['paths'][x] for x in ('env','prt','mod')}
  for k,p in sources.items():ch('source_sha256',0 if sha(p)==entry['sha256'][k] else 1,0,entry['stem']+'/'+k)
  prt=do_prt(sources['prt']);mod=do_mod(sources['mod'])
  ch('mode_count',abs(len(mod)-entry['modes']),0,entry['stem'])
  for mode,val,att in prt:
   ch('PRT_float32_matches_MOD',0 if np.float32(val).view(np.uint32)==np.float32(mod[mode-1].real).view(np.uint32) else 1,0,entry['stem']+'/'+str(mode))
   ch('attenuation_sign_reproduced',0 if np.sign(att)==np.sign(mod[mode-1].imag) else 1,0,entry['stem']+'/'+str(mode))
 # Verify saved main values using independent PRT input and MOD bits.
 main=readcsv(out/'PREC2_PRT_MOD_ALL_SAMPLES.csv')
 for row in main:
  entry=next(e for e in config['archive_sets'] if e['stem']==row['stem'])
  mod=do_mod(repo/entry['paths']['mod']);j=int(row['mode'])-1
  ch('main_k_MOD_reconstructed',abs(float(mod[j].real)-float(row['k_MOD_real'])),1e-14,row['stem']+'/'+row['mode'])
  ch('main_delta_rounded_reconstructed',abs(float(row['k_PRT_real'])-float(mod[j].real)-float(row['delta_k_PRT_minus_MOD'])),1e-14,row['stem']+'/'+row['mode'])
 # Independently shoot the continuous C-linear lossless Sturm-Liouville ODE.
 grid=readcsv(out/'PREC2_LOSSLESS_SUBSET_FDM.csv')
 selected=[r for r in grid if int(r['freq_hz'])==201]
 if len(selected)!=5:raise RuntimeError('missing 201Hz five subset modes')
 chosen=[selected[i] for i in config['cold_ivp_selected_rank_indices']]
 source=next(e for e in config['archive_sets'] if e['freq']==201 and e['mesh']==160001)
 z,c=env_speed(repo/source['paths']['env'])
 indep=[]
 for r in chosen:
  mode=int(r['mode']);guess=float(r['K_real_lossless_R6']);
  loose,n1=prufer_root(201,z,c,mode,guess,1.,1e-11,config['cold_ivp_wall_s'],began)
  tight,n2=prufer_root(201,z,c,mode,guess,.5,1e-12,config['cold_ivp_wall_s'],began)
  delta=abs(tight-loose);dr=abs(tight-guess);limit=config['phase_target_radians']/config['phase_reference_range_m']
  indep.append({'freq_hz':201,'mode':mode,'K_FDM_R6':guess,'K_PRUFER_loose':loose,'K_PRUFER_tight':tight,
                'prufer_step_delta':delta,'PRUFER_vs_FDM_delta':dr,'target_delta_k':limit,
                'bracket_expansions_loose':n1,'bracket_expansions_tight':n2,
                'IVP_precision_pass':delta<=limit,'independent_agreement_pass':dr<=limit})
  ch('IVP_two_tolerances_agree',delta,limit,str(mode))
  ch('IVP_vs_FDM_R6',dr,limit,str(mode))
 fields=list(indep[0])
 with (out/'PREC2_PRUFER_COLD_ROOTS.csv').open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(indep)
 with (out/'PREC2_COLD_CHECKS.csv').open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(checks[0]),lineterminator='\n');w.writeheader();w.writerows(checks)
 fail=sum(not c['PASS'] for c in checks)
 result={'PASS':fail==0,'checks':len(checks),'FAIL':fail,'cold_ivp_roots':len(indep),
         'scope':'PRT/MOD identity and subset lossless C-linear real-k IVP vs FDM; not full complex propagation',
         'max_IVP_vs_FDM_delta':max(x['PRUFER_vs_FDM_delta'] for x in indep),
         'new_acoustic_solvers_or_MC':0,'wall_s':time.monotonic()-began,'R4_percent':0,'state':'STOP'}
 (out/'PREC2_COLD_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print(json.dumps(result,indent=2));
 if fail:raise RuntimeError('COLD_NOT_ADMITTED; preserve all outputs')
if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--repo',type=Path,required=True);pa.add_argument('--config',type=Path,required=True);pa.add_argument('--out',type=Path,required=True)
 a=pa.parse_args();audit(a.repo.resolve(),json.loads(a.config.read_text(encoding='utf-8')),a.out.resolve())
