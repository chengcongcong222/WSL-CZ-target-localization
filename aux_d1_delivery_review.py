"""Post-seal read-only initializer and cost supplement; no optimizer/performance runs."""
import time,math,json
import numpy as np
from aux_d1_common import *

def run():
 verify();seal=read('INFERENCE_SEAL.json');tic=time.perf_counter()
 for b in seal['bindings']:
  if sha(ROOT/b['path'])!=b['sha256']:raise RuntimeError('Inference bytes changed')
 # Cache public arrays only after inference ends; never supply private values to initialization.
 raw=np.load(OUT/'public_inputs/test.npz');pub={k:raw[k] for k in raw.files};raw.close();used=np.load(OUT/'STARTS_USED.npy',mmap_mode='r');a=np.load(OUT/'ALL_TERMINALS.npy',mmap_mode='r');outputs=[json.loads(l) for l in (OUT/'CONFIGURATION_OUTPUTS.jsonl').read_text().splitlines()];low=np.array([45,-5,1,-15]);high=np.array([60,5,3,15]);checks=[];maxdiff=0
 for o in outputs:
  c=o['configuration'];i=o['record_index'];expected=np.full(4,.5)
  if o['resource']=='U5':
   t=pub['times'];p0=pub['reported_main'][i];p1=pub['reported_aux'][i];alpha=pub['bearing_main'][i][:,0];beta=pub['bearing_aux'][i][:,0];m0=np.tan(alpha);m1=np.tan(beta);det=np.cos(alpha)*np.sin(beta)-np.sin(alpha)*np.cos(beta)
   with np.errstate(divide='ignore',invalid='ignore'):
    xx=(p1[:,1]-p0[:,1]+m0*p0[:,0]-m1*p1[:,0])/(m0-m1);yy=p0[:,1]+m0*(xx-p0[:,0]);l0=(xx-p0[:,0])/np.cos(alpha);l1=(xx-p1[:,0])/np.cos(beta);xy=np.column_stack([xx,yy])
   valid=(abs(det)>1e-12)&(l0>0)&(l1>0)&np.isfinite(xy).all(axis=1)
   if valid.sum()>=4 and np.ptp(t[valid])>0:
    T=np.column_stack([np.ones(valid.sum()),t[valid]/1200]);points=xy[valid];coef=np.linalg.lstsq(T,points,rcond=None)[0]
    for _ in range(5):
     norm=np.hypot(*(points-T@coef).T);s=max(10.,float(np.median(norm)));weight=np.minimum(1.,1.5*s/np.maximum(norm,1e-15));W=np.sqrt(weight);coef=np.linalg.lstsq(T*W[:,None],points*W[:,None],rcond=None)[0]
    p=coef[0]-p0[0];v=coef[1]/1200;h=np.array([math.hypot(*p)/1000,math.degrees(math.atan2(p[1],p[0])),math.hypot(*v),math.degrees(math.atan2(v[1],v[0]))]);expected=np.clip((np.clip(h,low,high)-low)/(high-low),1e-12,1-1e-12) if np.isfinite(h).all() else expected
  delta=float(np.max(abs(used[c,1]-expected)));maxdiff=max(maxdiff,delta);checks.append({'configuration':c,'PASS':delta<=1e-10,'normalized_start_max_difference':delta})
 costs=[]
 for r,R in enumerate(['U0','U5']):
  bank=a[r::2];executed=np.isfinite(bank[:,:,4]);costs.append({'resource':R,'executed_terminals':int(executed.sum()),'actual_nfev':int(np.nansum(bank[:,:,8])),'actual_njev':int(np.nansum(bank[:,:,9])),'geometry_calls':int(np.nansum(bank[:,:,12])),'sum_fit_wall_seconds':float(np.nansum(bank[:,:,15])),'bearing_samples_per_residual':121*(r+1),'CPU_seconds':'NOT_SEPARATELY_INSTRUMENTED','comparison':'Equalstarts/nfevlimitnot equalCPU or observationcost'})
 csvwrite('COMPUTATIONAL_COST_BY_RESOURCE.csv',costs);record_clock('initializer_and_cost_cold_supplement',time.perf_counter()-tic,optimizer_calls=0)
 dump('INITIALIZER_COLD_SUPPLEMENT.json',{'PASS':all(c['PASS'] for c in checks),'checks':checks,'maximum_normalized_difference':maxdiff,'method':'Independent tangent-slope intersections;same frozen directed-ray validity and5IRLS;no truth','new_optimizer_calls':0,'post_A_scope':'Read-only saved-array supplement,not changed scientific code or inference'})
 print('INITIALIZER_SUPPLEMENT',all(c['PASS'] for c in checks),maxdiff,flush=True)
 if not all(c['PASS'] for c in checks):raise SystemExit(1)
if __name__=='__main__':run()