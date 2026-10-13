"""AUX-D1 private generator and zero-optimizer independent calibration."""
import argparse,time,math
import numpy as np
from aux_d1_common import *
from aux_d1_model import LOW,HIGH,SPAN,SIGMA,score,polar_state,wrap

def nominal(t):
 post=np.maximum(t-600.,0.);return np.column_stack([2*np.minimum(t,600)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)])
def make_record(h,error_seed,g,rep,side):
 t=np.arange(121,dtype=float)*10;main=nominal(t)
 rng=lambda k:np.random.default_rng(np.random.SeedSequence([error_seed,g,rep,k]))
 c,d,beta=rng(0).uniform([-.05,-.025,-10],[.05,.025,10]);an=np.column_stack([rng(1).standard_normal(121),rng(2).standard_normal(121)]);nav=np.stack([rng(3).standard_normal((121,2)),rng(4).standard_normal((121,2))],axis=1)
 off=5000*np.array([-np.sin(np.deg2rad(beta)),side*np.cos(np.deg2rad(beta))]);trueP=np.stack([main,main+off],axis=1)
 target=1000*h[0]*np.array([np.cos(np.deg2rad(h[1])),np.sin(np.deg2rad(h[1]))])+t[:,None]*h[2]*np.array([np.cos(np.deg2rad(h[3])),np.sin(np.deg2rad(h[3]))])
 dif=target[:,None,:]-trueP;observed=wrap(np.arctan2(dif[:,:,1],dif[:,:,0])+np.deg2rad([c-d,c+d])+SIGMA*an);reported=trueP+25*nav
 return {'times':t,'bearing_main':observed[:,:1].copy(),'bearing_aux':observed[:,1:].copy(),'reported_main':reported[:,0].copy(),'reported_aux':reported[:,1].copy()}, {'h_truth':h,'systematics':np.array([c,d,beta]),'side':side,'angle_standard_normals':an,'navigation_standard_normals':nav,'true_node_positions':trueP}

def bundle(records):
 public,private=zip(*records);pub={'times':public[0]['times']}
 for k in ['bearing_main','bearing_aux','reported_main','reported_aux']:pub[k]=np.stack([r[k] for r in public])
 priv={k:np.array([r[k] for r in private]) for k in private[0]};return pub,priv

def record_hashes(pub):
 data=[]
 for i in range(len(pub['bearing_main'])):
  def digest(keys):
   h=hashlib.sha256();h.update(pub['times'].tobytes())
   for k in keys:h.update(np.ascontiguousarray(pub[k][i]).tobytes())
   return h.hexdigest()
  data.append({'record_index':i,'MAIN_input_sha256':digest(['bearing_main','reported_main']),'AUX_input_sha256':digest(['bearing_aux','reported_aux'])})
 return data

def save(name,pub,priv):
 pubpath=OUT/'public_inputs'/f'{name}.npz';privpath=OUT/'evaluation_only'/f'{name}_private.npz';pubpath.parent.mkdir(exist_ok=True);privpath.parent.mkdir(exist_ok=True)
 np.savez_compressed(pubpath,**pub);np.savez_compressed(privpath,**priv)
 return pubpath,privpath

def run(phase):
 verify();spec=read('AUX_D1_TASK_SPEC_20261012_v1.json');tic=time.perf_counter()
 sentinel=OUT/(phase.upper()+'_START.json')
 if sentinel.exists():raise RuntimeError('Once-only generator already started; no redraw authorized')
 dump(sentinel.name,{'utc':utc(),'phase':phase,'task_SHA256':TASKSHA})
 if phase=='calibration':
  h=np.random.default_rng(2026101201).uniform(LOW,HIGH,size=(512,4));pub,priv=bundle([make_record(h[i],2026101202,i,0,-1 if i%2==0 else 1) for i in range(512)]);pubpath,privpath=save('calibration',pub,priv)
  scores=[]
  for i in range(512):
   tt=pub['times'];hh=priv['h_truth'][i];target=1000*hh[0]*np.array([np.cos(np.deg2rad(hh[1])),np.sin(np.deg2rad(hh[1]))])+tt[:,None]*hh[2]*np.array([np.cos(np.deg2rad(hh[3])),np.sin(np.deg2rad(hh[3]))])
   P=np.stack([pub['reported_main'][i],pub['reported_aux'][i]],axis=1);y=np.concatenate([pub['bearing_main'][i],pub['bearing_aux'][i]],axis=1);d=target[:,None,:]-P;e=wrap(np.arctan2(d[:,:,1],d[:,:,0])-y)/SIGMA
   hr=polar_state(np.r_[target[0],hh[2]*np.cos(np.deg2rad(hh[3])),hh[2]*np.sin(np.deg2rad(hh[3]))],P[0,0]);inside=bool(np.all((hr>=LOW)&(hr<=HIGH)))
   for R,n in [('U0',1),('U5',2)]:scores.append({'calibration_index':i,'resource':R,'Q':format_float(score(e[:,:n].ravel())),'truth_report_reference_in_domain':inside})
  csvwrite('CALIBRATION_SCORES.csv',scores);th={};sortids={}
  for R in ['U0','U5']:
   rec=[r for r in scores if r['resource']==R];order=sorted(rec,key=lambda r:(float(r['Q']),int(r['calibration_index'])));th[R]=float(order[499]['Q']);sortids[R]=[int(r['calibration_index']) for r in order]
  dump('CALIBRATION_THRESHOLDS.json',{'thresholds':th,'nearest_rank_1based':500,'samples_per_resource':512,'quantile':.975,'sorted_calibration_indices':sortids,'optimizer_calls':0,'scope':'Empirical screening threshold;not pergeometry97.5%coverage','reference':'Physicaltruth target versus reported positions;no clipping'})
  dump('THRESHOLD_FREEZE.json',{'utc':utc(),'thresholds':th,'bindings':[filemeta(pubpath),filemeta(privpath),filemeta(OUT/'CALIBRATION_SCORES.csv'),filemeta(OUT/'CALIBRATION_THRESHOLDS.json')],'test_generation_started':False})
  dump('PUBLIC_OBSERVATION_MANIFEST.json',{'coordinate':'globalx,y;angleatan2counterclockwisefrom+x;rad','position_unit':'m','time_unit':'s','public_only':'reportedpositions andobservedbearings','arrays':[filemeta(pubpath)],'calibration_records':record_hashes(pub)})
  dump('PRIVATE_GENERATOR_MANIFEST.json',{'truth_access':'Generator/calibrator andsealedindependentevaluation/auditonly;notestimator','arrays':[filemeta(privpath)],'calibration_geometry_seed':2026101201,'calibration_error_seed':2026101202,'streams':'SeedSequence([errorseed,geometry,rep,stream0..4]);uniform3,mainangle121,auxangle121,mainnav121x2,auxnav121x2'})
 else:
  if not (OUT/'THRESHOLD_FREEZE.json').exists():raise RuntimeError('thresholdfreeze required')
  for p in read('THRESHOLD_FREEZE.json')['bindings']:
   if sha(ROOT/p['path'])!=p['sha256']:raise RuntimeError('calibrationbindingchanged')
  import csv
  panel=list(csv.DictReader((ROOT/'results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_TRUTH_PANEL.csv').read_text(encoding='utf-8-sig').splitlines()));hh=np.array([[float(r[k]) for k in ['r0_km','theta0_deg','v_mps','psi_deg']] for r in panel]);assert len(hh)==12 and np.all((hh>=LOW)&(hh<=HIGH))
  records=[];plan=[]
  for rep in range(128):
   for g in range(12):records.append(make_record(hh[g],2026101203,g,rep,-1 if rep%2==0 else 1));plan.append({'record_index':len(records)-1,'geometry_index':g,'case_id':panel[g]['case_id'],'replicate':rep,'side':-1 if rep%2==0 else 1})
  pub,priv=bundle(records);pubpath,privpath=save('test',pub,priv);csvwrite('PAIRED_RECORD_INDEX.csv',plan)
  for name,path,key in [('PUBLIC_OBSERVATION_MANIFEST.json',pubpath,'test_records'),('PRIVATE_GENERATOR_MANIFEST.json',privpath,None)]:
   x=read(name);x['arrays'].append(filemeta(path));
   if key:x[key]=record_hashes(pub)
   dump(name,x)
  dump('TEST_GENERATION_SEAL.json',{'utc':utc(),'public_test':filemeta(pubpath),'private_test':filemeta(privpath),'threshold_freeze_sha256':sha(OUT/'THRESHOLD_FREEZE.json'),'plan':filemeta(OUT/'PAIRED_RECORD_INDEX.csv'),'records':1536,'test_error_seed':2026101203})
 duration=time.perf_counter()-tic;record_clock(phase,duration,optimizer_runs=0,records=512 if phase=='calibration' else 1536)
 print(phase,'COMPLETE',duration,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('phase',choices=['calibration','test']);run(p.parse_args().phase)