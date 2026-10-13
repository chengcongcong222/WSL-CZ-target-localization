"""Independent archived-input/Cartesian/selection/summary audit. No optimizer calls."""
import math,time,json,hashlib
import numpy as np
from aux_d1_common import *

def scalar_error(h,true):
 if h is None:return [math.inf]*4
 return [abs(h[0]-true[0])/true[0],abs((h[1]-true[1]+180)%360-180),abs(h[2]-true[2])/true[2],abs((h[3]-true[3]+180)%360-180)]
def nr(v,p):return sorted(v)[math.ceil(p*len(v))-1] if v else None

def audit():
 verify();tic=time.perf_counter();checks=[];unresolved=[];max_q_difference=0.;ranking_sensitive=0
 def check(n,yes,detail=''):checks.append({'check':n,'PASS':bool(yes),'detail':str(detail)})
 def near(x,y):return x==y if x is None or y is None or math.isinf(float(x)) or math.isinf(float(y)) else math.isclose(float(x),float(y),rel_tol=1e-10,abs_tol=1e-10)
 def independent_q(h,t,P,y):
  h=np.atleast_2d(h);th=np.deg2rad(h[:,1]);ps=np.deg2rad(h[:,3]);origin=P[0,0]
  xy=origin+np.stack([1000*h[:,0]*np.cos(th),1000*h[:,0]*np.sin(th)],axis=1);vel=np.stack([h[:,2]*np.cos(ps),h[:,2]*np.sin(ps)],axis=1);target=xy[:,None,:]+t[None,:,None]*vel[:,None,:];delta=target[:,:,None,:]-P[None,:,:,:];angles=np.arctan2(delta[...,1],delta[...,0]);err=np.arctan2(np.sin(angles-y),np.cos(angles-y))/math.radians(.05);f=err.reshape(len(h),-1)
  return np.sum(2*f*f/(np.sqrt(1+(f/1.5)**2)+1),axis=1,dtype=np.float64)/f.shape[1]
 design=read('DESIGN_FREEZE.json');threshold=read('CALIBRATION_THRESHOLDS.json')['thresholds'];manifest=read('PUBLIC_OBSERVATION_MANIFEST.json');private_manifest=read('PRIVATE_GENERATOR_MANIFEST.json')
 for b in manifest['arrays']+private_manifest['arrays']+read('THRESHOLD_FREEZE.json')['bindings']+read('INFERENCE_SEAL.json')['bindings']:check('hash:'+b['path'],sha(ROOT/b['path'])==b['sha256'])
 check('frozen_before_test',read('THRESHOLD_FREEZE.json')['utc']<read('TEST_GENERATION_SEAL.json')['utc']<read('INFERENCE_START.json')['utc']<read('INFERENCE_SEAL.json')['utc']<read('EVALUATION_START.json')['utc'])
 for phase,n,seed in [('calibration',512,2026101202),('test',1536,2026101203)]:
  pub=np.load(OUT/f'public_inputs/{phase}.npz');priv=np.load(OUT/f'evaluation_only/{phase}_private.npz');tt=pub['times'];check(phase+'times',np.array_equal(tt,np.arange(121)*10.))
  if phase=='calibration':check('calibration_geometry',np.array_equal(priv['h_truth'],np.random.default_rng(2026101201).uniform([45,-5,1,-15],[60,5,3,15],size=(512,4))))
  else:
   import csv
   truth=list(csv.DictReader((ROOT/'results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_TRUTH_PANEL.csv').read_text(encoding='utf-8-sig').splitlines()));hp=np.array([[float(r[k]) for k in ['r0_km','theta0_deg','v_mps','psi_deg']] for r in truth]);check('test_alloriginal12',np.array_equal(priv['h_truth'],np.tile(hp,(128,1))))
  for i in range(n):
   g,rep=(i,0) if phase=='calibration' else (i%12,i//12);side=-1 if (i if phase=='calibration' else rep)%2==0 else 1
   rng=lambda k:np.random.default_rng(np.random.SeedSequence([seed,g,rep,k]));c,d,beta=rng(0).uniform([-.05,-.025,-10],[.05,.025,10]);an=np.stack([rng(1).standard_normal(121),rng(2).standard_normal(121)],axis=1);nav=np.stack([rng(3).standard_normal((121,2)),rng(4).standard_normal((121,2))],axis=1)
   check(phase+str(i)+'innovations',np.array_equal(priv['systematics'][i],[c,d,beta]) and np.array_equal(priv['angle_standard_normals'][i],an) and np.array_equal(priv['navigation_standard_normals'][i],nav) and priv['side'][i]==side)
   main=np.array([[2*float(t),0] if t<=600 else [1200+2*(t-600)*math.cos(math.pi/12),2*(t-600)*math.sin(math.pi/12)] for t in tt]);offset=np.array([-5000*math.sin(math.radians(beta)),side*5000*math.cos(math.radians(beta))]);nodes=np.stack([main,main+offset],axis=1)
   check(phase+str(i)+'nodeformation',np.max(abs(priv['true_node_positions'][i]-nodes))<1e-10)
   P=np.stack([pub['reported_main'][i],pub['reported_aux'][i]],axis=1);check(phase+str(i)+'reported_nav',np.max(abs(P-(nodes+25*nav)))<1e-10)
   h=priv['h_truth'][i];target=np.array([1000*h[0]*math.cos(math.radians(h[1])),1000*h[0]*math.sin(math.radians(h[1]))])+tt[:,None]*np.array([h[2]*math.cos(math.radians(h[3])),h[2]*math.sin(math.radians(h[3]))]);delta=target[:,None,:]-nodes;orig=np.arctan2(delta[:,:,1],delta[:,:,0])+np.deg2rad([c-d,c+d])+math.radians(.05)*an;observed=np.concatenate([pub['bearing_main'][i],pub['bearing_aux'][i]],axis=1);diff=np.arctan2(np.sin(observed-orig),np.cos(observed-orig));check(phase+str(i)+'bearing_component',np.max(abs(diff))<1e-12)
   hasher=hashlib.sha256();hasher.update(tt.tobytes());hasher.update(np.ascontiguousarray(pub['bearing_main'][i]).tobytes());hasher.update(np.ascontiguousarray(pub['reported_main'][i]).tobytes());saved=manifest[phase+'_records'][i];check(phase+str(i)+'public_MAIN_identity',hasher.hexdigest()==saved['MAIN_input_sha256'])
 # Calibration score replay uses physicaltruth, not clipped reported coordinate state.
 cal=np.load(OUT/'public_inputs/calibration.npz');cp=np.load(OUT/'evaluation_only/calibration_private.npz');crows=csvread('CALIBRATION_SCORES.csv')
 for i,h in enumerate(cp['h_truth']):
  P=np.stack([cal['reported_main'][i],cal['reported_aux'][i]],axis=1);y=np.concatenate([cal['bearing_main'][i],cal['bearing_aux'][i]],axis=1);t=cal['times'];pos=np.array([1000*h[0]*math.cos(math.radians(h[1])),1000*h[0]*math.sin(math.radians(h[1]))])+t[:,None]*np.array([h[2]*math.cos(math.radians(h[3])),h[2]*math.sin(math.radians(h[3]))]);d=pos[:,None,:]-P;error=np.arctan2(np.sin(np.arctan2(d[:,:,1],d[:,:,0])-y),np.cos(np.arctan2(d[:,:,1],d[:,:,0])-y))/math.radians(.05)
  for ri,R in enumerate(['U0','U5']):
   v=error[:,:ri+1].ravel();q=math.fsum(2*float(e)**2/(math.sqrt(1+(float(e)/1.5)**2)+1) for e in v)/len(v);check('calQ'+str(i)+R,abs(q-float(crows[2*i+ri]['Q']))<=1e-10*max(1,abs(q)))
 for R in ['U0','U5']:check('cal500_'+R,sorted(float(r['Q']) for r in crows if r['resource']==R)[499]==threshold[R])
 # Every saved terminal score and strict original-score selection.
 pub=np.load(OUT/'public_inputs/test.npz');priv=np.load(OUT/'evaluation_only/test_private.npz');t=pub['times'];a=np.load(OUT/'ALL_TERMINALS.npy',mmap_mode='r');starts=np.load(OUT/'STARTS_USED.npy',mmap_mode='r');fixed=np.load(OUT/'SHARED_FIXED_126.npy');outputs=[json.loads(l) for l in (OUT/'CONFIGURATION_OUTPUTS.jsonl').read_text().splitlines()];est={int(r['configuration']):r for r in csvread('ESTIMATES.csv')};evalrows={int(r['configuration']):r for r in csvread('TRUTH_COMPATIBILITY_EVAL_ONLY.csv')};inderrors={};shared_errors={};sideerrors={};selection_risks={};low=np.array([45,-5,1,-15]);high=np.array([60,5,3,15]);scale=np.array([.05,.01,.02,.2]);limits=np.array([.1,1,.1,5.])
 for o in outputs:
  c=o['configuration'];i=o['record_index'];ri=0 if o['resource']=='U0' else 1;R=o['resource'];bank=a[c];P=np.stack([pub['reported_main'][i],pub['reported_aux'][i]],axis=1)[:,:ri+1];y=np.concatenate([pub['bearing_main'][i],pub['bearing_aux'][i]],axis=1)[:,:ri+1];finite=np.isfinite(bank[:,:6]).all(axis=1)
  qr=independent_q(bank[finite,:4],t,P,y);stored=bank[finite,4];diff=np.abs(qr-stored);max_q_difference=max(max_q_difference,float(diff.max(initial=0)));check('terminalQ'+str(c),np.all(diff<=1e-10*np.maximum(1,abs(stored))))
  valid=(bank[:,6]==1)&finite&np.all((bank[:,:4]>=low)&(bank[:,:4]<=high),axis=1)&(bank[:,4]<=threshold[R]);ids=list(np.flatnonzero(valid));selected=min(ids,key=lambda k:(float(bank[k,4]),int(k))) if ids else -1;sharedids=[k for k in ids if k!=1];shared=min(sharedids,key=lambda k:(float(bank[k,4]),int(k))) if sharedids else -1
  check('selection'+str(c),selected==o['selected_start'] and shared==o['shared_selected_start'] and ids==o['accepted_indices'])
  check('residualbudget'+str(c),np.all(bank[:,8]<=100)&np.all(bank[:,10]<=100)&np.all(bank[:,8]==bank[:,10]))
  check('scorecost'+str(c),np.all(abs(bank[finite,4]-2*bank[finite,5]/(121*(ri+1)))<1e-10*np.maximum(1,abs(bank[finite,4]))))
  if ri==1:check('sharedstartbytes'+str(i),starts[c,[0,*range(2,128)]].tobytes()==starts[c-1,[0,*range(2,128)]].tobytes())
  h=priv['h_truth'][i];ee=scalar_error(bank[selected,:4] if selected>=0 else None,h);se=scalar_error(bank[shared,:4] if shared>=0 else None,h);inderrors[c]=ee;shared_errors[c]=se;sideerrors[c]=int(priv['side'][i]);row=est[c]
  check('fixedstartarray'+str(c),np.array_equal(starts[c,2:],np.clip((fixed-low)/(high-low),1e-12,1-1e-12)))
  nm=np.stack([-np.sin(y[:,0]),np.cos(y[:,0])],axis=1);matrix=np.column_stack([nm,nm*t[:,None]]);ss=np.array([50000.,50000.,2.,2.]);coef,_,rank,_=np.linalg.lstsq(matrix*ss,np.sum(nm*P[:,0],axis=1),rcond=1e-12);cart=coef*ss;dd=cart[:2]-P[0,0];hh=np.array([math.hypot(*dd)/1000,math.degrees(math.atan2(dd[1],dd[0])),math.hypot(*cart[2:]),math.degrees(math.atan2(cart[3],cart[2]))]);base=hh if rank==4 and np.isfinite(hh).all() else (low+high)/2;uu=np.clip((np.clip(base,low,high)-low)/(high-low),1e-12,1-1e-12);check('maininitializer'+str(c),np.max(abs(uu-starts[c,0]))<1e-10)
  xy=np.array([1000*h[0]*math.cos(math.radians(h[1])),1000*h[0]*math.sin(math.radians(h[1]))]);vxy=np.array([h[2]*math.cos(math.radians(h[3])),h[2]*math.sin(math.radians(h[3]))]);pt=xy+t[:,None]*vxy;origin=xy-P[0,0];hr=np.array([math.hypot(*origin)/1000,math.degrees(math.atan2(origin[1],origin[0])),h[2],h[3]]);inside=bool(np.all((hr>=low)&(hr<=high)));delta=pt[:,None,:]-P;ang=np.arctan2(delta[:,:,1],delta[:,:,0]);er=np.arctan2(np.sin(ang-y),np.cos(ang-y))/math.radians(.05);qt=math.fsum(2*float(e)**2/(math.sqrt(1+(float(e)/1.5)**2)+1) for e in er.ravel())/er.size;tr=evalrows[c];check('truthQ'+str(c),abs(qt-float(tr['truth_Q']))<=1e-10*max(1,abs(qt)));check('truthdomain'+str(c),inside==(tr['truth_report_reference_in_domain']=='True'));check('truthdeliveryflag'+str(c),(inside and float(tr['truth_Q'])<=threshold[R] and selected<0)==(tr['TRUE_COMPATIBLE_BUT_NO_OUTPUT']=='True'))
  check('point_error'+str(c),all(near(float(row[p+'_error']),ee[j]) for j,p in enumerate(['r','theta','v','psi'])))
  if selected>=0:check('CSVfloat64literal'+str(c),all(float(row[p+'_hat'])==bank[selected,j] for j,p in enumerate(['r','theta','v','psi'])))
  # Greedy clustering independently, no truth centers.
  reps=[];labels=[]
  for k in sorted(ids,key=lambda k:(float(bank[k,4]),int(k))):
   label=-1
   for n,j in enumerate(reps):
    d=abs(bank[k,:4]-bank[j,:4]);d[[1,3]]=abs((bank[k,[1,3]]-bank[j,[1,3]]+180)%360-180)
    if max(d/scale)<=1:label=n;break
   if label<0:label=len(reps);reps.append(int(k))
   labels.append([int(k),label])
  check('cluster'+str(c),reps==o['envelope']['representative_start_ids'] and labels==o['envelope']['assignment'])
  risk=selected>=0 and any(ee[j]>limits[j] for j in range(4)) and any(all(x<=y for x,y in zip(scalar_error(bank[k,:4],h),limits)) for k in ids);selection_risks[c]=risk;check('selectionrisk'+str(c),risk==(row['better_terminal_not_selected']=='True'))
  if selected>=0 and len(ids)>1:
   ordered=sorted(ids,key=lambda k:(float(bank[k,4]),int(k)));margin=float(bank[ordered[1],4]-bank[ordered[0],4]);ranking_sensitive+=margin<=2e-10*max(1,float(bank[selected,4]))
  changed=np.flatnonzero((qr<=threshold[R])!=(stored<=threshold[R]));
  if len(changed):unresolved.append({'configuration':c,'terminal_indices':np.flatnonzero(finite)[changed].tolist(),'scope':'THRESHOLD_ARITHMETIC_ONLY;originalacceptancekept'})
 # Every saved conditional/unconditional quantile reconstructed from the original selected terminal.
 case=csvread('PAIRED_RESULTS_BY_CASE.csv');side=csvread('PAIRED_RESULTS_BY_SIDE.csv')
 for r in case+side:
  ids=[c for c in inderrors.keys() if c%2==(0 if r['resource']=='U0' else 1) and (c//2)%12==int(r['geometry_index']) and ('side' not in r or sideerrors[c]==int(r['side']))];collection=shared_errors if r.get('variant')=='SHARED' else inderrors;accepted=[c for c in ids if math.isfinite(collection[c][0])]
  check('denominator'+str(r.get('variant',r.get('side')))+r['resource']+r['geometry_index'],len(ids)==int(r['evaluated']) and len(accepted)==int(r['accepted']))
  for j,p in enumerate(['r','theta','v','psi']):
   for label,group in [('unconditional',ids),('conditional',accepted)]:
    for suffix,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99),('max',1.)]:
     expected=nr([collection[c][j] for c in group],q);token=r[p+'_'+label+'_'+suffix];check('quantile'+str(r.get('variant',r.get('side')))+r['resource']+r['geometry_index']+p+label+suffix,near(float(token) if token else None,expected))
 summary=csvread('SUMMARY_BY_RESOURCE.csv')
 for r in summary:
  cr=[x for x in case if x['variant']=='PRIMARY' and x['resource']==r['resource']]
  for p in ['r','theta','v','psi']:
   for cond in ['unconditional','conditional']:
    vals=[float(x[p+'_'+cond+'_median']) for x in cr if x[p+'_'+cond+'_median']];expected=sum(vals)/12 if len(vals)==12 else None;check('meanmedian'+r['resource']+p+cond,near(float(r[p+'_'+cond+'_mean_geometry_median']) if r[p+'_'+cond+'_mean_geometry_median'] else None,expected));vals=[float(x[p+'_'+cond+'_P95']) for x in cr if x[p+'_'+cond+'_P95']];check('worstP95'+r['resource']+p+cond,near(float(r[p+'_'+cond+'_worst_geometry_P95']) if r[p+'_'+cond+'_worst_geometry_P95'] else None,max(vals) if len(vals)==12 else None))
 # Main signal independently from saved primary case metrics.
 signal={p:[] for p in ['r','v','psi']};qualified=set()
 for g in range(12):
  aa=next(r for r in case if r['variant']=='PRIMARY' and r['resource']=='U0' and int(r['geometry_index'])==g);bb=next(r for r in case if r['variant']=='PRIMARY' and r['resource']=='U5' and int(r['geometry_index'])==g);valid=int(aa['evaluated'])==128 and int(bb['evaluated'])==128 and float(bb['output_rate'])>=.9 and float(bb['output_rate'])>=float(aa['output_rate'])-.02
  for p in signal:
   x=float(aa[p+'_unconditional_median']) if aa[p+'_unconditional_median'] else math.inf;y=float(bb[p+'_unconditional_median']) if bb[p+'_unconditional_median'] else math.inf
   if valid and math.isfinite(x) and x>0 and math.isfinite(y) and 1-y/x>=.3:signal[p].append(g);qualified.add(g)
 decision=read('DECISION.json');check('primarysignals',signal==decision['primary']['signal_geometries_by_parameter'] and sorted(qualified)==decision['primary']['qualifying_geometry_indices'])
 status=read('INFERENCE_SEAL.json');check('terminalcounts',sum(np.sum(np.isfinite(a[c,:,4])) for c in range(3072))==status['terminals'])
 check('science_boundaries',decision['new_propagation_calls']==0 and decision['R4_percent']==0 and decision['next_stage_execution_authorized'] is False)
 check('storage_RSS',storage_bytes()<=2*1024**3 and rss()<=8*1024**3)
 changed=__import__('subprocess').check_output(['git','diff','--name-only',PARENT],cwd=ROOT,text=True,encoding='utf-8').splitlines();check('historicalfilesunchanged',all(p.startswith('aux_d1_') or p.startswith('results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/') or p=='results/R4_MASTER/R4_EVIDENCE_LEDGER.csv' for p in changed))
 elapsed=time.perf_counter()-tic;record_clock('cold_audit',elapsed,optimizer_calls=0)
 result={'stage':'AUX-D1','PASS':all(r['PASS'] for r in checks),'checks':checks,'FAIL':[r for r in checks if not r['PASS']],'maximum_Q_absolute_difference':max_q_difference,'threshold_arithmetic_unresolved':unresolved,'near_tied_ranking_configurations':ranking_sensitive,'selection_rule':'Original float64 Q/start_id;independent geometric ties do not replacepoints','scientific_scope':'Cold arithmetic/input/summaryrebuild;noglobalcoverageorrealhardwarecertification','optimizer_calls':0,'elapsed_s':elapsed}
 dump('VALIDATION.json',result);print('COLD_AUDIT',result['PASS'],'FAIL',len(result['FAIL']),'QDIFF',max_q_difference,'UNRESOLVED',len(unresolved),flush=True)
 if not result['PASS']:raise SystemExit(1)

if __name__=='__main__':audit()