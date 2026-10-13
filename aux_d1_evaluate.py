"""Post-seal evaluation only: truth cannot affect inference selection."""
import math,time,json
import numpy as np
from aux_d1_common import *
from aux_d1_model import LOW,HIGH,SIGMA,wrap,score
PARAM=['r','theta','v','psi'];REFERENCE=np.array([.1,1.,.1,5.])
def errors(h,true):
 if h is None:return np.full(4,np.inf)
 h=np.asarray(h);true=np.asarray(true);d=np.abs(h-true);d[[1,3]]=np.abs(np.rad2deg(wrap(np.deg2rad(h[[1,3]]-true[[1,3]]))));d[0]/=true[0];d[2]/=true[2];return d

def quantile(v,p):
 v=sorted(float(x) for x in v);return v[math.ceil(p*len(v))-1] if v else None

def summarize(records):
 rows=[]
 for variant in ['PRIMARY','SHARED']:
  for R in ['U0','U5']:
   for g in range(12):
    rr=[r for r in records if r['resource']==R and r['geometry_index']==g];evaluated=[r for r in rr if r['evaluated']];acc=[r for r in evaluated if r['accepted' if variant=='PRIMARY' else 'shared_accepted']]
    key='error' if variant=='PRIMARY' else 'shared_error';row={'variant':variant,'resource':R,'geometry_index':g,'case_id':rr[0]['case_id'],'planned':128,'evaluated':len(evaluated),'accepted':len(acc),'unevaluated':128-len(evaluated),'output_rate':len(acc)/len(evaluated) if evaluated else None,'residual_abstentions':sum(r['reason']=='ABSTAIN_RESIDUAL' for r in evaluated) if variant=='PRIMARY' else None,'numerical_failures':sum(r['reason']=='NUMERICAL_OR_OPTIMIZER_FAILURE' for r in evaluated) if variant=='PRIMARY' else None,'large_any_count':sum(bool(np.any(r[key]>REFERENCE)) for r in acc),'better_terminal_not_selected':sum(r['better_terminal_not_selected'] for r in evaluated) if variant=='PRIMARY' else None}
    for j,p in enumerate(PARAM):
     for condition,source in [('unconditional',evaluated),('conditional',acc)]:
      vals=[r[key][j] for r in source]
      for suffix,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99),('max',1.)]:row[p+'_'+condition+'_'+suffix]=quantile(vals,q)
     row[p+'_large_count']=sum(r[key][j]>REFERENCE[j] for r in acc)
    rows.append(row)
 return rows

def classify(case,full):
 signals={p:[] for p in ['r','v','psi']};qualified=[];details=[]
 for g in range(12):
  a=next(r for r in case if r['resource']=='U0' and r['geometry_index']==g);b=next(r for r in case if r['resource']=='U5' and r['geometry_index']==g);valid=b['evaluated']==128 and a['evaluated']==128 and b['output_rate']>=.9 and b['output_rate']>=a['output_rate']-.02
  flags={};improvements={}
  for p in ['r','v','psi']:
   x=a[p+'_unconditional_median'];y=b[p+'_unconditional_median'];delta=1-y/x if x is not None and y is not None and math.isfinite(x) and x>0 and math.isfinite(y) else None
   flags[p]=bool(valid and delta is not None and delta>=.3);improvements[p]=delta
   if flags[p]:signals[p].append(g)
  if any(flags.values()):qualified.append(g)
  details.append({'geometry_index':g,'signal_by_parameter':flags,'improvements':improvements,'output_condition':valid})
 primary='AUX_D1_PAIRED_HORIZONTAL_GAIN_SUPPORTED_CONDITIONALLY' if len(qualified)>=8 else 'AUX_D1_INCREMENT_LIMITED_UNDER_REGISTERED_PROCESSOR'
 return {'primary_classification':primary if full else 'NOT_EVALUATED_INCOMPLETE_PANEL','qualifying_geometry_indices':qualified,'signal_geometries_by_parameter':signals,'by_geometry':details}

def json_safe(x):
 if isinstance(x,float) and not math.isfinite(x):return 'INF' if x>0 else '-INF'
 if isinstance(x,np.generic):return json_safe(x.item())
 if isinstance(x,np.ndarray):return json_safe(x.tolist())
 if isinstance(x,dict):return {k:json_safe(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [json_safe(v) for v in x]
 return x

def run():
 verify();seal=read('INFERENCE_SEAL.json')
 for b in seal['bindings']:
  if sha(ROOT/b['path'])!=b['sha256']:raise RuntimeError('sealed inference changed')
 if (OUT/'EVALUATION_START.json').exists():raise RuntimeError('Evaluation already started; report-only recovery separate')
 dump('EVALUATION_START.json',{'utc':utc(),'inference_seal_sha256':sha(OUT/'INFERENCE_SEAL.json'),'truth_read_after_inference':True});tic=time.perf_counter()
 pub=np.load(OUT/'public_inputs/test.npz');priv=np.load(OUT/'evaluation_only/test_private.npz');plan=csvread('PAIRED_RECORD_INDEX.csv');a=np.load(OUT/'ALL_TERMINALS.npy',mmap_mode='r');out={r['configuration']:r for r in (json.loads(l) for l in (OUT/'CONFIGURATION_OUTPUTS.jsonl').read_text().splitlines())};threshold=read('CALIBRATION_THRESHOLDS.json')['thresholds']
 records=[];estimates=[];truthrows=[];envrows=[];sharedrows=[]
 for i in range(1536):
  hh=priv['h_truth'][i];target=1000*hh[0]*np.array([np.cos(np.deg2rad(hh[1])),np.sin(np.deg2rad(hh[1]))])+pub['times'][:,None]*hh[2]*np.array([np.cos(np.deg2rad(hh[3])),np.sin(np.deg2rad(hh[3]))]);P=np.stack([pub['reported_main'][i],pub['reported_aux'][i]],axis=1);y=np.concatenate([pub['bearing_main'][i],pub['bearing_aux'][i]],axis=1);delta=target[:,None,:]-P;e=wrap(np.arctan2(delta[:,:,1],delta[:,:,0])-y)/SIGMA
  p0=target[0]-P[0,0];reported_h=np.array([np.hypot(*p0)/1000,np.rad2deg(np.arctan2(p0[1],p0[0])),hh[2],hh[3]]);inside=bool(np.all((reported_h>=LOW)&(reported_h<=HIGH)))
  for ri,R in enumerate(['U0','U5']):
   c=2*i+ri;o=out.get(c);evaluated=o is not None;accepted=evaluated and o['selected_start']>=0;sharedacc=evaluated and o['shared_selected_start']>=0;err=errors(o['selected_h'] if accepted else None,hh);shared_err=errors(o['shared_h'] if sharedacc else None,hh)
   risk=False
   if accepted and np.any(err>REFERENCE):risk=any(np.all(errors(a[c,int(k),:4],hh)<=REFERENCE) for k in o['accepted_indices'])
   record={'configuration':c,'record_index':i,'resource':R,'geometry_index':int(plan[i]['geometry_index']),'case_id':plan[i]['case_id'],'replicate':int(plan[i]['replicate']),'side':int(plan[i]['side']),'evaluated':evaluated,'accepted':accepted,'shared_accepted':sharedacc,'error':err,'shared_error':shared_err,'reason':o['reason'] if evaluated else 'NOT_EVALUATED','better_terminal_not_selected':risk};records.append(record)
   qr=score(e[:,:ri+1].ravel());truecompatible=inside and qr<=threshold[R];truthrows.append({'configuration':c,'resource':R,'truth_report_reference_in_domain':inside,'truth_Q':format_float(qr),'truth_Q_pass':qr<=threshold[R],'TRUE_COMPATIBLE_BUT_NO_OUTPUT':(truecompatible and not accepted) if evaluated else 'NOT_EVALUATED','evaluation_only':True})
   row={k:record[k] for k in ['configuration','record_index','resource','geometry_index','case_id','replicate','side','evaluated','accepted','reason','better_terminal_not_selected']};row['selected_start']=o['selected_start'] if evaluated else 'NA';row['selected_Q']=format_float(o['selected_Q']) if accepted else 'NA'
   for j,p in enumerate(PARAM):row[p+'_hat']=format_float(o['selected_h'][j]) if accepted else 'NA';row[p+'_error']=format_float(err[j]) if evaluated else 'NOT_EVALUATED'
   estimates.append(row)
   env=o['envelope'] if evaluated else None;envrows.append({'configuration':c,'resource':R,'accepted_count':env['accepted_count'] if env else 'NOT_EVALUATED','cluster_count':env['cluster_count'] if env else 'NOT_EVALUATED','minimum_json':json.dumps(env['minimum']) if env else 'NA','maximum_json':json.dumps(env['maximum']) if env else 'NA','representative_start_ids_json':json.dumps(env['representative_start_ids']) if env else 'NA','scope':'FINITE_TERMINAL_ENVELOPE_ONLY'})
   sr={'configuration':c,'resource':R,'evaluated':evaluated,'main_accepted':accepted,'shared_accepted':sharedacc,'main_start':o['selected_start'] if evaluated else 'NA','shared_start':o['shared_selected_start'] if evaluated else 'NA','changed_selected_start':evaluated and o['selected_start']!=o['shared_selected_start']}
   for j,p in enumerate(PARAM):sr[p+'_main_error']=format_float(err[j]) if evaluated else 'NOT_EVALUATED';sr[p+'_shared_error']=format_float(shared_err[j]) if evaluated else 'NOT_EVALUATED'
   sharedrows.append(sr)
 case=summarize(records);primary=[r for r in case if r['variant']=='PRIMARY'];shared=[r for r in case if r['variant']=='SHARED'];full=seal['status']=='COMPLETE';decision=classify(primary,full);shareddec=classify(shared,full);classes=[decision['primary_classification']] if full else ['PARTIAL_EXECUTION'];speedtail=[r['geometry_index'] for r in primary if r['resource']=='U5' and r['v_unconditional_P95'] is not None and r['v_unconditional_P95']>.1]
 if full and (len(decision['signal_geometries_by_parameter']['r'])>=8 or len(decision['signal_geometries_by_parameter']['psi'])>=8) and speedtail:classes.append('AUX_D1_POSITION_OR_HEADING_GAIN_WITH_SPEED_TAIL_LIMIT')
 pairrows=[]
 for g in [*range(12),'ALL']:
  pairs=[(records[2*i],records[2*i+1]) for i in range(1536) if g=='ALL' or int(plan[i]['geometry_index'])==g];both=[(x,y) for x,y in pairs if x['evaluated'] and y['evaluated']];row={'geometry':g,'planned_pairs':len(pairs),'evaluated_pairs':len(both)}
  for ka,kb,name in [(True,True,'both_accepted'),(True,False,'U0_only'),(False,True,'U5_only'),(False,False,'neither')]:row[name]=sum(x['accepted']==ka and y['accepted']==kb for x,y in both)
  for j,p in enumerate(PARAM):
   for relation in ['improved','equal','worse']:row[p+'_'+relation]=sum((y['error'][j]<x['error'][j] if relation=='improved' else y['error'][j]>x['error'][j] if relation=='worse' else y['error'][j]==x['error'][j]) for x,y in both)
   dd=[y['error'][j]-x['error'][j] for x,y in both if x['accepted'] and y['accepted']];row[p+'_both_accepted_delta_median']=quantile(dd,.5)
  for R,j in [('U0',0),('U5',1)]:row[R+'_large_any_count']=sum(r[j]['accepted'] and np.any(r[j]['error']>REFERENCE) for r in both);row[R+'_output_count']=sum(r[j]['accepted'] for r in both)
  pairrows.append(row)
 allpair=pairrows[-1]
 if allpair['U5_large_any_count']>allpair['U0_large_any_count'] or allpair['U5_output_count']<allpair['U0_output_count'] or any(summary[1][p+'_large_output_count']>summary[0][p+'_large_output_count'] for p in PARAM):classes.append('PAIRED_TRADEOFF_PRESENT')
 if any(r['TRUE_COMPATIBLE_BUT_NO_OUTPUT'] is True for r in truthrows):classes.append('FINITE_SEARCH_DELIVERY_LIMIT')
 if decision['primary_classification']!=shareddec['primary_classification']:classes.append('INITIALIZATION_DEPENDENCE_OBSERVED')
 siderows=[]
 for side in [-1,1]:
  for R in ['U0','U5']:
   for g in range(12):
    rr=[r for r in records if r['side']==side and r['resource']==R and r['geometry_index']==g];executed=[r for r in rr if r['evaluated']];ac=[r for r in executed if r['accepted']];row={'side':side,'resource':R,'geometry_index':g,'planned':64,'evaluated':len(executed),'accepted':len(ac),'unevaluated':64-len(executed)}
    for j,p in enumerate(PARAM):
     for cond,group in [('unconditional',executed),('conditional',ac)]:
      for label,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99),('max',1.)]:row[p+'_'+cond+'_'+label]=quantile([r['error'][j] for r in group],q)
    siderows.append(row)
 summary=[]
 for R in ['U0','U5']:
  rr=[r for r in primary if r['resource']==R];recordR=[r for r in records if r['resource']==R and r['evaluated']];row={'resource':R,'evaluated':len(recordR),'accepted':sum(r['accepted'] for r in recordR),'planned':1536}
  for p in PARAM:
   for cond in ['unconditional','conditional']:
    vals=[r[p+'_'+cond+'_median'] for r in rr];row[p+'_'+cond+'_mean_geometry_median']=sum(vals)/12 if all(v is not None for v in vals) else None
    vals=[r[p+'_'+cond+'_P95'] for r in rr];row[p+'_'+cond+'_worst_geometry_P95']=max(vals) if all(v is not None for v in vals) else None
   j=PARAM.index(p)
   for label,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99),('max',1.)]:row[p+'_pooled_unconditional_'+label]=quantile([r['error'][j] for r in recordR],q)
   row[p+'_large_output_count']=sum(r['accepted'] and r['error'][j]>REFERENCE[j] for r in recordR)
  row['large_any_output_count']=sum(r['accepted'] and np.any(r['error']>REFERENCE) for r in recordR);row['better_terminal_not_selected']=sum(r['better_terminal_not_selected'] for r in recordR);summary.append(row)
 for name,data in [('ESTIMATES.csv',estimates),('PAIRED_RESULTS_BY_CASE.csv',case),('PAIRED_RESULTS_BY_SIDE.csv',siderows),('PAIRED_OUTCOME_COUNTS.csv',pairrows),('FINITE_TERMINAL_ENVELOPES.csv',envrows),('SHARED_START_SUBSET_DIAGNOSTIC.csv',sharedrows),('TRUTH_COMPATIBILITY_EVAL_ONLY.csv',truthrows),('SUMMARY_BY_RESOURCE.csv',summary)]:csvwrite(name,data)
 dump('DECISION.json',json_safe({'stage':'AUX-D1','task_SHA256':TASKSHA,'parent_SHA':PARENT,'design_SHA':read('DESIGN_COMMIT_VERIFICATION.json')['design_SHA'],'execution_SHA':'CONTAINING_COMMIT','classifications':classes,'primary':decision,'shared_start_subset':shareddec,'U5_speed_P95_above10percent_geometries':speedtail,'summary':summary,'execution_status':seal['status'],'new_paired_test_records':1536,'new_calibration_paired_records':512,'executed_optimization_terminals':seal['terminals'],'no_radial_depth_acoustics':True,'new_propagation_calls':0,'signed_bearing':'ABSOLUTE_BEARING_SIGN_CONDITION_ASSUMED','confidence_scope':'FINITE_TERMINAL_ENVELOPE_ONLY','R4_percent':0,'next_stage_execution_authorized':False,'STOP':True}))
 record_clock('evaluation',time.perf_counter()-tic,optimizer_runs=0)
 dump('EVALUATION_SEAL.json',{'utc':utc(),'inference_seal_sha256':sha(OUT/'INFERENCE_SEAL.json'),'decision':filemeta(OUT/'DECISION.json'),'frozen_selection_unchanged':True})
 print('EVALUATION',classes,'SIGNALS',decision['signal_geometries_by_parameter'],flush=True)
if __name__=='__main__':run()