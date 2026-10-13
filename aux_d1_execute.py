"""Once-only single-process AUX-D1 inference; never opens generator truth archives."""
import time,json,os
import numpy as np
from aux_d1_common import *
from aux_d1_model import *
class BudgetStop(RuntimeError):pass

def run():
 design=verify();spec=read('AUX_D1_TASK_SPEC_20261012_v1.json')
 if (OUT/'INFERENCE_START.json').exists():raise RuntimeError('Already started: no second scientific run authorized')
 seal=read('TEST_GENERATION_SEAL.json')
 if sha(ROOT/seal['public_test']['path'])!=seal['public_test']['sha256']:raise RuntimeError('public inputidentity')
 threshold=read('CALIBRATION_THRESHOLDS.json')['thresholds'];tic=time.perf_counter();prior=elapsed_completed();dump('INFERENCE_START.json',{'utc':utc(),'design_SHA':read('DESIGN_COMMIT_VERIFICATION.json')['design_SHA'],'public_test_sha256':seal['public_test']['sha256'],'truth_inputs':False,'resource_views':3072})
 data=np.load(ROOT/seal['public_test']['path'],allow_pickle=False);t=data['times'];fixed=np.load(OUT/'SHARED_FIXED_126.npy',allow_pickle=False)
 terminals=np.lib.format.open_memmap(OUT/'ALL_TERMINALS.npy',mode='w+',dtype=np.float64,shape=(3072,128,len(FIELDS)));terminals[:]=np.nan
 used=np.lib.format.open_memmap(OUT/'STARTS_USED.npy',mode='w+',dtype=np.float64,shape=(3072,128,4));used[:]=np.nan
 initlog=[];completed=[];stopped='';peak=0;actual_nfev=0;terminal_count=0;next_config=0
 def guard():
  el=time.perf_counter()-tic
  if el>=18000 or prior+el>=21600:raise BudgetStop('TIME_HARD_CAP')
 with (OUT/'CONFIGURATION_OUTPUTS.jsonl').open('w',encoding='utf-8',newline='\n') as outputs,(OUT/'EXCEPTIONS.jsonl').open('w',encoding='utf-8',newline='\n') as exc:
  try:
   for i in range(1536):
    guard();peak=max(peak,rss())
    if peak>8*1024**3:raise BudgetStop('RSS_HARD_CAP')
    if storage_bytes()>2*1024**3:raise BudgetStop('STORAGE_HARD_CAP')
    Pm=np.ascontiguousarray(data['reported_main'][i][:,None,:]);ym=np.ascontiguousarray(data['bearing_main'][i]);first=None
    for ri,R in enumerate(['U0','U5']):
     config=2*i+ri;next_config=config
     P=Pm.copy() if ri==0 else np.ascontiguousarray(np.concatenate([Pm,data['reported_aux'][i][:,None,:]],axis=1));y=ym.copy() if ri==0 else np.ascontiguousarray(np.concatenate([ym,data['bearing_aux'][i]],axis=1))
     assert P[:,0].tobytes()==Pm[:,0].tobytes() and y[:,:1].tobytes()==ym.tobytes()
     starts_u,diag=starts(t,P,y,fixed);used[config]=starts_u
     if first is None:first=starts_u.copy()
     else:assert starts_u[[0,*range(2,128)]].tobytes()==first[[0,*range(2,128)]].tobytes()
     initlog.append({'configuration':config,'record_index':i,'resource':R,'diagnostics_json':json.dumps(diag,ensure_ascii=False,separators=(',',':'))})
     for start_id,u in enumerate(starts_u):
      guard()
      try:result=fit_one(t,P,y,u,guard)
      except Exception as err:
       exc.write(json.dumps({'configuration':config,'start_id':start_id,'utc':utc(),'error':repr(err)})+'\n');exc.flush();raise
      if result[8]>100 or result[10]>100:raise RuntimeError('Residual budget exceeded')
      terminals[config,start_id]=result;terminal_count+=1;actual_nfev+=int(result[8])
     a=terminals[config];sel,ids,reason=choose(a,threshold[R]);shared,_,sharedreason=choose(a,threshold[R],True);env=envelope(a,ids)
     record={'configuration':config,'record_index':i,'resource':R,'selected_start':sel,'selected_h':a[sel,:4].tolist() if sel>=0 else None,'selected_Q':float(a[sel,4]) if sel>=0 else None,'reason':reason,'shared_selected_start':shared,'shared_h':a[shared,:4].tolist() if shared>=0 else None,'shared_reason':sharedreason,'envelope':env,'accepted_indices':ids.tolist(),'solver_success_count':int(np.sum(a[:,6]==1)),'nonconverged_lower_Q':bool(sel>=0 and np.any((a[:,6]!=1)&(a[:,4]<a[sel,4])))}
     outputs.write(json.dumps(record,ensure_ascii=False,allow_nan=False,separators=(',',':'))+'\n');outputs.flush();completed.append(config);next_config=config+1
    terminals.flush();used.flush()
    dump('EXECUTION_CHECKPOINT.json',{'utc':utc(),'complete_pairs':i+1,'complete_configurations':len(completed),'terminals':terminal_count,'main_elapsed_s':time.perf_counter()-tic,'actual_nfev':actual_nfev,'peak_RSS_bytes':peak,'next_configuration':next_config})
    if (i+1)%12==0:print('PROGRESS',i+1,'/1536 pairs;',terminal_count,'terminals;',round(time.perf_counter()-tic,2),'s',flush=True)
  except Exception as err:stopped=repr(err);print('STOPPED',stopped,flush=True)
 terminals.flush();used.flush();csvwrite('INITIALIZER_AND_SEARCH_DIAGNOSTICS.csv',initlog or [{'configuration':'NOT_EVALUATED','record_index':'NA','resource':'NA','diagnostics_json':'{}'}])
 duration=time.perf_counter()-tic;record_clock('main_optimization',duration,complete_configurations=len(completed),terminals=terminal_count,nfev=actual_nfev)
 dump('TERMINAL_SCHEMA.json',{'file':'ALL_TERMINALS.npy','dtype':'float64','shape':[3072,128,len(FIELDS)],'fields':FIELDS,'unexecuted':'allNaN,never treated as zero or failure','score_order':'epochmajor,nodewithin epoch;np.sumfloat64stableterms','start_id_axis':True})
 bindings=[filemeta(OUT/n) for n in ['ALL_TERMINALS.npy','STARTS_USED.npy','CONFIGURATION_OUTPUTS.jsonl','INITIALIZER_AND_SEARCH_DIAGNOSTICS.csv','TERMINAL_SCHEMA.json']]
 dump('INFERENCE_SEAL.json',{'utc':utc(),'complete_configurations':len(completed),'complete_pairs':len(completed)//2,'terminals':terminal_count,'expected_terminals':393216,'status':'COMPLETE' if len(completed)==3072 else 'PARTIAL_EXECUTION','stop_reason':stopped,'unexecuted_configurations':[i for i in range(3072) if i not in completed],'bindings':bindings,'main_elapsed_s':duration,'peak_RSS_bytes':peak,'actual_nfev':actual_nfev,'new_propagation':0,'truth_read':False})
 print('INFERENCE_SEALED',len(completed),terminal_count,duration,flush=True)
if __name__=='__main__':run()
