"""Read-only numerical evaluation of sealed N2 full fields and all planned sentinel outcomes."""
from t2_n2_common import *

def main():
 seal=read(OUT/'SCIENCE_OUTPUT_FREEZE.json');g=geometries();ss=sstates();sns=sentinels();kready=seal['K_complete_variants'];internal=[];interp=[];cross=[];raw=[];cache={}
 for v in FVAR:
  for step in ([.05,.025] if v=='F2' else [.05]):
   for res in RES:cache[v,step,res]=np.load(LOCAL/'FFP_FIELDS'/f'{v}_{res}_step{step:.3f}.npz')
 def compare_pair(name,left,right,stepL=.05,stepR=.05,threshold=.002,dest=None,crossK=False):
  arrays={};valid={}
  for res in RES:
   a=cache[left,stepL,res]['field'];ma=cache[left,stepL,res]['completed'];b=cache[right,stepR,res]['field'] if not crossK else np.full_like(a,np.nan+1j*np.nan);mb=cache[right,stepR,res]['completed'] if not crossK else np.zeros_like(ma)
   if crossK and 'K3' in kready:
    full=np.load(LOCAL/'K_FIELDS'/f'K3_{res}_fine.npy',mmap_mode='r');cat=states('fine')
    for i,s in enumerate(ss):cid=int(np.flatnonzero(np.all(cat==s,axis=1))[0]);b[i]=full[cid];mb[i]=True
   arrays[res]=(a,b);valid[res]=ma&mb
  for fi,f in enumerate(FREQ):
   for z in [180,190,200,210,220]:
    indices=np.flatnonzero(ss[:,2]==z);complete=all(np.all(valid[res][indices,fi]) for res in RES)
    if complete:
     av=np.concatenate([arrays[res][0][indices,fi].ravel() for res in RES]);bv=np.concatenate([arrays[res][1][indices,fi].ravel() for res in RES]);c=np.vdot(av,bv)/np.vdot(av,av);e=float(np.linalg.norm(bv-c*av)/np.linalg.norm(bv));rr=float(np.linalg.norm(bv-av)/np.linalg.norm(bv));raw.append({'comparison':name,'frequency_hz':f,'source_z_m':z,'sentinels':len(indices),'resources':4,'entries':len(av),'global_factor_real':float(c.real),'global_factor_imag':float(c.imag),'raw_global_profiled_relative':e,'raw_unaligned_relative':rr,'status':'COMPLETED','factor_scope':'one_per_frequency_source_depth_all_sentinels_and4resources'})
    else:raw.append({'comparison':name,'frequency_hz':f,'source_z_m':z,'status':'NOT_EVALUATED','raw_global_profiled_relative':'NOT_EVALUATED','reason':'one_or_more_registered_inputs_missing'})
  for res in RES:
   a,b=arrays[res]
   for i,s in enumerate(sns):
    for fi,f in enumerate(FREQ):
     prefix={'comparison':name,'resource':res,'sentinel_id':s['sentinel_id'],'frequency_hz':f,'r_km':float(s['r_km']),'abs_theta_deg':float(s['abs_theta_deg']),'source_z_m':float(s['source_z_m']),'shape_limit':threshold}
     if valid[res][i,fi]:
      e=float(shape(a[i,fi][None],b[i,fi][None])[0]);rr=float(np.linalg.norm(b[i,fi]-a[i,fi])/np.linalg.norm(b[i,fi]));row=prefix|{'shape_error':e,'raw_unaligned_relative':rr,'norm_left':float(np.linalg.norm(a[i,fi])),'norm_right':float(np.linalg.norm(b[i,fi])),'PASS':e<=threshold,'status':'COMPLETED'}
     else:row=prefix|{'shape_error':'NOT_EVALUATED','raw_unaligned_relative':'NOT_EVALUATED','norm_left':'NOT_EVALUATED','norm_right':'NOT_EVALUATED','PASS':'NOT_EVALUATED','status':'NOT_EVALUATED'}
     dest.append(row)
 compare_pair('F1_VS_F2','F1','F2',dest=internal);compare_pair('F2_VS_F3','F2','F3',dest=internal);compare_pair('F2_STEP050_VS025','F2','F2',stepR=.025,threshold=.0005,dest=interp);compare_pair('F2_VS_K3','F2','F2',dest=cross,crossK=True)
 table(OUT/'FFP_INTERNAL_STABILITY.csv',internal);table(OUT/'FFP_RANGE_INTERPOLATION.csv',interp);table(OUT/'KRAKEN_FFP_CROSS_METHOD.csv',cross);table(OUT/'RAW_GLOBAL_COMPLEX_FACTOR_DIAGNOSTIC.csv',raw)
 full=rows(OUT/'FULL_CATALOGUE_KRAKEN_STABILITY.csv')
 def gate(items,expected,limit,yes,no,field='shape_error'):
  complete=[r for r in items if r.get('status','COMPLETED')=='COMPLETED'];values=[float(r[field]) for r in complete];n=len(complete);passed=n==expected and all(v<=limit for v in values)
  return {'status':yes if passed else no if n==expected else 'NOT_EVALUATED_PARTIAL_EXECUTION','PASS':passed,'planned':expected,'completed':n,'max_difference':max(values) if values else 'NOT_EVALUATED','threshold':limit,'failures':sum(v>limit for v in values)}
 win=gate([r for r in full if r['comparison'] in ['1800_TO_2000','2000_TO_2400']],48,.002,'WINDOW_BRACKET_STABLE_1800_TO_2400','WINDOW_CHANGE_MATERIAL','max');mesh=gate([r for r in full if r['comparison']=='MESH160001_TO240001'],24,.002,'K_MESH_STABLE_WITHIN_160001_240001','K_MESH_DIFFERENCE_MATERIAL','max');fm=gate(internal,288,.002,'FFP_FINITE_BRACKET_STABLE','FFP_NUMERIC_INCOMPLETE');ip=gate(interp,144,.0005,'FFP_INTERPOLATION_STABLE','FFP_NUMERIC_INCOMPLETE');ffpself={'status':'FFP_FINITE_BRACKET_STABLE' if fm['PASS'] and ip['PASS'] else 'FFP_NUMERIC_INCOMPLETE' if fm['completed'] or ip['completed'] else 'INDEPENDENT_SOLVER_NOT_READY','PASS':bool(fm['PASS'] and ip['PASS']),'mesh_and_spectral':fm,'interpolation':ip};cr=gate(cross,144,.002,'KRAKEN_FFP_SENTINEL_AGREEMENT','DIFFERENT_NUMERICAL_PREDICTIONS');cr['diagnostic_predictions_saved_even_if_not_admitted']=True
 if not ffpself['PASS']:cr['status']='CROSS_METHOD_NOT_EVALUATED';cr['PASS']=False;cr['reason']='FFP_SELF_BRACKET_NOT_ESTABLISHED; raw_differences_not_a_truth_error_bound'
 allpass=win['PASS'] and mesh['PASS'] and ffpself['PASS'] and cr['PASS'];calls=rows(OUT/'PROVIDER_CALLS.csv');decision={'stage':'T2-N2','overall':'T2_N2_SELECTED_MODEL_NUMERIC_BRACKETING_SUPPORTED' if allpass else 'T2_N2_NUMERIC_BRACKETING_NOT_ESTABLISHED','K_WINDOW':win,'K_MESH':mesh,'FFP_SELF':ffpself,'CROSS_METHOD':cr,'new_calls':{tool:sum(r['tool']==tool for r in calls) for tool in ['kraken','scooter','fields']},'max_call_budget':{'kraken':9,'scooter':9,'fields':144},'all_K_plans_complete':len(kready)==4,'all_FFP_plans_complete':seal['FFP_runtime_completed'],'historical_G0':'FORWARD_NUMERIC_NOT_ADMITTED_UNCHANGED','N1_template_advantage':'UNCHANGED_NOT_REEVALUATED','A_performance':'NOT_EVALUATED','B_performance':'NOT_EVALUATED','CHIGH_global':'UNCERTIFIED_ABOVE_TESTED_WINDOWS_AND_OUTSIDE_SENTINELS','SSP_locality_and_source_derivative':'UNCHANGED_NOT_EVALUATED','R4_percent':0,'physical_truth_certification':False,'research_lead_audit':'PENDING','state':'STOP','next_stage_authorized':False};dump(OUT/'NUMERICAL_DECISION.json',decision)
 dump(OUT/'NUMERICAL_EVALUATION_FREEZE.json',{'utc':utc(),'scientific_result_sha256':{p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name in ['FFP_INTERNAL_STABILITY.csv','FFP_RANGE_INTERPOLATION.csv','KRAKEN_FFP_CROSS_METHOD.csv','RAW_GLOBAL_COMPLEX_FACTOR_DIAGNOSTIC.csv','NUMERICAL_DECISION.json']},'main_forward_sources_unchanged':True});print('N2 NUMERICAL DECISION',decision,flush=True)
if __name__=='__main__':main()