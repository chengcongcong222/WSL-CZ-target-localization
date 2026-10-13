"""Post-seal truth evaluation and preregistered resource comparisons; no reselection."""
from t2_tilt_g0_core import *

def q(x,p):
 a=np.sort(np.asarray(x,float));return float(np.median(a)) if p==.5 else float(a[max(0,int(np.ceil(p*len(a)))-1)])
def finite(v):return 'INF' if not np.isfinite(v) else float(v)
def main():
 assert sha(OUT/'INFERENCE_RESULTS_NO_TRUTH.csv')==read(OUT/'INFERENCE_FREEZE.json')['result_SHA256']
 inf=rows(OUT/'INFERENCE_RESULTS_NO_TRUTH.csv');pub={int(r['record_id']):r for r in rows(OUT/'TEST_DFT_INDEX.csv')};private=read(LOCAL/'test/private/MANIFEST.json');truth={int(r['base_id']):r['truth'] for r in private['private_records']};g=geometries();providers={m:Provider(m) for m in [80001,160001]};ev=[];diagnostic=[];cache={};scores_cache={};bits_cache={};row_lookup={}
 for t in inf:
  res=t['resource'];cat=t['catalogue'];mesh=int(t['mesh']);rid=int(t['record_id']);key=(res,cat,mesh)
  if key not in scores_cache:
   scores_cache[key]=np.load(LOCAL/'inference'/f'SCORES_{res}_{cat}_n{mesh}.npy',mmap_mode='r');bits_cache[key]=np.load(LOCAL/'inference'/f'ACCEPTED_{res}_{cat}_n{mesh}.npy',mmap_mode='r');ids=np.load(LOCAL/'inference'/f'RECORD_IDS_{res}_{cat}_n{mesh}.npy');row_lookup[key]={int(v):j for j,v in enumerate(ids)}
  i=row_lookup[key][rid];s=catalogue(cat=='fine');accepted=np.unpackbits(bits_cache[key][i],bitorder='little')[:len(s)].astype(bool);b=truth[int(t['base_id'])];axes=[np.unique(s[:,j]) for j in range(3)];near=[a[np.argmin(abs(a-v))] for a,v in zip(axes,[b['r_km'],abs(b['theta_deg']),b['z_m']])];hret=bool(np.any(accepted&(s[:,0]==near[0])&(s[:,1]==near[1])));jret=bool(np.any(accepted&np.all(s==near,axis=1)));hc=(res,mesh,b['r_km'],b['theta_deg'],b['z_m'])
  if hc not in cache:cache[hc]=providers[mesh].field([[b['r_km'],b['theta_deg'],b['z_m']]],g[res])[0]
  hh=cache[hc];yy=np.load(pub[rid]['DFT_path']);trueD=float(1-score(hh[None],yy)[0]);valid=t['accepted_output']=='True';idsr=json.loads(t['range_labels']);cell=json.loads(t['range_cell_intervals']);pt=abs(float(t['best_r_km'])-b['r_km']);theta=abs(float(t['best_abs_theta_deg'])-abs(b['theta_deg']));kept=s[accepted];worst=float(np.max(abs(kept[:,0]-b['r_km']))) if valid else float('inf');snr=int(t['reference_SNR_dB']);nu=read(OUT/'NOISE_REFERENCE_FREEZE.json')['nu'][str(snr)];actual=10*np.log10(abs(hh)**2/nu)
  ev.append(dict(t,case_id=b['case_id'],generated_z_m=b['z_m'],replicate=b['replicate'],truth_r_km=b['r_km'],truth_theta_deg=b['theta_deg'],horizontal_label_retained=hret,joint_label_retained=jret,continuous_truth_D=trueD,continuous_truth_passes=trueD<=float(t['tau_D']),point_range_error_km=pt,unconditional_range_error_km=pt if valid else 'INF',point_abs_theta_error_deg=theta,full_set_max_label_error_km=finite(worst),min_actual_unit_source_SNR_dB=float(actual.min()),median_actual_unit_source_SNR_dB=float(np.median(actual)),remaining_depth_labels=json.dumps(np.unique(kept[:,2]).tolist()),left_right_unresolved=True))
 table(OUT/'ACCEPTANCE_AND_POINT_BY_CASE.csv',ev)
 # All9 geometry/depth cases, resources, C/SNR and both numeric checks remain separate.
 group={}
 for t in ev:
  key=tuple(t[k] for k in ['case_id','generated_z_m','resource','calibration','reference_SNR_dB','catalogue','mesh']);group.setdefault(key,[]).append(t)
 summary=[]
 for key,a in group.items():
  assert len(a)==16;valid=[t for t in a if t['accepted_output']=='True'];width=[float(t['envelope_km']) for t in a];error=[float(t['unconditional_range_error_km']) for t in a];row=dict(zip(['case_id','generated_z_m','resource','calibration','reference_SNR_dB','catalogue','mesh'],key));row.update({'n':16,'accepted':len(valid),'horizontal_retained':sum(t['horizontal_label_retained'] for t in a),'joint_retained':sum(t['joint_label_retained'] for t in a),'continuous_truth_passes':sum(t['continuous_truth_passes'] for t in a),'median_unconditional_envelope_km':finite(q(width,.5)),'P90_unconditional_envelope_km':finite(q(width,.9)),'P95_unconditional_envelope_km':finite(q(width,.95)),'median_conditional_envelope_km':finite(q([float(t['envelope_km']) for t in valid],.5)) if valid else 'NOT_EVALUATED','median_unconditional_point_error_km':finite(q(error,.5)),'P95_unconditional_point_error_km':finite(q(error,.95)),'median_conditional_point_error_km':finite(q([float(t['point_range_error_km']) for t in valid],.5)) if valid else 'NOT_EVALUATED','median_intervals':q([t['disconnected_intervals'] for t in a],.5),'min_actual_unit_source_SNR_dB':min(t['min_actual_unit_source_SNR_dB'] for t in a)})
  summary.append(row)
 table(OUT/'COMPARISON_BY_RESOURCE.csv',summary);table(OUT/'TILT_SIGN_AND_ANGLE_SENSITIVITY.csv',[r for r in summary if r['resource'].startswith('A_')]);table(OUT/'WEAK_FIELD_AND_NUMERICAL_DIAGNOSTICS.csv',[{k:t[k] for k in ['record_id','resource','calibration','reference_SNR_dB','case_id','generated_z_m','catalogue','mesh','min_actual_unit_source_SNR_dB','median_actual_unit_source_SNR_dB','continuous_truth_D','tau_D']} for t in ev])
 signals=[]
 for cat in ['primary','fine']:
  for mesh in [80001,160001]:
   for label,treat,ctrl in [('A','A_P4M','A0_0'),('B','B1_V4','B2_DUP')]:
    for case in ['H01','H06','H12']:
     def find(res):return next(r for r in summary if r['case_id']==case and r['generated_z_m']==200 and r['resource']==res and r['calibration']=='C0' and int(r['reference_SNR_dB'])==20 and r['catalogue']==cat and int(r['mesh'])==mesh)
     t=find(treat);c=find(ctrl);wt=float(t['median_unconditional_envelope_km']);wc=float(c['median_unconditional_envelope_km']);red=1-wt/wc if np.isfinite(wc) and wc>0 else float('-inf');yes=t['accepted']>=15 and t['horizontal_retained']>=15 and wt<=10 and red>=.3;signals.append({'comparison':label,'case_id':case,'catalogue':cat,'mesh':mesh,'treated':treat,'control':ctrl,'treated_accepted':t['accepted'],'treated_horizontal_retained':t['horizontal_retained'],'treated_median_envelope_km':finite(wt),'control_median_envelope_km':finite(wc),'envelope_reduction_fraction':finite(red),'signal':bool(yes)})
 table(OUT/'PRIMARY_SIGNAL_BY_GEOMETRY.csv',signals)
 counts={f'{l}_{cat}_{mesh}':sum(r['signal'] for r in signals if r['comparison']==l and r['catalogue']==cat and r['mesh']==mesh) for l in ['A','B'] for cat in ['primary','fine'] for mesh in [80001,160001]};reversal=any((counts[f'{l}_primary_160001']>=2)!=(counts[f'{l}_{cat}_{mesh}']>=2) for l in ['A','B'] for cat in ['primary','fine'] for mesh in [80001,160001]);numeric='NUMERICAL_OR_GRID_LIMITED' if reversal else 'NO_MAIN_SIGNAL_REVERSAL'
 dump(OUT/'DECISION.json',{'stage':'T2-TILT-G0','A': 'A_RANGE_SUPPORT_SIGNAL' if counts['A_primary_160001']>=2 else 'NOT_SUPPORTED_IN_TESTED_STATIC_CATALOGUE','B':'B_RANGE_SUPPORT_SIGNAL' if counts['B_primary_160001']>=2 else 'NOT_SUPPORTED_IN_TESTED_STATIC_CATALOGUE','primary_geometry_counts':{'A':counts['A_primary_160001'],'B':counts['B_primary_160001'],'denominator':3},'all_numeric_signal_counts':counts,'numerical_status':numeric,'propagation_status':'CONDITIONAL_MODAL_WINDOW_LIMIT','statistical_scope':'STATIC_EXPLORATORY_FINITE_CATALOGUE_ONLY','test_base_innovations':144,'test_resource_configs':5184,'paired_numeric_outputs':len(ev),'independent_horizontal_geometries':3,'R4_percent':0,'original_H3_P0_P1_SP1_RX_AUX_gate_unchanged':True,'research_lead_audit':'PENDING','state':'STOP','next_stage_authorized':False})
 print('EVALUATED',counts,'reversal',reversal,flush=True)
if __name__=='__main__':main()