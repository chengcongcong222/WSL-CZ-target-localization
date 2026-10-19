from common import *
def main():
 records=list(csv.DictReader((O/'P2_RECORD_INDEX.csv').open(encoding='utf-8')));rows=[];calrows=[]
 for rec in records:
  token=rec['token'];truth=np.load(RAW/'p2/truth_private'/(token+'.npz'))
  for mode in ['CAL','UNCAL']:
   z=json.loads((RAW/'p2/products'/(token+'_'+mode+'.json')).read_text());beta_true=float(truth['beta'])
   calrows.append(dict(token=token,mode=mode,beta_error_deg=float(np.rad2deg(z['epochs'][0]['beta_est_rad']-beta_true)),cal_residual=z['epochs'][0]['calibration_residual'],nav_model='common20m+independent15m each axis; correlated nodes; 25m marginal',clock_model='unknown fixed Uniform ±0.25s',association='all spectral candidates retained; no target identifier',signed_bearing='NOT_PROVIDED'))
   for k,e in enumerate(z['epochs']):
    actual=float(truth['angle'][k]);h=float(axis([truth['true_times'][k]])[0])+e['beta_est_rad'];true_u=np.cos(actual-h);ps=e['products'];present=bool(truth['present'][k]);false=bool(truth['false'][k]);targetf=float(truth['frequency'][k])
    matching=[p for p in ps if abs(p['f_rx_hz']-targetf)<=.75]
    valid=bool(matching) and present
    err=min((abs(np.rad2deg(np.arccos(np.clip(true_u,-1,1)))-p['axis_angle_deg']) for p in matching),default=np.inf) if present else None
    # Orientation uncertainty is a design bound from separate calibration, not truth clipping.
    du=np.deg2rad(.1);lo=min(np.cos(actual-h-du),np.cos(actual-h+du),true_u);hi=max(np.cos(actual-h-du),np.cos(actual-h+du),true_u)
    covered=present and any(p['u_lo']<=hi and p['u_hi']>=lo for p in matching)
    rows.append(dict(token=token,split=rec['split'],state_id=rec['state_id'],repeat=rec['repeat'],node=rec['node'],snr=int(rec['snr']),stress=rec['stress'],mode=mode,epoch=k,present=present,false_source=false,detected=bool(ps),target_component_extracted=valid,target_u_covered=covered,axis_error_deg=('INF' if err is not None and not np.isfinite(err) else err),n_frequency_branches=len(ps),n_signed_angle_branches=2*len(ps),source_row_kind='evaluation_only_true_frequency_not_used_in_estimator'))
 tab('P2_RECEIVER_BEARING_CASES.csv',rows);tab('P2_BIAS_NAV_SYNC_ASSOCIATION.csv',calrows)
 summary=[]
 for split in ['development','holdout','STRESS_CONTROL']:
  for stress in sorted(set(r['stress'] for r in rows if r['split']==split)):
   for snr in [0,10,20]:
    for mode in ['CAL','UNCAL']:
     z=[r for r in rows if r['split']==split and r['stress']==stress and r['snr']==snr and r['mode']==mode];p=[r for r in z if r['present']];absent=[r for r in z if not r['present']]
     k=sum(r['target_component_extracted'] for r in p);cov=sum(r['target_u_covered'] for r in p);fa=sum(r['detected'] for r in absent)
     errors=[np.inf if r['axis_error_deg']=='INF' else r['axis_error_deg'] for r in p]
     summary.append(dict(split=split,stress=stress,snr=snr,mode=mode,epochs=len(z),present_epochs=len(p),extracted=k,covered=cov,absent_epochs=len(absent),false_alarm=fa,extract_fraction=k/len(p) if p else None,cover_fraction=cov/len(p) if p else None,axis_error_P50_deg=qtiles(errors)[0],axis_error_P95_deg=qtiles(errors)[2],epoch_Wilson95_lower=wilson(k,len(p))[0],epoch_Wilson95_upper=wilson(k,len(p))[1],independence_warning='epochs clustered by trajectory; Wilson descriptive, not independent-track coverage'))
 tab('P2_EXTRACTION_SUMMARY.csv',summary)
 chosen=[r for r in summary if r['split'] in ['development','holdout'] and r['stress']=='MAIN_TRACK' and r['snr']==20 and r['mode']=='CAL']
 noise=[r for r in summary if r['stress']=='NOISE_ONLY' and r['mode']=='CAL']
 cal=[r for r in calrows if r['mode']=='CAL'];calok=sum(abs(r['beta_error_deg'])<=.1 and r['cal_residual']<=.02 for r in cal)/len(cal)
 positive=all(r['extract_fraction']>=.9 and r['cover_fraction']>=.95 and r['axis_error_P50_deg']!='INF' and r['axis_error_P50_deg']<=10 for r in chosen)
 noisek=all(r['false_alarm']/max(1,r['absent_epochs'])<=.05 for r in noise)
 gate=dict(G2a_nonoracle=positive and noisek and calok>=.99,G2b_scope='all branches/associations retained; independent calibration required; uncalibrated orientation unbounded',G2c='ARRAY_FRONTEND_CONDITIONAL_NOT_CZ_MULTIPATH',PHYSICAL_CZ_LINK='UNVERIFIED',calibration_within_design_fraction=calok,positive_main_gate=positive,noise_false_alarm_gate=noisek,P3_eligible=['S0_CAL_20','S2_B0_CAL_20','S2_B3_CAL_20'] if positive and noisek and calok>=.99 else [],unqualified_stress_kept=True,z='NOT_EVALUATED',source_truth_read_by_frontend=False)
 write('P2_GATE.json',gate)
 tab('P2_AMBIGUITY_BRANCH_REPORT.csv',[dict(token=r['token'],mode=r['mode'],epoch=r['epoch'],frequency_branches=r['n_frequency_branches'],signed_branches=r['n_signed_angle_branches'],association='UNASSIGNED_ALL_RETAINED') for r in rows])
 tab('P2_HORIZONTAL_PERFORMANCE.csv',[dict(stage='P2',status='NOT_EVALUATED',reason='P3 conditional branch search pending; no oracle bearing substitution')])
 print(gate,flush=True)
if __name__=='__main__':main()
