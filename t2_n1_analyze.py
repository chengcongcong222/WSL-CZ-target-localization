"""Once-only archive analysis: error distribution, deterministic score bounds and frozen probes."""
from t2_n1_common import *
from scipy.stats import spearmanr

def main():
 if (OUT/'ANALYSIS_STARTED.json').exists():raise RuntimeError('ONCE_ONLY analysis already started')
 verify_identities();dump(OUT/'ANALYSIS_STARTED.json',{'utc':utc(),'once_only':True,'new_solvers_MC_spectra':0});start=time.monotonic();summaries=[];aggregate=[];deciles=[];strata=[];joint=[];tops=[];bounds=[];data_files=[];max_old_delta=0.;total_count=0;total_above=0;group_fails=0
 for cat in CATS:
  states=catalog(cat);n=len(states);candidate_ids=np.arange(n)
  for res in RES:
   budget(start);hi=archived(res,cat,160001);lo=archived(res,cat,80001);bad=[]
   for mesh,h in [(160001,hi),(80001,lo)]:
    norm=np.linalg.norm(h,axis=-1);indices=np.argwhere((norm==0)|(~np.isfinite(norm))|(~np.isfinite(h).all(axis=-1)))
    for cid,fi in indices:bad.append({'resource':res,'catalogue':cat,'mesh':mesh,'candidate_id':int(cid),'frequency_hz':FREQ[int(fi)],'r_km':float(states[cid,0]),'abs_theta_deg':float(states[cid,1]),'z_m':float(states[cid,2]),'reason':'zero_or_nonfinite_spatial_vector'})
   if bad:table(OUT/'INVALID_VECTOR_IDS.csv',bad);raise RuntimeError('ARCHIVED_INPUT_INVALID_VECTOR preserved IDs, STOP')
   uh,ul,nh,nl,e=residual_error(hi,lo);old=np.load(ARCHIVE/f'SHAPE_ERROR_{res}_{cat}.npy',mmap_mode='r');delta=float(np.max(abs(e-old)));max_old_delta=max(max_old_delta,delta)
   if delta>TOL:raise RuntimeError('STORED_E_NOT_REPRODUCED')
   E=e.mean(axis=1);oldE=old.mean(axis=1);p=LOCAL/f'FIELDS_{res}_{cat}.npz';np.savez(p,e=e,norm_H=nh,norm_L=nl,E=E,E_from_stored_errors=oldE);data_files.append(p)
   for fi,f in enumerate(FREQ):
    prefix={'resource':res,'catalogue':cat,'frequency_hz':f};q=shape_stats(e[:,fi]);rho=float(spearmanr(e[:,fi],np.log10(nh[:,fi])).statistic);summaries.append(prefix|q|{'norm_min':float(nh[:,fi].min()),'norm_median':float(np.median(nh[:,fi])),'norm_max':float(nh[:,fi].max()),'spearman_error_vs_log10_norm':rho,'old_e_max_absolute_discrepancy':float(np.max(abs(e[:,fi]-old[:,fi]))),'worst_case_gate_fail':bool(q['max']>.002),'no_independence_claim':True});total_count+=n;total_above+=q['above_historical_0p002'];group_fails+=int(q['max']>.002)
    order=np.lexsort((candidate_ids,nh[:,fi]));labels=np.empty(n,int);labels[order]=np.arange(n)*10//n
    for d in range(10):
     sel=labels==d;deciles.append(prefix|{'norm_decile_low_to_high':d+1}|shape_stats(e[sel,fi])|{'norm_min':float(nh[sel,fi].min()),'norm_median':float(np.median(nh[sel,fi])),'norm_max':float(nh[sel,fi].max()),'tie_break':'candidate_id','weak_candidates_removed':0})
    masks=[('range','[45,50)',(states[:,0]>=45)&(states[:,0]<50)),('range','[50,55)',(states[:,0]>=50)&(states[:,0]<55)),('range','[55,60]',(states[:,0]>=55)&(states[:,0]<=60)),('angle','theta=0',states[:,1]==0),('angle','0<theta<5',(states[:,1]>0)&(states[:,1]<5)),('angle','theta=5',states[:,1]==5)]+[('depth',str(int(z)),states[:,2]==z) for z in [180,190,200,210,220]]
    for axis,label,sel in masks:strata.append(prefix|{'axis':axis,'stratum':label}|shape_stats(e[sel,fi]))
    for _,rl,rm in masks[:3]:
     for _,tl,tm in masks[3:6]:
      for z in [180,190,200,210,220]:joint.append(prefix|{'range_bin':rl,'angle_bin':tl,'source_depth_m':z}|shape_stats(e[rm&tm&(states[:,2]==z),fi]))
    topids=np.lexsort((candidate_ids,-e[:,fi]))[:5]
    for rank,cid in enumerate(topids,1):tops.append(prefix|{'rank':rank,'candidate_id':int(cid),'r_km':float(states[cid,0]),'abs_theta_deg':float(states[cid,1]),'z_m':float(states[cid,2]),'norm_H':float(nh[cid,fi]),'norm_L':float(nl[cid,fi]),'shape_error':float(e[cid,fi])})
   aggregate.append({'resource':res,'catalogue':cat}|shape_stats(e.ravel())|{'candidate_states':n,'frequencies':3,'not_independent_target_cases':True})
   bounds.append({'resource':res,'catalogue':cat}|quantiles(E)|{'E_above_0p002_descriptive':int(np.sum(E>.002)),'E_from_old_errors_max_absolute_discrepancy':float(np.max(abs(E-oldE))),'scope':'same_fixed_observation_between_two_archived_meshes','K8_divisor_sqrt8_used':False})
   print('SHAPE',res,cat,'vectors',3*n,'oldDelta',delta,flush=True)
 table(OUT/'FIELD_SHAPE_DISTRIBUTION.csv',summaries);table(OUT/'FIELD_SHAPE_DISTRIBUTION_AGGREGATED.csv',aggregate);table(OUT/'FIELD_NORM_DECILE_STRESS.csv',deciles);table(OUT/'STRATIFIED_FAILURE_COUNTS.csv',strata);table(OUT/'JOINT_GEOMETRY_FAILURE_COUNTS.csv',joint);table(OUT/'TOP_ERROR_CANDIDATES.csv',tops);table(OUT/'SCORE_PERTURBATION_BOUNDS_BY_RESOURCE.csv',bounds)
 maxerr=max(r['max'] for r in summaries);ok=total_count==433620 and total_above==8012 and group_fails==54 and abs(maxerr-.0073441923068232945)<=TOL
 dump(OUT/'HISTORICAL_COUNT_REPRODUCTION.json',{'PASS':ok,'frequency_vectors':total_count,'individual_above_0p002':total_above,'groupwise_max_failures':group_fails,'groupwise_denominator':54,'shape_error_max':maxerr,'old_error_max_absolute_discrepancy':max_old_delta,'original_admission_changed':False})
 if not ok:raise RuntimeError('HISTORICAL_COUNTS_NOT_REPRODUCED; no probe ranking')
 probe_rows=[];distance_rows=[];lookup={};max_excess=-float('inf');bound_checks=0
 for cat in CATS:
  states=catalog(cat);rs=np.unique(states[:,0]);nth=len(np.unique(states[:,1]));n=len(states)
  for res in RES:
   budget(start);uh,_=normalize(archived(res,cat,160001));ul,_=normalize(archived(res,cat,80001));E=np.load(LOCAL/f'FIELDS_{res}_{cat}.npz')['E'];stored_scores=[]
   for pid,r,t,z in PROBES:
    pids=np.flatnonzero(np.all(states==[r,t,z],axis=1));assert len(pids)==1;p=int(pids[0]);w=uh[p].copy();SH=probe_scores(uh,w);SL=probe_scores(ul,w);difference=abs(SH-SL);excess=difference-E;ex=float(excess.max());max_excess=max(max_excess,ex);bound_checks+=n
    if ex>TOL:raise RuntimeError('ANALYTIC_CANDIDATE_BOUND_VIOLATION')
    wrong=np.flatnonzero(abs(states[:,0]-r)>=1.);rh=wrong[ranking(SH[wrong])];rl=wrong[ranking(SL[wrong])];wh=int(rh[0]);wl=int(rl[0]);GH=float(SH[p]-SH[wh]);GL=float(SL[p]-SL[wl]);B=float(E[p]+E[wrong].max());gap_error=abs(GH-GL)
    if gap_error>B+TOL:raise RuntimeError('GAP_BOUND_VIOLATION')
    ranksH=ranking(SH);ranksL=ranking(SL);invH=np.empty(n,int);invL=np.empty(n,int);invH[ranksH]=np.arange(n);invL[ranksL]=np.arange(n);maxdiff=int(np.argmax(difference));others=wrong[wrong!=wh];best_lower=SH[wh]-E[wh];other_upper=float(np.max(SH[others]+E[others]));profileH=SH.reshape(len(rs),nth,5).max(axis=(1,2));profileL=SL.reshape(len(rs),nth,5).max(axis=(1,2));rrH=ranking(profileH);rrL=ranking(profileL);irH=np.empty(len(rs),int);irL=np.empty(len(rs),int);irH[rrH]=np.arange(len(rs));irL[rrL]=np.arange(len(rs))
    row={'resource':res,'catalogue':cat,'probe':pid,'probe_candidate_id':p,'probe_r_km':r,'probe_abs_theta_deg':t,'probe_z_m':z,'fixed_w_from_mesh':160001,'probe_score_H':float(SH[p]),'probe_score_L':float(SL[p]),'top_candidate_H':int(ranksH[0]),'top_candidate_L':int(ranksL[0]),'top_score_H':float(SH[ranksH[0]]),'top_score_L':float(SL[ranksL[0]]),'wrong_candidates':len(wrong),'strongest_wrong_ID_H':wh,'strongest_wrong_ID_L':wl,'strongest_wrong_r_H':float(states[wh,0]),'strongest_wrong_r_L':float(states[wl,0]),'strongest_wrong_abs_theta_H':float(states[wh,1]),'strongest_wrong_abs_theta_L':float(states[wl,1]),'strongest_wrong_z_H':float(states[wh,2]),'strongest_wrong_z_L':float(states[wl,2]),'strongest_wrong_score_H':float(SH[wh]),'strongest_wrong_score_L':float(SL[wl]),'wrong_top2_margin_H':float(SH[wh]-SH[rh[1]]),'wrong_top2_margin_L':float(SL[wl]-SL[rl[1]]),'G_H':GH,'G_L':GL,'G_sign_flip':bool(np.sign(GH)!=np.sign(GL)),'G_H_near_roundoff':abs(GH)<=TOL,'G_L_near_roundoff':abs(GL)<=TOL,'strongest_wrong_ID_changed':wh!=wl,'strongest_wrong_range_changed':float(states[wh,0])!=float(states[wl,0]),'E_probe':float(E[p]),'E_max_wrong':float(E[wrong].max()),'B_gap':B,'G_H_minus_L_absolute':gap_error,'gap_bound_covers':gap_error<=B+TOL,'max_candidate_score_difference':float(difference[maxdiff]),'max_score_change_candidate_id':maxdiff,'max_score_change_E':float(E[maxdiff]),'max_score_bound_excess':ex,'bound_violation_count_tolerance':int(np.sum(excess>TOL)),'state_rank_positions_changed':int(np.sum(invH!=invL)),'top10_state_overlap':len(set(ranksH[:10])&set(ranksL[:10])),'distance_rank_positions_changed':int(np.sum(irH!=irL)),'wrong_argmax_sign_robust_by_candidate_intervals':bool(best_lower>other_upper),'wrong_argmax_interval_clearance':float(best_lower-other_upper),'source_template_NOT_received_data':True}
    probe_rows.append(row);lookup[res,cat,pid]=row;stored_scores.append(np.stack([SH,SL]))
    for i,rv in enumerate(rs):distance_rows.append({'resource':res,'catalogue':cat,'probe':pid,'r_km':float(rv),'profiled_score_H':float(profileH[i]),'profiled_score_L':float(profileL[i]),'rank_H':int(irH[i]+1),'rank_L':int(irL[i]+1),'registered_wrong_range':bool(abs(rv-r)>=1.)})
   p=LOCAL/f'PROBE_SCORES_{res}_{cat}.npy';np.save(p,np.array(stored_scores));data_files.append(p);print('PROBES',res,cat,'5 full-catalogue fixed-H comparisons',flush=True)
 table(OUT/'TEMPLATE_PROBE_RANKING.csv',probe_rows);table(OUT/'DISTANCE_PROFILE_RANKING.csv',distance_rows);pairs=[]
 for label,treat,control in PAIRS:
  for cat in CATS:
   for pid,_,_,_ in PROBES:
    t=lookup[treat,cat,pid];c=lookup[control,cat,pid];advH=t['G_H']-c['G_H'];advL=t['G_L']-c['G_L'];B=t['B_gap']+c['B_gap'];robust=abs(advH)>B;d=abs(advH-advL)
    if d>B+TOL:raise RuntimeError('PAIR_ADVANTAGE_BOUND_VIOLATION')
    pairs.append({'pair':label,'treated':treat,'control':control,'catalogue':cat,'probe':pid,'G_treated_H':t['G_H'],'G_control_H':c['G_H'],'G_treated_L':t['G_L'],'G_control_L':c['G_L'],'Adv_H':advH,'Adv_L':advL,'B_treated':t['B_gap'],'B_control':c['B_gap'],'B_advantage_sum':B,'Adv_L_lower_from_H':advH-B,'Adv_L_upper_from_H':advH+B,'Adv_H_minus_L_absolute':d,'bound_covers':d<=B+TOL,'actual_advantage_sign_flip':bool(np.sign(advH)!=np.sign(advL)),'Adv_H_near_roundoff':abs(advH)<=TOL,'Adv_L_near_roundoff':abs(advL)<=TOL,'sign_robust_by_strict_bound':bool(robust),'registered_certificate':'SIGN_ROBUST_ACROSS_THE_TWO_ARCHIVED_MESHES_IN_THIS_TEMPLATE_PROBE' if robust else 'SIGN_UNRESOLVED_BY_BOUND','H_advantage_direction':'treated_positive' if advH>0 else 'control_positive' if advH<0 else 'zero','low_reference_signal_changed':False})
 table(OUT/'PAIR_ADVANTAGE_NUMERIC_INTERVALS.csv',pairs);assert len(probe_rows)==90 and len(pairs)==20
 tags={};global_sensitive=any(r['strongest_wrong_ID_changed'] or r['G_sign_flip'] for r in probe_rows) or any(r['actual_advantage_sign_flip'] for r in pairs)
 for label,treat,control in PAIRS:
  pp=[r for r in pairs if r['pair']==label];rr=[r for r in probe_rows if r['resource'] in [treat,control]];allrobust=all(r['sign_robust_by_strict_bound'] for r in pp) and len(set(r['H_advantage_direction'] for r in pp))==1;sensitive=any(r['strongest_wrong_ID_changed'] or r['G_sign_flip'] for r in rr) or any(r['actual_advantage_sign_flip'] for r in pp);labels=[]
  if allrobust:labels.append('ARCHIVED_TEMPLATE_NUMERIC_COMPARISON_ROBUST_IN_RECORDED_PROBES')
  if sensitive:labels.append('ARCHIVED_TEMPLATE_RANKING_OR_ADVANTAGE_SENSITIVE')
  if not allrobust:labels.append('ARCHIVED_TEMPLATE_NUMERIC_IMPACT_UNRESOLVED')
  tags[label]={'analysis_tags':labels,'strict_probe_certificates':sum(r['sign_robust_by_strict_bound'] for r in pp),'certificate_denominator':10,'actual_advantage_sign_flips':sum(r['actual_advantage_sign_flip'] for r in pp),'most_wrong_ID_changes_in_pair_resources':sum(r['strongest_wrong_ID_changed'] for r in rr),'gap_sign_flips_in_pair_resources':sum(r['G_sign_flip'] for r in rr),'all10_robust_same_direction':bool(allrobust),'G0_array_performance':'NOT_EVALUATED'}
 dump(OUT/'DECISION.json',{'stage':'T2-N1','status':'ARCHIVED_NUMERIC_DIAGNOSTIC_COMPLETE','pair_analysis':tags,'global_ranking_sensitive':global_sensitive,'HISTORICAL_G0_FORWARD_NUMERIC_NOT_ADMITTED_UNCHANGED':True,'A_NOT_EVALUATED':True,'B_NOT_EVALUATED':True,'CHIGH_UNCERTIFIED':True,'R4_percent':0,'new_KRAKEN_FIELD_other_propagation':0,'new_MC_waveform_spectrum_optimizer_FIM':0,'probe_observation_scope':'archived high-mesh templates, NOT real or simulated receiver measurements','research_lead_audit':'PENDING','state':'STOP','next_scientific_stage_authorized':False})
 dump(OUT/'SCORE_BOUND_CHECKS.json',{'candidate_probe_checks':bound_checks,'max_signed_bound_excess':max_excess,'all_candidate_bounds_covered_with_tolerance':max_excess<=TOL,'gap_checks':90,'pair_advantage_checks':20,'no_sqrtK_reduction':True,'all90_fixed_H_observation_between_models':True,'inequality_tolerance':TOL})
 evidence=[p for p in OUT.iterdir() if p.is_file() and p.suffix in ['.csv','.json'] and p.name not in ['ANALYSIS_FREEZE.json','EXECUTION_STARTED.json']];dump(OUT/'ANALYSIS_FREEZE.json',{'utc':utc(),'data_files':{str(p):sha(p) for p in data_files},'result_files':{p.name:sha(p) for p in evidence},'science_budget':budget(start),'normalization_or_ranking_reruns':0});print('ARCHIVE ANALYSIS SEALED',tags,flush=True)
if __name__=='__main__':main()