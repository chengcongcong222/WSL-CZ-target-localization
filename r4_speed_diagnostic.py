"""Frozen saved-scene causal diagnostic. Never creates draws or changes old scenes."""
from pathlib import Path
import argparse,math
import numpy as np
from r4_a1_new_dynamic import sha,load,dump,rows,write,git,scene,evaluate,METRICS,quant,trajectories
from r4_a1_new_estimator import estimate
from r4_speed_bias_estimator import estimate_bias_aware
OUT=Path('results/R4_A1_NEW_SPEED_BOTTLENECK_DIAGNOSTIC')
OLD=Path('results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC')
PARENT='e39182ac82b988e71b79d73a06a97076ca76f6bf'
TITLE='R4 A1-new: freeze speed bottleneck diagnostic'
VARIANTS=['D0','D1','D2','D3','D4']
def summary(data):
 out=[]
 for a in ['A','B']:
  for variant in VARIANTS:
   for case in sorted({x['case_id'] for x in data if x['anchor']==a}):
    group=[x for x in data if (x['anchor'],x['variant'],x['case_id'])==(a,variant,case)];row=dict(anchor=a,case_id=case,variant=variant,n_runs=len(group),failure_rate=sum(x['failed'] for x in group)/len(group))
    for m in METRICS+['absolute_speed_error_mps']:
     for tag,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99)]:row[m+'_'+tag]=quant([x[m] for x in group],q)
    for key,m in [('speed_rel_P95','speed_error'),('speed_abs_P95_mps','absolute_speed_error_mps'),('range_P95','range_error'),('bearing_P95','bearing_error_deg'),('heading_P95','heading_error_deg')]:row[key]=row[m+'_P95']
    row['PROJECT_speed_pass']=float(row['speed_rel_P95'])<=.1;row['STRONG_speed_pass']=float(row['speed_rel_P95'])<=.05;row['PROJECT_all_metrics_pass']=all(float(row[m+'_P95'])<=lim for m,lim in zip(METRICS,[.1,1,.1,5]));row['route_failure_rate_pass']=row['failure_rate']<=.01
    if variant=='D4':
     for m in ['common_bias_error_deg','diff_HALF_bias_error_deg']:row[m+'_P95']=quant([x[m] for x in group],.95)
     valid=[x for x in group if not x['failed']];row['common_at_bound_fraction_valid_fits']=sum(x['common_at_bound'] for x in valid)/len(valid) if valid else 'NA';row['diff_at_bound_fraction_valid_fits']=sum(x['diff_at_bound'] for x in valid)/len(valid) if valid else 'NA'
    else:
     for k in ['common_bias_error_deg_P95','diff_HALF_bias_error_deg_P95','common_at_bound_fraction_valid_fits','diff_at_bound_fraction_valid_fits']:row[k]='NOT_APPLICABLE'
    out.append(row)
 return out

def decide(summaries):
 A={variant:[x for x in summaries if x['anchor']=='A' and x['variant']==variant] for variant in VARIANTS}
 speed={v:sum(x['PROJECT_speed_pass'] for x in group) for v,group in A.items()};route=all(x['PROJECT_all_metrics_pass'] and x['route_failure_rate_pass'] for x in A['D4']);strong=all(x['STRONG_speed_pass'] for x in A['D4'])
 findings=[]
 if speed['D1']==12:findings.append('STATIC_BEARING_BIAS_IS_SUFFICIENT_TO_EXPLAIN_PROJECT_SPEED_FAILURE')
 if speed['D1']<12 and (speed['D2']==12 or speed['D3']==12):findings.append('NAVIGATION_ERROR_MATERIALLY_LIMITS_SPEED')
 if speed['D3']<12:findings.append('RANDOM_BEARING_PLUS_CURRENT_GEOMETRY_AND_4_STATE_ESTIMATOR_LIMIT_SPEED')
 if route:verdict='BIAS_AWARE_SPEED_ROUTE_WORTH_FRESH_VALIDATION'
 elif speed['D3']<12:verdict='SPEED_INFORMATION_INSUFFICIENT_AT_CURRENT_1200S_APPLICATION_SCENARIO'
 elif speed['D1']==12:verdict='BIAS_INFORMATION_PRESENT_BUT_ESTIMATOR_NOT_YET_ESTABLISHED'
 else:verdict='SPEED_ROUTE_NOT_ESTABLISHED_UNDER_CURRENT_SAVED_SCENE_DIAGNOSTIC'
 return dict(parent_A1_new_SHA=PARENT,A1_new_independent_audit='ACCEPTED_SPEED_BOTTLENECK',stage='R4_A1_NEW_SPEED_BOTTLENECK_DIAGNOSTIC',nature='DEVELOPMENT_CAUSAL_DIAGNOSTIC_NOT_CONFIRMATION',primary_anchor='A',secondary_anchor='B',primary_attribution=findings or ['NO_SINGLE_SOURCE_CLOSURE_ESTABLISHED'],attribution_scope='Frozen panel and current estimators only; oracle intervention findings are not unique-cause proofs. D3 failure is not universal physical nonidentifiability.',speed_route_decision=verdict,primary_PROJECT_speed_cases=speed,primary_STRONG_speed_cases_D4=sum(x['STRONG_speed_pass'] for x in A['D4']),D4_route_admitted=route,D4_strong_speed_diagnostic='BIAS_AWARE_ROUTE_MEETS_STRONG_SPEED_ON_DEVELOPMENT_DATA' if strong else 'NOT_ESTABLISHED',new_random_seeds=0,new_truth_cases=0,new_noise_draws=0,new_Monte_Carlo_realizations=0,new_application_performance_claim=0,R4_A1_percent=0,R4_percent=0,range='ESTABLISHED_AT_STRONG_TIER_IN_ACCEPTED_A1_NEW',bearing='ESTABLISHED_AT_STRONG_TIER_IN_ACCEPTED_A1_NEW',heading='ESTABLISHED_AT_PROJECT_TIER_NEAR_STRONG_IN_ACCEPTED_A1_NEW',speed='BOTTLENECK_DIAGNOSED_NO_SCIENTIFIC_CREDIT',next_recommended_stage='R4_A1_NEW_BIAS_AWARE_FRESH_VALIDATION' if route else 'STOP',automatic_fresh_validation=False,automatic_FIX2=False,depth='NOT_OPENED',new_acoustic_propagation=0,new_depth_score=0,audit_status='PENDING_RESEARCH_LEAD_INDEPENDENT_AUDIT',stop_after_commit_B=True)

def execute():
 p=load(OUT/'SPEED_DIAGNOSTIC_POLICY.json')
 for name,h in {**p['bindings'],**p['historical_bindings']}.items():assert sha(name)==h,name
 head=git('rev-parse','HEAD');assert head==git('ls-remote','origin','refs/heads/main').split()[0] and git('log','-1','--format=%s')==TITLE and git('rev-parse','HEAD^')==PARENT
 assert not (OUT/'SPEED_EXECUTION_START.json').exists(),'One saved-scene pass only'
 dump(OUT/'SPEED_EXECUTION_START.json',dict(policy_SHA=head,remote_verified_before_resolving=True,policy_sha256=sha(OUT/'SPEED_DIAGNOSTIC_POLICY.json')))
 # Identity byte copy before any branch resolves; no rerun as D0.
 (OUT/'D0_BASELINE_IDENTITY.csv').write_bytes((OLD/'A1_NEW_RUN_RESULTS.csv').read_bytes())
 oldpolicy=load(OLD/'A1_NEW_DESIGN_FREEZE.json');panel=rows(OLD/'A1_NEW_TRUTH_PANEL.csv');base=rows(OLD/'A1_NEW_RUN_RESULTS.csv');results=[];diagdata=[];jacdata=[]
 for a in ['A','B']:
  anchor=oldpolicy['anchors'][a];draws=np.load(OLD/f'A1_NEW_DRAWS_{a}.npz');original=[x for x in base if x['anchor']==a];t,main=trajectories()
  for ix,old in enumerate(original):
   truth=panel[ix//500];side=int(old['mirror']);tt,positions,bearings,errors=scene(truth,anchor,draws['uniforms'][ix],draws['normals'][ix],side);bias=np.deg2rad([errors['b_common_deg']-errors['b_diff_HALF_deg'],errors['b_common_deg']+errors['b_diff_HALF_deg']]);beta=math.radians(errors['beta_deg']);offset=anchor['baseline_m']*np.array([-math.sin(beta),side*math.cos(beta)]);true_positions=np.stack([main,main+offset],axis=1)
   for variant in VARIANTS:
    nodes=positions if variant in ['D0','D1','D4'] else true_positions;angles=bearings-bias if variant in ['D1','D3'] else bearings;key=dict(anchor=a,case_id=old['case_id'],realization=int(old['realization']),draw_index=ix,mirror=side,variant=variant,scope='IDENTITY_CONTROL' if variant=='D0' else 'CANDIDATE_DEVELOPMENT_ONLY' if variant=='D4' else 'ORACLE_DIAGNOSTIC_ONLY')
    diag={};jac=[];state=np.full(4,np.nan)
    try:
     if variant=='D0':
      failed=old['failed']=='True';reason=old['failure_reason'];state=np.array([float(old[k]) for k in ['x_hat_m','y_hat_m','vx_hat_mps','vy_hat_mps']]) if not failed else state;metrics={k:float(old[k]) for k in METRICS};metrics.update({k:old[k] for k in ['r_hat_km','theta_hat_deg','v_hat_mps','psi_hat_deg']})
     elif variant=='D4':
      state,diag,jac=estimate_bias_aware(t,nodes,angles,anchor['sigma_deg'],anchor['common_deg'],anchor['half_diff_deg']);failed=not diag['solver_success'];reason='NONE' if not failed else 'BIAS_AWARE_NONCONVERGENCE';metrics=evaluate(state,truth,nodes[0,0])
     else:
      state,diag=estimate(t,nodes,angles,anchor['sigma_deg']);failed=False;reason='NONE';metrics=evaluate(state,truth,nodes[0,0])
    except (ValueError,np.linalg.LinAlgError,FloatingPointError) as ex:failed=True;reason=str(ex);metrics={k:'INF' for k in METRICS};metrics.update({k:'NA' for k in ['r_hat_km','theta_hat_deg','v_hat_mps','psi_hat_deg']})
    if failed:
     for m in METRICS:metrics[m]='INF'
    absolute='INF' if failed else abs(float(metrics['v_hat_mps'])-float(truth['v_mps']))
    row=dict(**key,failed=failed,failure_reason=reason,x_hat_m=float(state[0]) if np.isfinite(state[0]) else 'NA',y_hat_m=float(state[1]) if np.isfinite(state[1]) else 'NA',vx_hat_mps=float(state[2]) if np.isfinite(state[2]) else 'NA',vy_hat_mps=float(state[3]) if np.isfinite(state[3]) else 'NA',**metrics,absolute_speed_error_mps=absolute,common_hat_deg=diag.get('common_hat_deg','NA'),diff_HALF_hat_deg=diag.get('diff_HALF_hat_deg','NA'),common_true_deg=errors['b_common_deg'],diff_HALF_true_deg=errors['b_diff_HALF_deg'],common_bias_error_deg=abs(diag['common_hat_deg']-errors['b_common_deg']) if variant=='D4' and not failed else 'INF' if variant=='D4' else 'NA',diff_HALF_bias_error_deg=abs(diag['diff_HALF_hat_deg']-errors['b_diff_HALF_deg']) if variant=='D4' and not failed else 'INF' if variant=='D4' else 'NA',common_at_bound=diag.get('common_at_bound',False),diff_at_bound=diag.get('diff_at_bound',False));results.append(row)
    if variant=='D4':
     fields=['solver_success','solver_status','nfev','cost','optimality','residual_RMS_deg','valid_intersections','invalid_intersections','initial_x_m','initial_y_m','initial_vx_mps','initial_vy_mps','initial_common_bias_deg','initial_diff_bias_deg','common_hat_deg','diff_HALF_hat_deg','common_at_bound','diff_at_bound'];diagdata.append(dict(**key,**{k:diag.get(k,'NA') for k in fields}))
     if not jac:
      jac=[dict(segment=s,n_epochs=n,numerical_rank=0,full_rank=False,smallest_largest_ratio=0.,condition_number='INF',**{f'singular_value_{j+1}':'NA' for j in range(6)}) for s,n in [('STRAIGHT',61),('POST_TURN',60),('FULL',121)]]
     jacdata.extend(dict(**key,failed=failed,**x) for x in jac)
   if (ix+1)%250==0:print(f'Saved scenes {a} {ix+1}/6000; D0-D4 complete',flush=True)
 write(OUT/'SPEED_VARIANT_RUN_RESULTS.csv',results);write(OUT/'BIAS_AWARE_ESTIMATOR_DIAGNOSTICS.csv',diagdata);write(OUT/'BIAS_AWARE_JACOBIAN_RUN_RESULTS.csv',jacdata);summaries=summary(results);write(OUT/'SPEED_VARIANT_CASE_SUMMARY.csv',summaries)
 attribution=[]
 for x in summaries:
  if x['variant']=='D0':continue
  baseline=next(y for y in summaries if y['anchor']==x['anchor'] and y['case_id']==x['case_id'] and y['variant']=='D0');b=float(baseline['speed_rel_P95']);v=float(x['speed_rel_P95']);attribution.append(dict(anchor=x['anchor'],case_id=x['case_id'],variant=x['variant'],D0_speed_rel_P95=b,variant_speed_rel_P95=v,P95_reduction_relative_to_D0=(b-v)/b,P95_absolute_difference=b-v,paired_saved_draws=True,causal_scope='INTERVENTION_COMPARISON; NONADDITIVE; NOT_VARIANCE_DECOMPOSITION'))
 write(OUT/'SPEED_ATTRIBUTION_SUMMARY.csv',attribution)
 js=[]
 for a in ['A','B']:
  for case in [x['case_id'] for x in panel]:
   for segment in ['STRAIGHT','POST_TURN','FULL']:
    group=[x for x in jacdata if (x['anchor'],x['case_id'],x['segment'])==(a,case,segment)];ratios=[x['smallest_largest_ratio'] for x in group];js.append(dict(anchor=a,case_id=case,segment=segment,n_runs=len(group),full_rank_rate=sum(x['full_rank'] for x in group)/len(group),ratio_min=min(ratios),ratio_median=quant(ratios,.5),ratio_P95=quant(ratios,.95),rank_min=min(x['numerical_rank'] for x in group),rank_median=quant([x['numerical_rank'] for x in group],.5),scope='SCALED_RAW_MEASUREMENT_JACOBIAN; NOT_ACCURACY_GATE; CONDITION_DEPENDS_ON_SCALING'))
 write(OUT/'BIAS_AWARE_JACOBIAN_SUMMARY.csv',js);d=decide(summaries);d['policy_SHA']=head;dump(OUT/'SPEED_DIAGNOSTIC_DECISION.json',d);print(d['speed_route_decision'],flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');args=ap.parse_args()
 if args.execute:execute()
 else:ap.error('--execute required')
