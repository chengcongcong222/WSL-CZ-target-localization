"""Saved-scene Gaussian MLE and local four-state information audit."""
from pathlib import Path
from collections import defaultdict
import argparse, math
import numpy as np
from scipy.optimize import least_squares
from r4_a1_new_estimator import SCALE, initialize, residual, jacobian
from r4_a1_new_dynamic import sha,load,dump,rows,write,git,scene,trajectories,evaluate,METRICS,quant
OUT=Path('results/R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY')
OLD=Path('results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC')
DIAG=Path('results/R4_A1_NEW_SPEED_BOTTLENECK_DIAGNOSTIC')
PARENT='6180e27edadb69aec1f68959b60bdb2036e719e5'
TITLE='R4 A1-new: freeze speed information efficiency audit'
WINDOWS=['PREFIX_300','STRAIGHT_600','PREFIX_900','FULL_1200','POST_TURN']
SCOPE='LOCAL_LINEAR_GAUSSIAN_EQUIVALENT;NOT_EMPIRICAL_P95;NOT_GLOBAL_GUARANTEE'

def truth_state(truth):
 r=float(truth['r0_km'])*1000;theta=math.radians(float(truth['theta0_deg']));v=float(truth['v_mps']);psi=math.radians(float(truth['psi_deg']))
 return np.array([r*math.cos(theta),r*math.sin(theta),v*math.cos(psi),v*math.sin(psi)])

def saved_scene(truth,a,u,z,side):
 t,_,angles,e=scene(truth,a,u,z,side)
 bias=np.deg2rad([e['b_common_deg']-e['b_diff_HALF_deg'],e['b_common_deg']+e['b_diff_HALF_deg']])
 _,main=trajectories();beta=math.radians(e['beta_deg']);offset=a['baseline_m']*np.array([-math.sin(beta),side*math.cos(beta)])
 return t,np.stack([main,main+offset],axis=1),angles-bias,e

def _solve(t,nodes,angles,sigma,initial):
 fit=least_squares(residual,initial/SCALE,jac=jacobian,args=(t,nodes,angles,sigma),loss='linear',f_scale=1.5,method='trf',max_nfev=100,ftol=1e-10,xtol=1e-10,gtol=1e-10)
 state=fit.x*SCALE;rr=residual(fit.x,t,nodes,angles,sigma)
 diag=dict(solver_success=bool(fit.success and np.isfinite(state).all()),solver_status=int(fit.status),nfev=int(fit.nfev),cost=float(fit.cost),optimality=float(fit.optimality),residual_RMS_deg=float(np.sqrt(np.mean(rr**2))*sigma))
 diag.update({k:float(v) for k,v in zip(['initial_x_m','initial_y_m','initial_vx_mps','initial_vy_mps'],initial)})
 return state,diag

def estimate_gaussian(t,node_positions,bearings,sigma_deg):
 initial,_=initialize(np.asarray(t),np.asarray(node_positions),np.asarray(bearings))
 return _solve(t,node_positions,bearings,sigma_deg,initial)

def estimate_oracle(t,node_positions,bearings,sigma_deg,true_cartesian_state):
 """ORACLE_INITIALIZATION_DIAGNOSTIC_ONLY; never an estimator claim."""
 return _solve(t,node_positions,bearings,sigma_deg,np.asarray(true_cartesian_state))

def masks(t):
 return [t<=300,t<=600,t<=900,t<=1200,t>600]

def information(t,nodes,sigma,state):
 # Known bearing sigma; no fitted residual-variance rescaling.
 J=jacobian(state/SCALE,t,nodes,np.zeros((len(t),2)),sigma)
 _,s,Vt=np.linalg.svd(J,full_matrices=False);rank=int(np.sum(s/s[0]>1e-10))
 Fscaled=J.T@J;H=J/SCALE;F=H.T@H
 weights=np.where(s/s[0]>1e-10,1/(s*s),0.)
 C=(SCALE[:,None]*((Vt.T*weights)@Vt))*SCALE[None,:]
 fsv=np.linalg.svd(F,compute_uv=False)
 x,y,vx,vy=state;r=math.hypot(x,y);v=math.hypot(vx,vy)
 gradients=[np.array([x/r,y/r,0,0]),np.array([-y/(r*r),x/(r*r),0,0]),np.array([0,0,vx/v,vy/v]),np.array([0,0,-vy/(v*v),vx/(v*v)])]
 variances=[max(0.,float(g@C@g)) for g in gradients] if rank==4 else [math.inf]*4
 sr,st,sv,sp=np.sqrt(variances)
 # Velocity information after profiling the two initial-position coordinates.
 Schur=F[2:,2:]-F[2:,:2]@np.linalg.solve(F[:2,:2],F[:2,2:])
 row=dict(n_epochs=len(t),n_measurements=2*len(t),rank_scaled_H=rank,full_rank=rank==4,condition_scaled_H=float(s[0]/s[-1]),condition_physical_F=float(fsv[0]/fsv[-1]),smallest_largest_scaled_H=float(s[-1]/s[0]),
  speed_variance_mps2=variances[2],sigma_speed_mps=sv,sigma_speed_rel=sv/v,speed_rel_1645=1.645*sv/v,speed_rel_196=1.96*sv/v,
  sigma_heading_deg=math.degrees(sp),heading_deg_196=1.96*math.degrees(sp),sigma_range_m=sr,range_rel_196=1.96*sr/r,sigma_bearing_deg=math.degrees(st),bearing_deg_196=1.96*math.degrees(st),
  velocity_schur_trace=float(np.trace(Schur)),scope=SCOPE,covariance_kind='INVERSE_F' if rank==4 else 'PSEUDOINVERSE_NOT_FINITE_BOUND')
 row.update({f'scaled_H_singular_{i+1}':float(s[i]) for i in range(4)})
 row.update({f'physical_F_singular_{i+1}':float(fsv[i]) for i in range(4)})
 return row,F,C

def summarize(data,info):
 groups=defaultdict(list);ig=defaultdict(list)
 for x in data:groups[x['anchor'],x['case_id'],x['variant']].append(x)
 for x in info:ig[x['anchor'],x['case_id'],x['window']].append(x)
 summary=[];eff=[]
 for key,g in groups.items():
  a,c,v=key;row=dict(anchor=a,case_id=c,variant=v,n_runs=len(g),n_failed=sum(x['failed'] for x in g),failure_rate=sum(x['failed'] for x in g)/len(g))
  for m in METRICS+['absolute_speed_error_mps']:
   for tag,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99)]:row[m+'_'+tag]=quant([x[m] for x in g],q)
  row.update(speed_rel_P95=row['speed_error_P95'],speed_abs_P95_mps=row['absolute_speed_error_mps_P95'],PROJECT_speed_pass=float(row['speed_error_P95'])<=.1,STRONG_speed_pass=float(row['speed_error_P95'])<=.05,
   PROJECT_other_metrics_pass=all(float(row[m+'_P95'])<=lim for m,lim in [('range_error',.1),('bearing_error_deg',1),('heading_error_deg',5)]),
   route_failure_rate_pass=row['failure_rate']<=.01)
  good=[x for x in g if not x['failed']];signed=np.array([x['signed_speed_error_mps'] for x in good])
  stats=dict(n_valid=len(good),signed_stats_scope='FINITE_ONLY_DESCRIPTIVE;FAILURES_RETAINED_IN_GATES',mean_speed_bias_mps=float(signed.mean()) if len(signed) else 'NA',
    median_speed_bias_mps=float(np.median(signed)) if len(signed) else 'NA',speed_std_mps=float(signed.std(ddof=1)) if len(signed)>1 else 'NA',speed_RMSE_mps=float(np.sqrt(np.mean(signed*signed))) if len(signed) else 'NA')
  row.update(stats);summary.append(row)
  cr=ig[a,c,'FULL_1200'];variance=float(np.mean([x['speed_variance_mps2'] for x in cr]))
  eligible=len(good)==len(g) and len(g)>1 and math.isfinite(variance) and float(stats['speed_std_mps'])>0
  eff.append(dict(anchor=a,case_id=c,variant=v,**stats,mean_conditional_local_CRLB_variance_mps2=variance,
   eta_v=variance/(float(stats['speed_std_mps'])**2) if eligible else 'NA',efficiency_eligible=eligible,
   bias_to_std=abs(float(stats['mean_speed_bias_mps']))/float(stats['speed_std_mps']) if eligible else 'NA',
   scope='MEAN_CONDITIONAL_LOCAL_UNBIASED_VARIANCE_REFERENCE;NOT_TOTAL_ERROR_BOUND;BIAS_REPORTED_SEPARATELY'))
 cs=[]
 for (a,c,w),g in ig.items():
  row=dict(anchor=a,case_id=c,window=w,n_saved_geometries=len(g),n_epochs=g[0]['n_epochs'],full_rank_rate=sum(x['full_rank'] for x in g)/len(g),
   min_rank=min(x['rank_scaled_H'] for x in g),scope=SCOPE,case_route_statistic='MAX_OVER_SAVED_DEPLOYMENT_ANGLES_AND_BOTH_MIRRORS')
  for m in ['sigma_speed_mps','sigma_speed_rel','speed_rel_1645','speed_rel_196','heading_deg_196','range_rel_196','bearing_deg_196','speed_variance_mps2','velocity_schur_trace']:
   vals=[x[m] for x in g];row[m+'_min']=min(vals);row[m+'_median']=quant(vals,.5);row[m+'_P95']=quant(vals,.95);row[m+'_max']=max(vals)
  row['PROJECT_CRLB_equivalent_pass']=float(row['speed_rel_196_max'])<=.1;row['STRONG_CRLB_equivalent_pass']=float(row['speed_rel_196_max'])<=.05;cs.append(row)
 return summary,eff,cs

def decide(summary,cr):
 groups={v:[x for x in summary if x['anchor']=='A' and x['variant']==v] for v in ['E0','E1','E2']}
 full=[x for x in cr if x['anchor']=='A' and x['window']=='FULL_1200']
 p={v:sum(x['PROJECT_speed_pass'] for x in g) for v,g in groups.items()}
 gap=max(abs(float(x['speed_rel_P95'])-float(next(y for y in groups['E2'] if y['case_id']==x['case_id'])['speed_rel_P95'])) for x in groups['E1'])
 rank_ok=all(x['full_rank_rate']==1 for x in full);above10=sum(float(x['speed_rel_196_max'])>.1 for x in full);above5=sum(float(x['speed_rel_196_max'])>.05 for x in full)
 rescue=p['E1']==12 and all(x['PROJECT_other_metrics_pass'] and x['route_failure_rate_pass'] for x in groups['E1'])
 if rescue:decision='GAUSSIAN_MLE_SPEED_ROUTE_WORTH_FRESH_VALIDATION';next_stage='FRESH_VALIDATION'
 elif not rank_ok:decision='CRLB_NUMERICAL_RANK_NOT_CLOSED';next_stage='STOP'
 elif p['E1']<12 and gap<=.02 and above10>0:decision='CURRENT_1200S_GEOMETRY_INFORMATION_LIMIT_CONFIRMED_BY_CRLB';next_stage='OBSERVATION_DESIGN'
 elif p['E1']<12 and above10==0:decision='ESTIMATOR_OR_NONLINEAR_EFFICIENCY_GAP_REMAINS';next_stage='ESTIMATOR_RESEARCH'
 elif gap>.02:decision='OPTIMIZER_OR_INITIALIZER_INEFFICIENCY_PRESENT';next_stage='ESTIMATOR_RESEARCH'
 else:decision='INFORMATION_EFFICIENCY_AUDIT_INCONCLUSIVE';next_stage='STOP'
 return dict(stage='R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY_AUDIT',parent_SHA=PARENT,parent_audit='ACCEPTED_CURRENT_1200S_BOTTLENECK',primary_anchor='A',secondary_anchor='B',nature='DEVELOPMENT_INFORMATION_EFFICIENCY_AUDIT_NOT_CONFIRMATION',route_decision=decision,
  next_recommended=next_stage,observation_design_stage='R4_A1_NEW_SPEED_OBSERVATION_DESIGN_GATE' if next_stage=='OBSERVATION_DESIGN' else 'NOT_OPENED',
  PROJECT_speed_cases=p,E1_STRONG_speed_cases=sum(x['STRONG_speed_pass'] for x in groups['E1']),E1_E2_max_case_P95_gap=gap,CRLB_cases_above_10=above10,CRLB_cases_above_5=above5,
  CRLB_scope=SCOPE,information_finding_scope='At least one saved primary geometry has local Gaussian 1.96-sigma speed equivalent above target; not an empirical P95 lower bound or universal physical impossibility.',
  new_Monte_Carlo_realizations=0,new_truth_cases=0,new_bearing_noise_draws=0,R4_A1_percent=0,R4_percent=0,depth='NOT_OPENED',automatic_next_stage=False,stop_after_commit_B=True,audit_status='PENDING_RESEARCH_LEAD_INDEPENDENT_AUDIT')

def execute():
 p=load(OUT/'SPEED_INFO_AUDIT_POLICY.json')
 for name,h in {**p['bindings'],**p['historical_bindings']}.items():assert sha(name)==h,name
 head=git('rev-parse','HEAD');assert head==git('ls-remote','origin','refs/heads/main').split()[0] and git('log','-1','--format=%s')==TITLE and git('rev-parse','HEAD^')==PARENT
 assert not (OUT/'SPEED_INFO_EXECUTION_START.json').exists(),'One saved-scene pass only'
 dump(OUT/'SPEED_INFO_EXECUTION_START.json',dict(policy_SHA=head,remote_verified_before_resolves_and_CRLB=True,policy_sha256=sha(OUT/'SPEED_INFO_AUDIT_POLICY.json')))
 source=rows(DIAG/'SPEED_VARIANT_RUN_RESULTS.csv');d3=[x for x in source if x['variant']=='D3'];write(OUT/'E0_D3_IDENTITY.csv',d3)
 anchors=load(OLD/'A1_NEW_DESIGN_FREEZE.json')['anchors'];panel=rows(OLD/'A1_NEW_TRUTH_PANEL.csv');data=[];info=[];Fs=[];Cs=[];paired=[]
 for a in ['A','B']:
  draws=np.load(OLD/f'A1_NEW_DRAWS_{a}.npz');uniforms=draws['uniforms'];normals=draws['normals'];base=[x for x in d3 if x['anchor']==a];anchor=anchors[a]
  for ix,old in enumerate(base):
   truth=panel[ix//500];side=int(old['mirror']);t,nodes,angles,errors=saved_scene(truth,anchor,uniforms[ix],normals[ix],side);true=truth_state(truth)
   key=dict(anchor=a,case_id=truth['case_id'],realization=int(old['realization']),draw_index=ix,mirror=side)
   fits={}
   for variant in ['E0','E1','E2']:
    diag=dict(solver_success=False,solver_status='NA',nfev='NA',cost='NA',optimality='NA',residual_RMS_deg='NA',initial_x_m='NA',initial_y_m='NA',initial_vx_mps='NA',initial_vy_mps='NA');state=np.full(4,np.nan);failed=True;reason='UNAVAILABLE'
    try:
     if variant=='E0':
      failed=old['failed']=='True';reason=old['failure_reason'];state=np.array([float(old[k]) for k in ['x_hat_m','y_hat_m','vx_hat_mps','vy_hat_mps']]) if not failed else state;metrics={k:float(old[k]) for k in METRICS};metrics.update({k:float(old[k]) if old[k]!='NA' else 'NA' for k in ['r_hat_km','theta_hat_deg','v_hat_mps','psi_hat_deg']});diag['solver_success']=not failed
     else:
      state,diag=estimate_gaussian(t,nodes,angles,anchor['sigma_deg']) if variant=='E1' else estimate_oracle(t,nodes,angles,anchor['sigma_deg'],true)
      failed=not diag['solver_success'];reason='NONE' if not failed else 'NONCONVERGENCE';metrics=evaluate(state,truth,nodes[0,0])
    except (ValueError,np.linalg.LinAlgError,FloatingPointError) as ex:failed=True;reason=str(ex);metrics={k:'INF' for k in METRICS};metrics.update({k:'NA' for k in ['r_hat_km','theta_hat_deg','v_hat_mps','psi_hat_deg']})
    if failed:
     for m in METRICS:metrics[m]='INF'
    signed='INF' if failed else float(metrics['v_hat_mps'])-float(truth['v_mps'])
    row=dict(**key,variant=variant,scope='E0_D3_IDENTITY' if variant=='E0' else 'OBSERVATION_ONLY_DEVELOPMENT' if variant=='E1' else 'ORACLE_INITIALIZATION_DIAGNOSTIC_ONLY',failed=failed,failure_reason=reason,
     **{k:float(v) if np.isfinite(v) else 'NA' for k,v in zip(['x_hat_m','y_hat_m','vx_hat_mps','vy_hat_mps'],state)},**metrics,signed_speed_error_mps=signed,absolute_speed_error_mps='INF' if failed else abs(signed),**diag)
    data.append(row);fits[variant]=row
   for mask,window in zip(masks(t),WINDOWS):
    stats,F,C=information(t[mask],nodes[mask],anchor['sigma_deg'],true);info.append(dict(**key,window=window,matrix_index=len(Fs),deployment_beta_deg=errors['beta_deg'],**stats));Fs.append(F);Cs.append(C)
   x,y=fits['E1'],fits['E2'];valid=not x['failed'] and not y['failed']
   paired.append(dict(**key,both_success=valid,absolute_signed_speed_difference_mps=abs(float(x['v_hat_mps'])-float(y['v_hat_mps'])) if valid else 'INF',absolute_cost_difference=abs(x['cost']-y['cost']) if valid else 'INF',max_scaled_state_difference=max(abs(float(x[k])-float(y[k]))/scale for k,scale in zip(['x_hat_m','y_hat_m','vx_hat_mps','vy_hat_mps'],SCALE)) if valid else 'INF'))
   if (ix+1)%500==0:print(f'Saved information audit {a} {ix+1}/6000; E0/E1/E2 and FIM complete',flush=True)
 write(OUT/'MLE_RUN_RESULTS.csv',data);write(OUT/'MLE_PAIRED_BASIN_DIAGNOSTICS.csv',paired);write(OUT/'CRLB_GEOMETRY_RUN_RESULTS.csv',info)
 np.savez_compressed(OUT/'CRLB_MATRIX_RECORDS.npz',F=np.asarray(Fs),C=np.asarray(Cs))
 summary,eff,cs=summarize(data,info);write(OUT/'MLE_CASE_SUMMARY.csv',summary);write(OUT/'SPEED_EFFICIENCY_SUMMARY.csv',eff)
 write(OUT/'CRLB_CASE_SUMMARY.csv',[x for x in cs if x['window']=='FULL_1200']);write(OUT/'CRLB_TIME_PREFIX.csv',[x for x in cs if x['window']!='POST_TURN']);write(OUT/'CRLB_SEGMENT_SUMMARY.csv',[x for x in cs if x['window'] in ['STRAIGHT_600','POST_TURN','FULL_1200']])
 turns=[];index={(x['anchor'],x['draw_index'],x['window']):x for x in info}
 for a in ['A','B']:
  for ix in range(6000):
   straight=index[a,ix,'STRAIGHT_600'];full=index[a,ix,'FULL_1200'];post=index[a,ix,'POST_TURN'];fi=Fs[full['matrix_index']];fs=Fs[straight['matrix_index']];fp=Fs[post['matrix_index']]
   turns.append(dict(anchor=a,case_id=full['case_id'],draw_index=ix,mirror=full['mirror'],F_additivity_max_abs=float(np.max(np.abs(fi-fs-fp))),speed_variance_reduction=1-full['speed_variance_mps2']/straight['speed_variance_mps2'],effective_velocity_information_trace_gain=full['velocity_schur_trace']/straight['velocity_schur_trace'],scope='ADDED_EPOCHS_AND_EXISTING_CHANGED_GEOMETRY;NOT_ISOLATED_TURN_CAUSAL_EFFECT'))
 write(OUT/'CRLB_INFORMATION_ACCUMULATION.csv',turns)
 d=decide(summary,cs);d['policy_SHA']=head;dump(OUT/'SPEED_INFO_DECISION.json',d);print(d['route_decision'],flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');args=ap.parse_args()
 if args.execute:execute()
 else:ap.error('--execute required')
