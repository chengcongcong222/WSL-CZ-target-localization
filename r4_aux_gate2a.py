"""Gate2A: preregistered single-factor nonideal geometry, no acoustic imports."""
from pathlib import Path
import argparse,math
import numpy as np
from r4_aux_gate1_crosstrack import sha,load,dump,read_csv,write_csv,git,quantiles,number
OUT=Path('results/R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET')
PARENT='1fde9dac9d848e5c8b9987b010182411fdbc79c6'
TITLE='R4 AUX Gate2A: freeze nonideal measurement error budget'
KEYS=['anchor','family','b_common_deg','b_diff_deg','sigma_pos_m','beta_deg']
def intersection(p1,p2,t1,t2,tol=1e-12):
 n=len(t1);p1=np.broadcast_to(p1,(n,2));p2=np.broadcast_to(p2,(n,2));d=p2-p1
 u1=np.column_stack((np.cos(t1),np.sin(t1)));u2=np.column_stack((np.cos(t2),np.sin(t2)));det=u1[:,0]*u2[:,1]-u1[:,1]*u2[:,0]
 parallel=np.abs(det)<=tol;safe=np.where(parallel,np.nan,det)
 l1=(d[:,0]*u2[:,1]-d[:,1]*u2[:,0])/safe;l2=(d[:,0]*u1[:,1]-d[:,1]*u1[:,0])/safe
 xy=p1+l1[:,None]*u1;behind=(l1<=0)|(l2<=0);nonfinite=~np.isfinite(xy).all(axis=1)
 return xy,np.abs(l1),parallel|behind|nonfinite,parallel,behind,nonfinite,det

def evaluate(row,z,p):
 a=p['anchors'][row['anchor']];r=float(row['range_km']);side=int(row['side']);b=a['baseline_km'];s=math.radians(a['sigma_random_deg']);beta=math.radians(float(row['beta_deg']))
 main_true=np.array([0.,0.]);aux_true=np.array([-b*math.sin(beta),side*b*math.cos(beta)])
 main_est=np.broadcast_to(main_true,(len(z),2)).copy();aux_est=np.broadcast_to(aux_true,(len(z),2)).copy()
 pos=float(row['sigma_pos_m'])/1000
 main_est+=pos*z[:,2:4]*np.array([1.,side]);aux_est+=pos*z[:,4:6]*np.array([1.,side])
 true2=math.atan2(-aux_true[1],r-aux_true[0]);common=math.radians(float(row['b_common_deg']));diff=math.radians(float(row['b_diff_deg']))
 theta1=side*s*z[:,0]+common-diff;theta2=true2+side*s*z[:,1]+common+diff
 xy,rhat,failed,parallel,behind,nonfinite,det=intersection(main_est,aux_est,theta1,theta2,p['parallel_det_tolerance'])
 er=np.abs(rhat-r)/r;position=np.linalg.norm(xy-np.array([r,0.]),axis=1)
 for v in [er,position]:v[failed]=np.inf
 out=dict(row);out.update(baseline_km=b,sigma_random_deg=a['sigma_random_deg'],true_aux_x_km=float(aux_true[0]),true_aux_y_km=float(aux_true[1]),crossing_angle_deg=abs(math.degrees(true2)),n_trials=len(z),execution_valid=True,
  failure_rate=float(failed.mean()),behind_sensor_rate=float(behind.mean()),parallel_rate=float(parallel.mean()),nonfinite_rate=float(nonfinite.mean()),ill_conditioned_rate=float((1/np.maximum(np.abs(det),np.finfo(float).tiny)>p['ill_condition_proxy_limit']).mean()))
 for name,v in [('relative_range',er),('position_2d_km',position)]:
  for tag,value in quantiles(v).items():out[name+'_'+tag]=number(value)
 n=len(z);q=np.sort(er);zz=1.959963984540054;delta=zz*math.sqrt(n*.95*.05)
 out['p95_rank_ci_low']=number(float(q[max(0,math.floor(n*.95-delta)-1)]));out['p95_rank_ci_high']=number(float(q[min(n-1,math.ceil(n*.95+delta)-1)]))
 for tag,t in [('5pct',.05),('10pct',.1)]:
  f=float((er<=t).mean());den=1+zz*zz/n;c=(f+zz*zz/(2*n))/den;h=zz*math.sqrt(f*(1-f)/n+zz*zz/(4*n*n))/den
  out['fraction_le'+tag]=f;out['fraction_'+tag+'_wilson_low']=c-h;out['fraction_'+tag+'_wilson_high']=c+h
 return out

def group_key(row):return tuple(str(row[k]) for k in KEYS)
def aggregate(data):
 groups={}
 for x in data:groups.setdefault(group_key(x),[]).append(x)
 summaries=[]
 for group in groups.values():
  x=group[0];worst=max(group,key=lambda y:float(y['relative_range_p95']));value=float(worst['relative_range_p95']);valid=len(group)==22 and len({float(y['range_km']) for y in group})==11 and all(y['execution_valid'] is True for y in group);finite=all(math.isfinite(float(y['relative_range_p95'])) for y in group)
  summaries.append(dict({k:x[k] for k in KEYS},baseline_km=x['baseline_km'],sigma_random_deg=x['sigma_random_deg'],n_ranges=11 if valid else len({y['range_km'] for y in group}),n_mirror_cells=len(group),all_execution_valid=valid,all_p95_finite=finite,worst_p95_relative=number(value),worst_range_km=worst['range_km'],worst_side=worst['side'],worst_p99_relative=number(max(float(y['relative_range_p99']) for y in group)),worst_p95_position_km=number(max(float(y['position_2d_km_p95']) for y in group)),max_failure_rate=max(float(y['failure_rate']) for y in group),full_range_5pct_pass=valid and finite and value<=.05,full_range_10pct_pass=valid and finite and value<=.1))
 return summaries

def cap_axis(rows,axis,tag):
 caps=sorted({abs(float(x[axis])) for x in rows});passing=[]
 for cap in caps:
  tested=[x for x in rows if abs(float(x[axis]))<=cap]
  if tested and all(x['full_range_'+tag+'_pass'] for x in tested):passing.append(cap)
 return max(passing) if passing else 'NOT_REACHED_ON_TESTED_AXIS'

def requirements(summaries,p):
 bias=[];position=[];deployment=[]
 for anchor in p['anchors']:
  group=[x for x in summaries if x['anchor']==anchor];row=dict(anchor=anchor,baseline_km=p['anchors'][anchor]['baseline_km'],sigma_random_deg=p['anchors'][anchor]['sigma_random_deg'])
  b=dict(row);nav=dict(row);dep=dict(row)
  common=[x for x in group if x['family']=='BIAS' and float(x['b_diff_deg'])==0];diff=[x for x in group if x['family']=='BIAS' and float(x['b_common_deg'])==0];pos=[x for x in group if x['family']=='POSITION'];beta=[x for x in group if x['family']=='DEPLOYMENT']
  for tag in ['5pct','10pct']:
   b['max_tested_abs_common_'+tag+'_deg']=cap_axis(common,'b_common_deg',tag);b['max_tested_abs_diff_'+tag+'_deg']=cap_axis(diff,'b_diff_deg',tag)
   nav['max_tested_sigma_pos_'+tag+'_m']=cap_axis(pos,'sigma_pos_m',tag);dep['max_tested_abs_beta_'+tag+'_deg']=cap_axis(beta,'beta_deg',tag)
  bias.append(b);position.append(nav);deployment.append(dep)
 return bias,position,deployment

def decision(bias,position,deployment,p):
 budgets=[]
 for b,n,g in zip(bias,position,deployment):
  vals=[b['max_tested_abs_common_5pct_deg'],b['max_tested_abs_diff_5pct_deg'],n['max_tested_sigma_pos_5pct_m'],g['max_tested_abs_beta_5pct_deg']]
  budgets.append(dict(anchor=b['anchor'],nonzero_all_four_axes=all(isinstance(v,(int,float)) and v>0 for v in vals)))
 numerical=any(x['nonzero_all_four_axes'] for x in budgets);confirmed=p['hardware_all_relevant_requirements_confirmed']
 outcome='AUX1_STATIC_MEASUREMENT_CHAIN_ADMITTED' if numerical and confirmed else 'AUX1_NONIDEAL_REQUIREMENTS_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED' if numerical else 'AUX1_5PCT_ARCHITECTURE_TOO_FRAGILE_UNDER_NONIDEAL_ERRORS'
 return dict(stage=p['stage'],parent_sha=PARENT,Gate1_independent_audit='ACCEPTED_TRADEOFF_ESTABLISHED',Gate2A_decision=outcome,bias_requirements=bias,position_requirements=position,deployment_requirements=deployment,numerical_anchor_admission=budgets,
  requirement_scope='SINGLE_FACTOR_AXES_ONLY; ALL11_TESTED_RANGES; NO_JOINT_BUDGET_OR_CONTINUOUS_TOLERANCE_GUARANTEE',actual_auxiliary_array_bearing_capability='UNKNOWN',actual_navigation_attitude_capability='UNKNOWN',hardware_all_relevant_requirements_confirmed=confirmed,
  client_confirmations_required=['auxiliary bearing-producing array / aperture / directed-bearing capability','per-node random DOA / heading / calibration error budget','residual common and differential bearing biases','per-axis navigation position accuracy','deployable5/7km lateral baseline and deployment orientation','target-bearing association and front-back resolution','common time reference / clock synchronization','array calibration and relative reference alignment'],R4_A1_NEW='NOT_OPENED',AUX_GATE2B='NOT_OPENED',depth='CLOSED_PENDING_INDEPENDENT_HORIZONTAL_ACQUISITION',R4_progress_percent=0,A2_A3_A4_B2_B3_B4_C='PAUSED / NOT_OPENED',P5='NOT_OPENED',new_acoustic_propagation=0,new_depth_scores=0,joint_stress_run=False,time_skew_or_association_mc=False,audit_status='PENDING_RESEARCH_LEAD_AUDIT',next_stage='STOP_FOR_RESEARCH_LEAD_AUDIT_AND_CLIENT_CONFIRMATION')

def execute():
 p=load(OUT/'AUX_GATE2A_DESIGN_FREEZE.json')
 for name,d in {**p['bindings'],**p['historical_bindings'],**p['external_bindings']}.items():assert sha(name)==d,name
 head=git('rev-parse','HEAD');remote=git('ls-remote','origin','refs/heads/main').split()[0]
 assert head==remote and git('log','-1','--format=%s')==TITLE and git('rev-parse','HEAD^')==PARENT
 assert not (OUT/'AUX_GATE2A_EXECUTION_START.json').exists(),'One frozen primary execution only'
 dump(OUT/'AUX_GATE2A_EXECUTION_START.json',dict(design_sha=head,remote_verified_before_mc=True,policy_sha256=sha(OUT/'AUX_GATE2A_DESIGN_FREEZE.json')))
 z=np.random.Generator(np.random.PCG64(p['seed'])).standard_normal((p['n_trials'],6));np.savez_compressed(OUT/'AUX_FROZEN_DRAWS.npz',z=z)
 data=[]
 for i,row in enumerate(read_csv(OUT/'AUX_GATE2A_ERROR_GRIDS.csv'),1):
  data.append(evaluate(row,z,p))
  if i%440==0:print(f'nonideal cells {i}/4268',flush=True)
 for family,name in [('BIAS','AUX_BIAS_RESULTS.csv'),('POSITION','AUX_POSITION_RESULTS.csv'),('DEPLOYMENT','AUX_DEPLOYMENT_RESULTS.csv')]:write_csv(OUT/name,[x for x in data if x['family']==family])
 summary=aggregate(data);write_csv(OUT/'AUX_FULL_RANGE_SUMMARY.csv',summary);bias,position,dep=requirements(summary,p)
 for name,items in [('AUX_BIAS_REQUIREMENT_FRONTIER.csv',bias),('AUX_POSITION_REQUIREMENT_FRONTIER.csv',position),('AUX_DEPLOYMENT_AZIMUTH_REQUIREMENT.csv',dep)]:write_csv(OUT/name,items)
 d=decision(bias,position,dep,p);d['design_sha']=head;dump(OUT/'AUX_GATE2A_DECISION.json',d);print(d['Gate2A_decision'],flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');args=ap.parse_args()
 if args.execute:execute()
 else:ap.error('--execute required')
