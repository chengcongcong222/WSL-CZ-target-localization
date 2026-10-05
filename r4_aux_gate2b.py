"""Gate2B frozen joint static geometry; historical Gate2A numerical core reused."""
from pathlib import Path
import argparse,itertools,math
import numpy as np
from r4_aux_gate2a import evaluate
from r4_aux_gate1_crosstrack import sha,load,dump,read_csv,write_csv,git,number
OUT=Path('results/R4_AUX_GATE2B_JOINT_ERROR_BUDGET')
PARENT='5d8b9c9f2d8efdb74a95da4ea546aae11457a24f'
TITLE='R4 AUX Gate2B: freeze joint static error budget'
FAMILIES=['T5_BUDGET_FAMILY','T10_BUDGET_FAMILY']
def grids(p,family):
 data=[]
 for anchor in p['anchors']:
  vector=p['vectors'][family][anchor]
  for lam in p['lambdas']:
   signs=[(0,0,0)] if lam==0 else list(itertools.product([-1,1],repeat=3))
   for sc,sd,sb in signs:
    for r in p['range_km']:
     for side in p['sides']:
      data.append(dict(cell_id=f'{family}-{len(data)+1:04d}',anchor=anchor,family=family,lambda_value=lam,sign_common=sc,sign_diff=sd,sign_beta=sb,b_common_deg=lam*vector[0]*sc,b_diff_deg=lam*vector[1]*sd,sigma_pos_m=lam*vector[2],beta_deg=lam*vector[3]*sb,range_km=r,side=side))
 return data

def aggregate(data,p):
 out=[]
 for family in FAMILIES:
  for anchor in p['anchors']:
   for lam in p['lambdas']:
    group=[x for x in data if x['family']==family and x['anchor']==anchor and float(x['lambda_value'])==lam]
    expected=22 if lam==0 else 176
    worst=max(group,key=lambda x:float(x['relative_range_p95']));v=float(worst['relative_range_p95'])
    valid=len(group)==expected and all(x['execution_valid'] is True for x in group)
    out.append(dict(family=family,anchor=anchor,lambda_value=lam,n_cells=len(group),all_execution_valid=valid,worst_joint_P95=number(v),worst_joint_P99=number(max(float(x['relative_range_p99']) for x in group)),worst_position_P95_km=number(max(float(x['position_2d_km_p95']) for x in group)),max_failure_rate=max(float(x['failure_rate']) for x in group),worst_range=worst['range_km'],worst_sign_common=worst['sign_common'],worst_sign_diff=worst['sign_diff'],worst_sign_beta=worst['sign_beta'],worst_mirror=worst['side'],full_range_T5_pass=valid and math.isfinite(v) and v<=.05,full_range_T10_pass=valid and math.isfinite(v) and v<=.1))
 return out

def decision(front,p):
 half={f:{a:next(x for x in front if x['family']==f and x['anchor']==a and x['lambda_value']==.5)['full_range_'+('T5' if f==FAMILIES[0] else 'T10')+'_pass'] for a in p['anchors']} for f in FAMILIES}
 primary=any(half[FAMILIES[0]].values());secondary=any(half[FAMILIES[1]].values())
 maxima=[]
 for f in FAMILIES:
  for a in p['anchors']:
   group=[x for x in front if x['family']==f and x['anchor']==a]
   maxima.append(dict(family=f,anchor=a,maximum_tested_lambda_T5=max([x['lambda_value'] for x in group if x['full_range_T5_pass']],default=None),maximum_tested_lambda_T10=max([x['lambda_value'] for x in group if x['full_range_T10_pass']],default=None)))
 return dict(stage=p['stage'],parent_sha=PARENT,Gate2A_independent_audit='ACCEPTED_SINGLE_FACTOR_REQUIREMENTS',primary_target='HALF_T5_JOINT_ENGINEERING_POINT',primary_decision='AUX1_HALF_T5_JOINT_BUDGET_ESTABLISHED' if primary else 'AUX1_HALF_T5_JOINT_BUDGET_NOT_ESTABLISHED',Gate2B_decision='AUX1_JOINT_STATIC_REQUIREMENT_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED' if primary else 'AUX1_JOINT_SUPPORTS_10PCT_NOT_5PCT' if secondary else 'AUX1_STATIC_GEOMETRY_TOO_FRAGILE_FOR_PROJECT_TARGET',half_point_pass=half,maximum_tested_lambdas=maxima,actual_hardware_capability='UNKNOWN',actual_auxiliary_array_bearing_capability='UNKNOWN',actual_navigation_attitude_capability='UNKNOWN',time_synchronization='NOT_NUMERICALLY_VALIDATED',target_association='NOT_NUMERICALLY_VALIDATED',automatic_further_numerical_work='STOP',AUX_GATE2C='NOT_OPENED',next_action='CLIENT_REQUIREMENT_CONFIRMATION_OR_RESEARCH_LEAD_ARCHITECTURE_DECISION',R4_A1_NEW='NOT_OPENED',R4_progress_percent=0,depth='CLOSED',new_acoustic_propagation=0,new_depth_scores=0,client_confirmations_required=p['client_confirmations_required'],audit_status='PENDING_RESEARCH_LEAD_AUDIT',scope='STATIC_SNAPSHOT; TWO_FROZEN_FAMILIES; TESTED_SIGN_CORNERS_AND_11_DISCRETE_RANGES_ONLY; NO_CONTINUOUS_HYPERBOX_GUARANTEE')

def execute():
 p=load(OUT/'AUX_GATE2B_DESIGN_FREEZE.json')
 for name,h in {**p['bindings'],**p['historical_bindings'],**p['external_bindings']}.items():assert sha(name)==h,name
 head=git('rev-parse','HEAD');remote=git('ls-remote','origin','refs/heads/main').split()[0]
 assert head==remote and git('log','-1','--format=%s')==TITLE and git('rev-parse','HEAD^')==PARENT
 assert not (OUT/'AUX_GATE2B_EXECUTION_START.json').exists(),'One primary execution only'
 dump(OUT/'AUX_GATE2B_EXECUTION_START.json',dict(design_sha=head,remote_verified_before_mc=True,policy_sha256=sha(OUT/'AUX_GATE2B_DESIGN_FREEZE.json')))
 z=np.random.Generator(np.random.PCG64(p['seed'])).standard_normal((p['n_trials'],6))
 np.savez_compressed(OUT/'AUX_FROZEN_DRAWS.npz',z=z,bearing=z[:,:2],main_position=z[:,2:4],aux_position=z[:,4:6])
 data=[]
 for family,stem in zip(FAMILIES,['T5','T10']):
  group=[]
  for i,row in enumerate(read_csv(OUT/f'AUX_{stem}_JOINT_GRID.csv'),1):
   group.append(evaluate(row,z,p))
   if i%220==0:print(f'{stem} joint cells {i}/1452',flush=True)
  write_csv(OUT/f'AUX_{stem}_JOINT_RESULTS.csv',group);data+=group
 front=aggregate(data,p)
 for f,name in zip(FAMILIES,['AUX_JOINT_LAMBDA_FRONTIER.csv','AUX_JOINT_T10_LAMBDA_FRONTIER.csv']):write_csv(OUT/name,[x for x in front if x['family']==f])
 d=decision(front,p);d['design_sha']=head;dump(OUT/'AUX_GATE2B_DECISION.json',d);print(d['Gate2B_decision'],flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');args=ap.parse_args()
 if args.execute:execute()
 else:ap.error('--execute required')
