import itertools
import numpy as np
from r4_aux_gate2b import grids,aggregate,decision,FAMILIES
from r4_aux_gate2a import evaluate
P=dict(anchors={'A':dict(baseline_km=5.,sigma_random_deg=.05),'B':dict(baseline_km=7.,sigma_random_deg=.075)},vectors={FAMILIES[0]:{'A':[.1,.05,50,20],'B':[.1,.05,50,20]},FAMILIES[1]:{'A':[.1,.1,100,20],'B':[.1,.1,200,20]}},lambdas=[0,.25,.5,.75,1],range_km=list(range(50,61)),sides=[-1,1],parallel_det_tolerance=1e-12,ill_condition_proxy_limit=1000,stage='TEST',client_confirmations_required=[])
def test_full_grid_and_family_vectors():
 for f in FAMILIES:
  g=grids(P,f);assert len(g)==1452
  for a in ['A','B']:
   for lam in P['lambdas']:
    group=[x for x in g if x['anchor']==a and x['lambda_value']==lam]
    assert len(group)==(22 if lam==0 else 176)
    assert {tuple(x[k] for k in ['sign_common','sign_diff','sign_beta']) for x in group}==({(0,0,0)} if lam==0 else set(itertools.product([-1,1],repeat=3)))
 assert next(x for x in grids(P,FAMILIES[1]) if x['anchor']=='B' and x['lambda_value']==.5)['sigma_pos_m']==100

def test_simultaneous_geometry_vs_independent_slope():
 row=next(x for x in grids(P,FAMILIES[0]) if x['lambda_value']==.5 and x['side']==1)
 z=np.zeros((100,6));z[:,2]=1;z[:,5]=-1
 result=evaluate(row,z,P);r=row['range_km'];b=5;beta=np.deg2rad(row['beta_deg']);aux=np.array([-b*np.sin(beta),b*np.cos(beta)]);pos=row['sigma_pos_m']/1000
 p1=pos*z[:,2:4];p2=aux+pos*z[:,4:6];m1=np.tan(np.deg2rad(row['b_common_deg']-row['b_diff_deg']));m2=np.tan(np.arctan2(-aux[1],r-aux[0])+np.deg2rad(row['b_common_deg']+row['b_diff_deg']));delta=p2-p1;xx=(delta[:,1]-m2*delta[:,0])/(m1-m2);yy=m1*xx
 assert abs(result['relative_range_p95']-abs(np.hypot(xx[0],yy[0])-r)/r)<1e-12

def front(t5,t10):
 return [dict(family=f,anchor=a,lambda_value=lam,full_range_T5_pass=(lam<=t5 if f==FAMILIES[0] else False),full_range_T10_pass=(lam<=t10),worst_joint_P95=.03) for f in FAMILIES for a in ['A','B'] for lam in P['lambdas']]
def test_primary_never_replaced_by_quarter():
 d=decision(front(.25,.5),P)
 assert d['primary_decision']=='AUX1_HALF_T5_JOINT_BUDGET_NOT_ESTABLISHED'
 assert d['Gate2B_decision']=='AUX1_JOINT_SUPPORTS_10PCT_NOT_5PCT'
def test_half_pass_still_stops_unknown_hardware():
 d=decision(front(.5,.5),P)
 assert d['Gate2B_decision']=='AUX1_JOINT_STATIC_REQUIREMENT_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED'
 assert d['actual_hardware_capability']=='UNKNOWN' and d['automatic_further_numerical_work']=='STOP' and d['R4_progress_percent']==0
