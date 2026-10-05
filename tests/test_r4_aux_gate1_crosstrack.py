import math
import numpy as np
from r4_aux_gate0_geometry import triangulate
from r4_aux_gate1_crosstrack import quantiles,summarize,frontiers,decide

def test_noiseless_crosstrack_intersection():
 for side in [-1,1]:
  xy,failed,*_=triangulate(np.array([0.,side*5]),np.array([0.]),np.array([-side*math.atan2(5,60)]))
  np.testing.assert_allclose(xy,[[60.,0.]],atol=2e-12);assert not failed[0]

def test_mirror_error_invariance():
 noise=np.array([[.001,-.002],[-.003,.002],[.002,.004]])
 output=[]
 for side in [-1,1]:
  xy,failed,*_=triangulate(np.array([0.,side*5]),noise[:,0]*side,-side*math.atan2(5,60)+noise[:,1]*side)
  output.append((np.linalg.norm(xy,axis=1),failed))
 np.testing.assert_array_equal(output[0][0],output[1][0]);np.testing.assert_array_equal(output[0][1],output[1][1])

def test_linear_sigma_closed_form_matches_matrix():
 r,b,s=60.,5.,math.radians(.1)
 H=np.array([[0,1/r],[b/(r*r+b*b),r/(r*r+b*b)]])
 C=np.linalg.inv(H.T@H/s**2)
 assert math.isclose(math.sqrt(C[0,0]),s*math.sqrt(r**4+(r*r+b*b)**2)/b,rel_tol=1e-12)

def test_infinite_failures_stay_in_unconditional_quantile():
 q=quantiles(np.array([0.]*94+[np.inf]*6));assert q['p95']==np.inf

def make_results(b,s,ranges,value):
 return [dict(baseline_km=b,sigma_deg=s,range_km=r,side=side,execution_valid=True,relative_range_p95=value,alpha_deg=math.degrees(math.atan2(b,r)),failure_rate=0.,nonlinear_to_linear_p95_ratio=1.) for r in ranges for side in [-1,1]]

def test_missing_one_range_prevents_full_range_pass():
 p=dict(sigma_deg=[.1],baseline_km=[5],range_km=list(range(50,61)))
 x=summarize(make_results(5,.1,list(range(50,60)),.001),p)[0]
 assert not x['all_ranges_execution_valid'] and not x['full_range_5pct_pass']

def test_worst_range_not_pooled():
 p=dict(sigma_deg=[.1],baseline_km=[5],range_km=list(range(50,61)))
 data=make_results(5,.1,p['range_km'],.01)
 data[-1]['relative_range_p95']=.11
 x=summarize(data,p)[0];assert x['worst_p95_relative']==.11 and not x['full_range_10pct_pass']

def test_minimum_tested_baseline_and_unreached_frontier():
 p=dict(sigma_deg=[.05,.1],baseline_km=[2,3,5],range_km=list(range(50,61)),fixed_baseline_km=[5])
 data=[]
 for s in p['sigma_deg']:
  for b in p['baseline_km']:data+=make_results(b,s,p['range_km'],.04 if s==.05 and b>=3 else .12)
 f,fixed=frontiers(summarize(data,p),p)
 assert f[0]['minimum_tested_B_for_5pct_km']==3
 assert f[1]['minimum_tested_B_for_5pct_km']=='NOT_REACHED_WITHIN_B<=10KM'
 assert fixed[0]['largest_tested_sigma_for_5pct_deg']==.05
 assert fixed[0]['next_tested_failing_sigma_for_5pct_deg']==.1

def test_tradeoff_pass_does_not_open_engineering_or_depth():
 p=dict(stage='TEST',sigma_deg=[.05,.1],baseline_km=[5,6],range_km=list(range(50,61)),fixed_baseline_km=[5])
 data=[]
 for s in p['sigma_deg']:
  for b in p['baseline_km']:data+=make_results(b,s,p['range_km'],.04 if s==.05 or b==6 else .06)
 summary=summarize(data,p);f,fixed=frontiers(summary,p);d=decide(summary,f,fixed,p)
 assert d['primary_architecture_decision']=='AUX1_BASELINE_BEARING_TRADEOFF_ESTABLISHED'
 assert d['R4_A1_NEW']=='NOT_OPENED' and d['R4_progress_percent']==0
 assert d['actual_bearing_hardware_capability']=='UNKNOWN'
