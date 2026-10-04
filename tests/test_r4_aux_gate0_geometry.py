import math
import numpy as np
from r4_aux_gate0_geometry import geometry,triangulate,covariance,quantile

def test_impossible_geometry():
 assert geometry(50,10,90,'NEAR',1) is None
 assert geometry(50,1,2,'FAR',1) is None

def test_both_roots_and_mirrors_have_correct_baseline():
 for branch in ['NEAR','FAR']:
  for side in [-1,1]:
   p,d=geometry(50,10,5,branch,side)
   assert abs(np.linalg.norm(p)-10)<1e-12
   assert abs(np.linalg.norm(np.array([50,0])-p)-d)<1e-12

def test_noiseless_intersection_uses_measurements_only():
 for branch in ['NEAR','FAR']:
  p,d=geometry(50,10,5,branch,1)
  xy,failed,*_=triangulate(p,np.array([0.]),np.array([-math.radians(5)]))
  np.testing.assert_allclose(xy,[[50,0]],atol=1e-12);assert not failed[0]

def test_parallel_is_failure():
 _,failed,parallel,*_=triangulate(np.array([0.,1.]),np.array([0.]),np.array([0.]))
 assert failed[0] and parallel[0]

def test_behind_sensor_is_failure():
 _,failed,_,behind,*_=triangulate(np.array([0.,1.]),np.array([0.]),np.array([math.pi/4]))
 assert failed[0] and behind[0]

def test_fim_jacobian_finite_difference():
 p,d=geometry(50,10,5,'FAR',1);H,F,C=covariance(50,p,.1);q=np.array([50.,0.]);eps=1e-5
 for i in range(2):
  step=np.zeros(2);step[i]=eps
  def bear(point):
   return np.array([math.atan2(point[1],point[0]),math.atan2(point[1]-p[1],point[0]-p[0])])
  np.testing.assert_allclose((bear(q+step)-bear(q-step))/(2*eps),H[:,i],atol=1e-11)
 np.testing.assert_allclose(F@C,np.eye(2),atol=1e-12)

def test_noise_sigma_scaling():
 p,d=geometry(50,10,5,'NEAR',1)
 np.testing.assert_allclose(covariance(50,p,.2)[2],4*covariance(50,p,.1)[2],atol=1e-12)

def test_reflection_covariance_invariance():
 p,d=geometry(50,10,5,'NEAR',1);m,_=geometry(50,10,5,'NEAR',-1)
 np.testing.assert_allclose(np.diag(covariance(50,p,.1)[2]),np.diag(covariance(50,m,.1)[2]),atol=1e-12)

def test_failures_not_dropped_from_percentile():
 assert quantile(np.array([0.]*94+[np.inf]*6),.95)==np.inf

def test_quantile_nearest_rank_no_interpolation():
 assert quantile(np.array([1.,2.,3.,100.]),.95)==100
