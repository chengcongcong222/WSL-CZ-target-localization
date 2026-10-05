import inspect
import numpy as np
import pytest
from r4_a1_new_estimator import estimate,residual,jacobian,SCALE

def fixture():
 t=np.arange(121)*10.;post=np.maximum(t-600,0);main=np.column_stack([2*np.minimum(t,600)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)]);nodes=np.stack([main,main+[0,5000]],axis=1);x=np.array([53000*np.cos(.01),53000*np.sin(.01),1.8*np.cos(.12),1.8*np.sin(.12)]);d=(x[:2]+t[:,None]*x[2:])[:,None,:]-nodes;b=np.arctan2(d[:,:,1],d[:,:,0]);return t,nodes,b,x

def test_noiseless_offgrid_acquisition():
 t,n,b,x=fixture();est,diag=estimate(t,n,b,.05)
 assert np.allclose(est,x,rtol=0,atol=1e-5)
 assert diag['valid_intersections']==121 and diag['solver_success']
def test_truth_arguments_rejected():
 t,n,b,x=fixture()
 assert list(inspect.signature(estimate).parameters)==['times','node_positions','bearings','sigma_deg']
 with pytest.raises(TypeError):estimate(t,n,b,.05,truth_state=x)
def test_analytic_jacobian():
 t,n,b,x=fixture();q=x/SCALE;h=1e-6;fd=np.column_stack([(residual(q+np.eye(4)[j]*h,t,n,b,.05)-residual(q-np.eye(4)[j]*h,t,n,b,.05))/(2*h) for j in range(4)])
 assert np.allclose(jacobian(q,t,n,b,.05),fd,rtol=1e-6,atol=1e-6)
def test_invalid_intersections_fail_without_truth_filter():
 t,n,b,x=fixture();b[:]=0
 with pytest.raises(ValueError):estimate(t,n,b,.05)
def test_navigation_observations_actually_used():
 t,n,b,x=fixture();shift=np.array([50.,-25.]);est,_=estimate(t,n+shift,b,.05)
 assert np.allclose(est[:2],x[:2]+shift,atol=1e-5)
 assert np.allclose(est[2:],x[2:],atol=1e-7)
