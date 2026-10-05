import inspect
from types import SimpleNamespace
import numpy as np
import pytest
import r4_speed_bias_estimator as candidate
from r4_a1_new_dynamic import rows,load,scene,OUT as OLD

def saved_fixture():
 p=load(OLD/'A1_NEW_DESIGN_FREEZE.json');truth=rows(OLD/'A1_NEW_TRUTH_PANEL.csv')[0];draw=np.load(OLD/'A1_NEW_DRAWS_A.npz');t,n,b,e=scene(truth,p['anchors']['A'],draw['uniforms'][0],draw['normals'][0],-1);return t,n,b,p['anchors']['A']
def test_bias_jacobian_finite_difference():
 t,n,b,a=saved_fixture();x,_=candidate.initialize(t,n,b);q=np.r_[x/candidate.SCALE4,.1,-.2];args=(t,n,b,a['sigma_deg'],a['common_deg'],a['half_diff_deg']);h=1e-6
 finite=np.column_stack([(candidate.residual(q+np.eye(6)[j]*h,*args)-candidate.residual(q-np.eye(6)[j]*h,*args))/(2*h) for j in range(6)])
 assert np.allclose(candidate.jacobian(q,*args),finite,rtol=1e-6,atol=1e-6)
def test_truth_argument_rejection():
 assert list(inspect.signature(candidate.estimate_bias_aware).parameters)==['times','node_positions','bearings','sigma_deg','common_bound_deg','diff_bound_deg']
 t,n,b,a=saved_fixture()
 with pytest.raises(TypeError):candidate.estimate_bias_aware(t,n,b,.05,.05,.025,truth_bias=.01)
def test_zero_bias_initialization_and_no_regularization(monkeypatch):
 captured={}
 def fake(fun,q0,**kw):
  captured.update(q0=q0,kw=kw);return SimpleNamespace(x=q0,success=True,status=1,nfev=1,cost=0.,optimality=0.)
 monkeypatch.setattr(candidate,'least_squares',fake);t,n,b,a=saved_fixture();candidate.estimate_bias_aware(t,n,b,.05,.05,.025)
 assert np.array_equal(captured['q0'][4:],[0.,0.])
 assert captured['kw']['max_nfev']==300 and captured['kw']['f_scale']==1.5
 assert np.array_equal(captured['kw']['bounds'][0][4:],[-1,-1])
 assert np.array_equal(captured['kw']['bounds'][1][4:],[1,1])
 assert captured['kw']['args'][0].shape==(121,)
def test_segment_definition():
 t,n,b,a=saved_fixture();x,_=candidate.initialize(t,n,b);j=candidate.structure(np.r_[x/candidate.SCALE4,0,0],t,n,b,.05,.05,.025)
 assert [x['n_epochs'] for x in j]==[61,60,121]
 assert all(len([k for k in x if k.startswith('singular_value_')])==6 for x in j)
