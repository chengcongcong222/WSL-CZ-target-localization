import inspect
from types import SimpleNamespace
import numpy as np
import r4_speed_information as m

def saved():
 panel=m.rows(m.OLD/'A1_NEW_TRUTH_PANEL.csv');a=m.load(m.OLD/'A1_NEW_DESIGN_FREEZE.json')['anchors']['A'];d=np.load(m.OLD/'A1_NEW_DRAWS_A.npz')
 t,n,b,_=m.saved_scene(panel[0],a,d['uniforms'][0],d['normals'][0],-1)
 return t,n,b,a,m.truth_state(panel[0])

def test_L2_only_loss_change_and_same_initializer(monkeypatch):
 t,n,b,a,_=saved();initial,_=m.initialize(t,n,b);seen={}
 def stub(fun,q,**kw):
  seen.update(kw);assert fun is m.residual;np.testing.assert_array_equal(q,initial/m.SCALE)
  return SimpleNamespace(x=q,success=True,status=1,nfev=1,cost=0.,optimality=0.)
 monkeypatch.setattr(m,'least_squares',stub);m.estimate_gaussian(t,n,b,a['sigma_deg'])
 assert seen['loss']=='linear' and seen['f_scale']==1.5 and seen['method']=='trf' and seen['max_nfev']==100 and seen['jac'] is m.jacobian
 assert all(seen[k]==1e-10 for k in ['ftol','xtol','gtol']) and 'bounds' not in seen

def test_observation_only_API():
 assert list(inspect.signature(m.estimate_gaussian).parameters)==['t','node_positions','bearings','sigma_deg']

def test_analytic_measurement_Jacobian():
 t,n,b,a,s=saved();q=s/m.SCALE;h=1e-6
 fd=np.column_stack([(m.residual(q+np.eye(4)[j]*h,t,n,b,a['sigma_deg'])-m.residual(q-np.eye(4)[j]*h,t,n,b,a['sigma_deg']))/(2*h) for j in range(4)])
 np.testing.assert_allclose(m.jacobian(q,t,n,b,a['sigma_deg']),fd,rtol=1e-5,atol=1e-5)

def test_physical_F_inverse_and_delta_method():
 t,n,b,a,s=saved();stats,F,C=m.information(t,n,a['sigma_deg'],s)
 assert stats['rank_scaled_H']==4
 np.testing.assert_allclose(F@C,np.eye(4),atol=1e-6)
 g=np.r_[0.,0.,s[2:]/np.linalg.norm(s[2:])]
 np.testing.assert_allclose(stats['speed_variance_mps2'],g@np.linalg.inv(F)@g,rtol=1e-7)

def test_frozen_prefix_partition_additivity():
 t,n,b,a,s=saved();ms=m.masks(t);assert [int(x.sum()) for x in ms]==[31,61,91,121,60]
 _,Fs,_=m.information(t[ms[1]],n[ms[1]],a['sigma_deg'],s)
 _,Fp,_=m.information(t[ms[4]],n[ms[4]],a['sigma_deg'],s)
 _,Ff,_=m.information(t,n,a['sigma_deg'],s)
 np.testing.assert_allclose(Fs+Fp,Ff,rtol=1e-12,atol=1e-12)
