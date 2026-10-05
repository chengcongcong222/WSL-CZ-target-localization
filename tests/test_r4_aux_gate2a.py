import math
import numpy as np
from r4_aux_gate2a import intersection,cap_axis,decision

def test_noiseless_known_nodes():
 xy,rhat,failed,*_=intersection(np.zeros(2),np.array([0.,5.]),np.array([0.]),np.array([-math.atan2(5,60)]))
 np.testing.assert_allclose(xy,[[60,0]],atol=1e-12);assert abs(rhat[0]-60)<1e-12 and not failed[0]

def test_estimated_common_translation_changes_absolute_not_relative_output():
 shift=np.array([.1,.2]);angles=(np.array([0.]),np.array([-math.atan2(5,60)]))
 x,r,*_=intersection(np.zeros(2),np.array([0.,5.]),*angles)
 y,s,*_=intersection(shift,shift+np.array([0.,5.]),*angles)
 np.testing.assert_allclose(y-x,[shift],atol=1e-12);np.testing.assert_allclose(r,s,atol=1e-12)

def test_signed_bias_mirror():
 output=[];noise=np.array([[.001,-.002],[-.003,.002]])
 for side in [-1,1]:
  theta1=side*(noise[:,0]+.002-.001);theta2=side*(-math.atan2(5,60)+noise[:,1]+.002+.001)
  xy,rhat,failed,*_=intersection(np.zeros(2),np.array([0.,side*5]),theta1,theta2);output.append((xy,rhat,failed))
 np.testing.assert_allclose(output[0][0][:,0],output[1][0][:,0],atol=1e-12);np.testing.assert_array_equal(output[0][1],output[1][1])

def test_parallel_and_behind_count_as_failures():
 _,_,failed,parallel,*_=intersection(np.zeros(2),np.array([0.,1.]),np.array([0.]),np.array([0.]))
 assert failed[0] and parallel[0]
 _,_,failed,_,behind,*_=intersection(np.zeros(2),np.array([0.,1.]),np.array([0.]),np.array([math.pi/4]))
 assert failed[0] and behind[0]

def test_axis_caps_cannot_hide_bad_sign_or_inner_node():
 rows=[{'bias':x,'full_range_5pct_pass':x not in [-.02]} for x in [-.1,-.02,-.01,0,.01,.02,.1]]
 assert cap_axis(rows,'bias','5pct')==.01

def test_nonzero_budgets_with_unknown_hardware_do_not_open_A1():
 p={'stage':'TEST','anchors':{'A':{}},'hardware_all_relevant_requirements_confirmed':False}
 b=[{'anchor':'A','max_tested_abs_common_5pct_deg':.1,'max_tested_abs_diff_5pct_deg':.02}]
 n=[{'anchor':'A','max_tested_sigma_pos_5pct_m':50}];g=[{'anchor':'A','max_tested_abs_beta_5pct_deg':20}]
 d=decision(b,n,g,p)
 assert d['Gate2A_decision']=='AUX1_NONIDEAL_REQUIREMENTS_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED'
 assert d['R4_A1_NEW']=='NOT_OPENED' and d['AUX_GATE2B']=='NOT_OPENED' and d['R4_progress_percent']==0
