import inspect
import numpy as np
import pytest
from r4_b1_depth_profile import profile_metrics,depth_profile,legacy_functions
Z=np.arange(150,251,5)

def test_exact_unique_mechanism():
    j=np.abs(Z-200.)/5
    m=profile_metrics(Z,j,200.,1e-10)
    assert m['z_hat_grid_m']==200 and m['true_depth_rank']==1 and m['min_tie_count']==1
    assert m['DeltaJ_second_db']==1 and m['DeltaJ_neighbor_db']==1 and m['K_z_db']==2
    assert m['W50_z_m']==50 and m['strict_interior_local_minima']==1

def test_wrong_truth_negative_margin_and_rank():
    j=np.abs(Z-180.)/5
    m=profile_metrics(Z,j,200.,1e-10)
    assert m['DeltaJ_second_db']<0 and m['true_depth_rank']>3 and m['absolute_depth_error_m']==20

def test_tie_not_unique():
    j=np.ones(21); j[9:12]=0
    m=profile_metrics(Z,j,200.,1e-10)
    assert m['min_tie_count']==3 and m['strict_interior_local_minima']==0
    assert m['z_hat_grid_m']==195 and m['W50_z_m']==10

def test_near_tie_rejected_for_unique_gate():
    j=np.ones(21); j[10]=0; j[11]=1e-12
    assert profile_metrics(Z,j,200,1e-10)['min_tie_count']==2

def test_flat_width_undefined():
    m=profile_metrics(Z,np.ones(21),200,1e-10)
    assert m['W50_z_m'] is None and m['min_tie_count']==21

def test_disconnected_half_rise_component():
    j=np.ones(21); j[3:6]=.1; j[10]=0
    assert profile_metrics(Z,j,200,1e-10)['W50_z_m']==0

def test_boundary_is_separate():
    m=profile_metrics(Z,np.arange(21),150,1e-10)
    assert m['boundary_local_minima']==1 and m['strict_interior_local_minima']==0 and m['K_z_db'] is None

@pytest.mark.parametrize('bad',[np.zeros(20),np.r_[np.zeros(20),np.nan]])
def test_invalid_profile(bad):
    with pytest.raises(ValueError): profile_metrics(Z,bad,200,1e-10)

def test_profile_api_truth_separation():
    assert list(inspect.signature(depth_profile).parameters)==['horizontal','observed','models','depth_labels','frequencies']
    assert 'reference_label' not in inspect.getsource(depth_profile)

def test_legacy_adapter_no_top_level_writes():
    namespace=legacy_functions()
    assert 'main' not in namespace and 'OUT' not in namespace and 'WORK' not in namespace
