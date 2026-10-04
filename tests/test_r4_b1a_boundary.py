import copy
import inspect
import itertools
import pytest
from r4_b1a_boundary import practical,strict,summarize,rectangles,nonmonotonic_witnesses,scientific_decision
from r4_b1_depth_profile import depth_profile
R=[-1,-.5,-.25,-.125,0,.125,.25,.5,1]
V=[-.2,-.1,-.05,-.025,0,.025,.05,.1,.2]

def point(cid='A',r=0,v=0,error=0,rank=1):
    p={'case_id':cid,'delta_r_km':r,'delta_v_mps':v,'absolute_depth_error_m':error,'true_depth_rank':rank,'z_hat_grid_m':200+error,'z_true_m':200,'DeltaJ_second_db':1,'DeltaJ_neighbor_db':1,'strict_exact_pass':error==0 and rank==1}
    p['practical_two_bin_pass']=practical(p)
    return p

def lattice(): return [point(str(c),r,v) for c,r,v in itertools.product(range(6),R,V)]

@pytest.mark.parametrize('error,rank,expected',[(10,3,True),(15,1,False),(0,4,False),(5,2,True)])
def test_primary_requires_both(error,rank,expected): assert practical(point(error=error,rank=rank))==expected

def test_strict_is_separate():
    p=point(error=5,rank=2)
    assert practical(p) and not strict(p,1e-10)
    p=point();p['DeltaJ_second_db']=1e-12
    assert practical(p) and not strict(p,1e-10)

def test_target_count_and_all_16():
    rows=rectangles(lattice(),[.125,.25,.5,1],[.025,.05,.1,.2])
    assert len(rows)==16 and next(r for r in rows if r['is_primary_target'])['n_points']==150
    assert sum(r['pareto_maximal_robust_rectangle'] for r in rows)==1
    assert rows[-1]['n_points']==486

def test_signed_interior_failure_cannot_average_away():
    points=lattice(); index=next(i for i,p in enumerate(points) if p['case_id']=='0' and p['delta_r_km']==-.125 and p['delta_v_mps']==.025)
    points[index]=point('0',-.125,.025,15,1)
    target=next(r for r in rectangles(points,[.125,.25,.5,1],[.025,.05,.1,.2]) if r['is_primary_target'])
    assert target['n_failed']==1 and target['n_practical_pass']==149 and not target['all_points_pass']
    assert scientific_decision(target['all_points_pass'])=='B1_HORIZONTAL_CONDITIONING_ENVELOPE_BELOW_TARGET'

def test_smaller_rectangle_never_rescues_target():
    points=lattice();index=next(i for i,p in enumerate(points) if p['case_id']=='0' and p['delta_r_km']==.25 and p['delta_v_mps']==.05)
    points[index]=point('0',.25,.05,20,5)
    rows=rectangles(points,[.125,.25,.5,1],[.025,.05,.1,.2])
    assert any(r['all_points_pass'] for r in rows)
    assert not next(r for r in rows if r['is_primary_target'])['all_points_pass']

def test_pareto_tradeoff_preserved():
    points=lattice()
    for p in points:
        if abs(p['delta_r_km'])>.25 and abs(p['delta_v_mps'])>.05:
            p['absolute_depth_error_m']=30;p['practical_two_bin_pass']=False
    rows=rectangles(points,[.125,.25,.5,1],[.025,.05,.1,.2])
    assert {(r['R_km'],r['V_mps']) for r in rows if r['pareto_maximal_robust_rectangle']}=={(.25,.2),(1,.05)}

def test_no_robust_rectangle():
    points=lattice()
    for p in points:
        if p['delta_r_km']!=0 or p['delta_v_mps']!=0: p['practical_two_bin_pass']=False
    assert not any(r['all_points_pass'] for r in rectangles(points,[.125,.25,.5,1],[.025,.05,.1,.2]))

def test_nonmonotonic_same_case_sign_and_axis():
    small=point('A',.125,0,30,5);large=point('A',.25,0);opposite=point('A',-.25,0);other=point('B',.25,0)
    witnesses=nonmonotonic_witnesses([small,large,opposite,other])
    assert len(witnesses)==1 and witnesses[0]['fixed_other_axis']

def test_joint_compensation_witness():
    small=point('A',.125,.025,15,4);large=point('A',.25,.05)
    witnesses=nonmonotonic_witnesses([small,large])
    assert len(witnesses)==1 and not witnesses[0]['fixed_other_axis']

def test_profiler_no_truth_or_lattice_gate_inputs():
    assert set(inspect.signature(depth_profile).parameters)=={'horizontal','observed','models','depth_labels','frequencies'}

def test_full_rectangles_no_mutation():
    points=lattice();old=copy.deepcopy(points)
    rectangles(points,[.125,.25,.5,1],[.025,.05,.1,.2])
    assert points==old
