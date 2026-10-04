import pytest
from r4_b1b_diagnostic import select_record,curve_metrics,route_decision

def rec(value,r=0.,v=0.):return {'J_z_db':value,'delta_r_km':r,'delta_v_mps':v}

def test_empirical_q25_selects_existing_value():
    values=[rec(i,float(i),0.) for i in range(15)]
    chosen,n,index=select_record(values,'Q25')
    assert chosen is values[3] and n==15 and index==4

def test_median_odd_support():
    values=[rec(i,float(i),0.) for i in range(25)]
    chosen,n,index=select_record(values,'MEDIAN')
    assert chosen is values[12] and n==25 and index==13

def test_ties_deterministic_horizontal_order():
    values=[rec(1,.25,0),rec(1,-.25,.05),rec(1,-.25,-.05)]
    assert select_record(values,'MIN')[0] is values[2]

def test_leave_center_removes_only_center():
    values=[rec(0),rec(1,0,.025),rec(2,.125,0)]
    assert select_record(values,'MIN')[0] is values[0]
    chosen,n,index=select_record(values,'MIN_NO_CENTER')
    assert chosen is values[1] and n==2 and index==1

def test_no_center_empty_support_invalid():
    with pytest.raises(ValueError):select_record([rec(0)],'MIN_NO_CENTER')

def test_unregistered_quantile_invalid():
    with pytest.raises(ValueError):select_record([rec(0)],'Q10')

def test_unique_depth_metrics_and_width():
    labels=list(range(150,251,5));scores=[abs(z-200)/5 for z in labels]
    m=curve_metrics(labels,scores,200,1e-10)
    assert m['z_hat_grid_m']==200 and m['true_depth_rank']==1 and m['DeltaJ_second_db']==1 and m['K_z_db']==2 and m['W50_z_m']==50

def test_flat_curve_ties_no_confidence_width():
    m=curve_metrics(list(range(150,251,5)),[1]*21,200,1e-10)
    assert m['min_tie_count']==21 and m['W50_z_m'] is None and m['strict_local_minima_count']==0

def test_wrong_depth_margin_negative():
    labels=list(range(150,251,5));m=curve_metrics(labels,[abs(z-180) for z in labels],200,1e-10)
    assert m['DeltaJ_second_db']<0 and m['absolute_depth_error_m']==20

def test_primary_min_center_success_cannot_alone_admit():
    points=[]
    for a in ['MIN','Q25']:
        for i in range(6):points.append({'aggregation':a,'frequency_configuration':'S7' if i<3 else 'S8','absolute_depth_error_m':0 if a=='MIN' else 30,'true_depth_rank':1})
    assert route_decision(points)['route_decision']=='JOINT_RVZ_PROFILE_ROUTE_NOT_SUPPORTED_BY_EXISTING_EVIDENCE'

def test_fixed_five_of_six_q25_admission():
    points=[]
    for a in ['MIN','Q25']:
        for i in range(6):points.append({'aggregation':a,'frequency_configuration':'S7' if i<3 else 'S8','absolute_depth_error_m':30 if a=='Q25' and i==0 else 0,'true_depth_rank':1})
    assert route_decision(points)['route_decision']=='JOINT_RVZ_PROFILE_ROUTE_WORTH_FRESH_DESIGN'

def test_all_six_min_rank_required():
    points=[]
    for a in ['MIN','Q25']:
        for i in range(6):points.append({'aggregation':a,'frequency_configuration':'S7' if i<3 else 'S8','absolute_depth_error_m':0,'true_depth_rank':4 if a=='MIN' and i==0 else 1})
    assert route_decision(points)['route_decision']=='JOINT_RVZ_PROFILE_ROUTE_NOT_SUPPORTED_BY_EXISTING_EVIDENCE'
