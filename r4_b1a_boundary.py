"""Frozen B1A aggregation only; no propagation or parameter tuning."""
import itertools

def practical(point):
    return point['absolute_depth_error_m']<=10 and point['true_depth_rank']<=3

def strict(point,tolerance):
    return point['z_hat_grid_m']==point['z_true_m'] and point['true_depth_rank']==1 and point['DeltaJ_second_db']>tolerance

def summarize(points):
    if not points: raise ValueError('Empty point set')
    n=len(points)
    return dict(n_points=n,n_exact=sum(p['absolute_depth_error_m']==0 for p in points),
                n_within_5m=sum(p['absolute_depth_error_m']<=5 for p in points),
                n_within_10m=sum(p['absolute_depth_error_m']<=10 for p in points),
                n_rank_le3=sum(p['true_depth_rank']<=3 for p in points),
                n_practical_pass=sum(p['practical_two_bin_pass'] for p in points),
                n_strict_pass=sum(p['strict_exact_pass'] for p in points),
                n_failed=sum(not p['practical_two_bin_pass'] for p in points),
                all_points_pass=all(p['practical_two_bin_pass'] for p in points),
                worst_depth_error_m=max(p['absolute_depth_error_m'] for p in points),
                worst_true_depth_rank=max(p['true_depth_rank'] for p in points),
                min_DeltaJ_second_db=min(p['DeltaJ_second_db'] for p in points),
                min_DeltaJ_neighbor_db=min(p['DeltaJ_neighbor_db'] for p in points))

def rectangles(points,range_limits,speed_limits):
    rows=[]
    for r,v in itertools.product(range_limits,speed_limits):
        group=[p for p in points if abs(p['delta_r_km'])<=r and abs(p['delta_v_mps'])<=v]
        rows.append({'R_km':r,'V_mps':v,**summarize(group)})
    passing=[p for p in rows if p['all_points_pass']]
    for p in rows:
        p['pareto_maximal_robust_rectangle']=p['all_points_pass'] and not any(q['R_km']>=p['R_km'] and q['V_mps']>=p['V_mps'] and (q['R_km']>p['R_km'] or q['V_mps']>p['V_mps']) for q in passing)
        p['is_primary_target']=p['R_km']==.25 and p['V_mps']==.05
    return rows

def nonmonotonic_witnesses(points):
    """All same-case, same closed orthant coordinatewise smaller-fail/larger-pass pairs."""
    rows=[]
    def compatible(a,b): return a==0 or b==0 or (a>0)==(b>0)
    for small in points:
        if small['practical_two_bin_pass']: continue
        for large in points:
            if small['case_id']!=large['case_id'] or not large['practical_two_bin_pass']: continue
            sr,sv=small['delta_r_km'],small['delta_v_mps']; lr,lv=large['delta_r_km'],large['delta_v_mps']
            if not (compatible(sr,lr) and compatible(sv,lv) and abs(sr)<=abs(lr) and abs(sv)<=abs(lv) and (abs(sr)<abs(lr) or abs(sv)<abs(lv))): continue
            rows.append({'case_id':small['case_id'],'smaller_dr_km':sr,'smaller_dv_mps':sv,'smaller_error_m':small['absolute_depth_error_m'],'smaller_true_rank':small['true_depth_rank'],
                         'larger_dr_km':lr,'larger_dv_mps':lv,'larger_error_m':large['absolute_depth_error_m'],'larger_true_rank':large['true_depth_rank'],
                         'fixed_other_axis':sr==lr or sv==lv})
    return rows

def scientific_decision(target_pass):
    return 'B1_CONDITIONAL_DEPTH_IDENTIFIABILITY_ESTABLISHED_WITH_RANGE_SPEED_ENVELOPE' if target_pass else 'B1_HORIZONTAL_CONDITIONING_ENVELOPE_BELOW_TARGET'
