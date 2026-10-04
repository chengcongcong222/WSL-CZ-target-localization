"""Saved-score-only B1B statistics; no scientific forward imports."""
import math

def select_record(records,aggregation):
    eligible=[r for r in records if aggregation!='MIN_NO_CENTER' or (r['delta_r_km'],r['delta_v_mps'])!=(0.,0.)]
    if not eligible: raise ValueError('No frozen support nodes')
    ordered=sorted(eligible,key=lambda r:(r['J_z_db'],r['delta_r_km'],r['delta_v_mps']))
    if aggregation in ('MIN','MIN_NO_CENTER'): index=0
    elif aggregation=='Q25': index=math.ceil(len(ordered)/4)-1
    elif aggregation=='MEDIAN': index=math.ceil(len(ordered)/2)-1
    else: raise ValueError('Unregistered aggregation')
    return ordered[index],len(ordered),index+1

def curve_metrics(labels,scores,reference,tolerance):
    """B1 frozen metric definitions, computed only from saved aggregate values."""
    if labels!=list(range(150,251,5)) or len(scores)!=21 or not all(math.isfinite(x) for x in scores): raise ValueError('Incomplete curve')
    i=labels.index(reference); best=min(range(21),key=lambda k:scores[k]); low=scores[best]
    second=min(scores[k] for k in range(21) if k!=i)-scores[i]
    neighbor=min(scores[k] for k in (i-1,i+1) if 0<=k<21)-scores[i]
    curvature=scores[i-1]-2*scores[i]+scores[i+1] if 0<i<20 else None
    span=max(scores)-low; width=None; component=[]
    if span>tolerance:
        mask=[(x-low)/span<=.5 for x in scores]; left=right=best
        while left>0 and mask[left-1]: left-=1
        while right<20 and mask[right+1]: right+=1
        component=labels[left:right+1]; width=labels[right]-labels[left]
    return dict(z_hat_grid_m=labels[best],true_depth_rank=1+sum(x<scores[i] for x in scores),
                tolerance_aware_true_depth_rank=1+sum(x<scores[i]-tolerance for x in scores),
                J_true_db=scores[i],J_min_db=low,min_tie_count=sum(x<=low+tolerance for x in scores),
                DeltaJ_second_db=second,DeltaJ_neighbor_db=neighbor,K_z_db=curvature,
                strict_local_minima_count=sum(scores[k]<scores[k-1] and scores[k]<scores[k+1] for k in range(1,20)),
                boundary_local_minima_count=int(scores[0]<scores[1])+int(scores[-1]<scores[-2]),
                W50_z_m=width,W50_component_labels=component,dynamic_range_db=span,
                absolute_depth_error_m=abs(labels[best]-reference),
                practical_two_bin_pass=abs(labels[best]-reference)<=10 and 1+sum(x<scores[i] for x in scores)<=3)

def summarize(metrics):
    return {'n_cases':len(metrics),'exact_depth':sum(m['absolute_depth_error_m']==0 for m in metrics),
            'within_10m':sum(m['absolute_depth_error_m']<=10 for m in metrics),'rank_le3':sum(m['true_depth_rank']<=3 for m in metrics),
            'practical_pass':sum(m['practical_two_bin_pass'] for m in metrics),'worst_error_m':max(m['absolute_depth_error_m'] for m in metrics),
            'worst_true_rank':max(m['true_depth_rank'] for m in metrics),'minimum_second_margin_db':min(m['DeltaJ_second_db'] for m in metrics)}

def route_decision(target_metrics):
    mins=[m for m in target_metrics if m['aggregation']=='MIN']; q25=[m for m in target_metrics if m['aggregation']=='Q25']
    if len(mins)!=6 or len(q25)!=6: raise ValueError('All six cases are mandatory')
    d1=all(m['true_depth_rank']<=3 for m in mins); d2=all(m['absolute_depth_error_m']<=10 for m in mins)
    d3=sum(m['absolute_depth_error_m']<=10 for m in q25)>=5 and sum(m['absolute_depth_error_m']<=10 for m in q25 if m['frequency_configuration']=='S7')>=2
    label='JOINT_RVZ_PROFILE_ROUTE_WORTH_FRESH_DESIGN' if d1 and d2 and d3 else 'JOINT_RVZ_PROFILE_ROUTE_NOT_SUPPORTED_BY_EXISTING_EVIDENCE'
    return {'D1_MIN_all_six_rank_le3':d1,'D2_MIN_all_six_within10m':d2,'D3_Q25_at_least_five_within10m_and_S7_at_least_two':d3,'route_decision':label}
