"""Frozen B1 score adapter. Loads only pure accepted R3 functions via AST."""
from pathlib import Path
import ast
import math
import numpy as np
ROOT=Path(__file__).resolve().parent
HISTORY=ROOT/'r3_closedloop_1d_fix.py'
FUNCTIONS=('platform_states_turn','target_states','range_traj','parse_mod','F_matrix','L_profile','demean')

def legacy_functions():
    tree=ast.parse(HISTORY.read_text(encoding='utf-8-sig'))
    nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in FUNCTIONS]
    assert {n.name for n in nodes}==set(FUNCTIONS)
    namespace={'np':np,'math':math,'T_TURN':600.,'U_PLAT':2.,'ZR':200.}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(HISTORY),'exec'),namespace)
    return namespace

def load_models():
    h=legacy_functions()
    return {f:h['parse_mod'](ROOT/'results/R3_C2_Yang_SA_depth/R3_C2_1/_kraken_zgrid'/f'zgrid_f{f}.mod') for f in (201,235,283,338)}

def raw_levels(horizontal,models,depth_labels):
    h=legacy_functions()
    ranges=h['range_traj'](np.arange(0.,1200.1,10.),horizontal[0]*1000.,*horizontal[1:],15.)
    values={f:h['L_profile'](mod,h['F_matrix'](mod,ranges),depth_labels) for f,mod in models.items()}
    return values

def centered(raw):
    result={}
    for f,levels in raw.items():
        out=np.asarray(levels,float).copy()
        for row in out:
            for s in (slice(0,61),slice(61,121)):
                row[s]-=float(np.mean(row[s]))
        result[f]=out
    return result

def depth_profile(horizontal,observed,models,depth_labels,frequencies):
    """No case ID, true depth, expected argmin, Gate or evaluator input."""
    values=centered(raw_levels(horizontal,models,depth_labels))
    scores=np.zeros(len(depth_labels))
    for f in frequencies:
        for iz in range(len(depth_labels)):
            e=values[f][iz]-observed[f]
            scores[iz]+=np.sum(e[:61]**2)+np.sum(e[61:]**2)
    return np.sqrt(scores/(121*len(frequencies)))

def profile_metrics(depth_labels,scores,reference_label,tolerance):
    """Evaluation only; ascending-label deterministic argmin and tolerance ties."""
    z=np.asarray(depth_labels,float); j=np.asarray(scores,float)
    if len(z)!=21 or len(j)!=21 or not np.isfinite(j).all():
        raise ValueError('Incomplete or invalid depth profile')
    i=int(np.flatnonzero(z==reference_label)[0]); best=int(np.argmin(j)); low=float(j[best])
    second=float(np.min(np.delete(j,i))-j[i])
    neighbors=[k for k in (i-1,i+1) if 0<=k<len(j)]
    neighbor=float(min(j[k] for k in neighbors)-j[i])
    curvature=float(j[i-1]-2*j[i]+j[i+1]) if 0<i<len(j)-1 else None
    span=float(j.max()-low)
    width=None; component=[]
    if span>tolerance:
        mask=(j-low)/span<=0.5
        left=right=best
        while left>0 and mask[left-1]: left-=1
        while right<len(j)-1 and mask[right+1]: right+=1
        component=z[left:right+1].tolist(); width=float(z[right]-z[left])
    return dict(z_hat_grid_m=float(z[best]),true_depth_rank=1+int(np.sum(j<j[i])),
                tolerance_aware_true_depth_rank=1+int(np.sum(j<j[i]-tolerance)),
                min_tie_count=int(np.sum(j<=low+tolerance)),J_true_db=float(j[i]),J_min_db=low,
                DeltaJ_second_db=second,DeltaJ_neighbor_db=neighbor,K_z_db=curvature,
                strict_interior_local_minima=int(np.sum((j[1:-1]<j[:-2])&(j[1:-1]<j[2:]))),
                boundary_local_minima=int(j[0]<j[1])+int(j[-1]<j[-2]),
                dynamic_range_db=span,W50_z_m=width,W50_component_labels=component,
                absolute_depth_error_m=float(abs(z[best]-reference_label)))
