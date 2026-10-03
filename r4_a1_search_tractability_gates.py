"""Frozen discrete threshold-hit and witness-cluster comparisons."""
from __future__ import annotations
import numpy as np
import pandas as pd
import r4_a1_search_tractability as t

KIND='EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM'
THRESHOLD=.001
STATE_COLUMNS=list(t.L.AXES)
WITNESS_COLUMNS=['case_id','solver','budget','depth_label_m','candidate_kind','J_exact',*STATE_COLUMNS,'bearing_cost']


def state_delta(left,right):
    delta=np.abs(np.asarray(left,float)-np.asarray(right,float))
    delta[...,1]=np.abs((np.asarray(left)[...,1]-np.asarray(right)[...,1]+180)%360-180)
    delta[...,3]=np.abs((np.asarray(left)[...,3]-np.asarray(right)[...,3]+180)%360-180)
    return delta


def state_agrees(left,right,tolerance=t.STATE_AGREEMENT):
    a=np.asarray(left,float);b=np.asarray(right,float)
    return bool(a.shape==(4,) and b.shape==(4,) and np.isfinite(a).all() and np.isfinite(b).all() and np.all(state_delta(a,b)<=tolerance))


def dedup_witnesses(frame):
    """Greedy exact-score-first representatives, separately in each depth."""
    if frame.empty:return frame.copy()
    if not (frame.candidate_kind==KIND).all():raise ValueError('Unrecognized witness provenance')
    ordered=frame.sort_values(['J_exact','depth_label_m',*STATE_COLUMNS],kind='stable').reset_index(drop=True)
    states=ordered[STATE_COLUMNS].to_numpy(float);labels=ordered.depth_label_m.to_numpy(float)
    if not np.isfinite(states).all() or not np.isfinite(ordered.J_exact).all():raise ValueError('Nonfinite witness')
    removed=np.zeros(len(ordered),bool);keep=[]
    for i in range(len(ordered)):
        if removed[i]:continue
        keep.append(i)
        same=(labels==labels[i]) & np.all(state_delta(states,states[i])<=t.DEDUP,axis=1)
        removed[same]=True
    return ordered.iloc[keep].reset_index(drop=True)


def subthreshold_clusters(frame):
    return dedup_witnesses(frame[frame.J_exact<THRESHOLD])


def clusters_contained(left,right):
    """Every left representative has a same-depth componentwise right match."""
    a=subthreshold_clusters(left);b=subthreshold_clusters(right)
    states=b[STATE_COLUMNS].to_numpy(float)
    labels=b.depth_label_m.to_numpy(float)
    for row in a.itertuples():
        state=np.array([getattr(row,col) for col in STATE_COLUMNS])
        matches=(labels==row.depth_label_m) & np.all(state_delta(states,state)<=t.DEDUP,axis=1)
        if not matches.any():return False
    return True


def result_valid(row):
    return bool(row['execution_valid'] and row['all_branches_execution_valid'] and row['n_depth_branches_completed']==21 and row['n_valid_branches']==21)


def result_recovered(row):
    j=float(row['best_exact_J'])
    return bool(np.isfinite(j) and j<THRESHOLD and row['recovered'])


def top_agrees(left,right):
    return bool(left['z_star_label_m']==right['z_star_label_m'] and state_agrees([left[c] for c in STATE_COLUMNS],[right[c] for c in STATE_COLUMNS]))


def raw_case_convergence(t2,t3,catalog_t2,catalog_t3):
    flags=dict(T2_execution_valid=result_valid(t2),T3_execution_valid=result_valid(t3),
        T2_recovered=result_recovered(t2),T3_recovered=result_recovered(t3),
        top_state_depth_agreement=top_agrees(t2,t3),T3_clusters_contained_in_T2=clusters_contained(catalog_t3,catalog_t2))
    return dict(**flags,raw_T2_T3_converged=all(flags.values()))


def dual_case_agreement(shgo,direct,shgo_catalog,direct_catalog):
    flags=dict(SHGO_execution_valid=result_valid(shgo),DIRECT_execution_valid=result_valid(direct),
        SHGO_recovered=result_recovered(shgo),DIRECT_recovered=result_recovered(direct),
        top_state_depth_agreement=top_agrees(shgo,direct),
        SHGO_clusters_contained_in_DIRECT=clusters_contained(shgo_catalog,direct_catalog),
        DIRECT_clusters_contained_in_SHGO=clusters_contained(direct_catalog,shgo_catalog))
    return dict(**flags,dual_solver_subthreshold_witness_cluster_agreement=all(flags.values()))


def development_gates(results,catalogs,case_ids):
    ids=list(case_ids)
    if len(ids)!=21 or len(set(ids))!=21:raise ValueError('Exactly 21 unique development identifiers required')
    comparisons=[];indexed={};all_execution=True
    for solver in ['SHGO','DIRECT']:
        frame=results[solver]
        if len(frame)!=63 or frame.duplicated(['case_id','budget']).any():raise ValueError('Incomplete raw results')
        if set(frame.solver)!={solver} or set(frame.case_id)!=set(ids) or set(frame.budget)!={'T1','T2','T3'}:raise ValueError('Unexpected raw identifiers')
        indexed[solver]=frame.set_index(['case_id','budget'])
        all_execution &= all(result_valid(row) for row in frame.to_dict('records'))
        for cid in ids:
            cat=catalogs[solver]
            c2=cat[(cat.case_id==cid)&(cat.budget=='T2')];c3=cat[(cat.case_id==cid)&(cat.budget=='T3')]
            flags=raw_case_convergence(indexed[solver].loc[(cid,'T2')],indexed[solver].loc[(cid,'T3')],c2,c3)
            comparisons.append(dict(case_id=cid,solver=solver,**flags))
    convergence=pd.DataFrame(comparisons)
    raw=bool(all_execution and convergence.raw_T2_T3_converged.all())
    dual=[]
    if raw:
        for cid in ids:
            cats=[catalogs[s][(catalogs[s].case_id==cid)&(catalogs[s].budget=='T3')] for s in ['SHGO','DIRECT']]
            flags=dual_case_agreement(indexed['SHGO'].loc[(cid,'T3')],indexed['DIRECT'].loc[(cid,'T3')],*cats)
            dual.append(dict(case_id=cid,**flags))
    dual_frame=pd.DataFrame(dual)
    agreement=bool(raw and len(dual_frame)==21 and dual_frame.dual_solver_subthreshold_witness_cluster_agreement.all())
    decision=dict(DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED=raw,DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED=agreement,
        dual_status='EVALUATED' if raw else 'NOT_REACHED_DUE_TO_DEVELOPMENT_GATE',all_execution_valid=bool(all_execution),
        scientific_comparison='threshold-hit / discrete witness-cluster convergence; no stationarity certificate',
        decision='DEVELOPMENT_GATES_PASS' if agreement else 'A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED')
    return decision,convergence,dual_frame


def exact_alias_pairs(catalog,reconstruct):
    """Strict independently reconstructed joint matches; no uniqueness claim."""
    verified=[]
    for row in catalog.itertuples():
        if row.J_exact>=1e-6:continue
        state=np.array([getattr(row,c) for c in STATE_COLUMNS])
        bearing_rms_deg,acoustic_rms_db=reconstruct(state,row.depth_label_m,row.case_id)
        if not np.isfinite([bearing_rms_deg,acoustic_rms_db]).all():raise ValueError('Invalid alias reconstruction')
        if bearing_rms_deg<1e-8 and acoustic_rms_db<1e-6:verified.append((row,state,bearing_rms_deg,acoustic_rms_db))
    pairs=[]
    for i,(left,ls,lb,lj) in enumerate(verified):
        for right,rs,rb,rj in verified[i+1:]:
            if left.case_id==right.case_id and not state_agrees(ls,rs,t.DEDUP):
                pairs.append(dict(case_id=left.case_id,left_solver=left.solver,right_solver=right.solver,
                    left_depth_label_m=left.depth_label_m,right_depth_label_m=right.depth_label_m,
                    left_bearing_rms_deg=lb,right_bearing_rms_deg=rb,left_acoustic_rms_db=lj,right_acoustic_rms_db=rj,
                    decision='A1_CONTINUOUS_ALIAS_FINDING_STOP'))
    return pd.DataFrame(pairs)
