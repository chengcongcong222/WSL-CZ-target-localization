"""Unreleased development harness. Real entry points stop before loading data."""
from __future__ import annotations
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import time
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
import r4_a1_search_tractability as t
import r4_a1_search_tractability_gates as g

DEVELOPMENT_RELEASED=False
FRESH_CONFIRMATION_RELEASED=False
PARENT='b4839b295777127ec0c8ade56b76c776db99148b'
BUDGET_SHA='40cd9c0b2cbc7a6fced83faf97009c25d652223bf7cbb77e7ebc6a981089672d'
CORE_SHA='0890fd3c5df680a2d56e13d04709dfe727d09d3b44aabfe4b539cab21de92e2b'
MANIFEST_COLUMNS=['case_id','panel_id','sigma_deg','seed','bearing_cutoff','observation_index','origin_group']
OBSERVATION_PATH=Path('results/R4_A1_FIX_CONTINUOUS_SEARCH/CASE_OBSERVATIONS.npz')


def require_release():
    if not DEVELOPMENT_RELEASED:raise SystemExit('DEVELOPMENT_NOT_RELEASED: research-lead audit required')


def verify_infrastructure():
    for path,expected in [(t.ROOT/'r4_a1_search_tractability.py',CORE_SHA),(t.OUT/'TRACTABILITY_BUDGET_FREEZE.json',BUDGET_SHA)]:
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise RuntimeError('FROZEN_INFRASTRUCTURE_CHANGED: stop before development')
    return json.loads((t.OUT/'TRACTABILITY_BUDGET_FREEZE.json').read_text())


def validate_manifest(frame):
    if set(frame.columns)!=set(MANIFEST_COLUMNS):raise ValueError('Manifest column whitelist violated')
    if len(frame)!=21 or frame.case_id.nunique()!=21 or frame.case_id.duplicated().any():raise ValueError('Exactly 21 unique cases required')
    if frame.case_id.isna().any() or not (frame.sigma_deg==.1).all():raise ValueError('Nominal noisy cases only')
    if not np.isfinite(frame[['sigma_deg','seed','bearing_cutoff','observation_index']].to_numpy(float)).all():raise ValueError('Nonfinite manifest')
    if (frame.bearing_cutoff<0).any() or (frame.seed<0).any() or (frame.observation_index<0).any():raise ValueError('Invalid manifest values')
    for col in ['seed','observation_index']:
        if not np.equal(frame[col],np.floor(frame[col])).all():raise ValueError('Noninteger manifest identifier')
    if frame.observation_index.duplicated().any():raise ValueError('Duplicate observation mapping')
    return frame.copy()


@dataclass(frozen=True)
class Observation:
    bearing_rad: np.ndarray
    relative_tl: np.ndarray


def map_observations(manifest,archive):
    frame=validate_manifest(manifest);ids=np.asarray(archive['case_ids'])
    if ids.ndim!=1 or ids.dtype.kind not in 'US' or len(set(ids.tolist()))!=len(ids):raise ValueError('Observation identifiers must be unique strings')
    positions={str(cid):i for i,cid in enumerate(ids.astype(str))}
    bearing=np.asarray(archive['bearing_rad'],float);tl=np.asarray(archive['relative_tl'],float)
    if bearing.shape!=(len(ids),121) or tl.shape!=(len(ids),3,121):raise ValueError('Corrupted observation shapes')
    mapped={}
    for row in frame.itertuples():
        if row.case_id not in positions or positions[row.case_id]!=row.observation_index:raise ValueError('Explicit case-id observation mapping mismatch')
        index=positions[row.case_id];b=bearing[index].copy();a=tl[index].copy()
        if not np.isfinite(b).all() or not np.isfinite(a).all():raise ValueError('Corrupted observation values: '+row.case_id)
        b.setflags(write=False);a.setflags(write=False);mapped[row.case_id]=Observation(b,a)
    return mapped


def load_execution_inputs():
    require_release();budget=verify_infrastructure()
    freeze=json.loads((t.OUT/'DEVELOPMENT_EXECUTION_FREEZE.json').read_text())
    for path,digest in freeze['sha256'].items():
        if hashlib.sha256((t.ROOT/path).read_bytes()).hexdigest()!=digest:raise RuntimeError('EXECUTION_FREEZE_CHANGED')
    manifest=validate_manifest(pd.read_csv(t.OUT/'DEVELOPMENT_CASE_MANIFEST.csv',float_precision='round_trip'))
    with np.load(t.ROOT/OBSERVATION_PATH,allow_pickle=False) as archive:mapped=map_observations(manifest,archive)
    return budget,manifest,mapped


def solver_options(budget,solver,level):
    config=deepcopy(budget[solver.lower()][level])
    cap=int(config.pop('hard_budget_limit'))
    options=config if solver=='SHGO' else config['options']
    return cap,options


def counters_consistent(wrapper,record):
    accepted=[row for row in wrapper.ledger if row['admitted']];blocked=[row for row in wrapper.ledger if not row['admitted']]
    counts=wrapper.counters();unique={(row['depth_label_m'],row['state_sha256']) for row in accepted}
    good=(all(record.get(key)==value for key,value in counts.items()) and
        counts['n_objective_requests']<=counts['hard_budget_limit'] and
        counts['n_objective_attempts']==len(wrapper.ledger)==counts['n_objective_requests']+counts['n_blocked_requests'] and
        len(accepted)==counts['n_objective_requests'] and len(blocked)==counts['n_blocked_requests'] and
        [row['request'] for row in accepted]==list(range(1,len(accepted)+1)) and
        [row['attempt'] for row in wrapper.ledger]==list(range(1,len(wrapper.ledger)+1)) and
        sum(row['exact_forward_calls'] for row in accepted)==counts['n_exact_forward_evaluations'] and
        sum(row['cache_hit'] for row in accepted)==counts['n_cache_hits'] and len(unique)==counts['n_unique_states_evaluated'] and
        all(not row['cache_hit'] or row['exact_forward_calls']==0 for row in accepted) and
        all(row['exact_forward_calls']==0 and not row['cache_hit'] for row in blocked))
    if counts['hard_budget_triggered']:
        good &= bool(counts['n_objective_requests']==counts['hard_budget_limit'] and counts['n_blocked_requests']==1 and
            wrapper.ledger[-1]['admitted']==False and record['termination_reason']=='HARD_BUDGET_EXCEEDED' and
            'ObjectiveBudgetExceeded' in json.loads(record['exception_chain']))
    else:good &= bool(counts['n_blocked_requests']==0 and record['termination_reason']=='SOLVER_RETURN')
    return bool(good and record.get('hard_budget_enforcement_supported',False))


def extract_witnesses(wrapper,case_id,solver,level,depth_label):
    rows=[]
    for e in wrapper.cache.values():
        if not e.feasible or e.J_exact is None:continue
        if e.depth_label_m!=depth_label or not np.isfinite(e.J_exact) or e.J_exact<0:raise ValueError('Invalid evaluated witness')
        if not np.isfinite(e.state).all() or np.any(e.state<t.inherited.prior.LOW) or np.any(e.state>t.inherited.prior.HIGH):raise ValueError('Invalid witness state')
        rows.append(dict(case_id=case_id,solver=solver,budget=level,depth_label_m=e.depth_label_m,candidate_kind=g.KIND,J_exact=e.J_exact,
            bearing_cost=e.bearing_cost,**dict(zip(g.STATE_COLUMNS,map(float,e.state)))))
    return pd.DataFrame(rows,columns=g.WITNESS_COLUMNS)


def rebuild_witnesses(models,observation,cutoff,frame):
    if frame.empty:return frame.copy()
    converted=frame.rename(columns={'depth_label_m':'z_star_label_m'})
    verified=t.verify_catalog(models,observation.bearing_rad,observation.relative_tl,cutoff,converted)
    result=frame.copy();result['J_exact']=verified.recomputed_J_exact.to_numpy();result['bearing_cost']=verified.recomputed_bearing_cost.to_numpy()
    return result.sort_values(['J_exact','depth_label_m',*g.STATE_COLUMNS],kind='stable').reset_index(drop=True)


def run_branch(models,observation,cutoff,case_id,solver,level,depth_label,budget,*,objective_factory=None,runner=None,verifier=None):
    if objective_factory is None:require_release()
    factory=t.ExactObjective if objective_factory is None else objective_factory
    run=t.run_solver if runner is None else runner;verify=rebuild_witnesses if verifier is None else verifier
    cap,options=solver_options(budget,solver,level);started=time.perf_counter()
    wrapper=None;record={};catalog=pd.DataFrame(columns=g.WITNESS_COLUMNS);invalid='';validation_count=0
    try:
        engine=factory(models,observation.bearing_rad,observation.relative_tl,cutoff,depth_label)
        wrapper=t.BudgetObjective(engine,cap)
        _,record=run(solver,wrapper,options)
        if not counters_consistent(wrapper,record):invalid='COUNTER_OR_ENFORCEMENT_INVALID'
        catalog=extract_witnesses(wrapper,case_id,solver,level,depth_label)
        validation_count=len(catalog);catalog=verify(models,observation,cutoff,catalog)
    except Exception as error:
        invalid=type(error).__name__+': '+str(error)
        if not record:record=dict(termination_reason='UNRELATED_EXCEPTION',solver_success=False,solver_message=repr(error),exception_chain=json.dumps(t.exception_chain(error)))
        catalog=pd.DataFrame(columns=g.WITNESS_COLUMNS)
    counts=wrapper.counters() if wrapper is not None else dict(n_objective_requests=0,n_objective_attempts=0,n_blocked_requests=0,n_exact_forward_evaluations=0,n_unique_states_evaluated=0,n_cache_hits=0,hard_budget_limit=cap,hard_budget_triggered=False)
    branch={**record,**counts}
    branch.update(case_id=case_id,solver=solver,budget=level,depth_label_m=depth_label,
        execution_valid=not bool(invalid),execution_invalid_reason=invalid,n_validation_exact_forward_evaluations=validation_count,
        elapsed_seconds=time.perf_counter()-started,best_branch_J_exact=float(catalog.J_exact.min()) if len(catalog) else float('nan'))
    for name in g.STATE_COLUMNS:branch['best_branch_'+name]=float(catalog.iloc[0][name]) if len(catalog) else float('nan')
    branch['best_branch_bearing_cost']=float(catalog.iloc[0].bearing_cost) if len(catalog) else float('nan')
    return branch,catalog,list(wrapper.ledger) if wrapper is not None else []


def aggregate_case(case_id,solver,level,cutoff,branches,witnesses):
    labels=[row['depth_label_m'] for row in branches];nvalid=sum(row['execution_valid'] for row in branches)
    valid=len(branches)==21 and len(set(labels))==21 and set(labels)==set(t.L.PROFILE) and nvalid==21
    ranked=witnesses.sort_values(['J_exact','depth_label_m',*g.STATE_COLUMNS],kind='stable').reset_index(drop=True)
    top=ranked.iloc[0] if len(ranked) else None
    j=float(top.J_exact) if top is not None else float('nan')
    result=dict(case_id=case_id,solver=solver,budget=level,best_exact_J=j,recovered=bool(np.isfinite(j) and j<g.THRESHOLD),
        bearing_cutoff=float(cutoff),n_depth_branches_completed=len(branches),n_valid_branches=nvalid,
        execution_valid=bool(valid),all_branches_execution_valid=bool(valid),
        z_star_label_m=float(top.depth_label_m) if top is not None else float('nan'),
        bearing_cost=float(top.bearing_cost) if top is not None else float('nan'))
    for name in g.STATE_COLUMNS:result[name]=float(top[name]) if top is not None else float('nan')
    for output,field in [('total_objective_requests','n_objective_requests'),('total_exact_forward_evaluations','n_exact_forward_evaluations'),('total_unique_states','n_unique_states_evaluated'),('total_cache_hits','n_cache_hits')]:
        result[output]=sum(row[field] for row in branches)
    return result


def run_case(models,observation,case_id,cutoff,solver,level,budget,**test_dependencies):
    if not test_dependencies:require_release()
    branches=[];catalogs=[];requests={}
    for z in t.L.PROFILE:
        record,catalog,ledger=run_branch(models,observation,cutoff,case_id,solver,level,float(z),budget,**test_dependencies)
        branches.append(record);catalogs.append(catalog);requests[float(z)]=ledger
    witnesses=pd.concat(catalogs,ignore_index=True)
    result=aggregate_case(case_id,solver,level,cutoff,branches,witnesses)
    return result,g.dedup_witnesses(witnesses),pd.DataFrame(branches),requests,witnesses


def run_family(models,manifest,observations,solver,budget):
    require_release();results=[];catalogs=[];logs=[];all_witnesses=[]
    for row in manifest.itertuples():
        for level in ['T1','T2','T3']:
            verify_infrastructure()
            result,catalog,log,requests,raw=run_case(models,observations[row.case_id],row.case_id,row.bearing_cutoff,solver,level,budget)
            results.append(result);catalogs.append(catalog);logs.append(log);all_witnesses.append(raw)
            for depth,ledger in requests.items():
                t.save(f'DEVELOPMENT_REQUEST_LOGS/{solver}/{row.case_id}/{level}_z{int(depth)}.csv',pd.DataFrame(ledger))
            t.save(f'DEVELOPMENT_{solver}_RESULTS.csv',pd.DataFrame(results))
            t.save(f'{solver}_WITNESS_CATALOG.csv',pd.concat(catalogs,ignore_index=True))
            t.save(f'{solver}_EVALUATED_WITNESSES.csv',pd.concat(all_witnesses,ignore_index=True))
            t.save(f'{solver}_BUDGET_LOG.csv',pd.concat(logs,ignore_index=True))
    return pd.DataFrame(results),pd.concat(catalogs,ignore_index=True)


def run_fresh_confirmation(decision,generator):
    if not (decision['DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED'] and decision['DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED']):
        return 'NOT_REACHED_DUE_TO_DEVELOPMENT_GATE'
    if not FRESH_CONFIRMATION_RELEASED:raise SystemExit('FRESH_CONFIRMATION_NOT_RELEASED: separate frozen design and approval required')
    return generator()


def _execute_released_development():
    require_release();budget,manifest,observations=load_execution_inputs()
    models=t.inherited.inputs();results={};catalogs={}
    with threadpool_limits(limits=2):
        for solver in ['SHGO','DIRECT']:results[solver],catalogs[solver]=run_family(models,manifest,observations,solver,budget)
        decision,raw,dual=g.development_gates(results,catalogs,manifest.case_id.tolist())
        t.save('DEVELOPMENT_RAW_WITNESS_CONVERGENCE.csv',raw)
        if not dual.empty:t.save('DUAL_SOLVER_AGREEMENT.csv',dual)
        def reconstruct(state,depth,case_id):
            observation=observations[case_id];b,r=t.L.geometry(state)
            residual=(b[0]-observation.bearing_rad+np.pi)%(2*np.pi)-np.pi
            feature=t.L.direct_features(models,r,[depth])[0,0]
            return float(np.degrees(np.sqrt(np.mean(residual**2)))),float(np.sqrt(np.mean((feature-observation.relative_tl)**2)))
        aliases=g.exact_alias_pairs(pd.concat(list(catalogs.values()),ignore_index=True).query("budget == 'T3'"),reconstruct)
        t.save('ALIAS_DIAGNOSTIC.csv',aliases)
        if len(aliases):decision['decision']='A1_CONTINUOUS_ALIAS_FINDING_STOP'
        fresh='NOT_REACHED_DUE_TO_ALIAS_STOP' if len(aliases) else ('NOT_RELEASED' if decision['DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED'] else 'NOT_REACHED_DUE_TO_DEVELOPMENT_GATE')
        decision.update(noisy_development_runs=21,n_raw_case_runs=126,R4_progress_percent=0,fresh_confirmation=fresh,
            NO_EXACT_ALIAS_FINDING_IN_TESTED_CATALOGS=not bool(len(aliases)))
        t.json_write('DEVELOPMENT_EXECUTION_DECISION.json',decision)
    return decision


def execute_development():
    require_release()
    try:return _execute_released_development()
    except Exception as error:
        t.json_write('DEVELOPMENT_EXECUTION_INVALID.json',dict(execution_valid=False,
            termination_reason=type(error).__name__,message=str(error),exception_chain=t.exception_chain(error),
            DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED=False,DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED=False,
            fresh_confirmation='NOT_REACHED_DUE_TO_EXECUTION_INVALID',R4_progress_percent=0,
            run_count='Use retained branch/request ledgers; an aborted run is never labeled zero'))
        raise


if __name__=='__main__':execute_development()
