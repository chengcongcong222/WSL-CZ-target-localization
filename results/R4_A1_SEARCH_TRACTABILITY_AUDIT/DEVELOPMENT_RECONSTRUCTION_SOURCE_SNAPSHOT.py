"""Post-run mechanical reconstruction; never invokes an optimizer or changes Gates."""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import time
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
import r4_a1_search_tractability as t
import r4_a1_search_tractability_development as d
import r4_a1_search_tractability_audit as a

OUT=t.OUT
AXES=list(t.L.AXES)
TOL_J=1e-7
TOL_BEARING=1e-12


def read(name):return pd.read_csv(OUT/name,float_precision='round_trip')


def delta(left,right):
    x=np.asarray(left,float);y=np.asarray(right,float);v=np.abs(x-y)
    for i in [1,3]:v[...,i]=np.abs((x[...,i]-y[...,i]+180)%360-180)
    return v


def representatives(frame):
    f=frame.sort_values(['J_exact','depth_label_m',*AXES],kind='stable').reset_index(drop=True)
    states=f[AXES].to_numpy(float);depth=f.depth_label_m.to_numpy(float);keep=[];active=np.ones(len(f),bool)
    for i in range(len(f)):
        if not active[i]:continue
        keep.append(i)
        active[(depth==depth[i])&np.all(delta(states,states[i])<=np.array(t.DEDUP),axis=1)]=False
    return f.iloc[keep].reset_index(drop=True)


def contains(left,right):
    x=representatives(left[left.J_exact<.001]);y=representatives(right[right.J_exact<.001])
    ys=y[AXES].to_numpy(float);yd=y.depth_label_m.to_numpy(float)
    return all(np.any((yd==r.depth_label_m)&np.all(delta(ys,np.array([getattr(r,k) for k in AXES]))<=t.DEDUP,axis=1)) for r in x.itertuples())


def agrees(left,right):
    return bool(left.z_star_label_m==right.z_star_label_m and np.all(delta(left[AXES].to_numpy(float),right[AXES].to_numpy(float))<=t.STATE_AGREEMENT))


def audit():
    checks=[];reconstructed=[];started=time.perf_counter()
    def check(label,ok,detail=''):
        checks.append(dict(check=label,pass_check=bool(ok),detail=str(detail)))
        if not ok:print('FAILED',label,detail,flush=True)
    # Accepted files remain immutable even after actual scientific outputs exist.
    for manifest_name in ['METHOD_FREEZE.json','DEVELOPMENT_EXECUTION_FREEZE.json']:
        for path,h in json.loads((OUT/manifest_name).read_text())['sha256'].items():
            check('Frozen '+path,hashlib.sha256((t.ROOT/path).read_bytes()).hexdigest()==h)
    check('Budget and core SHA',hashlib.sha256((OUT/'TRACTABILITY_BUDGET_FREEZE.json').read_bytes()).hexdigest()==d.BUDGET_SHA and hashlib.sha256((t.ROOT/'r4_a1_search_tractability.py').read_bytes()).hexdigest()==d.CORE_SHA)
    check('Protected historical identities',len(a.frozen_inputs())==2989)
    authorization=json.loads((OUT/'DEVELOPMENT_RELEASE_AUTHORIZATION.json').read_text())
    process=json.loads((OUT/'DEVELOPMENT_PROCESS_METADATA.json').read_text())
    check('Formal release and fresh prohibition',authorization['development_released'] and not authorization['fresh_confirmation_released'] and process['runtime_development_released'] and not process['fresh_confirmation_released'])
    check('Successful complete process',process['returncode']==0)
    manifest=read('DEVELOPMENT_CASE_MANIFEST.csv');d.validate_manifest(manifest);ids=manifest.case_id.tolist()
    with np.load(t.ROOT/d.OBSERVATION_PATH,allow_pickle=False) as archive:observations=d.map_observations(manifest,archive)
    budget=json.loads((OUT/'TRACTABILITY_BUDGET_FREEZE.json').read_text());results={};catalogs={};raws={};logs={};branch_valid={}
    for solver in ['SHGO','DIRECT']:
        results[solver]=read(f'DEVELOPMENT_{solver}_RESULTS.csv');catalogs[solver]=read(f'{solver}_WITNESS_CATALOG.csv');raws[solver]=read(f'{solver}_EVALUATED_WITNESSES.csv');logs[solver]=read(f'{solver}_BUDGET_LOG.csv')
        f=results[solver];log=logs[solver];raw=raws[solver]
        check('Complete raw case-budget set '+solver,len(f)==63 and not f.duplicated(['case_id','budget']).any() and set(f.case_id)==set(ids) and set(f.budget)=={'T1','T2','T3'})
        check('Complete branch set '+solver,len(log)==1323 and not log.duplicated(['case_id','budget','depth_label_m']).any())
        check('Witness provenance '+solver,(raw.candidate_kind=='EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM').all() and (catalogs[solver].candidate_kind=='EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM').all())
        for row in log.itertuples():
            tag=f'{solver}/{row.case_id}/{row.budget}/z{int(row.depth_label_m)}'
            limit=budget[solver.lower()][row.budget]['hard_budget_limit']
            ledger=read(f'DEVELOPMENT_REQUEST_LOGS/{solver}/{row.case_id}/{row.budget}_z{int(row.depth_label_m)}.csv')
            admitted=ledger[ledger.admitted];refused=ledger[~ledger.admitted]
            unique=admitted[['depth_label_m','state_sha256']].drop_duplicates()
            counts=(len(admitted)==row.n_objective_requests<=limit==row.hard_budget_limit and len(ledger)==row.n_objective_attempts and len(refused)==row.n_blocked_requests and len(unique)==row.n_unique_states_evaluated and int(admitted.cache_hit.sum())==row.n_cache_hits and int(ledger.exact_forward_calls.sum())==row.n_exact_forward_evaluations)
            sequence=(admitted.request.tolist()==list(range(1,len(admitted)+1)) and ledger.attempt.tolist()==list(range(1,len(ledger)+1)) and not (admitted.cache_hit&(admitted.exact_forward_calls>0)).any() and not (refused.exact_forward_calls>0).any())
            hard=(row.n_objective_requests==limit and row.n_blocked_requests==1 and not ledger.iloc[-1].admitted and row.termination_reason=='HARD_BUDGET_EXCEEDED' and not row.solver_success and 'ObjectiveBudgetExceeded' in json.loads(row.exception_chain)) if row.hard_budget_triggered else (len(refused)==0 and row.termination_reason=='SOLVER_RETURN')
            check('All counters '+tag,counts);check('Request sequence/cache/refusal '+tag,sequence);check('Hard cap semantics '+tag,hard)
            keys=[]
            for r in admitted.itertuples():
                state=np.array([getattr(r,k) for k in AXES],dtype='<f8')
                keys.append(hashlib.sha256(state.tobytes()).hexdigest()==r.state_sha256)
            check('Roundtrip physical-state keys '+tag,all(keys))
            feasible=admitted[admitted.feasible].drop_duplicates(['depth_label_m','state_sha256'])
            saved=raw[(raw.case_id==row.case_id)&(raw.budget==row.budget)&(raw.depth_label_m==row.depth_label_m)]
            expected_keys={(r.depth_label_m,r.state_sha256) for r in feasible.itertuples()}
            saved_keys={(r.depth_label_m,hashlib.sha256(np.array([getattr(r,k) for k in AXES],dtype='<f8').tobytes()).hexdigest()) for r in saved.itertuples()}
            check('All feasible distinct witnesses retained '+tag,len(saved)==len(saved_keys) and saved_keys==expected_keys and row.n_validation_exact_forward_evaluations==len(saved))
            # Native declaration does not certify validity; ledger consistency and score rebuild do.
            branch_valid[(solver,row.case_id,row.budget,row.depth_label_m)]=bool(counts and sequence and hard and row.execution_valid)
        for cid in ids:
            for level in ['T1','T2','T3']:
                frame=raw[(raw.case_id==cid)&(raw.budget==level)]
                expected=representatives(frame)
                exported=catalogs[solver][(catalogs[solver].case_id==cid)&(catalogs[solver].budget==level)].reset_index(drop=True)
                check('Frozen dedup representatives '+solver+'/'+cid+'/'+level,expected.to_dict('records')==exported.to_dict('records'))
    models=t.inherited.inputs();memo={};all_raw=pd.concat(list(raws.values()),ignore_index=True)
    # Cold independent geometry/RMS for every unique exported case/depth/state.
    # This memo is audit-only, created after searches finish and never supplied to a solver.
    pending=defaultdict(dict)
    for row in all_raw.itertuples():
        state=np.array([getattr(row,k) for k in AXES],dtype='<f8');key=(row.case_id,float(row.depth_label_m),state.tobytes())
        pending[(row.case_id,float(row.depth_label_m))][key]=state
    calls=0;unique_states=0
    with threadpool_limits(limits=2):
        for (cid,depth),entries in pending.items():
            items=list(entries.items());obs=observations[cid]
            for begin in range(0,len(items),16):
                part=items[begin:begin+16];geometry=[a.independent_geometry(state) for _,state in part]
                bearings=np.array([x[0] for x in geometry]);ranges=np.array([x[1] for x in geometry])
                features=t.L.direct_features(models,ranges,[depth])[:,0]
                scores=np.sqrt(np.mean((features-obs.relative_tl)**2,axis=(1,2)))
                residual=(bearings-obs.bearing_rad+np.pi)%(2*np.pi)-np.pi
                costs=np.sum(residual**2,axis=1);rms_deg=np.degrees(np.sqrt(np.mean(residual**2,axis=1)))
                calls+=1;unique_states+=len(part)
                for i,(key,state) in enumerate(part):memo[key]=(float(scores[i]),float(costs[i]),float(rms_deg[i]))
            print('RECONSTRUCT',cid,'depth',depth,'unique states',len(items),'total',unique_states,flush=True)
    cutoffs=manifest.set_index('case_id').bearing_cutoff.to_dict();bad_score_cases=set();max_score_delta=0.;max_cost_delta=0.
    for solver,frame in raws.items():
        for (cid,level),group in frame.groupby(['case_id','budget'],sort=False):
            differences=[];bearing_differences=[];valid=[]
            for row in group.itertuples():
                state=np.array([getattr(row,k) for k in AXES],dtype='<f8');j,cost,_=memo[(cid,float(row.depth_label_m),state.tobytes())]
                dj=abs(j-row.J_exact);dc=abs(cost-row.bearing_cost);differences.append(dj);bearing_differences.append(dc)
                valid.append(dj<=TOL_J and dc<=TOL_BEARING and cost<=cutoffs[cid]+1e-14 and np.all(state>=t.inherited.prior.LOW) and np.all(state<=t.inherited.prior.HIGH))
            mj=max(differences,default=0.);mc=max(bearing_differences,default=0.);max_score_delta=max(max_score_delta,mj);max_cost_delta=max(max_cost_delta,mc)
            check('Independent geometry/direct score '+solver+'/'+cid+'/'+level,all(valid),(len(group),mj,mc))
            if not all(valid):bad_score_cases.add((solver,cid,level))
            reconstructed.append(dict(solver=solver,case_id=cid,budget=level,n_witnesses=len(group),maximum_J_delta=mj,maximum_bearing_cost_delta=mc,pass_check=all(valid)))
    for solver in ['SHGO','DIRECT']:
        raw=raws[solver];log=logs[solver]
        for row in results[solver].itertuples():
            tag=solver+'/'+row.case_id+'/'+row.budget;branches=log[(log.case_id==row.case_id)&(log.budget==row.budget)]
            witnesses=raw[(raw.case_id==row.case_id)&(raw.budget==row.budget)].sort_values(['J_exact','depth_label_m',*AXES],kind='stable')
            complete=len(branches)==21 and set(branches.depth_label_m)==set(t.L.PROFILE)
            valid_count=sum(branch_valid[(solver,row.case_id,row.budget,float(z))] for z in branches.depth_label_m)
            valid=complete and valid_count==21 and (solver,row.case_id,row.budget) not in bad_score_cases
            check('Case branch completion/validity '+tag,row.n_depth_branches_completed==len(branches) and row.n_valid_branches==valid_count and row.execution_valid==valid and row.all_branches_execution_valid==valid)
            for aggregate,branch in [('total_objective_requests','n_objective_requests'),('total_exact_forward_evaluations','n_exact_forward_evaluations'),('total_unique_states','n_unique_states_evaluated'),('total_cache_hits','n_cache_hits')]:check('Case aggregate '+aggregate+'/'+tag,getattr(row,aggregate)==int(branches[branch].sum()))
            if len(witnesses):
                top=witnesses.iloc[0];state=np.array([getattr(row,k) for k in AXES],dtype='<f8');j,cost,_=memo[(row.case_id,float(row.z_star_label_m),state.tobytes())]
                check('Raw top exact/state/depth '+tag,row.best_exact_J==top.J_exact and row.z_star_label_m==top.depth_label_m and all(getattr(row,k)==top[k] for k in AXES) and abs(j-row.best_exact_J)<=TOL_J and abs(cost-row.bearing_cost)<=TOL_BEARING)
                check('Strict reconstructed recovery '+tag,row.recovered==bool(j<.001))
            else:check('No feasible witness recovery '+tag,not row.recovered and np.isnan(row.best_exact_J))
    raw_export=read('DEVELOPMENT_RAW_WITNESS_CONVERGENCE.csv');raw_expected=[];indexed={s:results[s].set_index(['case_id','budget']) for s in results}
    all_valid=all(results[s].execution_valid.all() for s in results)
    for solver in ['SHGO','DIRECT']:
        for cid in ids:
            left=indexed[solver].loc[(cid,'T2')];right=indexed[solver].loc[(cid,'T3')]
            cats=[catalogs[solver][(catalogs[solver].case_id==cid)&(catalogs[solver].budget==b)] for b in ['T2','T3']]
            flags=dict(T2_execution_valid=bool(left.execution_valid),T3_execution_valid=bool(right.execution_valid),T2_recovered=bool(left.recovered and left.best_exact_J<.001),T3_recovered=bool(right.recovered and right.best_exact_J<.001),top_state_depth_agreement=agrees(left,right),T3_clusters_contained_in_T2=contains(cats[1],cats[0]))
            value=all(flags.values());raw_expected.append(dict(case_id=cid,solver=solver,**flags,raw_T2_T3_converged=value))
            exported=raw_export[(raw_export.case_id==cid)&(raw_export.solver==solver)]
            check('Independent raw comparison '+solver+'/'+cid,len(exported)==1 and all(bool(exported.iloc[0][key])==bool(v) for key,v in dict(flags,raw_T2_T3_converged=value).items()))
    raw_gate=bool(all_valid and all(r['raw_T2_T3_converged'] for r in raw_expected));dual_expected=[]
    if raw_gate:
        dual_export=read('DUAL_SOLVER_AGREEMENT.csv')
        for cid in ids:
            left=indexed['SHGO'].loc[(cid,'T3')];right=indexed['DIRECT'].loc[(cid,'T3')]
            cats=[catalogs[s][(catalogs[s].case_id==cid)&(catalogs[s].budget=='T3')] for s in ['SHGO','DIRECT']]
            value=bool(left.execution_valid and right.execution_valid and left.recovered and right.recovered and agrees(left,right) and contains(cats[0],cats[1]) and contains(cats[1],cats[0]))
            dual_expected.append(value)
            check('Independent dual comparison '+cid,len(dual_export[dual_export.case_id==cid])==1 and bool(dual_export[dual_export.case_id==cid].iloc[0].dual_solver_subthreshold_witness_cluster_agreement)==value)
    else:check('Dual comparison not reached',not (OUT/'DUAL_SOLVER_AGREEMENT.csv').exists())
    dual_gate=bool(raw_gate and len(dual_expected)==21 and all(dual_expected))
    joint=defaultdict(list)
    for solver,frame in catalogs.items():
        for row in frame[frame.budget=='T3'].itertuples():
            state=np.array([getattr(row,k) for k in AXES],dtype='<f8');j,_,bearing_rms=memo[(row.case_id,float(row.depth_label_m),state.tobytes())]
            if j<1e-6 and bearing_rms<1e-8:joint[row.case_id].append(state)
    alias=any(np.any(delta(left,right)>t.DEDUP) for states in joint.values() for i,left in enumerate(states) for right in states[i+1:])
    try:alias_frame=read('ALIAS_DIAGNOSTIC.csv');exported_alias=bool(len(alias_frame))
    except pd.errors.EmptyDataError:exported_alias=False
    check('Strict independent alias finding',alias==exported_alias,'tested catalogs only')
    decision=json.loads((OUT/'DEVELOPMENT_EXECUTION_DECISION.json').read_text())
    check('Development Gate aggregation',decision['DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED']==raw_gate and decision['DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED']==dual_gate)
    expected_decision='A1_CONTINUOUS_ALIAS_FINDING_STOP' if alias else ('DEVELOPMENT_GATES_PASS' if dual_gate else 'A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED')
    check('Scientific decision reconstruction',decision['decision']==expected_decision)
    check('Fresh never generated',not d.FRESH_CONFIRMATION_RELEASED and not list(OUT.glob('FRESH_*')))
    check('All raw executions retained',sum(len(x) for x in results.values())==126 and sum(len(x) for x in logs.values())==2646)
    totals={key:int(sum(f[key].sum() for f in logs.values())) for key in ['n_objective_requests','n_exact_forward_evaluations','n_unique_states_evaluated','n_cache_hits','n_blocked_requests','n_validation_exact_forward_evaluations']}
    failed_checks=[r for r in checks if not r['pass_check']]
    summary=dict(reconstruction_valid=not bool(failed_checks),integrity_checks_passed=len(checks)-len(failed_checks),integrity_checks_failed=len(failed_checks),failed_checks=failed_checks,
        independent_modal_batch_calls=calls,independent_unique_case_depth_states_reconstructed=unique_states,maximum_J_delta=max_score_delta,maximum_bearing_cost_delta=max_cost_delta,optimizer_and_driver_totals=totals,
        raw_convergence={s:sum(r['raw_T2_T3_converged'] for r in raw_expected if r['solver']==s) for s in ['SHGO','DIRECT']},
        recovery_counts={s:{level:int(results[s][results[s].budget==level].recovered.sum()) for level in ['T1','T2','T3']} for s in ['SHGO','DIRECT']},
        dual_solver_agreement=sum(dual_expected) if raw_gate else 'NOT_REACHED',exact_alias_finding=alias,scientific_decision=expected_decision if not failed_checks else 'EXECUTION_INVALID',elapsed_seconds=time.perf_counter()-started,
        execution_invalid_cases=[dict(solver=s,case_id=r.case_id,budget=r.budget) for s,f in results.items() for r in f[~f.execution_valid].itertuples()],
        failed_recovery_cases=[dict(solver=s,case_id=r.case_id,budget=r.budget) for s,f in results.items() for r in f[~f.recovered].itertuples()],
        failed_convergence_cases=[dict(solver=r['solver'],case_id=r['case_id']) for r in raw_expected if not r['raw_T2_T3_converged']],
        hard_budget_enforcement_failures=int(sum((f.termination_reason=='HARD_BUDGET_ENFORCEMENT_FAILED').sum() for f in logs.values())),fresh_confirmation='NOT_RELEASED' if dual_gate else 'NOT_REACHED',R4_progress_percent=0,A2_depth_SSP_P5='UNOPENED')
    t.save('DEVELOPMENT_RECONSTRUCTION_AUDIT.csv',pd.DataFrame(checks));t.save('DEVELOPMENT_SCORE_RECONSTRUCTION.csv',pd.DataFrame(reconstructed));t.json_write('DEVELOPMENT_RECONSTRUCTION_SUMMARY.json',summary)
    print('RECONSTRUCTION COMPLETE',json.dumps({k:v for k,v in summary.items() if k not in ['failed_recovery_cases','failed_convergence_cases','failed_checks']},indent=2),flush=True)
    if failed_checks:raise SystemExit('EXECUTION_INVALID: reconstruction failed; no algorithm repair or rerun authorized')


if __name__=='__main__':audit()
