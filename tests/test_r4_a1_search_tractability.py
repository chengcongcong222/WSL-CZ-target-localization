import inspect
import numpy as np
import pandas as pd
import pytest
from unittest.mock import patch
from threadpoolctl import threadpool_limits
import r4_a1_search_tractability as t
import r4_a1_search_tractability_audit as a
import hashlib


class Quadratic:
    def __init__(self):self.calls=0;self.preparations=0
    def prepare(self,q):
        self.preparations+=1
        return t.Prepared(np.r_[q,0.,0.],np.zeros((1,121)),0.,True)
    def key(self,p):return (200.,p.state.astype('<f8').tobytes())
    def evaluate(self,p):
        self.calls+=1
        score=float(np.sum((p.state[:2]-.37)**2))
        return t.Evaluation(score,score,p.state,0.,200.,True,1)


@pytest.fixture(scope='module')
def physical():
    with threadpool_limits(limits=2):
        models=t.inherited.inputs();state=np.array([51.37,.32,2.11,3.42])
        bearing,ranges=t.L.geometry(state)
        tl=t.L.direct_features(models,ranges,[200.])[0,0]
        yield models,state,bearing[0],tl


@pytest.mark.parametrize('limit',[1,2,7])
def test_wrapper_boundary_and_cache_requests(limit):
    raw=Quadratic();wrapper=t.BudgetObjective(raw,limit)
    values=[wrapper(np.array([.31,.62])) for _ in range(limit)]
    assert len(set(values))==1
    assert wrapper.n_objective_requests==limit and wrapper.n_exact_forward_evaluations==1
    assert wrapper.n_cache_hits==limit-1 and wrapper.n_unique_states_evaluated==1
    assert raw.preparations==limit
    with pytest.raises(t.ObjectiveBudgetExceeded):wrapper(np.array([.8,.9]))
    assert wrapper.n_objective_requests==limit and wrapper.n_blocked_requests==1
    assert wrapper.n_objective_attempts==limit+1 and raw.preparations==limit and raw.calls==1
    assert sum(row['admitted'] for row in wrapper.ledger)==limit


@pytest.mark.parametrize('solver',['SHGO','DIRECT'])
@pytest.mark.parametrize('limit',[1,3,7,20,64,128])
def test_native_exception_propagation_and_no_late_callback(solver,limit):
    wrapper=t.BudgetObjective(Quadratic(),limit)
    options={'n':256,'iters':1} if solver=='SHGO' else {'maxfun':4096,'maxiter':100}
    native,record=t.run_solver(solver,wrapper,options)
    assert native is None and not record['solver_success']
    assert record['termination_reason']=='HARD_BUDGET_EXCEEDED'
    assert record['hard_budget_enforcement_supported']
    assert wrapper.n_objective_requests==limit and wrapper.n_blocked_requests==1
    assert wrapper.n_objective_attempts==limit+1
    expected='["ObjectiveBudgetExceeded"]' if solver=='SHGO' else '["SystemError", "ObjectiveBudgetExceeded"]'
    assert record['exception_chain']==expected


def test_shgo_local_refinement_and_cache_hits_consume_budget():
    wrapper=t.BudgetObjective(Quadratic(),12)
    _,record=t.run_solver('SHGO',wrapper,{'n':8,'iters':1,'minimizer_kwargs':{'method':'SLSQP','options':{'maxiter':100,'ftol':1e-14}}})
    assert record['termination_reason']=='HARD_BUDGET_EXCEEDED'
    assert wrapper.n_objective_requests==12 and wrapper.n_blocked_requests==1
    assert wrapper.n_cache_hits>0 and wrapper.n_exact_forward_evaluations<12


def test_unrelated_system_error_is_not_relabelled():
    wrapper=t.BudgetObjective(Quadratic(),20)
    with patch.object(t,'direct',side_effect=SystemError('unrelated foreign boundary error')):
        with pytest.raises(SystemError,match='unrelated'):t.run_solver('DIRECT',wrapper,{})


def test_repeated_callbacks_after_exhaustion_are_rejected():
    wrapper=t.BudgetObjective(Quadratic(),1)
    def broken(fn,bounds,**options):
        fn(np.array([.2,.3]))
        for _ in range(2):
            try:fn(np.array([.4,.5]))
            except t.ObjectiveBudgetExceeded:pass
        fn(np.array([.4,.5]))
    with patch.object(t,'shgo',side_effect=broken):
        _,record=t.run_solver('SHGO',wrapper,{})
    assert record['termination_reason']=='HARD_BUDGET_ENFORCEMENT_FAILED'
    assert not record['hard_budget_enforcement_supported'] and wrapper.n_objective_requests==1


def test_cached_and_uncached_exact_identity(physical):
    models,state,bearing,tl=physical
    engine=t.ExactObjective(models,bearing,tl,13.3*np.radians(.1)**2,200.)
    q=np.array([.41,.63]);cached=t.BudgetObjective(engine,3);uncached=t.BudgetObjective(engine,3,cache_enabled=False)
    with threadpool_limits(limits=2):
        a=[cached(q) for _ in range(3)];b=[uncached(q) for _ in range(3)]
    assert a==b and len(set(a))==1
    assert cached.n_exact_forward_evaluations==1 and cached.n_cache_hits==2
    assert uncached.n_exact_forward_evaluations==3 and uncached.n_cache_hits==0
    assert cached.n_unique_states_evaluated==uncached.n_unique_states_evaluated==1


def test_exact_self_match_and_radial_transform(physical):
    models,state,bearing,tl=physical
    u=state[2]*np.cos(np.radians(state[3]-state[1]))
    q=np.array([(1000*state[0]-45000)/15000,(u-.94)/2.06])
    engine=t.ExactObjective(models,bearing,tl,1e-20,200.)
    p=engine.prepare(q)
    np.testing.assert_allclose(p.state,state,rtol=0,atol=1e-10)
    assert p.feasible and p.bearing_cost<1e-20
    assert engine.evaluate(p).J_exact<1e-7


def test_geometry_reconstruction_and_fixed_domains(physical):
    models,state,bearing,tl=physical
    engine=t.ExactObjective(models,bearing,tl,13.3*np.radians(.1)**2,200.)
    p=engine.prepare(np.array([.41,.63]));r,theta,v,psi=p.state
    times=np.arange(121)*10.;post=np.maximum(times-600,0.)
    xp=2*np.minimum(times,600)+2*post*np.cos(np.pi/12);yp=2*post*np.sin(np.pi/12)
    x=1000*r*np.cos(np.radians(theta))+v*times*np.cos(np.radians(psi))-xp
    y=1000*r*np.sin(np.radians(theta))+v*times*np.sin(np.radians(psi))-yp
    np.testing.assert_allclose(p.ranges[0],np.sqrt(x*x+y*y),rtol=0,atol=2e-11)
    assert t.R_BOUNDS==(45000.,60000.) and t.U_BOUNDS==(.94,3.)


def test_bearing_and_state_rejection_precedes_exact_forward(physical):
    models,state,bearing,tl=physical
    restrictive=t.ExactObjective(models,bearing,tl,0.,200.)
    wrapper=t.BudgetObjective(restrictive,2)
    for q in [np.array([.05,.01]),np.array([1.1,.5])]:assert wrapper(q)>=1e6
    assert wrapper.n_exact_forward_evaluations==0 and not len(t.candidate_catalog(wrapper))


def test_inherited_depth_mapping(physical):
    models,_,_,_=physical
    for model in models.values():
        _,actual=t.L.effective_depths(model,[155.,180.,200.,220.])
        np.testing.assert_array_equal(actual,[154.,180.,200.,220.])


def test_catalog_preserves_separated_states_and_exact_order():
    wrapper=t.BudgetObjective(Quadratic(),4)
    for score,r in [(3.,50.006),(2.,50.0049),(1.,50.)]:
        state=np.array([r,.3,2.,3.])
        wrapper.cache[(200.,state.astype('<f8').tobytes())]=t.Evaluation(score,score,state,0.,200.,True,1)
    catalog=t.candidate_catalog(wrapper)
    assert list(catalog.J_exact)==[1.,3.]
    assert list(catalog.r_km)==[50.,50.006]
    assert (catalog.candidate_kind=='EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM').all()


def test_static_search_path_exclusion():
    names=['truth','errors','neighborhood','oracle','basin_width','basin_direction']
    functions=[t.ExactObjective.__init__,t.ExactObjective.prepare,t.ExactObjective.key,t.ExactObjective.evaluate,t.BudgetObjective.__call__,t.run_solver,t.candidate_catalog,t.verify_catalog,t.inherited.radial_profile,t.L.geometry,t.L.direct_features]
    for fun in functions:
        source=inspect.getsource(fun).lower()
        assert not any(name in source or name in str(inspect.signature(fun)).lower() for name in names),(fun.__name__,source)


@pytest.mark.parametrize('solver',['SHGO','DIRECT'])
def test_runtime_sentinel_solver_extraction_and_verification(physical,solver,tmp_path):
    models,state,bearing,_=physical
    # Arbitrary zero relative TL makes this a compatibility fixture, not recovery evidence.
    tl=np.zeros((3,121));cut=13.3*np.radians(.1)**2
    wrapper=t.BudgetObjective(t.ExactObjective(models,bearing,tl,cut,200.),16)
    options={'n':8,'iters':1} if solver=='SHGO' else {'maxfun':128,'maxiter':20}
    with patch.object(t.L,'errors',side_effect=AssertionError('evaluation access')),patch.object(t.L,'neighborhood',side_effect=AssertionError('evaluation access')),patch.object(t.L,'generate_observation',side_effect=AssertionError('generation access')),patch.object(t.inherited,'basin_geometry',side_effect=AssertionError('diagnostic access')):
        with threadpool_limits(limits=2):
            _,record=t.run_solver(solver,wrapper,options)
            catalog=t.candidate_catalog(wrapper)
            assert len(catalog)>0
            path=tmp_path/'candidates.csv';catalog.to_csv(path,index=False,float_format='%.17g')
            loaded=pd.read_csv(path,float_precision='round_trip')
            checks=t.verify_catalog(models,bearing,tl,cut,loaded)
    assert checks.pass_check.all() and record['hard_budget_enforcement_supported']


@pytest.mark.parametrize('bad',[0,-1,True,2.5])
def test_invalid_budget_rejected(bad):
    with pytest.raises(ValueError):t.BudgetObjective(Quadratic(),bad)


@pytest.mark.parametrize('checkout',[b'x,y\n1,2\n',b'x,y\r\n1,2\r\n'])
def test_historical_text_checkout_line_endings_are_portable(tmp_path,checkout):
    canonical=b'x,y\n1,2\n';original=b'x,y\r\n1,2\r\n'
    p=tmp_path/'evidence.csv';p.write_bytes(checkout)
    oid=hashlib.sha1(b'blob 8\0'+canonical).hexdigest()
    assert a.verify_frozen_input(p,hashlib.sha256(original).hexdigest(),hashlib.sha256(canonical).hexdigest(),oid,True) in ['RAW_BYTES_IDENTICAL','GIT_TEXT_CRLF_LF_ONLY']
    p.write_bytes(checkout.replace(b'1,2',b'1,3'))
    with pytest.raises(AssertionError):a.verify_frozen_input(p,hashlib.sha256(original).hexdigest(),hashlib.sha256(canonical).hexdigest(),oid,True)


def test_binary_or_explicit_no_text_file_never_uses_line_ending_fallback(tmp_path):
    canonical=b'x,y\n1,2\n';original=b'x,y\r\n1,2\r\n'
    p=tmp_path/'opaque.bin';p.write_bytes(canonical)
    oid=hashlib.sha1(b'blob 8\0'+canonical).hexdigest()
    with pytest.raises(AssertionError):a.verify_frozen_input(p,hashlib.sha256(original).hexdigest(),hashlib.sha256(canonical).hexdigest(),oid,False)


# Harness tests use synthetic observations and non-acoustic wiring substitutes.
import json
import subprocess
import sys
import r4_a1_search_tractability_development as d
import r4_a1_search_tractability_gates as g


def synthetic_manifest():
    return pd.DataFrame([dict(case_id=f'SYNTH_{i:02}',panel_id='FIXTURE',sigma_deg=.1,seed=i+1,bearing_cutoff=.001,observation_index=i,origin_group='SYNTHETIC_METADATA') for i in range(21)])[d.MANIFEST_COLUMNS]


def synthetic_archive(manifest):
    return dict(case_ids=manifest.case_id.to_numpy(dtype=str),bearing_rad=np.repeat(np.arange(21)[:,None],121,axis=1),relative_tl=np.repeat(np.arange(21)[:,None,None],3*121,axis=1).reshape(21,3,121))


def witness(score=.0005,r=50.,depth=200.,case_id='SYNTH',solver='SHGO',budget='T3',theta=.3,psi=3.):
    return dict(case_id=case_id,solver=solver,budget=budget,depth_label_m=depth,candidate_kind=g.KIND,J_exact=score,r_km=r,theta_deg=theta,v_mps=2.,psi_deg=psi,bearing_cost=0.)


def case_result(score=.0005,r=50.,depth=200.,**updates):
    value=dict(best_exact_J=score,recovered=score<.001,r_km=r,theta_deg=.3,v_mps=2.,psi_deg=3.,z_star_label_m=depth,
        execution_valid=True,all_branches_execution_valid=True,n_depth_branches_completed=21,n_valid_branches=21)
    value.update(updates);return value


class SemanticObjective:
    """Non-acoustic synthetic substitute used only for harness wiring tests."""
    def __init__(self,models,bearing,relative_tl,cutoff,depth):self.depth=depth
    def prepare(self,q):return t.Prepared(np.array([50.+float(q[0])*.01,.3,2.,3.]),np.zeros((1,121)),0.,True)
    def key(self,p):return (self.depth,p.state.astype('<f8').tobytes())
    def evaluate(self,p):
        score=.0002 if self.depth==200. else .0009
        return t.Evaluation(score,score,p.state,0.,self.depth,True,0)


def semantic_verifier(models,observation,cutoff,catalog):return catalog.sort_values('J_exact',kind='stable').reset_index(drop=True)


def synthetic_record(wrapper,solver,reason='SOLVER_RETURN',supported=True):
    return dict(solver=solver,termination_reason=reason,solver_success=reason=='SOLVER_RETURN',solver_message='SYNTHETIC_WIRING_TEST',
        exception_chain=json.dumps(['ObjectiveBudgetExceeded'] if wrapper.hard_budget_triggered else []),hard_budget_enforcement_supported=supported,**wrapper.counters())


def semantic_runner(solver,wrapper,options):
    wrapper(np.array([.2,.4]));wrapper(np.array([.2,.4]))
    return None,synthetic_record(wrapper,solver)


@pytest.fixture
def harness_budget():return d.verify_infrastructure()


@pytest.fixture
def semantic_observation():return d.Observation(np.zeros(121),np.zeros((3,121)))


def test_real_manifest_metadata_only_and_21_case_projection():
    frame=pd.read_csv(t.OUT/'DEVELOPMENT_CASE_MANIFEST.csv',float_precision='round_trip');d.validate_manifest(frame)
    assert len(frame)==frame.case_id.nunique()==21 and set(frame.columns)==set(d.MANIFEST_COLUMNS)
    with np.load(t.ROOT/d.OBSERVATION_PATH,allow_pickle=False) as archive:mapped=d.map_observations(frame,archive)
    assert len(mapped)==21 and all(not item.bearing_rad.flags.writeable for item in mapped.values())


@pytest.mark.parametrize('column',['truth_r','z_true','top1_error','old_recovered','old_J','r_km','oracle_center'])
def test_manifest_rejects_privileged_columns(column):
    frame=synthetic_manifest();frame[column]=0
    with pytest.raises(ValueError,match='whitelist'):d.validate_manifest(frame)


@pytest.mark.parametrize('bad',['missing','duplicate','noiseless','fractional_index','duplicate_index'])
def test_manifest_rejects_incomplete_or_non_nominal_panel(bad):
    frame=synthetic_manifest()
    if bad=='missing':frame=frame.iloc[:-1]
    if bad=='duplicate':frame.loc[0,'case_id']=frame.loc[1,'case_id']
    if bad=='noiseless':frame.loc[0,'sigma_deg']=0.
    if bad=='fractional_index':frame['observation_index']=frame.observation_index.astype(float);frame.loc[0,'observation_index']=.5
    if bad=='duplicate_index':frame.loc[0,'observation_index']=1
    with pytest.raises(ValueError):d.validate_manifest(frame)


def test_observation_mapping_uses_ids_not_manifest_row_order():
    frame=synthetic_manifest();archive=synthetic_archive(frame);shuffled=frame.sample(frac=1,random_state=19)
    mapped=d.map_observations(shuffled,archive)
    for row in shuffled.itertuples():assert np.all(mapped[row.case_id].bearing_rad==row.observation_index)
    bad=shuffled.copy();bad.loc[0,'observation_index']=1;bad.loc[1,'observation_index']=0
    with pytest.raises(ValueError,match='mapping'):d.map_observations(bad,archive)


@pytest.mark.parametrize('bad',['duplicate_id','missing_id','nonfinite','shape'])
def test_corrupted_observations_rejected(bad):
    frame=synthetic_manifest();archive=synthetic_archive(frame)
    if bad=='duplicate_id':archive['case_ids'][0]=archive['case_ids'][1]
    if bad=='missing_id':archive['case_ids'][0]='ABSENT'
    if bad=='nonfinite':archive['relative_tl']=archive['relative_tl'].astype(float);archive['relative_tl'][0,0,0]=np.nan
    if bad=='shape':archive['bearing_rad']=archive['bearing_rad'][:,:120]
    with pytest.raises(ValueError):d.map_observations(frame,archive)


def test_all_solver_budget_branch_objects_are_fresh(harness_budget,semantic_observation):
    engines=[];wrappers=[];caches=[]
    def factory(*args):
        engine=SemanticObjective(*args);engines.append(engine);return engine
    def runner(solver,wrapper,options):
        assert wrapper.n_objective_requests==0 and not wrapper.cache and not wrapper.ledger
        wrappers.append(wrapper);caches.append(wrapper.cache);return semantic_runner(solver,wrapper,options)
    for solver in ['SHGO','DIRECT']:
        for level in ['T1','T2','T3']:
            result,catalog,logs,requests,raw=d.run_case(None,semantic_observation,'SYNTH_ONLY',.001,solver,level,harness_budget,objective_factory=factory,runner=runner,verifier=semantic_verifier)
            assert result['n_depth_branches_completed']==result['n_valid_branches']==21
            assert result['execution_valid'] and result['z_star_label_m']==200.
            assert set(logs.depth_label_m)==set(t.L.PROFILE) and len(requests)==21
            assert result['total_objective_requests']==42 and result['total_cache_hits']==21
            assert result['total_exact_forward_evaluations']==0 and len(raw)==21
    assert len(engines)==126 and len({id(x) for x in engines})==126
    assert len({id(x) for x in wrappers})==len({id(x) for x in caches})==126


def test_hard_cap_exhaustion_can_recover(harness_budget,semantic_observation):
    def runner(solver,wrapper,options):
        for _ in range(wrapper.limit):wrapper(np.array([.2,.4]))
        with pytest.raises(t.ObjectiveBudgetExceeded):wrapper(np.array([.8,.9]))
        return None,synthetic_record(wrapper,solver,'HARD_BUDGET_EXCEEDED')
    result,_,logs,_,_=d.run_case(None,semantic_observation,'SYNTH_CAP',.001,'SHGO','T1',harness_budget,objective_factory=SemanticObjective,runner=runner,verifier=semantic_verifier)
    assert result['recovered'] and result['execution_valid'] and not logs.solver_success.any()
    assert (logs.n_objective_requests==16).all() and (logs.n_blocked_requests==1).all()


@pytest.mark.parametrize('bad',['enforcement','counter','uncaught','reconstruction'])
def test_execution_invalid_never_passes_gate(bad,harness_budget,semantic_observation):
    def runner(solver,wrapper,options):
        _,record=semantic_runner(solver,wrapper,options)
        if bad=='uncaught':raise SystemError('unrelated')
        if bad=='enforcement':record['termination_reason']='HARD_BUDGET_ENFORCEMENT_FAILED';record['hard_budget_enforcement_supported']=False
        if bad=='counter':record['n_cache_hits']+=1
        return None,record
    def verify(*args):
        if bad=='reconstruction':raise AssertionError('saved exact score mismatch')
        return semantic_verifier(*args)
    result,_,logs,_,_=d.run_case(None,semantic_observation,'SYNTH_INVALID',.001,'SHGO','T1',harness_budget,objective_factory=SemanticObjective,runner=runner,verifier=verify)
    assert not result['execution_valid'] and result['n_valid_branches']==0 and not logs.execution_valid.any()
    assert not g.result_valid(result)


def test_incomplete_branch_union_invalid(harness_budget,semantic_observation):
    result,_,logs,_,raw=d.run_case(None,semantic_observation,'SYNTH_COUNT',.001,'SHGO','T1',harness_budget,objective_factory=SemanticObjective,runner=semantic_runner,verifier=semantic_verifier)
    incomplete=d.aggregate_case('SYNTH_COUNT','SHGO','T1',.001,logs.to_dict('records')[:-1],raw)
    assert result['best_exact_J']==.0002 and result['z_star_label_m']==200.
    assert incomplete['recovered'] and not incomplete['execution_valid']


def test_state_circular_angles_and_exact_depth_agreement():
    assert g.state_agrees([50.,359.9995,2.,179.9999],[50.,.0001,2.,-179.9999])
    assert not g.state_agrees([50.,0.,2.,3.],[50.,.002,2.,3.])
    assert not g.top_agrees(case_result(),case_result(depth=200.+1e-10))


def test_depth_separate_score_first_witness_dedup():
    frame=pd.DataFrame([witness(.0009,r=50.0049),witness(.0001),witness(.0002,depth=205.),witness(.0003,r=50.006)])
    reduced=g.dedup_witnesses(frame)
    assert list(reduced.J_exact)==[.0001,.0002,.0003]
    with pytest.raises(ValueError):g.dedup_witnesses(frame.assign(candidate_kind='LOCAL_MINIMUM'))


def test_subthreshold_cluster_containment_is_directional():
    small=pd.DataFrame([witness()]);extra=pd.DataFrame([witness(),witness(r=50.02)])
    assert g.clusters_contained(small,extra) and not g.clusters_contained(extra,small)
    assert not g.raw_case_convergence(case_result(),case_result(),small,extra)['raw_T2_T3_converged']
    assert g.raw_case_convergence(case_result(),case_result(),extra,small)['raw_T2_T3_converged']
    assert not g.clusters_contained(small,small.assign(J_exact=.001))
    assert not g.clusters_contained(small,small.assign(depth_label_m=205.))


def test_dual_matching_is_bidirectional_and_threshold_is_strict():
    small=pd.DataFrame([witness()]);extra=pd.DataFrame([witness(),witness(r=50.02)])
    assert g.dual_case_agreement(case_result(),case_result(),small,small)['dual_solver_subthreshold_witness_cluster_agreement']
    assert not g.dual_case_agreement(case_result(),case_result(),small,extra)['dual_solver_subthreshold_witness_cluster_agreement']
    assert not g.result_recovered(case_result(score=.001))


def synthetic_development_tables():
    ids=[f'SYNTH_{i:02}' for i in range(21)];results={};catalogs={}
    for solver in ['SHGO','DIRECT']:
        results[solver]=pd.DataFrame([dict(case_result(),case_id=cid,solver=solver,budget=level) for cid in ids for level in ['T1','T2','T3']])
        catalogs[solver]=pd.DataFrame([witness(case_id=cid,solver=solver,budget=level) for cid in ids for level in ['T1','T2','T3']])
    return ids,results,catalogs


def test_full_21_synthetic_gate_logic_and_fail_stops_dual(monkeypatch):
    ids,results,catalogs=synthetic_development_tables();decision,raw,dual=g.development_gates(results,catalogs,ids)
    assert decision['DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED'] and decision['DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED']
    assert len(raw)==42 and len(dual)==21
    results['SHGO'].loc[1,'recovered']=False
    monkeypatch.setattr(g,'dual_case_agreement',lambda *args:pytest.fail('dual comparison reached after raw failure'))
    decision,raw,dual=g.development_gates(results,catalogs,ids)
    assert not decision['DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED'] and dual.empty
    assert d.run_fresh_confirmation(decision,lambda:pytest.fail('fresh generator called'))=='NOT_REACHED_DUE_TO_DEVELOPMENT_GATE'


def test_fresh_confirmation_not_released_even_if_development_passes():
    with pytest.raises(SystemExit,match='FRESH_CONFIRMATION_NOT_RELEASED'):
        d.run_fresh_confirmation(dict(DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED=True,DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED=True),lambda:pytest.fail('unreleased generator called'))


@pytest.mark.parametrize('solver',['SHGO','DIRECT'])
def test_harness_runtime_sentinels_native_solver_and_extraction(solver,harness_budget,semantic_observation):
    with patch.object(t.L,'errors',side_effect=AssertionError('privileged access')),patch.object(t.L,'neighborhood',side_effect=AssertionError('privileged access')),patch.object(t.L,'generate_observation',side_effect=AssertionError('privileged access')),patch.object(t.inherited,'basin_geometry',side_effect=AssertionError('privileged access')):
        branch,catalog,ledger=d.run_branch(None,semantic_observation,.001,'SYNTH_NATIVE',solver,'T1',200.,harness_budget,objective_factory=SemanticObjective,verifier=semantic_verifier)
    assert branch['execution_valid'] and len(catalog)>0 and branch['n_objective_requests']<=16


def test_unreleased_real_driver_stops_before_input_loading(monkeypatch):
    monkeypatch.setattr(d,'load_execution_inputs',lambda:pytest.fail('scientific input loading reached'))
    with pytest.raises(SystemExit,match='DEVELOPMENT_NOT_RELEASED'):d.execute_development()
    with pytest.raises(SystemExit,match='DEVELOPMENT_NOT_RELEASED'):d.run_branch(None,None,.001,'ANY','SHGO','T1',200.,{})
    r=subprocess.run([sys.executable,str(t.ROOT/'r4_a1_search_tractability_development.py')],capture_output=True,text=True)
    assert r.returncode!=0 and 'DEVELOPMENT_NOT_RELEASED' in r.stderr


def test_frozen_options_are_copied_without_mutating_budget(harness_budget):
    before=json.dumps(harness_budget,sort_keys=True)
    for solver in ['SHGO','DIRECT']:
        cap,options=d.solver_options(harness_budget,solver,'T1');assert cap==16
        options.clear()
    assert json.dumps(harness_budget,sort_keys=True)==before


def test_budget_tamper_stops_before_execution(tmp_path,monkeypatch):
    original_core=(t.ROOT/'r4_a1_search_tractability.py').read_bytes();original_budget=(t.OUT/'TRACTABILITY_BUDGET_FREEZE.json').read_bytes()
    (tmp_path/'r4_a1_search_tractability.py').write_bytes(original_core);(tmp_path/'TRACTABILITY_BUDGET_FREEZE.json').write_bytes(original_budget+b' ')
    monkeypatch.setattr(t,'ROOT',tmp_path);monkeypatch.setattr(t,'OUT',tmp_path)
    with pytest.raises(RuntimeError,match='FROZEN_INFRASTRUCTURE_CHANGED'):d.verify_infrastructure()


def test_alias_uses_strict_independent_joint_reconstruction():
    frame=pd.DataFrame([witness(score=1e-7),witness(score=1e-7,r=50.02,solver='DIRECT')])
    assert len(g.exact_alias_pairs(frame,lambda *args:(0.,0.)))==1
    assert g.exact_alias_pairs(frame,lambda *args:(1e-8,0.)).empty
    assert g.exact_alias_pairs(frame,lambda *args:(0.,1e-6)).empty


@pytest.mark.parametrize('solver',['SHGO','DIRECT'])
def test_harness_sentinel_through_direct_modal_verification(physical,solver,harness_budget):
    models,_,bearing,_=physical
    observation=d.Observation(bearing,np.zeros((3,121)))
    with patch.object(t.L,'errors',side_effect=AssertionError('privileged access')),patch.object(t.L,'neighborhood',side_effect=AssertionError('privileged access')),patch.object(t.L,'generate_observation',side_effect=AssertionError('privileged access')),patch.object(t.inherited,'basin_geometry',side_effect=AssertionError('privileged access')):
        with threadpool_limits(limits=2):
            record,catalog,_=d.run_branch(models,observation,13.3*np.radians(.1)**2,'SYNTH_EXACT_ONLY',solver,'T1',200.,harness_budget,objective_factory=t.ExactObjective)
    assert record['execution_valid'] and len(catalog)>0
    assert record['n_validation_exact_forward_evaluations']==len(catalog)


def test_entry_records_execution_invalid_for_loading_exception(monkeypatch,tmp_path):
    monkeypatch.setattr(d,'DEVELOPMENT_RELEASED',True)
    monkeypatch.setattr(d,'_execute_released_development',lambda:(_ for _ in ()).throw(ValueError('synthetic corrupted observation')))
    monkeypatch.setattr(t,'OUT',tmp_path)
    with pytest.raises(ValueError):d.execute_development()
    record=json.loads((tmp_path/'DEVELOPMENT_EXECUTION_INVALID.json').read_text())
    assert record['execution_valid']==False and record['DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED']==False
