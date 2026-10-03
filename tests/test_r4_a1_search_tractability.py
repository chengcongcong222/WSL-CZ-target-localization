import inspect
import numpy as np
import pandas as pd
import pytest
from unittest.mock import patch
from threadpoolctl import threadpool_limits
import r4_a1_search_tractability as t


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
