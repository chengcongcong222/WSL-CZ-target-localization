"""Analytic logic tests and fixed modal enclosure fixtures only."""
import ast,inspect,itertools,json,math
from pathlib import Path
import numpy as np
import pytest
from flint import arb
from r4_cz_envelope_representation import basis,feature,score
from r4_cz_envelope_support import (Cell,Bounds,Certificate,Observation,SupportEngine,HardBudget,BudgetExceeded,
    COUNTERS,DOMAIN_LOW,DOMAIN_HIGH,DEPTH_LABELS,COARSE_WIDTHS,FINE_WIDTHS,search_budget,exact,adjacent,components,raw_pair)
from r4_cz_envelope_model import (ModalModel,precision,geometry,project_balls,certified_basis,lower_bound_bearing_cost)
from r4_cz_envelope_guard import estimator_scope
from r4_cz_envelope_evaluation import retention

ROOT=Path(__file__).resolve().parents[1]
POLICY=json.loads((ROOT/'results/R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN/CZ_NUMERICAL_CALIBRATION_POLICY.json').read_text(encoding='utf-8-sig'))
OBS=Observation(np.zeros(121),np.zeros(12),1.)
ROOT_CELL=Cell(DOMAIN_LOW,DOMAIN_HIGH)
B0=Bounds(arb(0),arb(0))

def ample():return HardBudget(dict.fromkeys(COUNTERS,10000),60000)

def cert(cell,lo=0,hi=0,feature_lo=0,feature_hi=0,method='ANALYTIC_TOY_INTERVAL_V1',depths=DEPTH_LABELS):
    return Certificate(cell,Bounds(exact(lo),exact(hi)),Bounds(exact(feature_lo),exact(feature_hi)),depths,method)

def compatible(cell,observation,budget,baseline):
    budget.charge(n_bound_evaluations=1)
    return cert(cell)

def unknown(cell,observation,budget,baseline):
    budget.charge(n_bound_evaluations=1)
    return cert(cell,0,10,0,10)

@pytest.mark.parametrize('n',[61,60])
def test_q2_basis(n):
    b=basis(n)
    assert np.max(abs(b.sum(axis=0)))<1e-14
    assert np.max(abs(b.T@b-np.eye(2)))<1e-14
    with precision():
        rows=certified_basis(n)
        for k in (0,1):
            assert sum((r[k] for r in rows),arb(0)).contains(0)
            assert sum((r[k]**2 for r in rows),arb(0)).contains(1)
        assert sum((r[0]*r[1] for r in rows),arb(0)).contains(0)

def test_static_offset_invariance():
    x=np.tile(np.linspace(-1,1,121)**2,(3,1));offsets=np.zeros((3,121))
    for i,(a,b) in enumerate(POLICY['analytic_tests']['offsets_line_window']):offsets[i,:61]=a;offsets[i,61:]=b
    assert score(feature(x),feature(x+offsets))<1e-12

def test_exact_repeat_determinism():
    x=np.tile(np.linspace(-1,1,121),(3,1));assert np.array_equal(feature(x),feature(x))

@pytest.mark.parametrize('x',[np.zeros((3,120)),np.zeros((2,121)),np.zeros((3,121),dtype=complex),np.full((3,121),np.nan)])
def test_m0_rejects_invalid_shape_or_phase(x):
    with pytest.raises(ValueError):feature(x)

def test_feature_dimension_and_equal_score_norm():
    assert feature(np.zeros((3,121))).shape==(12,)
    assert score(np.ones(12),np.zeros(12))==pytest.approx(math.sqrt(2))

@pytest.mark.parametrize('field',['truth','range_truth','truth_motion','truth_depth','phase','oracle','source_spectrum'])
def test_prohibited_inputs(field):
    with pytest.raises(TypeError):feature(np.zeros((3,121)),**{field:0})
    with pytest.raises(TypeError):Observation(np.zeros(121),np.zeros(12),1.,**{field:0})

def test_n0_constant_polynomial_control():
    assert np.array_equal(feature(np.ones((3,121))),np.zeros(12))
    levels=[]
    with precision():
        for _ in range(3):levels.append([arb(7)]*61+[arb(-19)]*60)
        assert all(v.contains(0) for v in project_balls(levels))

def test_candidate_trajectory_domain_and_containment():
    with precision():
        angles,ranges=geometry(ROOT_CELL)
        for state in POLICY['fixed_states_r_km_theta_deg_v_mps_psi_deg']:
            a,r=geometry(Cell(tuple(state),tuple(state)))
            assert all(x.contains(y) for x,y in zip(angles,a))
            assert all(x.contains(y) for x,y in zip(ranges,r))
        assert all(x>0 for x in ranges)
    with pytest.raises(ValueError):Cell((44,-5,1,-15),DOMAIN_HIGH)

@pytest.fixture(scope='module')
def model():return ModalModel()

def test_finite_depth_enumeration(model):
    assert DEPTH_LABELS==tuple(range(150,251,5))
    assert len(model.mapping)==63
    for f in (201,235,283):assert [r['profile_label_m'] for r in model.mapping if r['frequency_hz']==f]==list(DEPTH_LABELS)
    with pytest.raises(ValueError):model.levels_arb(ROOT_CELL,201,ample())

def toy(cell,observation,budget,baseline):
    budget.charge(n_bound_evaluations=1)
    delta=exact(cell.low[0]).union(exact(cell.high[0]))-50
    return Certificate(cell,B0,None if baseline else Bounds(delta.abs_lower(),delta.abs_upper()),DEPTH_LABELS,'ANALYTIC_TOY_INTERVAL_V1')

def test_certified_bound_preserves_analytic_feasible_states():
    # Analytically known feasible interval |r-50| <= .25; no sampled rejection.
    result=SupportEngine(toy,OBS,.25,(.25,10,2,30),search_budget(200)).run()
    for r in [49.75,50.,50.25]:
        assert any(c.low[0]<=r<=c.high[0] for c in result.retained)
    for rec in result.records:
        if rec.status=='REJECTED_BY_CERTIFIED_BOUND':assert rec.certificate.feature.lower>.25

def test_unresolved_and_budget_cells_retained():
    result=SupportEngine(unknown,OBS,0,COARSE_WIDTHS,search_budget(1)).run()
    assert result.termination_reason=='HARD_CAP' and not result.coverage_closed
    assert result.metrics()['hull_width_km']==15
    volume=sum(math.prod(c.widths) for c in result.retained)
    assert volume==math.prod(ROOT_CELL.widths)
    assert len(result.retained)==2

def test_terminal_loose_bound_retained():
    result=SupportEngine(unknown,OBS,0,ROOT_CELL.widths,search_budget(1)).run()
    assert result.records[0].status=='RETAINED_UNRESOLVED'
    assert result.termination_reason=='TERMINAL_UNRESOLVED'
    assert not result.coverage_closed

def test_empty_support_is_failure():
    def reject(cell,observation,budget,baseline):
        budget.charge(n_bound_evaluations=1);return cert(cell,2,2)
    result=SupportEngine(reject,OBS,0,COARSE_WIDTHS,search_budget(1)).run()
    assert result.metrics()['empty_support'] and result.metrics()['hull_width_km'] is None
    assert retention(result,[50,0,2,0])=={'joint_retained':False,'range_retained':False}

def test_rejection_requires_valid_certificate():
    assert cert(ROOT_CELL,2,2,method='CENTER_SAMPLE').verdict(OBS,0,False)=='RETAINED_UNRESOLVED'
    assert cert(ROOT_CELL,0,0,2,2,depths=(200,)).verdict(OBS,0,False)=='RETAINED_UNRESOLVED'
    bad=Certificate(ROOT_CELL,Bounds(arb('nan'),arb(1)),None,(),'ARB_WHOLE_CELL_V1')
    assert bad.verdict(OBS,0,True)=='RETAINED_UNRESOLVED'
    def invalid(cell,observation,budget,baseline):return 5
    result=SupportEngine(invalid,OBS,0,COARSE_WIDTHS,search_budget(1)).run()
    assert not result.execution_valid and result.metrics()['hull_width_km']==15

def test_coarse_fine_fresh_independent():
    engines=[SupportEngine(compatible,OBS,0,w,search_budget(1)) for w in [COARSE_WIDTHS,FINE_WIDTHS]]
    assert engines[0].queue is not engines[1].queue and engines[0].cache is not engines[1].cache
    assert engines[0].budget is not engines[1].budget
    assert tuple(w/2 for w in COARSE_WIDTHS)==FINE_WIDTHS
    a,b=[e.run() for e in engines]
    assert a.initial_domain==b.initial_domain==ROOT_CELL
    assert a.accounting==b.accounting
    with pytest.raises(RuntimeError):engines[0].run()
    pair=raw_pair(compatible,OBS,0,{'COARSE':1,'FINE':1})
    assert all(r.initial_domain==ROOT_CELL for r in pair)

@pytest.mark.parametrize('key',COUNTERS)
def test_hard_cap_n_n_plus_one(key):
    b=HardBudget(dict.fromkeys(COUNTERS,3),100)
    b.charge(**{key:3})
    with pytest.raises(BudgetExceeded):b.charge(**{key:1})
    assert b.counts[key]==3 and b.attempts[key]==4 and b.blocked_transactions==1

def test_aggregate_cap_and_atomic_transaction():
    b=HardBudget(dict.fromkeys(COUNTERS,10),3)
    b.charge(n_cell_requests=1,n_bound_evaluations=2)
    with pytest.raises(BudgetExceeded):b.charge(n_profile_evaluations=1,n_forward_state_evaluations=1)
    assert b.counts['n_forward_state_evaluations']==0 and sum(b.counts.values())==3

def test_cache_accounting():
    engine=SupportEngine(compatible,OBS,0,COARSE_WIDTHS,search_budget(2))
    a=engine.certificate(ROOT_CELL);b=engine.certificate(ROOT_CELL)
    assert a is b and engine.budget.counts['n_cell_requests']==2 and engine.budget.counts['n_cache_hits']==1
    assert engine.budget.counts['n_bound_evaluations']==1
    with pytest.raises(BudgetExceeded):engine.certificate(ROOT_CELL)

def test_evaluation_runtime_sentinel():
    with estimator_scope(),pytest.raises(RuntimeError):retention(None,[50,0,2,0])
    def leak(cell,observation,budget,baseline):retention(None,[50,0,2,0])
    with pytest.raises(RuntimeError):SupportEngine(leak,OBS,0,COARSE_WIDTHS,search_budget(1)).run()

def test_b0_cz_separation():
    a=SupportEngine(toy,OBS,.25,ROOT_CELL.widths,search_budget(1,True),baseline=True).run()
    b=SupportEngine(toy,OBS,.25,ROOT_CELL.widths,search_budget(1)).run()
    assert a.coverage_closed and a.records[0].status=='RETAINED_COMPATIBLE'
    assert not b.coverage_closed and b.records[0].status=='RETAINED_UNRESOLVED'

def test_adjacency_deterministic_and_no_corner_merge():
    a,b=ROOT_CELL.split(COARSE_WIDTHS)
    assert adjacent(a,b) and components([a,b])==components([b,a])==1
    c=Cell((45,-5,1,-15),(46,-4,2,-14));d=Cell((46,-4,2,-14),(47,-3,3,-13))
    assert not adjacent(c,d) and components([c,d])==2

def test_boundary_censoring():
    result=SupportEngine(compatible,OBS,0,COARSE_WIDTHS,search_budget(1)).run()
    assert result.records[0].flags==('BOUNDARY_CENSORED',)
    assert result.metrics()['boundary_censored']
    c=Cell((46,-4,1.1,-14),(47,-3,1.2,-13));assert not c.boundary_censored

def test_unreleased_development_hard_stop(monkeypatch):
    import r4_cz_envelope_development as d
    def fail(*args,**kwargs):pytest.fail('No filesystem access allowed before hard stop')
    monkeypatch.setattr(Path,'read_bytes',fail)
    assert d.DEVELOPMENT_RELEASED is False
    with pytest.raises(RuntimeError,match='not been released'):d.main()

def test_observation_immutable_and_metadata_validation():
    with pytest.raises(ValueError):OBS.coefficients[0]=1
    with pytest.raises(ValueError):OBS.coefficients.setflags(write=True)
    with pytest.raises(ValueError):feature(np.zeros((3,121)),frequencies=(201,235,284))
    with pytest.raises(ValueError):feature(np.zeros((3,121)),times=np.arange(121))

def test_estimator_source_whitelist():
    banned={'truth','truth_range','truth_state','errors','old_recovered','old_best','oracle','basin'}
    for name in ['r4_cz_envelope_support.py','r4_cz_envelope_model.py','r4_cz_envelope_representation.py']:
        tree=ast.parse((ROOT/name).read_text(encoding='utf-8-sig'))
        identifiers={n.id for n in ast.walk(tree) if isinstance(n,ast.Name)}|{n.arg for n in ast.walk(tree) if isinstance(n,ast.arg)}
        assert not banned&identifiers
        imports=[n.module or '' for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
        assert not any('evaluation' in m or 'r4_a1' in m for m in imports)
        assert not any(isinstance(n,ast.Constant) and isinstance(n.value,str) and ('.npz' in n.value or 'CASE_' in n.value) for n in ast.walk(tree))

@pytest.mark.parametrize('kind',['ROOT','TINY'])
def test_real_modal_enclosure_contains_declared_fixture(model,kind):
    state=tuple(POLICY['fixed_states_r_km_theta_deg_v_mps_psi_deg'][0]);point=Cell(state,state)
    eps=POLICY['non_development_cell_bound_tests']['tiny_width_each_axis']
    cell=ROOT_CELL if kind=='ROOT' else Cell(state,tuple(a+eps for a in state))
    b=ample()
    with precision():
        outer=project_balls(model.levels_arb(cell,150,b))
    with precision(256):inner=project_balls(model.levels_arb(point,150,b))
    assert all(p.contains(q) for p,q in zip(outer,inner))
    # Enclosure tests use endpoint/midpoint checks only as regression, not proof.
    with precision(256):
        for candidate in [cell.low,cell.high,tuple((a+d)/2 for a,d in zip(cell.low,cell.high))]:
            predicted=project_balls(model.levels_arb(Cell(candidate,candidate),150,b))
            assert all(p.contains(q) for p,q in zip(outer,predicted))


def test_bearing_wrap_uncertainty_retains():
    o=Observation(np.full(121,np.pi),np.zeros(12),1.)
    with precision():b=lower_bound_bearing_cost(ROOT_CELL,o,ample())
    assert b.valid() and b.lower==0 and b.upper>0

def test_bound_cap_prevents_geometry_work(model,monkeypatch):
    import r4_cz_envelope_model as module
    def fail(*args,**kwargs):pytest.fail('Geometry must not start before bound reservation')
    monkeypatch.setattr(module,'geometry',fail)
    with pytest.raises(BudgetExceeded):model.certificate(ROOT_CELL,OBS,ample_zero_bound())

def ample_zero_bound():
    b=ample();b.limits['n_bound_evaluations']=0;return b

def test_direct_certificate_entry_sentinel():
    def leak(cell,observation,budget,baseline):retention(None,[50,0,2,0])
    engine=SupportEngine(leak,OBS,0,COARSE_WIDTHS,search_budget(1))
    with pytest.raises(RuntimeError):engine.certificate(ROOT_CELL)

def test_baseline_model_has_no_acoustic_evaluation(model):
    b=search_budget(1,True)
    c=model.certificate(ROOT_CELL,OBS,b,True)
    assert c.feature is None and b.counts['n_forward_state_evaluations']==0 and b.counts['n_profile_evaluations']==0
    assert b.counts['n_bound_evaluations']==1

def test_profile_reservation_is_atomic_before_prediction(model,monkeypatch):
    b=ample();b.limits['n_profile_evaluations']=20
    def fail(*args,**kwargs):pytest.fail('No profile calculation before atomic admission')
    monkeypatch.setattr(model,'_levels_arb_batch',fail)
    with pytest.raises(BudgetExceeded):model.certificate(ROOT_CELL,OBS,b)
    assert b.counts['n_bound_evaluations']==1 and b.counts['n_forward_state_evaluations']==0 and b.counts['n_profile_evaluations']==0
