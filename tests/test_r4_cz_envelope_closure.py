"""Non-development closure semantics and structure tests; original 47 tests retained."""
import ast,json,math
from pathlib import Path
import numpy as np
import pytest
from flint import arb
from r4_cz_envelope_support import (Cell,Bounds,Certificate,Observation,DEPTH_LABELS,DOMAIN_LOW,DOMAIN_HIGH,COARSE_WIDTHS,FINE_WIDTHS,COUNTERS,HardBudget,search_budget,BudgetExceeded)
from r4_cz_envelope_support_corrected import CorrectedEngine,structural_analysis,raw_pair
from r4_cz_envelope_evaluation import retention

ROOT=Path(__file__).resolve().parents[1]
POLICY=json.loads((ROOT/'results/R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN/CZ_CLOSURE_CORRECTION_POLICY.json').read_text(encoding='utf-8-sig'))
OBS=Observation(np.zeros(121),np.zeros(12),1.)
ROOT_CELL=Cell(DOMAIN_LOW,DOMAIN_HIGH)
ZERO=Bounds(arb(0),arb(0));LOOSE=Bounds(arb(0),arb(2))

def unknown(cell,observation,budget,baseline):
    budget.charge(n_bound_evaluations=1)
    return Certificate(cell,ZERO,LOOSE,DEPTH_LABELS,'ANALYTIC_TOY_INTERVAL_V1') if not baseline else Certificate(cell,LOOSE,None,(),'ANALYTIC_TOY_INTERVAL_V1')

def single_high_path(cell,observation,budget,baseline):
    # Exact endpoint feasibility x=DOMAIN_HIGH, not a sampled rejection.
    # If any upper face is strictly below that endpoint, its distance is positive.
    budget.charge(n_bound_evaluations=1)
    rejected=any(b<endpoint for b,endpoint in zip(cell.high,DOMAIN_HIGH))
    bearing=Bounds(arb(2),arb(2)) if rejected else LOOSE
    return Certificate(cell,bearing,LOOSE,DEPTH_LABELS,'ANALYTIC_TOY_INTERVAL_V1')

@pytest.mark.parametrize('label,widths,splits,minimum',[('COARSE',COARSE_WIDTHS,(6,4,4,4),37),('FINE',FINE_WIDTHS,(7,5,5,5),45)])
def test_structural_derivation_and_frozen_candidates(label,widths,splits,minimum):
    a=structural_analysis(widths)
    assert a['axis_splits']==splits and a['structural_minimum_cell_requests']==minimum
    assert a['path_split_depth']==sum(splits)
    assert POLICY['structural_budget_analysis'][label]['structural_minimum_cell_requests']==minimum
    assert all(n>=minimum for n in POLICY['structural_budget_analysis'][label]['candidate_cell_cap_ladder'])

@pytest.mark.parametrize('widths',[COARSE_WIDTHS,FINE_WIDTHS])
def test_t1_below_structural_minimum_cannot_close(widths):
    minimum=structural_analysis(widths)['structural_minimum_cell_requests']
    result=CorrectedEngine(single_high_path,OBS,1e-7,widths,search_budget(minimum-1)).run()
    assert result.execution_valid and not result.domain_partition_closed and not result.search_budget_closed
    assert result.termination_reason=='HARD_CAP'
    assert result.accounting['n_cell_requests']==minimum-1
    assert retention(result,DOMAIN_HIGH)['joint_retained']

@pytest.mark.parametrize('widths',[COARSE_WIDTHS,FINE_WIDTHS])
def test_t2_exact_minimum_finishes_single_path(widths):
    minimum=structural_analysis(widths)['structural_minimum_cell_requests']
    result=CorrectedEngine(single_high_path,OBS,1e-7,widths,search_budget(minimum)).run()
    assert result.domain_partition_closed and result.search_budget_closed and result.execution_valid
    assert result.accounting['n_cell_requests']==minimum and not result.compatibility_certified
    assert result.maximum_depth==structural_analysis(widths)['path_split_depth']
    assert retention(result,DOMAIN_HIGH)=={'joint_retained':True,'range_retained':True}
    assert result.diagnostics()['n_terminal_possible']==1
    assert all(all(w<=t for w,t in zip(c.widths,widths)) for c in result.scientific_cells)

def test_t3_terminal_possible_closes_without_upper_certificate():
    result=CorrectedEngine(unknown,OBS,1e-7,ROOT_CELL.widths,search_budget(1)).run()
    assert result.domain_partition_closed and result.search_budget_closed and not result.compatibility_certified
    assert result.records[0].status=='RETAINED_POSSIBLE_AT_TERMINAL'
    assert result.scientific_metrics()['hull_width_km']==15

def test_t4_budget_unresolved_blocks_gate_and_retains_whole_domain():
    result=CorrectedEngine(unknown,OBS,1e-7,COARSE_WIDTHS,search_budget(1)).run()
    assert not result.domain_partition_closed and not result.search_budget_closed
    assert all(r.status=='RETAINED_UNRESOLVED_BUDGET' for r in result.records)
    assert sum(math.prod(c.widths) for c in result.retained)==math.prod(ROOT_CELL.widths)
    with pytest.raises(RuntimeError):result.scientific_metrics()

@pytest.mark.parametrize('baseline',[False,True])
def test_b0_and_cz_use_same_outer_semantics(baseline):
    result=CorrectedEngine(unknown,OBS,1e-7,ROOT_CELL.widths,search_budget(1,baseline),baseline=baseline).run()
    assert result.domain_partition_closed and result.records[0].status=='RETAINED_POSSIBLE_AT_TERMINAL'

def test_invalid_bound_keeps_domain_and_invalidates_execution():
    def invalid(cell,observation,budget,baseline):
        return Certificate(cell,Bounds(arb('nan'),arb(1)),LOOSE,DEPTH_LABELS,'ARB_WHOLE_CELL_V1')
    result=CorrectedEngine(invalid,OBS,1e-7,COARSE_WIDTHS,search_budget(1)).run()
    assert not result.execution_valid and not result.domain_partition_closed and not result.search_budget_closed
    assert result.records[0].status=='RETAINED_INVALID_BOUND'
    assert result.diagnostics()['retained_range_hull']==(45,60)

def test_certificate_sampling_never_allows_rejection():
    def sample(cell,observation,budget,baseline):
        return Certificate(cell,Bounds(arb(2),arb(2)),LOOSE,DEPTH_LABELS,'CENTER_SAMPLE')
    result=CorrectedEngine(sample,OBS,1e-7,COARSE_WIDTHS,search_budget(1)).run()
    assert not result.execution_valid and result.retained==(ROOT_CELL,)

def test_certified_compatible_diagnostic_is_separate():
    def full(cell,observation,budget,baseline):
        budget.charge(n_bound_evaluations=1);return Certificate(cell,ZERO,ZERO,DEPTH_LABELS,'ANALYTIC_TOY_INTERVAL_V1')
    result=CorrectedEngine(full,OBS,1e-7,ROOT_CELL.widths,search_budget(1)).run()
    assert result.compatibility_certified and result.domain_partition_closed
    assert result.records[0].status=='RETAINED_CERTIFIED_COMPATIBLE'

def test_empty_outer_partition_is_not_scientific_success():
    def reject(cell,observation,budget,baseline):
        budget.charge(n_bound_evaluations=1);return Certificate(cell,Bounds(arb(2),arb(2)),LOOSE,DEPTH_LABELS,'ANALYTIC_TOY_INTERVAL_V1')
    result=CorrectedEngine(reject,OBS,1e-7,COARSE_WIDTHS,search_budget(1)).run()
    assert result.domain_partition_closed and result.scientific_metrics()['empty_support']
    assert result.scientific_metrics()['hull_width_km'] is None

def test_corrected_runtime_sentinel_and_input_exclusion():
    def leak(cell,observation,budget,baseline):retention(None,DOMAIN_HIGH)
    engine=CorrectedEngine(leak,OBS,1e-7,COARSE_WIDTHS,search_budget(1))
    with pytest.raises(RuntimeError):engine.run()
    tree=ast.parse((ROOT/'r4_cz_envelope_support_corrected.py').read_text(encoding='utf-8-sig'))
    banned={'truth','truth_range','truth_state','errors','old_recovered','old_best','oracle','basin'}
    names={n.id for n in ast.walk(tree) if isinstance(n,ast.Name)}|{n.arg for n in ast.walk(tree) if isinstance(n,ast.arg)}
    assert not names&banned
    assert not any(isinstance(n,ast.ImportFrom) and ('evaluation' in (n.module or '') or 'r4_a1' in (n.module or '')) for n in ast.walk(tree))

def test_corrected_raw_independence():
    engines=[CorrectedEngine(single_high_path,OBS,1e-7,w,search_budget(structural_analysis(w)['structural_minimum_cell_requests'])) for w in [COARSE_WIDTHS,FINE_WIDTHS]]
    a,b=[e.run() for e in engines]
    assert engines[0].cache is not engines[1].cache and engines[0].queue is not engines[1].queue
    assert a.initial_domain==b.initial_domain==ROOT_CELL
    assert a.domain_partition_closed and b.domain_partition_closed
    with pytest.raises(RuntimeError):engines[0].run()

@pytest.mark.parametrize('key',COUNTERS)
def test_inherited_strict_counter_limits(key):
    b=HardBudget(dict.fromkeys(COUNTERS,2),100);b.charge(**{key:2})
    with pytest.raises(BudgetExceeded):b.charge(**{key:1})
    assert b.counts[key]==2 and b.blocked_transactions==1
