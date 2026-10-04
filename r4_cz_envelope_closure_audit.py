"""Fixed, preregistered non-development physical closure ladder only."""
from __future__ import annotations
import csv,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from r4_cz_envelope_model import ModalModel,point_geometry
from r4_cz_envelope_representation import feature
from r4_cz_envelope_support import Observation,HardBudget,COUNTERS,COARSE_WIDTHS,FINE_WIDTHS,search_budget
from r4_cz_envelope_support_corrected import CorrectedEngine,structural_analysis
from r4_cz_envelope_evaluation import retention

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN'
SCRATCH=Path('D:/ProjectStorage/WSL-CZ/gate3a-closure-20261004')
POLICY_SHA='98f4e3184752025be73da2ba2bad568bdc0e9992'

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def main():
    policy_path=OUT/'CZ_CLOSURE_CORRECTION_POLICY.json';policy=load(policy_path)
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()!=POLICY_SHA or subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],cwd=ROOT,text=True).split()[0]!=POLICY_SHA:
        raise RuntimeError('Correction policy must be committed and pushed before fixtures')
    old=subprocess.check_output(['git','show',POLICY_SHA+':'+policy_path.relative_to(ROOT).as_posix()],cwd=ROOT)
    if policy_path.read_bytes().replace(b'\r\n',b'\n')!=old.replace(b'\r\n',b'\n'):raise RuntimeError('Correction policy changed')
    for rel,expected in {**policy['input_sha256'],**policy['enclosure_source_sha256']}.items():
        if sha(ROOT/rel)!=expected:raise RuntimeError('Frozen input/source changed: '+rel)
    if (OUT/'CZ_CLOSURE_FIXTURE_SUMMARY.json').exists():raise RuntimeError('Fixture ladder already sealed; no rerun or extra ladder')
    generator=HardBudget(dict.fromkeys(COUNTERS,32),192);model=ModalModel();fixtures=[]
    for f in policy['fixtures']:
        generator.charge(n_validation_evaluations=2)
        levels=model.levels_numpy(f['state'],f['depth_label_m'],generator)
        bearing,_=point_geometry(f['state'])
        observation=Observation(bearing,feature(levels),policy['fixture_bearing_cutoff']['value_radian_SSE'])
        fixtures.append((f['id'],observation,tuple(f['state'])))
    rows=[];details=[];selected={};valid=True
    print(json.dumps({'fixture_ladder_started':True,'fixtures':[f[0] for f in fixtures],'policy':POLICY_SHA}),flush=True)
    for resolution,widths in [('COARSE',COARSE_WIDTHS),('FINE',FINE_WIDTHS)]:
        a=policy['structural_budget_analysis'][resolution]
        if structural_analysis(widths)['structural_minimum_cell_requests']!=a['structural_minimum_cell_requests']:raise RuntimeError('Structural derivation changed')
        for tier,cap in enumerate(a['candidate_cell_cap_ladder'],1):
            passing=True
            for fixture_id,observation,reference in fixtures:
                for mode in ['B0','CZ']:
                    name=f'{resolution}_L{tier}_{fixture_id}_{mode}';start=time.perf_counter();baseline=mode=='B0'
                    print(json.dumps({'fixture_run_start':name,'cell_cap':cap}),flush=True)
                    def progress(count,depth,queue,records):
                        if count%20==0:
                            print(json.dumps({'fixture_run':name,'requests':count,'depth':depth,'queue':queue,'records':records,'elapsed_seconds':round(time.perf_counter()-start,3)}),flush=True)
                    engine=CorrectedEngine(model.certificate,observation,1e-7,widths,search_budget(cap,baseline),baseline=baseline,progress=progress,
                        storage_envelope_bytes=a['storage_record_envelope_bytes_per_raw_run'])
                    result=engine.run()
                    # Immutable result is sealed before the isolated evaluator sees reference.
                    generator.charge(n_validation_evaluations=1);kept=retention(result,reference)
                    diagnostic=result.diagnostics();elapsed=time.perf_counter()-start
                    row={'fixture_id':fixture_id,'resolution':resolution,'mode':mode,'ladder_tier':tier,'hard_cell_cap':cap,
                        'n_cell_requests':diagnostic['n_cell_requests'],'n_forward_state_evaluations':diagnostic['n_forward_state_evaluations'],
                        'n_profile_evaluations':diagnostic['n_profile_evaluations'],'n_bound_evaluations':diagnostic['n_bound_evaluations'],
                        'n_cache_hits':diagnostic['n_cache_hits'],'n_validation_evaluations':diagnostic['n_validation_evaluations'],
                        'n_rejected':diagnostic['n_rejected'],'n_terminal_possible':diagnostic['n_terminal_possible'],'n_certified_compatible':diagnostic['n_certified_compatible'],
                        'n_budget_unresolved':diagnostic['n_budget_unresolved'],'n_invalid_bound':diagnostic['n_invalid_bound'],
                        'maximum_depth':diagnostic['maximum_depth'],'queue_peak':diagnostic['queue_peak'],'estimated_storage_peak_bytes':diagnostic['estimated_storage_peak_bytes'],
                        'retained_range_hull':json.dumps(diagnostic['retained_range_hull']),'domain_partition_closed':result.domain_partition_closed,
                        'search_budget_closed':result.search_budget_closed,'compatibility_certified':result.compatibility_certified,'execution_valid':result.execution_valid,
                        'reference_joint_retained':kept['joint_retained'],'reference_range_retained':kept['range_retained'],'elapsed_seconds':elapsed,'termination_reason':result.termination_reason}
                    rows.append(row)
                    with (OUT/'CZ_CLOSURE_FIXTURE_RESULTS.csv').open('w',encoding='utf-8',newline='') as handle:
                        w=csv.DictWriter(handle,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
                    # Full box partition/certificates are infrastructure evidence; never case data.
                    boxes=[{'low':r.cell.low,'high':r.cell.high,'depth':r.depth,'status':r.status,'reason':r.reason,'flags':r.flags,
                        'bearing_lower':str(r.certificate.bearing.lower) if r.certificate else None,
                        'bearing_upper':str(r.certificate.bearing.upper) if r.certificate else None,
                        'feature_lower':str(r.certificate.feature.lower) if r.certificate and r.certificate.feature else None,
                        'feature_upper':str(r.certificate.feature.upper) if r.certificate and r.certificate.feature else None,
                        'profile_labels':r.certificate.depth_labels if r.certificate else ()} for r in result.records]
                    details.append({'run_id':name,'diagnostics':diagnostic,'sealed_partition':boxes,'reference_evaluation':kept})
                    print(json.dumps({'fixture_run_done':name,'cap':cap,'closed':result.domain_partition_closed,'valid':result.execution_valid,'reference_retained':kept['joint_retained'],'elapsed_seconds':round(elapsed,3),'requests':row['n_cell_requests'],'rejected':row['n_rejected'],'terminal_possible':row['n_terminal_possible']}),flush=True)
                    invalid=not result.execution_valid or not kept['joint_retained'] or not kept['range_retained']
                    if invalid:valid=False;passing=False;break
                    passing=passing and result.domain_partition_closed and result.search_budget_closed
                if not valid:break
            if not valid:break
            if passing:selected[resolution]=cap;break
        if not valid:break
        if resolution not in selected:selected[resolution]=None
    for resolution in ['COARSE','FINE']:selected.setdefault(resolution,None)
    decision=('CZ_ENVELOPE_PRE_RUN_IMPLEMENTATION_INVALID' if not valid else
        'CZ_ENVELOPE_CORRECTED_PRE_RUN_READY_FOR_INDEPENDENT_AUDIT' if all(x is not None for x in selected.values()) else
        'CZ_ENVELOPE_PRE_RUN_SEARCH_CLOSURE_NOT_ESTABLISHED')
    final={}
    for resolution in ['COARSE','FINE']:
        relevant=[r for r in rows if r['resolution']==resolution]
        tier=max((r['ladder_tier'] for r in relevant),default=0)
        latest=[r for r in relevant if r['ladder_tier']==tier]
        final[resolution]={'tier':tier,'cell_cap':latest[0]['hard_cell_cap'] if latest else None,
            'CZ_fixtures_closed':sum(r['domain_partition_closed'] for r in latest if r['mode']=='CZ'),'B0_fixtures_closed':sum(r['domain_partition_closed'] for r in latest if r['mode']=='B0'),
            'fixture_count':len(fixtures),'references_retained':all(r['reference_joint_retained'] and r['reference_range_retained'] for r in latest),'all_run_rows':latest}
    write(OUT/'CZ_CLOSURE_FIXTURE_PARTITIONS.json',{'policy_commit':POLICY_SHA,'scope':'preregistered infrastructure fixtures only, not scientific case support','runs':details})
    write(OUT/'CZ_CLOSURE_FIXTURE_SUMMARY.json',{'stage':policy['stage'],'policy_commit':POLICY_SHA,'policy_sha256':sha(policy_path),'decision':decision,'execution_valid':valid,
        'selected_cell_caps':selected,'final_ladder_diagnostics':final,'reference_retained_all_runs':all(r['reference_joint_retained'] and r['reference_range_retained'] for r in rows),
        'raw_run_count':len(rows),'generator_evaluator_budget':generator.report(),'forbidden_input_used':False,'development_numerical_runs':0,'development_released':False,'fresh':'NOT_GENERATED',
        'source_sha256':{name:sha(ROOT/name) for name in ['r4_cz_envelope_support_corrected.py','r4_cz_envelope_closure_audit.py']},
        'original_enclosure_source_sha256':policy['enclosure_source_sha256'],'measurement_scope':'noise-free matched non-development infrastructure fixtures; not CZ acquisition performance',
        'ladder_selection_applied':policy['ladder_selection'],'no_fourth_budget':True})
    print(json.dumps({'closure_fixture_audit_done':True,'decision':decision,'selected_caps':selected,'development_runs':0}),flush=True)

if __name__=='__main__':main()
