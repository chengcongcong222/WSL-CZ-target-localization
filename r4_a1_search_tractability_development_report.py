'''Render development evidence without modifying the accepted harness or its artifacts.'''
from pathlib import Path
import csv
import hashlib
import json
import subprocess

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_A1_SEARCH_TRACTABILITY_AUDIT'
PARENT='c5c9ba2e1258e73486e62df26499e03cce204f27'
RELEASE='d731cf51d72b41dea7fd4c7a48fcf21a623f7bc4'

def load(name):return json.loads((OUT/name).read_text(encoding='utf-8'))
def write(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
def failure_csv(name,rows,fields):
    with (OUT/name).open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(rows)

def main():
    audit=load('DEVELOPMENT_RECONSTRUCTION_SUMMARY.json')
    native=load('DEVELOPMENT_EXECUTION_DECISION.json')
    tests=load('PYTEST_DEVELOPMENT_POST_RUN.json')
    assert tests['returncode']==0
    for freeze in ['METHOD_FREEZE.json','DEVELOPMENT_EXECUTION_FREEZE.json','CHECKPOINT_EVIDENCE_FREEZE.json']:
        for path,digest in load(freeze)['sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==RELEASE
    alias_forwards=0
    for solver in ['SHGO','DIRECT']:
        with (OUT/(solver+'_WITNESS_CATALOG.csv')).open(newline='',encoding='utf-8') as handle:
            alias_forwards+=sum(row['budget']=='T3' and float(row['J_exact'])<1e-6 for row in csv.DictReader(handle))
    raw=native['DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED']
    dual=native['DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED'] if raw else 'NOT_REACHED'
    decision={
        'parent_harness_commit':PARENT,'release_authorization_commit':RELEASE,
        'scientific_decision':audit['scientific_decision'],
        'DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED':bool(raw and audit['reconstruction_valid']),
        'DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED':dual,
        'dual_solver_agreement':audit['dual_solver_agreement'],
        'recovery_counts':audit['recovery_counts'],'raw_convergence':audit['raw_convergence'],
        'exact_alias':'YES' if audit['exact_alias_finding'] else 'NO_IN_TESTED_CATALOGS',
        'execution_invalid_cases':audit['execution_invalid_cases'],
        'failed_recovery_cases':audit['failed_recovery_cases'],
        'failed_convergence_cases':audit['failed_convergence_cases'],
        'optimizer_and_driver_totals':audit['optimizer_and_driver_totals'],
        'driver_alias_validation_exact_forward_evaluations':alias_forwards,
        'post_run_independent_forward_state_evaluations':audit['independent_unique_case_depth_states_reconstructed'],
        'post_run_independent_forward_batch_calls':audit['independent_modal_batch_calls'],
        'hard_budget_enforcement_failures':audit['hard_budget_enforcement_failures'],
        'mechanical_reconstruction_valid':audit['reconstruction_valid'],
        'mechanical_reconstruction_checks_passed':audit['integrity_checks_passed'],
        'mechanical_reconstruction_checks_failed':audit['integrity_checks_failed'],
        'post_run_tests':tests['output'].strip(),
        'fresh_confirmation':audit['fresh_confirmation'],'fresh_confirmation_released':False,
        'R4_progress_percent':0,'A2_depth_SSP_P5':'UNOPENED',
        'next_action':'STOP; await third independent audit; no further experiment authorized',
        'execution_workspace':str(ROOT),'original_C_workspace_at_execution':'c5c9ba2e1258e73486e62df26499e03cce204f27; C drive full; final synchronization reported separately',
    }
    write('DEVELOPMENT_FINAL_DECISION.json',decision)
    failure_csv('DEVELOPMENT_FAILED_RECOVERY_CASES.csv',audit['failed_recovery_cases'],['solver','case_id','budget'])
    failure_csv('DEVELOPMENT_FAILED_CONVERGENCE_CASES.csv',audit['failed_convergence_cases'],['solver','case_id'])
    failure_csv('DEVELOPMENT_EXECUTION_INVALID_CASES.csv',audit['execution_invalid_cases'],['solver','case_id','budget'])
    totals=audit['optimizer_and_driver_totals']
    rows='\n'.join('| '+solver+' | '+' | '.join(str(audit['recovery_counts'][solver][b])+'/21' for b in ['T1','T2','T3'])+' | '+str(audit['raw_convergence'][solver])+'/21 |' for solver in ['SHGO','DIRECT'])
    report=f'''# Frozen 21-case development result

Decision: `{audit['scientific_decision']}`. R4 remains 0%. Fresh confirmation was not authorized and did not run. A2, depth development, SSP and P5 remain unopened. This stage stops after submission.

Accepted parent harness: `{PARENT}`. Runtime release authorization: `{RELEASE}`. The source release flags remain false. Algorithm, objective, budgets, solver choices, manifest, Gate logic and tolerances remain byte-identical to the accepted harness. No failed case was retried or removed.

The single frozen execution completed 21 noisy cases x two solvers x three independent raw budgets = 126 case-budget runs, each with 21 inherited depth branches = 2646 branches. All SHGO runs preceded all DIRECT runs. Every branch used a fresh engine, controller and cache with admitted-request caps 16/64/256. Recovery means only that the best feasible exact evaluated witness has J < 0.001 dB. These are inherited synthetic matched E0 controls with perfect relative-TL observations and noisy bearings; the findings have that model and panel scope.

| Solver | T1 threshold hits | T2 threshold hits | T3 threshold hits | Raw T2/T3 convergence |
|---|---:|---:|---:|---:|
{rows}

Raw budget convergence Gate: `{decision['DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED']}`. Dual solver agreement: `{audit['dual_solver_agreement']}`. Dual Gate: `{dual}`. The dual comparison is conditional on the raw Gate and cannot receive a scientific PASS from informal partial comparisons.

Exact alias: `{decision['exact_alias']}`. This finding applies only to exported tested catalogs and does not establish global uniqueness. All candidates retain `EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM`; neither optimizer capture basins nor local/global minima are certified. Failure of these capped searches does not establish physical non-identifiability or universal computational intractability.

Optimizer admitted objective requests: {totals['n_objective_requests']}; exact forward evaluations: {totals['n_exact_forward_evaluations']}; unique physical states: {totals['n_unique_states_evaluated']}; cache hits: {totals['n_cache_hits']}; refused requests: {totals['n_blocked_requests']}. Mandatory driver cold validation made {totals['n_validation_exact_forward_evaluations']} further exact forward evaluations outside the optimizer cap. The frozen alias prefilter selected {alias_forwards} witnesses for additional driver exact forward reconstruction. Hard-budget enforcement failures: {audit['hard_budget_enforcement_failures']}. Cap exhaustion is an expected accounting event and does not itself declare search or execution failure.

Mechanical reconstruction passed {audit['integrity_checks_passed']} checks and failed {audit['integrity_checks_failed']}. It reads all 2646 request ledgers, independently reconstructs counters, physical-state hashes, every retained feasible witness, representatives, branch/case aggregates, raw convergence, conditional dual evaluation and strict joint-alias status. Independent geometry and cold direct-modal RMS reconstruction covered {audit['independent_unique_case_depth_states_reconstructed']} unique case/depth/state combinations in {audit['independent_modal_batch_calls']} batched forward calls; maximum J difference = {audit['maximum_J_delta']:.12g} dB and maximum bearing cost difference = {audit['maximum_bearing_cost_delta']:.12g} rad^2. This audit-only memo was created after all searches finished and was never supplied to a solver.

Execution-invalid raw runs: {len(audit['execution_invalid_cases'])}. Failed recovery raw runs: {len(audit['failed_recovery_cases'])}. Failed convergence comparisons: {len(audit['failed_convergence_cases'])}. Their full case lists are preserved in DEVELOPMENT_EXECUTION_INVALID_CASES.csv, DEVELOPMENT_FAILED_RECOVERY_CASES.csv and DEVELOPMENT_FAILED_CONVERGENCE_CASES.csv. Raw results, request ledgers, full evaluated witnesses, representative catalogs, exception chains and numerical failures remain preserved.

Post-run unchanged harness tests: {tests['output'].strip()}

Actual execution workspace: `{ROOT}`. The original C drive had no free space, so execution used a verified independent D-drive checkout with all frozen identities preserved. The original C checkout requires a separate disk-space check before synchronization. Historical accepted reports and manifests remain intact; this report and DEVELOPMENT_FINAL_DECISION.json are new evidence for the third independent audit.
'''
    (OUT/'DEVELOPMENT_RUN_REPORT.md').write_text(report,encoding='utf-8')
    (OUT/'DEVELOPMENT_LOCAL_VALIDATION.md').write_text(f'''# Development local validation

Single frozen experiment: 126 raw runs / 2646 branches. No rerun and no algorithm repair.

Mechanical reconstruction: {audit['integrity_checks_passed']} passed / {audit['integrity_checks_failed']} failed. See DEVELOPMENT_RECONSTRUCTION_AUDIT.csv and DEVELOPMENT_SCORE_RECONSTRUCTION.csv.

Pre-release tests: 74 passed after creating the missing runtime_tmp parent directory. The earlier setup attempt had 67 passed / 7 WinError 3 setup errors before any noisy run; no frozen algorithm was changed. See DEVELOPMENT_ENVIRONMENT_NOTES.json. Post-run unchanged tests:\n\n{tests['output']}

Frozen METHOD, EXECUTION and accepted CHECKPOINT hashes revalidated after execution. New evidence receives its own manifest; historical evidence remains unchanged. No fresh experiment was released.
''',encoding='utf-8')
    write('DEVELOPMENT_GPT_SYNC.json',decision)
    master=ROOT/'results/R4_MASTER';progress=json.loads((master/'R4_PROGRESS.json').read_text(encoding='utf-8'))
    progress.update(current_stage='A1_SEARCH_TRACTABILITY_DEVELOPMENT_COMPLETE',independent_audit_status='FIX2_ACCEPTED_AT_7ab2484; PRE_RUN_ACCEPTED_AT_b4839b2; HARNESS_ACCEPTED_AT_c5c9ba2; DEVELOPMENT_PENDING_THIRD_AUDIT',next_recommended_stage=decision['next_action'],tractability_checkpoint='DEVELOPMENT_EXECUTION_HARNESS_FREEZE_ACCEPTED',tractability_scientific_decision=audit['scientific_decision'],tractability_noisy_development_runs=21,tractability_raw_case_budget_runs=126,tractability_depth_branch_runs=2646,development_released=True,fresh_confirmation_released=False,tractability_development_gates={'raw_budget_convergence':decision['DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED'],'dual_solver_agreement':dual},overall_progress_percent=0)
    (master/'R4_PROGRESS.json').write_text(json.dumps(progress,indent=2)+'\n',encoding='utf-8')
    plan=(master/'R4_PLAN.md').read_text(encoding='utf-8')
    old='Only the A1 search-tractability DEVELOPMENT_EXECUTION_HARNESS_FREEZE checkpoint is currently authorized. A2/B/SSP/P5 are not opened. R4 remains at 0%.'
    plan=plan.replace(old,'The authorized frozen 21-case noisy development stage is complete. No further experiment is authorized. A2/B/SSP/P5 are not opened. R4 remains at 0%.')
    plan=plan.replace('The new harness remains development_released=false; noisy_development_runs=0. It awaits a second independent audit.','At the harness-freeze checkpoint, development_released=false and noisy_development_runs=0; the second audit and separate development release were pending.')
    plan+=f'''\n## Accepted harness and completed frozen development\n\nThe lead accepted `{PARENT}` and authorized only the 21-case noisy development, recorded before search at `{RELEASE}`. The unchanged runtime-released harness executed all 126 raw runs / 2646 depth branches. Its decision is `{audit['scientific_decision']}`; raw Gate = {raw}; dual Gate = {dual}. All failures and raw request/witness evidence are preserved under ../R4_A1_SEARCH_TRACTABILITY_AUDIT/. See DEVELOPMENT_RUN_REPORT.md and DEVELOPMENT_FINAL_DECISION.json. This stage stops pending the third independent audit; no T4, budget/tolerance change, solver substitution, FIX3/FIX4 or fresh confirmation is authorized. R4 remains 0%; A2/depth/SSP/P5 remain unopened.\n'''
    (master/'R4_PLAN.md').write_text(plan,encoding='utf-8')
    ledger=master/'R4_EVIDENCE_LEDGER.csv'
    with ledger.open(newline='',encoding='utf-8') as handle:records=list(csv.DictReader(handle))
    for row in records:
        if row['stage']=='A1-SEARCH-TRACTABILITY-EXECUTION-HARNESS':row.update(status='DEVELOPMENT_EXECUTION_HARNESS_FREEZE_ACCEPTED',independent_audit_status='ACCEPTED_AT_c5c9ba2; DEVELOPMENT_RELEASE_AT_d731cf5')
    records.append(dict(stage='A1-SEARCH-TRACTABILITY-DEVELOPMENT',artifact='../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_RUN_REPORT.md',status=audit['scientific_decision'],evidence_scope='21 noisy cases; 126 independent raw runs; 2646 branches; full request/witness evidence; mechanical reconstruction; fresh never released',independent_audit_status='PENDING_THIRD_INDEPENDENT_AUDIT; STOP; NO_SCIENTIFIC_GATE_CREDIT',credited_weight_percent=0))
    with ledger.open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=['stage','artifact','status','evidence_scope','independent_audit_status','credited_weight_percent']);writer.writeheader();writer.writerows(records)
    accepted_paths={p for name in ['METHOD_FREEZE.json','DEVELOPMENT_EXECUTION_FREEZE.json','CHECKPOINT_EVIDENCE_FREEZE.json'] for p in load(name)['sha256']}
    files=[p for p in OUT.rglob('*') if p.is_file() and p.relative_to(ROOT).as_posix() not in accepted_paths and p.name!='DEVELOPMENT_RESULTS_EVIDENCE_FREEZE.json']
    portable=[ROOT/'r4_a1_search_tractability_development_reconstruction.py',Path(__file__),master/'R4_PROGRESS.json',master/'R4_PLAN.md',ledger]
    for source,name in [(portable[0],'DEVELOPMENT_RECONSTRUCTION_SOURCE_SNAPSHOT.py'),(portable[1],'DEVELOPMENT_REPORT_SOURCE_SNAPSHOT.py')]:
        snapshot=OUT/name;snapshot.write_bytes(source.read_bytes());files.append(snapshot)
    text_artifacts={}
    for source in portable:
        data=source.read_bytes();canonical=data.replace(b'\r\n',b'\n')
        text_artifacts[source.relative_to(ROOT).as_posix()]={'sha256_at_execution':hashlib.sha256(data).hexdigest(),'canonical_lf_sha256':hashlib.sha256(canonical).hexdigest(),'canonical_git_blob_oid':hashlib.sha1(b'blob '+str(len(canonical)).encode()+b'\0'+canonical).hexdigest(),'allowed_transport_change':'CRLF_LF_ONLY'}
    freeze={'parent_harness_commit':PARENT,'release_authorization_commit':RELEASE,'stage':'FROZEN_21_CASE_DEVELOPMENT','fresh_confirmation_released':False,'text_artifacts':text_artifacts,'sha256':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(files))}}
    write('DEVELOPMENT_RESULTS_EVIDENCE_FREEZE.json',freeze)
    print(json.dumps({k:v for k,v in decision.items() if k not in ['failed_recovery_cases','failed_convergence_cases']},indent=2))

if __name__=='__main__':main()
