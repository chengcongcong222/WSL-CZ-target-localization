"""Read-only historical verification and pre-run sealing, no numeric observation reads."""
from __future__ import annotations
import argparse,ast,concurrent.futures,csv,hashlib,json,math,re,subprocess
from decimal import Decimal,ROUND_CEILING
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN'
HIST=ROOT/'results/R4_A1_SEARCH_TRACTABILITY_AUDIT'
SCRATCH=Path('D:/ProjectStorage/WSL-CZ/gate3a-20261004')
PARENT='e56ba746ebe7370d1dee3d0a181bc348b5fecc5f'
POLICY_SHA='a393c7b45c92c42fd148d9c162faeb5ea35160c7'
CORE=['r4_cz_envelope_representation.py','r4_cz_envelope_support.py','r4_cz_envelope_model.py']
SOURCES=CORE+['r4_cz_envelope_guard.py','r4_cz_envelope_evaluation.py','r4_cz_envelope_development.py','r4_cz_envelope_prerun_audit.py','r4_cz_envelope_prerun_integrity.py','tests/test_r4_cz_envelope.py']
MASTER={'results/R4_MASTER/R4_PLAN.md','results/R4_MASTER/R4_EVIDENCE_LEDGER.csv'}
DECISION='CZ_ENVELOPE_PRE_RUN_FREEZE_READY_FOR_INDEPENDENT_AUDIT'

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def rows(p):return list(csv.DictReader(p.open(encoding='utf-8-sig',newline='')))
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):h.update(b)
    return h.hexdigest()
def digest(b):return hashlib.sha256(b).hexdigest()
def canonical(b):return b.replace(b'\r\n',b'\n')
def oid(b):return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)
def text_identity(p):
    b=p.read_bytes();c=canonical(b)
    return {'raw_sha256':digest(b),'canonical_lf_sha256':digest(c),'canonical_git_blob_oid':oid(c),'bytes':len(b),'transport':'CRLF_LF_ONLY'}
def write(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def seal():
    checks=[]
    def check(ok,label,expected='',observed='',scope='infrastructure'):
        checks.append({'check_id':label,'status':'PASS' if ok else 'FAIL','expected':str(expected),'observed':str(observed),'scope':scope})
        if not ok:raise AssertionError(label)
    before=load(SCRATCH/'protected-before.json')['files']
    check(git('rev-parse','HEAD').decode().strip()==POLICY_SHA,'policy_commit_is_local_HEAD',POLICY_SHA)
    check(git('ls-remote','origin','refs/heads/main').decode().split()[0]==POLICY_SHA,'policy_commit_was_pushed',POLICY_SHA)
    check(len(before)==7194,'protected_inventory_count',7194,len(before))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        actual=dict(zip(before,pool.map(lambda rel:{'sha256':sha(ROOT/rel),'bytes':(ROOT/rel).stat().st_size},before)))
    for rel,expected in before.items():check(actual[rel]==expected,'protected:'+rel,expected['sha256'],actual[rel]['sha256'],'preexisting raw bytes')
    write(SCRATCH/'protected-after.json',{'files':actual})
    total_raw=0
    for name in ['METHOD_FREEZE.json','DEVELOPMENT_EXECUTION_FREEZE.json','CHECKPOINT_EVIDENCE_FREEZE.json','DEVELOPMENT_RESULTS_EVIDENCE_FREEZE.json']:
        bindings=load(HIST/name)['sha256']
        for rel,expected in bindings.items():
            check(actual[rel]['sha256']==expected,'historical_raw:'+name+':'+rel,expected,actual[rel]['sha256'],'accepted scientific evidence')
        total_raw+=len(bindings)
    check(total_raw==2796,'raw_historical_binding_count',2796,total_raw)
    old=rows(HIST/'FROZEN_INPUT_HASHES.csv');raw=lf=0
    for r in old:
        rel=r['path'];ok=actual[rel]['sha256']==r['sha256']
        if ok:raw+=1
        else:
            b=canonical((ROOT/rel).read_bytes());ok=r['allow_line_endings']=='True' and digest(b)==r['canonical_lf_sha256'] and oid(b)==r['baseline_git_blob_oid'];lf+=int(ok)
        check(ok,'historical_input:'+rel,r['sha256'],actual[rel]['sha256'],'raw or prior approved LF transport')
    check(raw+lf==2989,'historical_input_count',2989,raw+lf)
    prior='9cfa99c14bb80d08e6f16210ebf89fd9c541b850'
    for rel,identity in load(HIST/'DEVELOPMENT_RESULTS_EVIDENCE_FREEZE.json')['text_artifacts'].items():
        b=canonical(git('show',prior+':'+rel) if rel in MASTER else (ROOT/rel).read_bytes())
        check(digest(b)==identity['canonical_lf_sha256'] and oid(b)==identity['canonical_git_blob_oid'],'portable_prior:'+rel,identity['canonical_lf_sha256'],digest(b))
    for rel,oldname in [('results/R4_MASTER/R4_PLAN.md','plan-before.md'),('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv','ledger-before.csv')]:
        check((ROOT/rel).read_bytes().startswith((SCRATCH/oldname).read_bytes()),'master_append_only:'+rel)
        check(canonical((SCRATCH/oldname).read_bytes())==canonical(git('show',PARENT+':'+rel)),'master_old_prefix_accepted_parent:'+rel)
    lr=rows(ROOT/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv')
    check(len(rows(SCRATCH/'ledger-before.csv'))==29 and len(lr)==33 and all(float(r['credited_weight_percent'])==0 for r in lr),'ledger_29_old_plus_4_zero_credit')
    progress=load(ROOT/'results/R4_MASTER/R4_PROGRESS.json');check(progress['overall_progress_percent']==0 and sum(progress['weights_percent'].values())==100,'R4_progress_and_weights_unchanged')
    changed=set(git('diff','--name-only',POLICY_SHA).decode().splitlines());check(changed==MASTER,'only_two_preexisting_management_files_changed',sorted(MASTER),sorted(changed))
    check(not git('diff','--cached','--name-only').strip(),'staging_clean_before_seal')
    policy=load(OUT/'CZ_NUMERICAL_CALIBRATION_POLICY.json');cal=load(OUT/'CZ_NUMERICAL_CALIBRATION_SUMMARY.json')
    rel=(OUT/'CZ_NUMERICAL_CALIBRATION_POLICY.json').relative_to(ROOT).as_posix()
    check(canonical((ROOT/rel).read_bytes())==canonical(git('show',POLICY_SHA+':'+rel)),'pushed_policy_bytes_immutable')
    for rel,expected in policy['input_sha256'].items():check(sha(ROOT/rel)==expected,'policy_input:'+rel,expected,sha(ROOT/rel))
    for rel,expected in cal['source_sha256'].items():check(sha(ROOT/rel)==expected,'calibration_final_code:'+rel,expected,sha(ROOT/rel))
    check(cal['policy_commit']==POLICY_SHA and cal['policy_raw_sha256']==sha(OUT/'CZ_NUMERICAL_CALIBRATION_POLICY.json'),'calibration_policy_binding')
    data=rows(OUT/'CZ_NUMERICAL_CALIBRATION.csv');check(len(data)==15,'fixed_calibration_fixture_count',15,len(data))
    state_set={json.dumps(x) for x in policy['fixed_states_r_km_theta_deg_v_mps_psi_deg']}
    check(all(r['state'] in state_set and int(r['depth_label_m']) in policy['fixed_calibration_depth_labels_m'] for r in data[:12]),'fixed_12_corner_depths_only')
    for i,(state,z) in enumerate((state,z) for state in policy['fixed_states_r_km_theta_deg_v_mps_psi_deg'] for z in policy['fixed_calibration_depth_labels_m']):
        check(json.loads(data[i]['state'])==state and int(data[i]['depth_label_m'])==z,'policy_fixture_order:'+str(i))
    delta=max(float(r[k]) for r in data for k in ['repeat_dB','numpy_vs_arb256_dB','arb160_vs_arb256_dB','offset_invariance_dB'])
    a=policy['numerical_allowance'];quantum=Decimal(str(a['round_up_quantum_dB']))
    tau=float((max(Decimal(str(a['minimum_floor_dB'])),Decimal(a['safety_factor'])*Decimal.from_float(delta))/quantum).to_integral_value(rounding=ROUND_CEILING)*quantum)
    check(delta==cal['delta_num_dB'] and tau==cal['tau_F_dB'] and delta<=a['maximum_allowed_delta_num_dB'],'preregistered_tau_formula',tau,cal['tau_F_dB'])
    rt=policy['runtime_calibration'];cost=max(1e-3,max(x['elapsed_seconds'] for x in cal['enclosure_fixtures']))
    caps={label:max(1,min(rt['structural_cell_ceiling'][label],math.floor(rt['safety_time_fraction']*rt[label.lower()+'_runtime_envelope_seconds']/cost))) for label in ['COARSE','FINE']}
    check(cal['runtime_cell_cost_seconds']==cost and cal['raw_cell_caps']==caps,'preregistered_runtime_caps',caps,cal['raw_cell_caps'])
    check(all(x['depth_profiles']==21 for x in cal['enclosure_fixtures']) and len(cal['depth_mapping'])==63,'calibration_full_profile_mapping')
    representation=load(OUT/'CZ_REPRESENTATION_FREEZE.json');threshold=load(OUT/'CZ_THRESHOLD_FREEZE.json');grid=load(OUT/'CZ_GRID_FREEZE.json');budget=load(OUT/'CZ_BUDGET_FREEZE.json');tests=load(OUT/'CZ_PRE_RUN_TESTS.json');metadata=load(OUT/'CZ_INPUT_METADATA.json')
    check(representation['id']=='WINDOW_SHAPE_Q2' and representation['dimension']==12 and representation['weights']=='equal three lines and two windows' and representation['tau_F_dB']==tau and representation['sensitivity']=='NONE','single_frozen_representation')
    check(representation['nuisance']=='N0_STATIC_PER_LINE_PER_WINDOW_OFFSET' and representation['N1']=='NOT_RUN' and representation['N2']=='SYMBOLIC_ONLY_NOT_RUN','nuisance_scope')
    check(threshold['tier']=='NOMINAL' and threshold['hull_width_max_km']==5 and threshold['C_r_max']==1/3 and threshold['range_component_count_max']==3 and threshold['joint_broad_component_count_max']==6 and threshold['coarse_fine_endpoint_change_max_km']==threshold['coarse_fine_width_change_max_km']==.5,'lead_nominal_threshold_values')
    check(threshold['joint_reference_retention_fraction']==threshold['range_reference_retention_fraction']==1 and threshold['all_case_pass_required'] and threshold['required_primary_case_count']==21 and threshold['coverage_closed'] and threshold['same_case_B0_additional_contraction_required'],'all_case_retention_coverage_and_B0')
    check(grid['COARSE_max_widths']==policy['grid_rule']['COARSE_terminal_max_widths'] and grid['FINE_max_widths']==policy['grid_rule']['FINE_terminal_max_widths'] and grid['FINE_max_widths']==[x/2 for x in grid['COARSE_max_widths']],'fixed_raw_two_level_resolution')
    check(grid['domain_low']==[45,-5,1,-15] and grid['domain_high']==[60,5,3,15] and grid['profile_depth_labels_m']==list(range(150,251,5)),'inherited_domain_and_profile')
    check(budget['cell_caps']==caps and budget['runtime_policy']==rt and not budget['cap_tuning_from_development'],'strict_budget_freeze_and_no_retuning')
    for mode in ['B0','CZ']:
        for label,n in caps.items():
            b=budget['strict_operation_budgets'][mode][label];l=b['limits']
            check(l['n_cell_requests']==l['n_cache_hits']==n and l['n_validation_evaluations']==0 and l['n_bound_evaluations']==n*(1 if mode=='B0' else 22) and l['n_forward_state_evaluations']==l['n_profile_evaluations']==(0 if mode=='B0' else 21*n) and b['hard_budget_limit']==n*(3 if mode=='B0' else 66),'hard_budget_units:'+mode+':'+label)
    check(tests['tests']==tests['passed']==47 and tests['failed']==tests['errors']==tests['skipped']==0 and len(tests['required_categories_passed'])==20,'47_tests_20_required_categories')
    check(tests['junit_sha256']==sha(OUT/'CZ_PRE_RUN_TESTS.xml'),'junit_identity')
    check(metadata['primary_count']==len(set(metadata['primary_case_ids']))==21 and metadata['case_numerical_feature_search_runs']==metadata['development_observation_payload_reads']==0 and all(not x['numeric_payload_read'] for x in metadata['archive_member_headers']),'21_ids_headers_only_no_numeric_observations')
    check(metadata['manifest_sha256']==sha(ROOT/'results/R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_CASE_MANIFEST.csv') and metadata['archive_sha256']==sha(ROOT/'results/R4_A1_FIX_CONTINUOUS_SEARCH/CASE_OBSERVATIONS.npz'),'archive_and_manifest_binding')
    check(cal['development_numerical_runs']==0 and not cal['development_observations_loaded'] and metadata['fresh']=='NOT_GENERATED' and not tests['development_data_used'],'no_development_or_fresh')
    forbidden={'truth','truth_range','truth_state','errors','old_recovered','old_best','oracle','basin'}
    for name in CORE:
        tree=ast.parse((ROOT/name).read_text(encoding='utf-8-sig'))
        ids={n.id for n in ast.walk(tree) if isinstance(n,ast.Name)}|{n.arg for n in ast.walk(tree) if isinstance(n,ast.arg)}
        imports=[n.module or '' for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
        constants=[n.value for n in ast.walk(tree) if isinstance(n,ast.Constant) and isinstance(n.value,str)]
        check(not forbidden&ids and not any('evaluation' in x or 'r4_a1' in x for x in imports) and not any('.npz' in x or 'CASE_' in x for x in constants),'estimator_AST_exclusion:'+name)
    check('DEVELOPMENT_RELEASED = False' in (ROOT/'r4_cz_envelope_development.py').read_text(encoding='utf-8-sig'),'development_hard_stop_source')
    for name in SOURCES:
        source=(ROOT/name).read_text(encoding='utf-8-sig');ast.parse(source);check(chr(65533) not in source,'source_utf8_syntax:'+name)
    for p in list(OUT.glob('*.md'))+[ROOT/'results/R4_MASTER/R4_PLAN.md']:
        txt=p.read_text(encoding='utf-8-sig');check(chr(65533) not in txt,'document_utf8:'+p.name)
        for target in re.findall(r'\]\(([^)]+)\)',txt):
            if target.startswith(('http:','https:','#')):continue
            dest=p.parent/target.split('#')[0]
            check(dest.exists() or dest.name in ['CZ_PRE_RUN_DECISION.json','CZ_PRE_RUN_DESIGN_FREEZE.json'],'relative_link:'+p.name+':'+target)
    subprocess.run(['git','diff','--check'],cwd=ROOT,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    check(True,'git_diff_check')
    inventory_digest=digest(json.dumps(before,sort_keys=True,separators=(',',':')).encode())
    csvpath=OUT/'CZ_PRE_RUN_INTEGRITY.csv'
    with csvpath.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['check_id','status','expected','observed','scope']);w.writeheader();w.writerows(checks)
    write(OUT/'CZ_PRE_RUN_INTEGRITY_SUMMARY.json',{'status':'PASS','checks':len(checks),'failed':0,'protected_files':7194,'protected_files_changed':[],
        'protected_inventory_sha256':inventory_digest,'historical_raw_bindings':total_raw,'historical_input_bindings':2989,'raw_input_matches':raw,'LF_only_input_matches':lf,
        'prior_portable_text_bindings':5,'master_old_ledger_rows':29,'master_new_rows':4,'progress_percent':0,'development_numerical_runs':0,
        'CSV_sha256':sha(csvpath),'metadata_scope':'ZIP/NPY header shapes only; hashes are byte identities, not numerical payload processing',
        'storage':'canonical D workspace; D operational files, dependencies, temp and caches','passed_check_ids_sha256':digest('\n'.join(x['check_id'] for x in checks).encode())})
    stop_rules=['Stop after PRE-RUN commit/push; wait independent audit','No 21-case feature/search before separate lead release','No sensitivity, N1, third optimizer, extra budget or third grid','Retain unresolved/loose/capped domain; closure false','If enclosure or implementation invalid, stop; no point-grid fallback','No fresh, A2/depth development/SSP/P5; R4=0%']
    write(OUT/'CZ_PRE_RUN_DECISION.json',{'stage':'R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN_FREEZE','parent_gate3_commit':PARENT,'calibration_policy_commit':POLICY_SHA,
        'decision':DECISION,'whole_cell_enclosure':'ESTABLISHED','enclosure_scope':'mathematically conservative fixed modal model, not useful contraction or development coverage/convergence',
        'truth_exclusion':'ESTABLISHED','primary_representation':'WINDOW_SHAPE_Q2','sensitivity':'NONE','primary_nuisance':'N0_STATIC_PER_LINE_PER_WINDOW_OFFSET','threshold_tier':'NOMINAL',
        'tau_F_dB':tau,'raw_cell_caps':caps,'unit_tests_passed':47,'integrity_checks_passed':len(checks),'development_numerical_feature_search_runs':0,'development_released':False,
        'fresh':'NOT_GENERATED','N1':'NOT_RUN','N2':'SYMBOLIC_ONLY_NOT_RUN','scientific_gate_result':'NOT_RUN','R4_progress_percent':0,'A2_depth_SSP_P5':'UNOPENED',
        'stop_rules':stop_rules,'stop_after_commit_push':True,'current_real_receiver_M1':'UNKNOWN_NOT_FROZEN','real_UUV_evidence':'NOT_VERIFIED'})
    bindings={name:text_identity(ROOT/name) for name in SOURCES}
    for p in sorted(OUT.iterdir()):
        if p.is_file() and p.name!='CZ_PRE_RUN_DESIGN_FREEZE.json':bindings[p.relative_to(ROOT).as_posix()]=text_identity(p)
    for rel in MASTER:bindings[rel]=text_identity(ROOT/rel)
    write(OUT/'CZ_PRE_RUN_DESIGN_FREEZE.json',{'stage':'R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN_FREEZE','parent_gate3_commit':PARENT,'calibration_policy_commit':POLICY_SHA,
        'immutable_after_commit':True,'sha256':bindings,'frozen_model_and_observation_input_sha256':policy['input_sha256'],
        'protected_preexisting_inventory_sha256':inventory_digest,'depth_labels_m':list(range(150,251,5)),
        'primary_representation':'WINDOW_SHAPE_Q2','primary_nuisance':'N0_STATIC_PER_LINE_PER_WINDOW_OFFSET','sensitivity':'NONE','threshold_tier':'NOMINAL','tau_F_dB':tau,
        'COARSE_max_widths':grid['COARSE_max_widths'],'FINE_max_widths':grid['FINE_max_widths'],'hard_cell_caps':caps,
        'backend':policy['backend'],'dependency':cal['dependency'],'stop_rules':stop_rules,'development_released':False,'development_numerical_runs':0,'fresh':'NOT_GENERATED',
        'R4_progress_percent':0,'A2_depth_SSP_P5':'UNOPENED','binding_note':'Freeze manifest bound by enclosing PRE-RUN Git commit. Other new sources/artifacts and append-only master texts bind raw bytes plus canonical-LF SHA256 and Git blob identity; only CRLF/LF transport change allowed. Frozen historical evidence checked under prior published raw/LF rules, never rewritten.'})
    print(json.dumps({'decision':DECISION,'checks':len(checks),'failed':0,'protected_files':len(before),'historical_raw_bindings':total_raw,'historical_inputs':raw+lf,'development_runs':0}))

def verify_sealed():
    manifest=load(OUT/'CZ_PRE_RUN_DESIGN_FREEZE.json')
    for rel,b in manifest['sha256'].items():
        actual=text_identity(ROOT/rel)
        if actual['raw_sha256']!=b['raw_sha256'] and (actual['canonical_lf_sha256']!=b['canonical_lf_sha256'] or actual['canonical_git_blob_oid']!=b['canonical_git_blob_oid']):raise RuntimeError('Sealed input changed: '+rel)
    for rel,expected in manifest['frozen_model_and_observation_input_sha256'].items():
        if sha(ROOT/rel)!=expected:raise RuntimeError('Frozen model/observation changed: '+rel)
    checks=rows(OUT/'CZ_PRE_RUN_INTEGRITY.csv')
    summary=load(OUT/'CZ_PRE_RUN_INTEGRITY_SUMMARY.json')
    if not all(x['status']=='PASS' for x in checks) or len(checks)!=summary['checks']:raise RuntimeError('Integrity evidence mismatch')
    if load(OUT/'CZ_PRE_RUN_DECISION.json')['decision']!=DECISION:raise RuntimeError('Wrong pre-run decision')
    print(json.dumps({'sealed_bindings_verified':len(manifest['sha256']),'integrity_checks':len(checks),'sealed_verification':'PASS','development_runs':0}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--seal',action='store_true');parser.add_argument('--verify-sealed',action='store_true');args=parser.parse_args()
    if args.seal and args.verify_sealed:parser.error('Choose one operation')
    if args.seal:
        if (OUT/'CZ_PRE_RUN_DESIGN_FREEZE.json').exists():raise RuntimeError('Already frozen; never reseal')
        seal()
    elif args.verify_sealed:verify_sealed()
    else:parser.error('Explicit --seal or --verify-sealed required')
