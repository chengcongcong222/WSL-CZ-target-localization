"""Seal infrastructure closure evidence; no development observations or new forward runs."""
from __future__ import annotations
import argparse,ast,concurrent.futures,csv,hashlib,json,math,re,shutil,subprocess,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN'
SCRATCH=Path('D:/ProjectStorage/WSL-CZ/gate3a-closure-20261004')
PARENT='e9316c031d40f11fc07772959d439b965d876cd5'
POLICY_SHA='98f4e3184752025be73da2ba2bad568bdc0e9992'
SOURCES=['r4_cz_envelope_support_corrected.py','r4_cz_envelope_closure_audit.py','r4_cz_envelope_closure_integrity.py','tests/test_r4_cz_envelope_closure.py']
MASTER={'results/R4_MASTER/R4_PLAN.md','results/R4_MASTER/R4_EVIDENCE_LEDGER.csv'}
NEW_NAMES={'CZ_FIXTURE_GEOMETRIC_BUDGET_ANALYSIS.json','CZ_CLOSURE_BUDGET_LIMITATION.md','CZ_PRE_RUN_INDEPENDENT_AUDIT_FINDING.md','CZ_CLOSURE_CORRECTION_POLICY.json','CZ_STRUCTURAL_BUDGET_ANALYSIS.json','CZ_CLOSURE_FIXTURE_RESULTS.csv','CZ_CLOSURE_FIXTURE_SUMMARY.json','CZ_CLOSURE_FIXTURE_PARTITIONS.json','CZ_COVERAGE_SEMANTICS.md','CZ_CORRECTED_BUDGET_FREEZE.json','CZ_CORRECTED_PRE_RUN_TESTS.json','CZ_CORRECTED_PRE_RUN_TESTS.xml','CZ_CORRECTED_PRE_RUN_INTEGRITY.csv','CZ_CORRECTED_PRE_RUN_INTEGRITY_SUMMARY.json','CZ_CORRECTED_PRE_RUN_DECISION.json','CZ_CORRECTED_PRE_RUN_FREEZE.json','GPT_SYNC_CORRECTED.md'}

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def rows(p):return list(csv.DictReader(p.open(encoding='utf-8-sig',newline='')))
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def digest(b):return hashlib.sha256(b).hexdigest()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(4194304),b''):h.update(block)
    return h.hexdigest()
def canonical(b):return b.replace(b'\r\n',b'\n')
def oid(b):return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
def git(*a):return subprocess.check_output(['git',*a],cwd=ROOT)
def identity(p):
    b=p.read_bytes();c=canonical(b)
    return {'raw_sha256':digest(b),'canonical_lf_sha256':digest(c),'git_blob_oid':oid(c),'transport':'CRLF_LF_ONLY'}
def widths(cell):return tuple(b-a for a,b in zip(cell[0],cell[1]))
def split(cell,terminal):
    axis=max(range(4),key=lambda i:(widths(cell)[i]/terminal[i],-i));mid=(cell[0][axis]+cell[1][axis])/2
    hi=list(cell[1]);hi[axis]=mid;lo=list(cell[0]);lo[axis]=mid
    return (cell[0],tuple(hi)),(tuple(lo),cell[1])

def geometric_ideal_minimum(state,terminal):
    stats={'internal_path_boxes':0,'ideal_rejected_siblings':0,'reference_terminal_boxes':0}
    def recurse(cell):
        if not all(a<=v<=b for a,v,b in zip(cell[0],state,cell[1])):
            stats['ideal_rejected_siblings']+=1;return 1
        if all(w<=t for w,t in zip(widths(cell),terminal)):
            stats['reference_terminal_boxes']+=1;return 1
        stats['internal_path_boxes']+=1
        return 1+sum(recurse(child) for child in split(cell,terminal))
    minimum=recurse(((45.,-5.,1.,-15.),(60.,5.,3.,15.)))
    return {'ideal_reference_geometry_minimum_requests':minimum,**stats}

def seal():
    checks=[]
    def check(ok,label,expected='',actual='',scope='infrastructure only'):
        checks.append({'check_id':label,'status':'PASS' if ok else 'FAIL','expected':str(expected),'actual':str(actual),'scope':scope})
        if not ok:raise AssertionError(label)
    policy=load(OUT/'CZ_CLOSURE_CORRECTION_POLICY.json');summary=load(OUT/'CZ_CLOSURE_FIXTURE_SUMMARY.json');data=rows(OUT/'CZ_CLOSURE_FIXTURE_RESULTS.csv');partitions=load(OUT/'CZ_CLOSURE_FIXTURE_PARTITIONS.json')['runs']
    check(git('rev-parse','HEAD').decode().strip()==POLICY_SHA,'policy_local_HEAD')
    check(git('ls-remote','origin','refs/heads/main').decode().split()[0]==POLICY_SHA,'policy_pushed_before_fixtures')
    old=load(SCRATCH/'protected-before.json')['files'];check(len(old)==7220,'protected_preexisting_file_count',7220,len(old))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        after=dict(zip(old,pool.map(lambda rel:{'sha256':sha(ROOT/rel),'bytes':(ROOT/rel).stat().st_size},old)))
    for rel,b in old.items():check(b==after[rel],'protected:'+rel,b['sha256'],after[rel]['sha256'],'historical raw bytes unchanged')
    write(SCRATCH/'protected-after.json',{'files':after})
    for rel,b in load(OUT/'CZ_PRE_RUN_DESIGN_FREEZE.json')['sha256'].items():
        raw=(ROOT/rel).read_bytes() if rel not in MASTER else git('show',PARENT+':'+rel)
        canon=canonical(raw)
        check(digest(raw)==b['raw_sha256'] or (digest(canon)==b['canonical_lf_sha256'] and oid(canon)==b['canonical_git_blob_oid']),'rejected_checkpoint_binding:'+rel)
    hist=ROOT/'results/R4_A1_SEARCH_TRACTABILITY_AUDIT';total=0
    for name in ['METHOD_FREEZE.json','DEVELOPMENT_EXECUTION_FREEZE.json','CHECKPOINT_EVIDENCE_FREEZE.json','DEVELOPMENT_RESULTS_EVIDENCE_FREEZE.json']:
        for rel,b in load(hist/name)['sha256'].items():check(after[rel]['sha256']==b,'historical_raw:'+name+':'+rel,b,after[rel]['sha256']);total+=1
    check(total==2796,'historical_raw_binding_count',2796,total)
    raw_ok=lf_ok=0
    for r in rows(hist/'FROZEN_INPUT_HASHES.csv'):
        rel=r['path'];ok=after[rel]['sha256']==r['sha256']
        if ok:raw_ok+=1
        else:
            c=canonical((ROOT/rel).read_bytes());ok=r['allow_line_endings']=='True' and digest(c)==r['canonical_lf_sha256'] and oid(c)==r['baseline_git_blob_oid'];lf_ok+=int(ok)
        check(ok,'historical_input:'+rel)
    check(raw_ok+lf_ok==2989,'historical_input_count',2989,raw_ok+lf_ok)
    for rel,b in load(hist/'DEVELOPMENT_RESULTS_EVIDENCE_FREEZE.json')['text_artifacts'].items():
        c=canonical(git('show','9cfa99c14bb80d08e6f16210ebf89fd9c541b850:'+rel) if rel in MASTER else (ROOT/rel).read_bytes())
        check(digest(c)==b['canonical_lf_sha256'] and oid(c)==b['canonical_git_blob_oid'],'historical_portable:'+rel)
    for name in ['CZ_CLOSURE_CORRECTION_POLICY.json','CZ_STRUCTURAL_BUDGET_ANALYSIS.json','CZ_PRE_RUN_INDEPENDENT_AUDIT_FINDING.md']:
        rel=(OUT/name).relative_to(ROOT).as_posix();check(canonical((OUT/name).read_bytes())==canonical(git('show',POLICY_SHA+':'+rel)),'policy_commit_bytes:'+name)
    for rel,expected in {**policy['input_sha256'],**policy['enclosure_source_sha256']}.items():check(sha(ROOT/rel)==expected,'unchanged_scientific_input:'+rel,expected,sha(ROOT/rel))
    for rel,expected in summary['source_sha256'].items():check(sha(ROOT/rel)==expected,'executed_corrected_source:'+rel,expected,sha(ROOT/rel))
    check(summary['policy_commit']==POLICY_SHA and summary['policy_sha256']==sha(OUT/'CZ_CLOSURE_CORRECTION_POLICY.json'),'fixture_policy_binding')
    scientific=policy['immutable_scientific_configuration']
    check(scientific['representation']=='WINDOW_SHAPE_Q2' and scientific['nuisance']=='N0_STATIC_PER_LINE_PER_WINDOW_OFFSET' and scientific['tau_F_dB']==1e-7 and scientific['hull_max_km']==5 and scientific['raw_endpoint_width_change_max_km']==.5,'unchanged_scientific_choices')
    for resolution,terminal in [('COARSE',(0.25,.625,.125,1.875)),('FINE',(.125,.3125,.0625,.9375))]:
        a=policy['structural_budget_analysis'][resolution];ks=[]
        for e,w in zip((15,10,2,30),terminal):
            k=0
            while e/2**k>w:k+=1
            ks.append(k)
        minimum=2*sum(ks)+1
        check(a['axis_splits']==ks and a['structural_minimum_cell_requests']==minimum,'automatic_structure:'+resolution,minimum,a['structural_minimum_cell_requests'])
        maxcap=min(a['structural_ceiling'],math.floor(.8*a['runtime_envelope_seconds_per_raw_run']/a['old_frozen_cell_cost_seconds']),a['storage_record_envelope_bytes_per_raw_run']//a['conservative_record_bytes'])
        check(a['candidate_cell_cap_ladder']==[minimum,4*minimum,maxcap] and all(x>=minimum for x in a['candidate_cell_cap_ladder']),'preregistered_ladder:'+resolution)
    geometric=load(OUT/'CZ_FIXTURE_GEOMETRIC_BUDGET_ANALYSIS.json')
    check(not geometric['budget_policy_changed'] and not geometric['fourth_ladder_allowed'],'geometric_diagnosis_no_retuning')
    for entry in geometric['results']:
        fixture=next(f for f in policy['fixtures'] if f['id']==entry['fixture_id'])
        analysis=policy['structural_budget_analysis'][entry['resolution']]
        independent=geometric_ideal_minimum(fixture['state'],analysis['terminal_max_widths'])
        check(all(entry[k]==v for k,v in independent.items()),'automatic_fixture_geometric_minimum:'+entry['fixture_id']+':'+entry['resolution'])
        check(entry['registered_maximum_cap']==analysis['candidate_cell_cap_ladder'][-1] and entry['maximum_satisfies_geometric_necessary_minimum']==(entry['registered_maximum_cap']>=independent['ideal_reference_geometry_minimum_requests']),'fixture_geometric_budget_feasibility_report:'+entry['fixture_id']+':'+entry['resolution'])
    check(len(data)==len(partitions)==summary['raw_run_count'],'run_rows_partition_count')
    mapping={p['run_id']:p for p in partitions};references={f['id']:f['state'] for f in policy['fixtures']}
    result_map={}
    for r in data:
        rid=f"{r['resolution']}_L{r['ladder_tier']}_{r['fixture_id']}_{r['mode']}";p=mapping[rid];d=p['diagnostics'];boxes=p['sealed_partition'];resolution=r['resolution'];mode=r['mode'];cap=int(r['hard_cell_cap']);tier=int(r['ladder_tier'])
        terminal=tuple(policy['structural_budget_analysis'][resolution]['terminal_max_widths'])
        check(cap==policy['structural_budget_analysis'][resolution]['candidate_cell_cap_ladder'][tier-1],'run_registered_cap:'+rid)
        check(d['n_cell_requests']<=cap and int(r['n_cell_requests'])==d['n_cell_requests'],'request_cap:'+rid)
        for name in ['n_forward_state_evaluations','n_profile_evaluations','n_bound_evaluations','n_cache_hits','n_validation_evaluations']:
            check(d[name]<=d['limits'][name] and int(r[name])==d[name],'category_cap:'+rid+':'+name)
        check(sum(d[k] for k in d['limits'])==d['charged_total']<=d['hard_budget_limit'],'aggregate_cap:'+rid)
        check(d['n_validation_evaluations']==0 and (mode!='B0' or d['n_forward_state_evaluations']==d['n_profile_evaluations']==0),'search_no_evaluation_and_B0_no_acoustics:'+rid)
        # Independently reconstruct the entire dyadic leaf partition, not just volume.
        leaves={(tuple(b['low']),tuple(b['high'])):b for b in boxes};used=set()
        check(len(leaves)==len(boxes),'no_duplicate_boxes:'+rid)
        def visit(cell,depth):
            if cell in leaves:
                b=leaves[cell];used.add(cell);check(b['depth']==depth,'box_depth:'+rid+':'+str(len(used)))
                if b['status'] in ('RETAINED_CERTIFIED_COMPATIBLE','RETAINED_POSSIBLE_AT_TERMINAL'):
                    check(all(w<=t for w,t in zip(widths(cell),terminal)),'retained_terminal:'+rid+':'+str(len(used)))
                return
            check(any(w>t for w,t in zip(widths(cell),terminal)),'no_partition_hole:'+rid+':'+str(depth))
            for child in split(cell,terminal):visit(child,depth+1)
        visit(((45.,-5.,1.,-15.),(60.,5.,3.,15.)),0)
        check(len(used)==len(boxes),'no_overlapping_nested_or_orphan_boxes:'+rid)
        leaf_depth=max(b['depth'] for b in boxes)
        check(d['maximum_depth']<=leaf_depth<=d['maximum_depth']+1,'dequeued_vs_frontier_depth:'+rid)
        counts={s:sum(b['status']==s for b in boxes) for s in ['REJECTED_BY_CERTIFIED_BOUND','RETAINED_CERTIFIED_COMPATIBLE','RETAINED_POSSIBLE_AT_TERMINAL','RETAINED_UNRESOLVED_BUDGET','RETAINED_INVALID_BOUND']}
        check(sum(counts.values())==len(boxes),'known_statuses_only:'+rid)
        for status,key in [('REJECTED_BY_CERTIFIED_BOUND','n_rejected'),('RETAINED_CERTIFIED_COMPATIBLE','n_certified_compatible'),('RETAINED_POSSIBLE_AT_TERMINAL','n_terminal_possible'),('RETAINED_UNRESOLVED_BUDGET','n_budget_unresolved'),('RETAINED_INVALID_BOUND','n_invalid_bound')]:
            check(counts[status]==d[key]==int(r[key]),'pruning_diagnostic:'+rid+':'+key)
        retained=[b for b in boxes if b['status']!='REJECTED_BY_CERTIFIED_BOUND'];x=references[r['fixture_id']]
        joint=any(all(a<=v<=b for a,v,b in zip(c['low'],x,c['high'])) for c in retained);range_kept=any(c['low'][0]<=x[0]<=c['high'][0] for c in retained)
        check(joint==p['reference_evaluation']['joint_retained'] and range_kept==p['reference_evaluation']['range_retained'],'independent_reference_retention:'+rid)
        closed=d['execution_valid'] and d['search_budget_closed'] and not counts['RETAINED_UNRESOLVED_BUDGET'] and not counts['RETAINED_INVALID_BOUND']
        check(d['domain_partition_closed']==closed and str(closed)==r['domain_partition_closed'],'closure_semantics:'+rid)
        check(d['estimated_storage_peak_bytes']<=policy['structural_budget_analysis'][resolution]['storage_record_envelope_bytes_per_raw_run'],'conservative_storage_envelope:'+rid)
        passed=d['execution_valid'] and closed and joint and range_kept
        result_map[(resolution,tier,r['fixture_id'],mode)]=passed
    selected={}
    for resolution in ['COARSE','FINE']:
        a=policy['structural_budget_analysis'][resolution];selected[resolution]=None
        for tier,cap in enumerate(a['candidate_cell_cap_ladder'],1):
            values=[result_map.get((resolution,tier,f['id'],mode)) for f in policy['fixtures'] for mode in ['B0','CZ']]
            if all(v is True for v in values):selected[resolution]=cap;break
        check(selected[resolution]==summary['selected_cell_caps'][resolution],'first_complete_tier_selection:'+resolution)
        if selected[resolution] is not None:
            check(selected[resolution]>=max(r['ideal_reference_geometry_minimum_requests'] for r in geometric['results'] if r['resolution']==resolution),'selected_cap_satisfies_all_fixture_geometric_minima:'+resolution)
        if summary['execution_valid'] and selected[resolution] is None:
            check(all((resolution,3,f['id'],mode) in result_map for f in policy['fixtures'] for mode in ['B0','CZ']),'max_failure_all_registered_fixtures_tested:'+resolution)
    expected=('CZ_ENVELOPE_PRE_RUN_IMPLEMENTATION_INVALID' if not summary['execution_valid'] else 'CZ_ENVELOPE_CORRECTED_PRE_RUN_READY_FOR_INDEPENDENT_AUDIT' if all(v is not None for v in selected.values()) else 'CZ_ENVELOPE_PRE_RUN_SEARCH_CLOSURE_NOT_ESTABLISHED')
    check(summary['decision']==expected,'decision_from_closed_partitions',expected,summary['decision'])
    check(summary['development_numerical_runs']==0 and summary['development_released'] is False and summary['fresh']=='NOT_GENERATED' and not summary['forbidden_input_used'],'zero_development_and_no_fresh')
    testtree=ET.parse(SCRATCH/'tests.xml');suites=testtree.findall('.//testsuite');count=sum(int(s.get('tests',0)) for s in suites)
    check(count==69 and sum(int(s.get(k,0)) for s in suites for k in ['failures','errors','skipped'])==0,'69_unit_tests_pass')
    shutil.copyfile(SCRATCH/'tests.xml',OUT/'CZ_CORRECTED_PRE_RUN_TESTS.xml')
    write(OUT/'CZ_CORRECTED_PRE_RUN_TESTS.json',{'tests':count,'passed':count,'failed':0,'errors':0,'skipped':0,'original_unchanged_tests':47,'new_closure_tests':22,
        'categories':['automatic structural derivation','T1 below-min cannot close','T2 exact-min single path closes','T3 terminal possible closes','T4 budget unresolved blocks Gate','B0 same semantics','invalid bound retained','no sampling rejection','compatibility diagnostic separated','empty support failure','sentinel/exclusion','raw independence','all six strict caps'],
        'junit_sha256':sha(OUT/'CZ_CORRECTED_PRE_RUN_TESTS.xml'),'development_data_used':False,'scientific_result':False})
    for name in SOURCES:
        source=(ROOT/name).read_text(encoding='utf-8-sig');tree=ast.parse(source);check(chr(65533) not in source,'source_syntax_utf8:'+name)
        if name=='r4_cz_envelope_support_corrected.py':
            forbidden={'truth','truth_range','truth_state','errors','old_recovered','old_best','oracle','basin'}
            ids={n.id for n in ast.walk(tree) if isinstance(n,ast.Name)}|{n.arg for n in ast.walk(tree) if isinstance(n,ast.arg)}
            check(not ids&forbidden and not any(isinstance(n,ast.ImportFrom) and ('evaluation' in (n.module or '') or 'r4_a1' in (n.module or '')) for n in ast.walk(tree)),'corrected_estimator_AST_exclusion')
    check('DEVELOPMENT_RELEASED = False' in (ROOT/'r4_cz_envelope_development.py').read_text(encoding='utf-8-sig'),'unchanged_development_hard_stop')
    for rel,previous in [('results/R4_MASTER/R4_PLAN.md','plan-before.md'),('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv','ledger-before.csv')]:
        check((ROOT/rel).read_bytes().startswith((SCRATCH/previous).read_bytes()),'master_append_only:'+rel)
        check(canonical((SCRATCH/previous).read_bytes())==canonical(git('show',PARENT+':'+rel)),'master_previous_prefix_accepted_parent:'+rel)
    ledger=rows(ROOT/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv');check(len(rows(SCRATCH/'ledger-before.csv'))==33 and len(ledger)==36 and all(float(r['credited_weight_percent'])==0 for r in ledger),'33_old_plus_3_zero_credit_ledger_rows')
    check(load(ROOT/'results/R4_MASTER/R4_PROGRESS.json')['overall_progress_percent']==0,'unchanged_R4_progress')
    check(set(git('diff','--name-only',POLICY_SHA).decode().splitlines())==MASTER,'only_master_existing_files_changed')
    check(not git('diff','--cached','--name-only').strip(),'no_preexisting_staging')
    for p in [*OUT.glob('*.md'),ROOT/'results/R4_MASTER/R4_PLAN.md']:
        check(chr(65533) not in p.read_text(encoding='utf-8-sig'),'document_utf8:'+p.name)
    subprocess.run(['git','diff','--check'],cwd=ROOT,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE);check(True,'git_diff_check')
    csvpath=OUT/'CZ_CORRECTED_PRE_RUN_INTEGRITY.csv'
    with csvpath.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['check_id','status','expected','actual','scope']);w.writeheader();w.writerows(checks)
    write(OUT/'CZ_CORRECTED_PRE_RUN_INTEGRITY_SUMMARY.json',{'status':'PASS','checks':len(checks),'failed':0,'protected_old_files':7220,'protected_files_changed':[],
        'rejected_checkpoint_bindings_verified':27,'historical_raw_bindings':2796,'historical_input_bindings':2989,'raw_matches':raw_ok,'LF_only_matches':lf_ok,
        'fixture_partitions_independently_verified':len(data),'method':'Raw/SHA/Git bindings; exact deterministic dyadic partition reconstruction; metadata accounting and post-seal reference membership; no new modal/feature computations and no development arrays',
        'CSV_sha256':sha(csvpath),'development_numerical_runs':0,'master_rows_appended':3,'R4_percent':0})
    write(OUT/'CZ_CORRECTED_BUDGET_FREEZE.json',{'policy_commit':POLICY_SHA,'selection_status':'SELECTED' if all(x is not None for x in selected.values()) else 'SEARCH_CLOSURE_NOT_ESTABLISHED_NOT_SELECTED',
        'selected_cell_caps':selected,'structural_minimum_cell_requests':{k:v['structural_minimum_cell_requests'] for k,v in policy['structural_budget_analysis'].items()},
        'candidate_ladders':{k:v['candidate_cell_cap_ladder'] for k,v in policy['structural_budget_analysis'].items()},'final_ladder_diagnostics':summary['final_ladder_diagnostics'],
        'selected_scientific_run_budgets':None if any(v is None for v in selected.values()) else {k:{'n_cell_requests':v,'CZ_forward_profile_limit_each':21*v,'CZ_bound_limit':22*v,'CZ_aggregate_limit':66*v,'B0_bound_limit':v,'B0_aggregate_limit':3*v,'cache_limit':v,'validation_limit':0} for k,v in selected.items()},
        'fixture_specific_geometric_budget_limitation':geometric,'no_fourth_budget':True,'runtime_storage_policy':policy['structural_budget_analysis'],'development_released':False})
    stop=['Stop after correction checkpoint commit/push','No development feature/search or fresh','No C4/F4, representation/noise/tau/grid/Gate changes','Keep rejected pre-run and accepted science evidence historical','R4=0%; N1/A2/depth/SSP/P5 unopened']
    write(OUT/'CZ_CORRECTED_PRE_RUN_DECISION.json',{'stage':policy['stage'],'parent_rejected_pre_run_commit':PARENT,'closure_policy_commit':POLICY_SHA,
        'lead_parent_audit':'CZ_ENVELOPE_PRE_RUN_FREEZE_NOT_ACCEPTED','decision':expected,'whole_cell_lower_bound_enclosure':'UNCHANGED','coverage_semantics':'CONSERVATIVE_OUTER_PARTITION',
        'structural_minimum_cell_requests':{k:v['structural_minimum_cell_requests'] for k,v in policy['structural_budget_analysis'].items()},'selected_cell_caps':selected,
        'closure_fixture_final_diagnostics':summary['final_ladder_diagnostics'],'reference_retained':summary['reference_retained_all_runs'],'tests_passed':69,'integrity_checks_passed':len(checks),
        'execution_valid':summary['execution_valid'],'development_numerical_feature_search_runs':0,'development_released':False,'fresh':'NOT_GENERATED','R4_percent':0,'A2_depth_SSP_P5':'UNOPENED',
        'additional_budget_limitation':geometric,'claim':'Infrastructure closure feasibility only; coarse cap below midpoint-specific ideal minimum; no CZ information/recovery or physics impossibility conclusion','stop_rules':stop,'stop_after_commit_push':True})
    bindings={name:identity(ROOT/name) for name in SOURCES}
    for name in sorted(NEW_NAMES-{'CZ_CORRECTED_PRE_RUN_FREEZE.json'}):bindings[(OUT/name).relative_to(ROOT).as_posix()]=identity(OUT/name)
    for name in MASTER:bindings[name]=identity(ROOT/name)
    write(OUT/'CZ_CORRECTED_PRE_RUN_FREEZE.json',{'stage':policy['stage'],'parent_rejected_pre_run_commit':PARENT,'closure_policy_commit':POLICY_SHA,'immutable_after_commit':True,
        'bindings':bindings,'original_enclosure_and_input_sha256':{**policy['input_sha256'],**policy['enclosure_source_sha256']},
        'old_protected_inventory_sha256':digest(json.dumps(old,sort_keys=True,separators=(',',':')).encode()),'coverage_semantics':'CONSERVATIVE_OUTER_PARTITION',
        'selected_cell_caps':selected,'structural_minimum_requests':{k:v['structural_minimum_cell_requests'] for k,v in policy['structural_budget_analysis'].items()},
        'unchanged_scientific_configuration':scientific,'decision':expected,'stop_rules':stop,'development_released':False,'development_numerical_runs':0,
        'binding_note':'This manifest bound by correction Git commit. Rejected e9316c sources/artifacts unchanged; its old master bindings checked at e9316c, current master append-only. New bindings permit CRLF/LF transport only.'})
    print(json.dumps({'decision':expected,'integrity_checks':len(checks),'failed':0,'protected_files':7220,'fixture_runs':len(data),'selected_caps':selected,'development_runs':0}))

def verify():
    m=load(OUT/'CZ_CORRECTED_PRE_RUN_FREEZE.json')
    for rel,b in m['bindings'].items():
        current=identity(ROOT/rel)
        if current['raw_sha256']!=b['raw_sha256'] and (current['canonical_lf_sha256']!=b['canonical_lf_sha256'] or current['git_blob_oid']!=b['git_blob_oid']):raise RuntimeError('Corrected freeze changed: '+rel)
    for rel,b in m['original_enclosure_and_input_sha256'].items():
        if sha(ROOT/rel)!=b:raise RuntimeError('Original science input changed: '+rel)
    data=rows(OUT/'CZ_CORRECTED_PRE_RUN_INTEGRITY.csv')
    if not all(r['status']=='PASS' for r in data):raise RuntimeError('Integrity failure')
    print(json.dumps({'verified_bindings':len(m['bindings']),'integrity_checks':len(data),'decision':m['decision'],'development_runs':0,'sealed_verification':'PASS'}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seal',action='store_true');p.add_argument('--verify',action='store_true');a=p.parse_args()
    if a.seal and a.verify:p.error('Choose one operation')
    if a.seal:
        if (OUT/'CZ_CORRECTED_PRE_RUN_FREEZE.json').exists():raise RuntimeError('Already frozen; never reseal')
        seal()
    elif a.verify:verify()
    else:p.error('Explicit --seal or --verify required')
