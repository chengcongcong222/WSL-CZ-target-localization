"""Final metadata/integrity validation from saved outputs; no scientific rerun."""
from pathlib import Path
import json,csv,hashlib,subprocess
O=Path(__file__).resolve().parents[1];R=O.parents[1];W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_NUM03')
def h(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,a):p.write_text(json.dumps(a,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
checks=[]
def ck(n,b,d=None):checks.append({'check':n,'PASS':bool(b),'detail':d})
f=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'));ck('all1043_frozen_identities',all(h(a['path'])==a['sha256'] for a in f['execution_inputs']))
seal=json.loads((O/'PRIMARY_OUTPUT_FREEZE.json').read_text(encoding='utf-8'));ck('all_primary_outputs_unchanged_after_cold',all(h(a['path'])==a['sha256'] and Path(a['path']).stat().st_size==a['bytes'] for a in seal['files']))
for name in ['GPT_SYNC_NUM03.md','ANALYTIC_REFERENCE_DERIVATION.md','NUM03_NINE_CASE_SENSITIVITY.csv','NUM03_HISTORICAL_FAILURE_RECHECK.csv','NUM03_INDEPENDENT_COLD_AUDIT.json','NUM03_RESEARCH_DECISION.md','NUM03_FINAL_DECISION.json','RAW_DATA_INDEX.csv','STOP.json']:
 ck('required_'+name,(O/name).is_file() and (O/name).stat().st_size>0)
final=json.loads((O/'NUM03_FINAL_DECISION.json').read_text(encoding='utf-8'));cold=json.loads((O/'NUM03_INDEPENDENT_COLD_AUDIT.json').read_text(encoding='utf-8'));primary=json.loads((O/'NUM03_PRIMARY_RESULTS.json').read_text(encoding='utf-8'));summary=read(O/'NUM03_NINE_CASE_SENSITIVITY.csv');domains=read(O/'NUM03_FULL_DOMAIN_SUMMARY.csv');mp=read(O/'NUM03_MP_REFERENCE.csv');hist=read(O/'NUM03_HISTORICAL_FAILURE_RECHECK.csv')
ck('cold73820_all_PASS',cold['PASS'] and cold['checks']==73820 and cold['FAIL_count']==0)
ck('all9_cases_one_family',len(summary)==9 and len({a['case_id'] for a in summary})==9 and {a['source_family'] for a in summary}=={'WOA23_17N118E'})
ck('all72018_reference_rows',len(read(O/'NUM03_FULL_DOMAIN_INDEPENDENT_REFERENCE.csv'))==72018 and len(read(O/'NUM03_RETURN_CURVES.csv'))==72018)
ck('all_coarse_metadata_matches_independent_saved_grid',len(domains)==18 and all(a['coarse_grid_reference_metadata_match']=='True' for a in domains))
ck('all30_MP_controls',len(mp)==30 and all(a['PASS']=='True' for a in mp))
ck('all84_ODE_outputs_and_limits',len(read(O/'NUM03_ODE_COMPARISON.csv'))==84 and all(a['PASS']=='True' and a['step_PASS']=='True' for a in read(O/'NUM03_ODE_COMPARISON.csv')))
ck('all18_first_candidate_checks',len(read(O/'NUM03_FIRST_CANDIDATE_INDEPENDENT.csv'))==18 and all(a['PASS']=='True' for a in read(O/'NUM03_FIRST_CANDIDATE_INDEPENDENT.csv')))
ck('all6_historical_failures_rechecked_not_relabelled',len(hist)==6 and all(a['old_failure_unchanged']=='True' for a in hist))
ck('geometry9_and_layered_limits',final['geometrically_supported_cases']==9 and final['labels']['NUMERICAL_REFERENCE_VALIDATED'] and final['labels']['FIRST_RETURN_GEOMETRY_SUPPORTED'] and final['labels']['REPRESENTATION_SENSITIVITY_CHARACTERIZED'] and final['labels']['ENERGY_ATTRIBUTION_UNRESOLVED'] and final['labels']['YANG_MMAR_PERFORMANCE_NOT_EVALUATED'])
ck('zero_new_calls_MC_and_no_followup',final['resource']['new_acoustic_solver_calls']==0 and final['resource']['new_MC']==0 and not final['automatic_coherent_provider_development'] and not final['automatic_Yang_MMAR_performance'])
# Preserve the development parser failure as well as accepted numerical tests.
ck('development_failure_v002_retained',not json.loads((O/'development/v002/RESULTS.json').read_text(encoding='utf-8'))['PASS'])
ck('development16113_PASS',json.loads((O/'DEVELOPMENT_GATE.json').read_text(encoding='utf-8'))['numerical_checks_passed']==16113)
for name in ['integrator.py','independent_reference.py']:ck(name+'_unchanged_since_v001',h(O/'scripts'/name)==h(O/'development/v001'/name)==h(O/'development/v003'/name)==h(O/'development/v004'/name))
cmd=['git','-c','safe.directory='+str(R)];design='a641df8e72adf7acec4a71e295c06c78b7458079'
for q in (O/'scripts').glob('*.py'):
 if q.name in ['build_delivery.py','validate_delivery.py']:continue
 b=subprocess.check_output(cmd+['show',design+':'+q.relative_to(R).as_posix()],cwd=R);ck('Commit_A_'+q.name,hashlib.sha256(b).hexdigest()==h(q))
ct=int(subprocess.check_output(cmd+['show','-s','--format=%ct',design],cwd=R,text=True));ck('Commit_A_before_science',ct<=json.loads((W/'FORMAL_STARTED.json').read_text(encoding='utf-8'))['time_unix'])
tracked=subprocess.check_output(cmd+['diff','--name-only','-z',f['parent_SHA']],cwd=R).decode('utf-8').split('\x00');ck('only_NUM03_tracked_changes',all(a.startswith('results/ENV_CAP_NUM03/') for a in tracked if a))
status=subprocess.check_output(cmd+['status','--porcelain=v1','-z'],cwd=R).decode('utf-8').split('\x00');kept={s for s in status if s.startswith('?? ') and not s.startswith('?? results/ENV_CAP_NUM03/')};ck('original16_untracked_preserved',kept==set(f['old_untracked_preserved']))
for phase in ['ANALYSIS','COLD']:
 a=json.loads((O/(phase+'_PROCESS_USAGE.json')).read_text(encoding='utf-8'));ck(phase+'_single_execution_and_budget',a['returncode']==0 and a['hard_limit_failure'] is None and a['tree_peak_RSS_bytes']<=f['process_tree_RSS_max_bytes'])
ck('science_wall_under_limit',final['resource']['science_wall_s']<=f['science_wall_max_s'])
size=sum(p.stat().st_size for p in O.rglob('*') if p.is_file());work=sum(p.stat().st_size for p in W.rglob('*') if p.is_file());ck('Git_under96MiB',size<f['new_git_max_bytes'],size);ck('D_under4GiB',size+work<f['new_D_max_bytes'],size+work)
result={'PASS':all(a['PASS'] for a in checks),'scope':'Delivery metadata/integrity only; scientific evidence comes from the one frozen analysis and one independent audit','checks':checks,'Git_results_bytes_before_manifest':size,'D_work_bytes':work,'scientific_cold_PASS':cold['PASS']};js(O/'DELIVERY_VALIDATION.json',result)
with (O/'FILE_HASHES_SHA256.csv').open('w',encoding='utf-8',newline='') as out:
 d=csv.DictWriter(out,['path','bytes','sha256']);d.writeheader()
 for p in sorted(O.rglob('*')):
  if p.is_file() and p.name!='FILE_HASHES_SHA256.csv':d.writerow({'path':p.relative_to(O).as_posix(),'bytes':p.stat().st_size,'sha256':h(p)})
print(json.dumps({'PASS':result['PASS'],'checks':len(checks),'failures':[a['check'] for a in checks if not a['PASS']],'Git_MiB':size/1024**2,'new_data_MiB':(size+work)/1024**2}));raise SystemExit(0 if result['PASS'] else 2)
