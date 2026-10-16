"""Final delivery integrity, not a second scientific audit and never a scientific PASS override."""
from pathlib import Path
import json,csv,hashlib,subprocess
O=Path(__file__).resolve().parents[1];R=O.parents[1];W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_CZROB02')
def h(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,a):p.write_text(json.dumps(a,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
checks=[]
def ck(n,b,d=None):checks.append({'check':n,'PASS':bool(b),'detail':d})
f=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'));ck('all931_frozen_input_and_code_identities',all(h(a['path'])==a['sha256'] for a in f['execution_inputs']))
seal=json.loads((O/'PRIMARY_OUTPUT_FREEZE.json').read_text(encoding='utf-8'));ck('primary_outputs_and_local_exit_records_unchanged_after_cold',all(h(a['path'])==a['sha256'] and Path(a['path']).stat().st_size==a['bytes'] for a in seal['files']))
required=['GPT_SYNC_CZROB02.md','CZROB02_SSP_TOPOLOGY.csv','CZROB02_PATH_ATTRIBUTION.csv','CZROB02_INDEPENDENT_COLD_AUDIT.json','CZROB02_RESEARCH_DECISION.md','CZROB02_FINAL_DECISION.json','RAW_DATA_INDEX.csv','STOP.json','PRIMARY_OUTPUT_FREEZE.json']
for name in required:ck('required_'+name,(O/name).is_file() and (O/name).stat().st_size>0)
res=json.loads((O/'CZROB02_PRIMARY_RESULTS.json').read_text(encoding='utf-8'));final=json.loads((O/'CZROB02_FINAL_DECISION.json').read_text(encoding='utf-8'));cold=json.loads((O/'CZROB02_INDEPENDENT_COLD_AUDIT.json').read_text(encoding='utf-8'));fail=[a for a in cold['checks_detail'] if not a['PASS']]
ck('all9_cases_retained',len(res['cases'])==9 and len({a['case_id'] for a in res['cases']})==9)
ck('single_source_family',{a['source_family'] for a in res['cases']}=={'WOA23_17N118E'})
ck('all18_failed_quadrature_gates_reported',all(not m['quadrature_PASS'] for a in res['cases'] for m in a['methods'].values()) and final['A']['numeric_admitted_cases']==0)
ck('all6_cold_failures_retained_no_relabel',len(fail)==6 and not cold['PASS'] and not final['science_gate_PASS'] and final['cold_audit']['FAIL_count']==6)
ck('all_ODE_cancelled_not_claimed_pass',all(m['ODE_cases']==0 and m['ODE_PASS'] is None for a in res['cases'] for m in a['methods'].values()) and final['A']['independent_ODE_runs']==0)
ck('energy_attribution_unresolved_and_no_provider_opened',final['B']['exclusive_energy_attribution']=='EXCLUSIVE_ENERGY_ATTRIBUTION_UNRESOLVED' and final['C']['invest_coherent_provider_now']=='DEFER')
ck('two_science_launches_only',len(list(W.glob('ANALYSIS_LAUNCH.json')))==1 and len(list(W.glob('COLD_LAUNCH.json')))==1)
cmd=['git','-c','safe.directory='+str(R)];design='3d2a56e1821c3c94949346bb5d9e80170dfc221b'
for name in ['physics.py','analysis.py','cold.py','preflight.py','run_once.py']:
 p=O/'scripts'/name;b=subprocess.check_output(cmd+['show',design+':'+p.relative_to(R).as_posix()],cwd=R);ck('Commit_A_code_'+name,h(p)==hashlib.sha256(b).hexdigest())
ct=int(subprocess.check_output(cmd+['show','-s','--format=%ct',design],cwd=R,text=True));launch=json.loads((W/'FORMAL_STARTED.json').read_text(encoding='utf-8'));ck('A_precedes_main',ct<=launch['time_unix'])
tracked=subprocess.check_output(cmd+['diff','--name-only','-z',f['parent_SHA']],cwd=R).decode('utf-8').split('\x00');ck('only_new_stage_tracked_paths',all(p.startswith('results/ENV_CAP_CZROB02/') for p in tracked if p))
status=subprocess.check_output(cmd+['status','--porcelain=v1','-z'],cwd=R).decode('utf-8').split('\x00');kept={s for s in status if s.startswith('?? ') and not s.startswith('?? results/ENV_CAP_CZROB02/')};ck('16_old_untracked_preserved',kept==set(f['precheck_untracked_preserved']))
for phase,rc in [('ANALYSIS',0),('COLD',2)]:
 u=json.loads((O/(phase+'_PROCESS_USAGE.json')).read_text(encoding='utf-8'));ck(phase+'_execution_and_budget',u['returncode']==rc and u['hard_limit_failure'] is None and u['tree_peak_RSS_bytes']<=f['process_tree_RSS_max_bytes'])
ck('science_wall_under_limit',final['budget']['science_wall_s']<=f['science_wall_max_s'])
size=sum(p.stat().st_size for p in O.rglob('*') if p.is_file());work=sum(p.stat().st_size for p in W.rglob('*') if p.is_file());ck('Git_under64MiB',size<f['new_git_max_bytes'],size);ck('D_under4GiB',size+work<f['new_D_max_bytes'],size+work)
output={'PASS':all(a['PASS'] for a in checks),'scope':'Delivery identities/completeness only, scientific admission remains FALSE','science_GATE_PASS':False,'cold_audit_PASS':False,'cold_audit_FAIL_count':6,'git_result_bytes_before_manifest':size,'local_work_bytes':work,'checks':checks};js(O/'DELIVERY_VALIDATION.json',output)
with (O/'FILE_HASHES_SHA256.csv').open('w',encoding='utf-8',newline='') as out:
 d=csv.DictWriter(out,['path','bytes','sha256']);d.writeheader()
 for p in sorted(O.rglob('*')):
  if p.is_file() and p.name!='FILE_HASHES_SHA256.csv':d.writerow({'path':p.relative_to(O).as_posix(),'bytes':p.stat().st_size,'sha256':h(p)})
print(json.dumps({'delivery_PASS':output['PASS'],'science_PASS':False,'cold_FAIL':6,'checks':len(checks),'delivery_failures':[a['check'] for a in checks if not a['PASS']],'Git_result_MiB':size/1024**2}));raise SystemExit(0 if output['PASS'] else 2)
