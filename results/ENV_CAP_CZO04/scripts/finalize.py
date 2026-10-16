"""Delivery checks only; no forward or information-function execution."""
from pathlib import Path
import json,csv,hashlib,subprocess,ast
O=Path(__file__).resolve().parents[1];R=O.parents[1];W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_CZO04')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def js(p,v):Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
checks=[]
def ck(n,v):checks.append(dict(check=n,PASS=bool(v)))
freeze=read(O/'DESIGN_FREEZE.json');decision=read(O/'CZO04_FINAL_DECISION.json');audit=read(O/'CZO04_INDEPENDENT_COLD_AUDIT.json')
ck('input_identity_1141',all(sha(x['path'])==x['sha256'] for x in freeze['execution_inputs']))
ck('primary_scientific_output_seal',all(sha(x['path'])==x['sha256'] for x in read(O/'PRIMARY_OUTPUT_FREEZE.json')['files']))
ck('cold_PASS',audit['PASS'] and audit['FAIL']==0)
ck('one_analysis_marker',len(list(W.glob('ANALYSIS_LAUNCH.json')))==1)
ck('one_cold_marker',len(list(W.glob('COLD_LAUNCH.json')))==1)
for phase in ['ANALYSIS','COLD']:
 a=read(O/(phase+'_PROCESS_USAGE.json'));ck(phase+'_valid',a['returncode']==0 and a['hard_limit_failure'] is None)
ck('wall_budget',decision['science_wall_s']<freeze['science_wall_max_s']);ck('memory_budget',decision['science_peak_RSS_bytes']<freeze['process_tree_RSS_max_bytes'])
ck('zeros_solver_MC_receiver',decision['new_external_solvers']==decision['new_MC']==decision['new_receive_records']==0)
ck('R4_STOP',decision['R4']=='0%' and decision['status']=='STOP')
for name in ['GPT_SYNC_CZO04.md','CZO04_RESEARCH_DECISION.md','ZERO_CONTROLS_AND_INFORMATION_DERIVATION.md','OBSERVABLE_AND_EXTRACTION_GAPS.md','CONDITION_AND_NOISE_SUMMARY.csv','LOCAL_EFFECTIVE_INFORMATION.csv','FINITE_CANDIDATE_SUPPORT.csv','RESOURCE_SEPARATION_SCALES.csv','RAW_DATA_INDEX.csv','CZO04_INDEPENDENT_COLD_AUDIT.json']:
 ck('deliverable/'+name,(O/name).is_file() and (O/name).stat().st_size>0)
for p in (O/'scripts').glob('*.py'):ast.parse(p.read_text(encoding='utf-8-sig'));ck('syntax/'+p.name,True)
for p in O.glob('FIG*.png'):ck('figure/'+p.name,p.stat().st_size>5000)
with (O/'NOISELESS_GRID_POINT_RETENTION_CONTROL.csv').open(encoding='utf-8-sig') as f:a=list(csv.DictReader(f))
ck('108_grid_controls_per_condition',len(a)==324 and all(x['noiseless_regression_PASS']=='True' for x in a))
ck('two_conditional_witnesses_are_one_geometry_mirror',len(decision['primary_stable_witnesses'])==2 and {x['scene'] for x in decision['primary_stable_witnesses']}=={'G02P','G02M'})
ck('four_offgrid_empty_scenes',len(decision['finite_grid_empty_bearing_scenes'])==4)
ck('36_remote_review_score_copies',len(list((O/'data').glob('*_SCORES.npz')))==36 and all(sha(p)==sha(W/p.name) for p in (O/'data').glob('*_SCORES.npz')))
ck('all_D_paths',all(str(p.resolve()).lower().startswith('d:') for root in [W,O] for p in root.rglob('*') if p.is_file()))
ck('new_D_budget',sum(p.stat().st_size for root in [O,W] for p in root.rglob('*') if p.is_file())<freeze['new_D_max_bytes'])
ck('new_Git_budget',sum(p.stat().st_size for p in O.rglob('*') if p.is_file())<freeze['new_Git_max_bytes'])
js(O/'DELIVERY_VALIDATION.json',dict(PASS=all(a['PASS'] for a in checks),checks=len(checks),FAIL=sum(not a['PASS'] for a in checks),scope='delivery/metadata checks; not a second scientific cold audit',items=checks))
js(O/'RUN_STATE.json',dict(stage='ENV_CAP_CZO04',state='STOP',analysis_completed_once=True,cold_audit_completed_once=True,scientific_retry_permitted=False,commit_A=decision['design_SHA'],next_action='Commit B push/verify only; do not relaunch analysis or cold',historical_R4='0%'))
js(O/'BUDGET_USAGE.json',dict(science_wall_s=decision['science_wall_s'],science_peak_tree_RSS_bytes=decision['science_peak_RSS_bytes'],new_D_snapshot_bytes=sum(p.stat().st_size for root in [O,W] for p in root.rglob('*') if p.is_file()),new_Git_snapshot_bytes=sum(p.stat().st_size for p in O.rglob('*') if p.is_file()),snapshot_scope='after reporting, before this budget/manifest finalization; exact final storage also recorded outside Git after push',new_solver_calls=0,new_MC=0,new_receiver_records=0,all_new_paths='D:',science_retry=0))
files=[dict(path=str(p.relative_to(O)).replace('\\','/'),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(O.rglob('*')) if p.is_file() and p.name!='FILE_HASHES_SHA256.csv']
with (O/'FILE_HASHES_SHA256.csv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['path','bytes','sha256']);w.writeheader();w.writerows(files)
print(json.dumps(dict(PASS=all(a['PASS'] for a in checks),checks=len(checks),manifest_files=len(files),new_Git_bytes=sum(p.stat().st_size for p in O.rglob('*') if p.is_file()),new_D_total_bytes=sum(p.stat().st_size for root in [O,W] for p in root.rglob('*') if p.is_file()))))
assert all(a['PASS'] for a in checks)
