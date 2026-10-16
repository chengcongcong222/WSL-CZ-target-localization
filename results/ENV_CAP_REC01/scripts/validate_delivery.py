"""Read-only final delivery integrity checks, not a second scientific audit."""
from pathlib import Path
import csv,json,hashlib,subprocess,time
O=Path(__file__).resolve().parents[1];R=O.parents[1];W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_REC01');OLD=O.parent/'ENV_CAP_FULL_CAMPAIGN'
def h(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,a):p.write_text(json.dumps(a,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
checks=[]
def ck(n,b,detail=None):checks.append({'check':n,'PASS':bool(b),'detail':detail})
f=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'));ck('all825_input_and_frozen_code_still_identical',all(h(a['path'])==a['sha256'] for a in f['execution_inputs']))
pre=json.loads((O/'REC01_CZ_PHYSICS_ADJUDICATION_PRE_COLD.json').read_text(encoding='utf-8'));post=json.loads((O/'REC01_CZ_PHYSICS_ADJUDICATION.json').read_text(encoding='utf-8'));ck('post_cold_metadata_only_no_scientific_case_changes',pre['cases']==post['cases'])
# Historical references used in the review: preserve recorded working bytes, do not normalize old files.
refs=[]
for a in json.loads((OLD/'E0_HISTORICAL_INPUT_IDENTITY.json').read_text(encoding='utf-8')):
 p=R/a['repository_relative_path'];actual=h(p) if p.is_file() else None;expected=a['sha256_working'];ck('historical_reference_'+a['evidence_alias'],actual==expected);refs.append({'source':'frozen_E0_history','path':str(p),'sha256':actual,'expected_from_old_identity':expected,'match':actual==expected})
for p in sorted(Path('D:/ProjectStorage/WSL-CZ/T2_PREC3/official_source').glob('*.f90')):refs.append({'source':'official_source_snapshot_build_link_UNKNOWN','path':str(p),'sha256':h(p),'expected_from_old_identity':'SUPPLEMENTAL_READ_ONLY_INDEX','match':'NOT_A_BINARY_IDENTITY_CERTIFICATE'})
for pat in ['R4_T2_N3_PREC1*/GPT_SYNC*','R4_T2_PREC2*/GPT_SYNC*','R4_T2_PREC3*/GPT_SYNC*']:
 for p in sorted((R/'results').glob(pat)):refs.append({'source':'reviewed_history_read_only','path':str(p),'sha256':h(p),'expected_from_old_identity':'unchanged Git tracked history','match':'SEE_GIT_PARENT_DIFF'})
with (O/'REC01_REFERENCE_SOURCE_INDEX.csv').open('w',encoding='utf-8',newline='') as out:
 d=csv.DictWriter(out,['source','path','sha256','expected_from_old_identity','match']);d.writeheader();d.writerows(refs)
for p in ['GPT_SYNC_ENV_CAP_REC01.md','REC01_E1_NUMERICAL_RECONSTRUCTION.csv','REC01_CZ_PHYSICS_ADJUDICATION.json','REC01_MIRROR_CONTROL.json','REC01_COHERENT_FORWARD_FEASIBILITY.md','REC01_RESOURCE_AND_RESEARCH_GAPS.md','REC01_INDEPENDENT_COLD_AUDIT.json','REC01_FINAL_DECISION.json','STOP.json','RAW_DATA_INDEX.csv']:
 ck('required_'+p,(O/p).is_file() and (O/p).stat().st_size>0)
rows=list(csv.DictReader((O/'REC01_FULL_RANGE_ENERGY.csv').open(encoding='utf-8-sig')));ck('all31716_energy_rows_present',len(rows)==31716)
for a in rows:
 ck_key=False
# One aggregate integrity check, not a second scientific reconstruction.
ck('all15_units_retained',len(post['cases'])==15 and len(set(c['case_id'] for c in post['cases']))==15)
cold=json.loads((O/'REC01_INDEPENDENT_COLD_AUDIT.json').read_text(encoding='utf-8'));ck('frozen_independent_audit_PASS',cold['PASS'] and cold['FAIL_count']==0 and cold['checks']==1137)
ck('all12_old_failures_mapped',len(list(csv.DictReader((O/'REC01_LEGACY_FAILURE_ATTRIBUTION.csv').open(encoding='utf-8-sig'))))==12)
ck('zero_windows_retained',sum(a['exact_zero_window']=='True' for a in rows)==565)
ck('zero_new_solver_MC',json.loads((O/'REC01_FINAL_DECISION.json').read_text(encoding='utf-8'))['R3']['new_solver_calls']==0)
oldstop=json.loads((OLD/'STOP.json').read_text(encoding='utf-8'));ck('old_STOP_file_identity_preserved',any(a['path']==str(OLD/'STOP.json') and h(OLD/'STOP.json')==a['sha256'] for a in f['execution_inputs']))
# Proof of Commit A ordering and scope, with exact blob identities for the frozen code.
cmd=['git','-c','safe.directory='+str(R)];design='ebb1f45fddf468c04bedfd027e8bce54d1cbec65'
for p in (O/'scripts').glob('*.py'):
 if p.name not in ('rec01_analysis.py','rec01_cold.py','preflight.py','run_once.py'):continue
 b=subprocess.check_output(cmd+['show',design+':'+p.relative_to(R).as_posix()],cwd=R);ck('Commit_A_code_'+p.name,hashlib.sha256(b).hexdigest()==h(p))
commit_unix=int(subprocess.check_output(cmd+['show','-s','--format=%ct',design],cwd=R,text=True));start=json.loads((W/'FORMAL_STARTED.json').read_text(encoding='utf-8'))['time_unix'];ck('Commit_A_before_formal_science',commit_unix<=start)
diff=subprocess.check_output(cmd+['diff','--name-only','-z',f['parent_SHA']],cwd=R).decode('utf-8').split('\x00');ck('only_new_REC01_tracked_paths_changed',all(x.startswith('results/ENV_CAP_REC01/') for x in diff if x))
sizes={'Git_results_bytes':sum(p.stat().st_size for p in O.rglob('*') if p.is_file()),'D_work_bytes':sum(p.stat().st_size for p in W.rglob('*') if p.is_file())};ck('Git_under64MiB',sizes['Git_results_bytes']<64*1024**2);ck('new_D_under4GiB',sum(sizes.values())<4*1024**3)
ck('one_analysis_and_one_cold_launch',len(list(W.glob('ANALYSIS_LAUNCH.json')))==1 and len(list(W.glob('COLD_LAUNCH.json')))==1)
result={'PASS':all(a['PASS'] for a in checks),'scope':'Delivery integrity only; no scientific retry or recomputation','sizes':sizes,'checks':checks};js(O/'DELIVERY_VALIDATION.json',result)
# Manifest excludes only itself; Git Commit B seals its bytes.
with (O/'FILE_HASHES_SHA256.csv').open('w',encoding='utf-8',newline='') as out:
 d=csv.DictWriter(out,['path','bytes','sha256']);d.writeheader()
 for p in sorted(O.rglob('*')):
  if p.is_file() and p.name!='FILE_HASHES_SHA256.csv':d.writerow({'path':p.relative_to(O).as_posix(),'bytes':p.stat().st_size,'sha256':h(p)})
print(json.dumps({'PASS':result['PASS'],'checks':len(checks),'failures':[a['check'] for a in checks if not a['PASS']],'sizes':sizes},ensure_ascii=False));raise SystemExit(0 if result['PASS'] else 2)
