"""Delivery consistency checks and manifest; no scientific execution."""
from pathlib import Path
import ast,json,csv,hashlib,subprocess
O=Path(__file__).resolve().parent;R=O.parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def js(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def rr(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def save(n,x):
 with (O/n).open('x',encoding='utf-8') as f:json.dump(x,f,indent=2,ensure_ascii=False);f.write('\n')
checks=[]
def ck(k,v):checks.append(dict(check=k,PASS=bool(v)))
for p in O.glob('*.py'):ast.parse(p.read_text(encoding='utf-8-sig'));ck('syntax:'+p.name,True)
f=js(O/'SOURCE_COMPLETENESS_FINAL.json');a=js(O/'AUDIT_FINAL/INDEPENDENT_SOURCE_COMPLETENESS_AUDIT.json');b=js(O/'INDEPENDENT_DISCOVERY_GAP_CLOSURE.json')
ck('full_audit_not_rewritten',a['PASS_count']==671 and a['FAIL_count']==1 and a['status']=='FAIL')
ck('supplement_all14_pass',b['status']=='PASS' and len(b['checks'])==14 and all(x['PASS'] for x in b['checks']))
ck('semantic_table_complete',len(rr(O/'SOURCE_SEMANTIC_ADJUDICATION_FINAL.csv'))==3802)
ck('registry_binding',sha(O/'HISTORICAL_STATE_EXCLUSIONS_V3.csv')==f['registry_sha256'] and len(rr(O/'HISTORICAL_STATE_EXCLUSIONS_V3.csv'))==762)
ck('formal_authorization_false',f['FORMAL_EXECUTION_AUTHORIZED'] is False and js(O/'REPAIR_02_STOP_RECORD.json')['FORMAL_EXECUTION_AUTHORIZED'] is False)
ck('original_tracked_files_unchanged',subprocess.check_output(['git','diff','--name-only','HEAD'],cwd=R).decode().strip()=='')
ck('parent',subprocess.check_output(['git','rev-parse','HEAD'],cwd=R).decode().strip()=='65f9cec61c05e43b3a628f91b078bfdcc4b5318b')
ck('only_authorized_untracked_directory',subprocess.check_output(['git','status','--porcelain'],cwd=R).decode().strip()=='?? results/CZ_D2/HISTORY_REPAIR_02/')
ck('no_binary_science_outputs',not any(p.suffix in ('.npz','.npy','.mod','.shd','.wav') for p in O.rglob('*') if p.is_file()))
assert all(x['PASS'] for x in checks),checks
save('DELIVERY_VALIDATION.json',dict(status='PASS',checks=checks,scope='Artifact consistency only; independent science audit and formal MC not run',formal_execution_authorized=False,remote_before_commit={'research/cz-d2':'65f9cec61c05e43b3a628f91b078bfdcc4b5318b','main':'a9a80191ee4d0f0cecd9db296198606b63a68080'}))
files=[dict(path=p.relative_to(O).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(O.rglob('*')) if p.is_file() and p.name!='DELIVERY_SHA256_MANIFEST.json']
save('DELIVERY_SHA256_MANIFEST.json',dict(scope='All REPAIR02 files except manifest; Git commit seals manifest',files=files,total_bytes=sum(x['bytes'] for x in files)))
print(json.dumps(dict(PASS=len(checks),files=len(files),bytes=sum(x['bytes'] for x in files))))
