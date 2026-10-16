from pathlib import Path
import json,csv,hashlib,subprocess,time
R=Path.cwd();O=R/'results/ENV_CAP_CZROB02';P=R/'results/ENV_CAP_REC01';start=time.time();inputs={}
def h(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for a in json.loads((P/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))['execution_inputs']:
 inputs[a['path']]={'path':a['path'],'sha256':a['sha256'],'source':'REC01_frozen_input'}
for a in csv.DictReader((P/'FILE_HASHES_SHA256.csv').open(encoding='utf-8-sig')):
 p=P/a['path'];assert p.stat().st_size==int(a['bytes']);inputs[str(p)]={'path':str(p),'sha256':a['sha256'],'source':'REC01_delivery'}
for a in json.loads((R/'results/ENV_CAP_FULL_CAMPAIGN/E0_HISTORICAL_INPUT_IDENTITY.json').read_text(encoding='utf-8')):
 p=R/a['repository_relative_path'];inputs[str(p)]={'path':str(p),'sha256':a['sha256_working'],'source':'historical_E0_reference'}
for p in [P/'FILE_HASHES_SHA256.csv',P/'DESIGN_FREEZE.json']:
 inputs[str(p)]={'path':str(p),'sha256':h(p),'source':'REC01_manifest_and_contract'}
fail=[]
for a in inputs.values():
 if h(a['path'])!=a['sha256']:fail.append(a['path'])
cmd=['git','-c','safe.directory='+str(R)];parent=subprocess.check_output(cmd+['rev-parse','HEAD'],text=True).strip();assert parent=='4e53bf035625aeddf4d8ab6ee009242c06726a9c'
status=subprocess.check_output(cmd+['status','--porcelain=v1','-z']).decode('utf-8').split('\x00');assert all(s.startswith('?? ') for s in status if s)
res={'PASS':not fail,'parent_SHA':parent,'inputs':list(inputs.values()),'input_count':len(inputs),'failures':fail,'untracked_preserved':[s for s in status if s and not s.startswith('?? results/ENV_CAP_CZROB02')],'wall_s':time.time()-start,'prior_reports_read':['GPT_SYNC_ENV_CAP_REC01.md','REC01_INDEPENDENT_COLD_AUDIT.json','REC01_CZ_PHYSICS_ADJUDICATION.json','REC01_COHERENT_FORWARD_FEASIBILITY.md']};(O/'PRECHECK_IDENTITY.json').write_text(json.dumps(res,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({k:v for k,v in res.items() if k not in ('inputs','untracked_preserved')},ensure_ascii=False));raise SystemExit(0 if not fail else 2)
