from pathlib import Path
import json,csv,hashlib,subprocess
R=Path.cwd();P=R/'results/ENV_CAP_CZROB02';O=R/'results/ENV_CAP_NUM03';inputs={}
def h(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for a in json.loads((P/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))['execution_inputs']:inputs[a['path']]={'path':a['path'],'sha256':a['sha256'],'source':'prior_frozen_input'}
for a in csv.DictReader((P/'FILE_HASHES_SHA256.csv').open(encoding='utf-8-sig')):
 p=P/a['path'];assert p.stat().st_size==int(a['bytes']);inputs[str(p)]={'path':str(p),'sha256':a['sha256'],'source':'CZROB02_delivery'}
for p in [P/'FILE_HASHES_SHA256.csv',P/'DESIGN_FREEZE.json']:inputs[str(p)]={'path':str(p),'sha256':h(p),'source':'prior_manifest'}
for a in inputs.values():assert h(a['path'])==a['sha256'],a['path']
cmd=['git','-c','safe.directory='+str(R)];parent=subprocess.check_output(cmd+['rev-parse','HEAD'],text=True).strip();assert parent=='392710ecaa69e02b6247068cb83d062a0e610995';status=subprocess.check_output(cmd+['status','--porcelain=v1','-z']).decode('utf-8').split('\x00');assert all(a.startswith('?? ') for a in status if a)
r={'PASS':True,'parent_SHA':parent,'input_count':len(inputs),'inputs':list(inputs.values()),'untracked_preserved':[a for a in status if a and not a.startswith('?? results/ENV_CAP_NUM03')],'old_failures':json.loads((P/'CZROB02_INDEPENDENT_COLD_AUDIT.json').read_text(encoding='utf-8'))['FAIL_count']};(O/'PRECHECK_IDENTITY.json').write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({k:v for k,v in r.items() if k not in ('inputs','untracked_preserved')}))
