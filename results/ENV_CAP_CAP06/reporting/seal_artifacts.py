from pathlib import Path
import json,csv,hashlib,time
O=Path(__file__).resolve().parents[1];W=Path(r'D:\ProjectStorage\WSL-CZ\ENV_CAP_CAP06');start=time.monotonic()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def tab(name,a):
 with (O/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(a[0]));w.writeheader();w.writerows(a)
checks=[]
for rel,h in json.loads((O/'CODE_FREEZE.json').read_text()).items():checks.append(dict(identity=rel,PASS=sha(O/rel)==h))
expected={}
for name in ['P1_RAW_DATA_INDEX.csv','P2_RAW_DATA_INDEX.csv','P2_PRODUCT_INDEX.csv']:
 for r in csv.DictReader((O/name).open(encoding='utf-8')):expected[str(Path(r['path']).resolve())]=r['sha256']
raw=[]
for i,p in enumerate(sorted((W/'raw').rglob('*'))):
 if not p.is_file():continue
 h=sha(p);old=expected.get(str(p.resolve()));raw.append(dict(path=str(p),bytes=p.stat().st_size,sha256=h,prior_index_verified=(old==h) if old else 'NEW_MATRIX_INDEX',scope='private label access restricted in frontend; data synthetic'))
 if old:checks.append(dict(identity=str(p),PASS=old==h))
 if i%2000==0:print('sealed raw',i,flush=True)
tab('RAW_ARCHIVE_INDEX.csv',raw)
fig=[]
for name,source in [('P1_TOY_NUISANCE.png','P1_NUISANCE_INFO.json'),('P2_RECEIVER_EXTRACTION.png','P2_EXTRACTION_SUMMARY.csv'),('P2_STRESS_BOUNDARY.png','P2_EXTRACTION_SUMMARY.csv')]:
 p=O/'figures'/name;fig.append(dict(figure='figures/'+name,sha256=sha(p),source=source,source_sha256=sha(O/source),scope='TOY or PLANE_WAVE_CONTROL; not actual CZ; clustered sample denominators'))
tab('FIGURE_PROVENANCE.csv',fig)
required=json.loads((O/'task_package/12_STAGE_SPEC.json').read_text(encoding='utf-8'))['required_final_files']
for name in required:checks.append(dict(identity='required:'+name,PASS=(O/name).is_file()))
assert all(x['PASS'] for x in checks)
size=sum(p.stat().st_size for p in O.rglob('*') if p.is_file());assert size<120*2**20
out=dict(PASS=True,scope='Administrative final identity/required-delivery seal, not a repeated scientific cold audit',identity_checks=len(checks),FAIL=0,raw_files=len(raw),raw_bytes=sum(r['bytes'] for r in raw),Git_result_bytes_before_index=size,wall_s=time.monotonic()-start,required_files=required,code_unchanged_from_A=True)
(O/'FINAL_DELIVERY_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
art=[dict(path=str(p.relative_to(O)).replace('\\','/'),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(O.rglob('*')) if p.is_file() and p.name!='ARTIFACT_INDEX.csv'];tab('ARTIFACT_INDEX.csv',art)
print(out)
