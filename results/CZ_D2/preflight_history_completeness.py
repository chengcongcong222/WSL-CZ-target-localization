"""Read-only historical-exclusion gate. Never generates panel/noise or fits.
Secondary NPY decoder below uses only zipfile/ast/struct, not NumPy's reader.
"""
from pathlib import Path
import ast,csv,hashlib,json,math,struct,zipfile,subprocess,time
import numpy as np
R=Path(__file__).resolve().parents[2]; O=Path(__file__).resolve().parent
A='dc91982e67921eeb7a6aac4e12c482a980d5ad89'
S=R/'results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT'
P=S/'evaluation_only/calibration_private.npz'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(name,x):
 with (O/name).open('x',encoding='utf-8') as f:json.dump(x,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def ang(a,b):return abs((a-b+180)%360-180)
def neighbor(h,o,w):return all((ang(h[j],o[j]) if j in [1,3] else abs(h[j]-o[j]))<=w[j] for j in range(4) if o[j] is not None)
start=time.perf_counter()
with (O/'PREFLIGHT_HISTORY_CHECK_STARTED.json').open('x',encoding='utf-8') as f:json.dump({'scope':'historical_input_preflight_only','formal_main_started':False,'new_panel_generated':False},f)
checks=[]
def ck(k,v,detail=''):
 checks.append({'check':k,'PASS':bool(v),'detail':detail})
 if not v:raise RuntimeError('Evidence reconstruction failure: '+k)
ck('HEAD_is_approved_A',subprocess.check_output(['git','rev-parse','HEAD'],cwd=R).decode().strip()==A)
manifest=read(O/'DESIGN_FILE_MANIFEST.json')
for row in manifest['files']:ck('A_identity:'+row['path'],sha(O/row['path'])==row['sha256'])
identity=read(O/'SOURCE_IDENTITIES.json')
for row in identity['sources']:ck('source_identity:'+row['path'],sha(R/row['path'])==row['sha256'])
source_paths=[P,S/'THRESHOLD_FREEZE.json',S/'PRIVATE_GENERATOR_MANIFEST.json',S/'CALIBRATION_SCORES.csv',S/'AUX_D1_TASK_SPEC_20261012_v1.json',R/'aux_d1_generate.py',S/'DESIGN_FREEZE.json']
source_records=[]
for p in source_paths:
 rel=p.relative_to(R).as_posix();b=subprocess.check_output(['git','show',A+':'+rel],cwd=R)
 ck('archive_in_A:'+rel,b.replace(b'\r\n',b'\n')==p.read_bytes().replace(b'\r\n',b'\n'))
 source_records.append({'path':rel,'sha256':sha(p),'bytes':p.stat().st_size,'parent_git_blob_sha256':hashlib.sha256(b).hexdigest()})
for name,key in [('THRESHOLD_FREEZE.json','bindings'),('PRIVATE_GENERATOR_MANIFEST.json','arrays')]:
 row=next(x for x in read(S/name)[key] if x['path'].replace('\\','/').endswith('/calibration_private.npz'))
 ck('historical_seal:'+name,sha(P)==row['sha256'])
with np.load(P,allow_pickle=False) as a:h=a['h_truth'].copy()
ck('archive_truth_shape_and_finite',h.shape==(512,4) and np.isfinite(h).all())
# Independent standard-library decode of the archived NPY member, no RNG/re-generation.
with zipfile.ZipFile(P) as z:raw=z.read('h_truth.npy')
ck('NPY_magic',raw[:6]==b'\x93NUMPY')
major,minor=raw[6:8];size=2 if major==1 else 4
n=struct.unpack('<H' if size==2 else '<I',raw[8:8+size])[0];offset=8+size+n
header=ast.literal_eval(raw[8+size:offset].decode('latin1'))
ck('NPY_layout',header['shape']==(512,4) and header['descr']=='<f8' and not header['fortran_order'])
flat=struct.unpack('<2048d',raw[offset:]);ind=[flat[i*4:i*4+4] for i in range(512)]
ck('independent_NPY_values',all(h[i,j]==ind[i][j] for i in range(512) for j in range(4)))
reg=list(csv.DictReader((O/'HISTORICAL_STATE_EXCLUSIONS.csv').open(encoding='utf-8')))
old=[[float(x[k]) if x[k] else None for k in ['r_km','theta_deg','v_mps','psi_deg']] for x in reg]
ck('registry_216',len(old)==216)
ck('source_not_in_registered_12',P.relative_to(R).as_posix() not in {x['source'] for x in reg})
expanded=[]
for x in old:
 expanded.extend([x,[x[0],-x[1],x[2],None if x[3] is None else -x[3]]])
origin=[45,-5,1,-15];step=[1,.5,.2,1]
out=[]
for i,x in enumerate(h):
 exact=any(all(v is not None and abs(float(x[j])-v)<1e-11 for j,v in enumerate(o)) for o in expanded)
 near=any(neighbor(x,o,[.1,.02,.01,.05] if o[2] is not None else [.01,.005,0,0]) for o in expanded)
 grid=any(abs((x[j]-origin[j])/step[j]-round((x[j]-origin[j])/step[j]))<=1e-6 for j in range(4))
 out.append(dict(archive_index=i,r_km=float(x[0]),theta_deg=float(x[1]),v_mps=float(x[2]),psi_deg=float(x[3]),exact_registered_or_mirror=exact,within_registered_exclusion=near,rejected_by_grid_guard=grid,source='AUX_D1_historical_calibration_NOT_NEW_CZ_D2_PANEL'))
# Independent plain-Python comparison, separate loops/units and no helper above.
ind_uncovered=0
for x in ind:
 found=False
 for o in old:
  for s in (1,-1):
   b=(o[0],s*o[1],o[2],None if o[3] is None else s*o[3]);width=(.1,.02,.01,.05) if o[2] is not None else (.01,.005,0,0)
   passed=True
   for j in range(4):
    if b[j] is None:continue
    d=abs(x[j]-b[j]);d=min(d%360,360-d%360) if j in (1,3) else d
    if d>width[j]:passed=False
   if passed:found=True
 ongrid=any(abs(x[j]-(origin[j]+round((x[j]-origin[j])/step[j])*step[j]))<=1e-6*step[j] for j in range(4))
 if not found and not ongrid:ind_uncovered+=1
uncovered=sum(not x['within_registered_exclusion'] and not x['rejected_by_grid_guard'] for x in out)
ck('independent_uncovered_count',uncovered==ind_uncovered)
scores=list(csv.DictReader((S/'CALIBRATION_SCORES.csv').open(encoding='utf-8-sig')))
ck('historical_calibration_used_512_both_resources',len(scores)==1024 and {int(x['calibration_index']) for x in scores}==set(range(512)))
with (O/'UNREGISTERED_HISTORICAL_CALIBRATION_STATES.csv').open('x',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),lineterminator='\n');w.writeheader();w.writerows(out)
result={'stage':'CZ-D2_COMMIT_B_PREFLIGHT','design_SHA':A,'status':'HISTORICAL_EXCLUSION_REGISTRY_INCOMPLETE','scientific_execution_admitted':False,'historical_registry_rows':216,'historical_registry_source_count':len({x['source'] for x in reg}),'discovered_archive_rows':512,'unique_discovered_states':len({tuple(x) for x in ind}),'exact_matches':sum(x['exact_registered_or_mirror'] for x in out),'covered_by_existing_neighborhood':sum(x['within_registered_exclusion'] for x in out),'covered_by_grid_guard':sum(x['rejected_by_grid_guard'] for x in out),'uncovered_by_both_rules':uncovered,'sources':source_records,'evidence_checks':checks,'read_only_elapsed_s':time.perf_counter()-start,'formal_main_calls':0,'formal_cold_calls':0,'new_panel_samples':0,'new_observations':0,'new_MC':0,'optimizer_calls':0,'propagation_calls':0,'scan_scope':'stopped_at_first_confirmed_omitted_source;not_exhaustive_historical_catalog_certification','no_registry_or_seed_change':True,'STOP':True}
dump('HISTORY_COMPLETENESS_PREFLIGHT.json',result)
dump('PREFLIGHT_EVIDENCE_REVIEW.json',{'scope':'historical_seals_and_independent_standard_library_NPY_decode;not_formal_scientific_cold_audit','checks':len(checks),'FAIL':0,'independent_uncovered_count':ind_uncovered,'status':'OMITTED_ARCHIVE_CONFIRMED','formal_cold_calls':0})
with (O/'CONFIGURATION_EXECUTION_STATUS.csv').open('x',encoding='utf-8',newline='') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['quality','geometry','planned_case_id','replicate','planned_aux_side','bias_support','status','reason'])
 for q in ['Q1','Q2']:
  for b in ['B0','B3']:
   for i in range(12):
    for rep in range(500):
     for support in ['P0','P1']:w.writerow([q,b,f'D{i+1:02d}',rep,-1 if rep<250 else 1,support,'NOT_EVALUATED','HISTORICAL_EXCLUSION_REGISTRY_INCOMPLETE'])
print(json.dumps({k:result[k] for k in ['status','discovered_archive_rows','exact_matches','covered_by_existing_neighborhood','covered_by_grid_guard','uncovered_by_both_rules','formal_main_calls','formal_cold_calls']},indent=2))
