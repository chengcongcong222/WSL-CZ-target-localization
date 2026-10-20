"""Independent read-only re-preflight. Does not import repair/discovery or science code."""
from pathlib import Path
import csv,json,hashlib,math,zipfile,struct,ast,subprocess,time
R=Path(__file__).resolve().parents[3];O=Path(__file__).resolve().parent;D=O.parent
STOP='3c7291bd6d4683796deccb352d1e8c105db63d87';SNAP='a580499475cc4148f1195e6256c4813c6e225726'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def emit(n,x):
 with (O/n).open('x',encoding='utf-8') as f:json.dump(x,f,indent=2,ensure_ascii=False,allow_nan=False);f.write('\n')
# Independent NPY decoder; no numpy.load and no regeneration of calibration RNG.
def decode(p,member):
 with zipfile.ZipFile(p) as z:b=z.read(member+'.npy')
 assert b[:6]==b'\x93NUMPY';nsize=2 if b[6]==1 else 4;n=struct.unpack('<H' if nsize==2 else '<I',b[8:8+nsize])[0];a=8+nsize+n;h=ast.literal_eval(b[8+nsize:a].decode('latin1'))
 assert h['descr']=='<f8' and not h['fortran_order'];vals=struct.unpack('<'+'d'*(len(b[a:])//8),b[a:]);return h['shape'],vals
checks=[]
def ck(k,v,detail=''):checks.append(dict(check=k,PASS=bool(v),detail=detail))
emit('INDEPENDENT_REPREFLIGHT_STARTED.json',dict(scope='historical identity/rule audit only; not formal science cold audit',formal_cold_audit_calls=0,new_panel=0))
t0=time.perf_counter();ck('HEAD_remains_STOP',subprocess.check_output(['git','rev-parse','HEAD'],cwd=R).decode().strip()==STOP)
for fn in ['DESIGN_FILE_MANIFEST.json','COMMIT_B_STOP_MANIFEST.json']:
 for x in read(D/fn)['files']:ck('preserved:'+x['path'],sha(D/x['path'])==x['sha256'])
ck('all_historical_tracked_content_unchanged',subprocess.run(['git','diff','--quiet',SNAP,'--','.',':(exclude)results/CZ_D2/**'],cwd=R).returncode==0)
old=rows(D/'HISTORICAL_STATE_EXCLUSIONS.csv');v2=rows(O/'HISTORICAL_STATE_EXCLUSIONS_V2.csv');ck('V1_preserved_all_columns',v2[:len(old)]==old)
fields=['r_km','theta_deg','v_mps','psi_deg'];csvcache={};jsoncache={}
for i,x in enumerate(old):
 p=R/x['source'];idx=x['record']
 if p.suffix=='.csv':
  if str(p) not in csvcache:csvcache[str(p)]=rows(p)
  rr=csvcache[str(p)][int(idx)]
  if 'truth_r0_km' in rr:keys=['truth_r0_km','truth_theta0_deg','truth_v','truth_psi_deg']
  elif 'truth_r_km' in rr:keys=['truth_r_km','truth_theta_deg','truth_v','truth_psi_deg']
  elif 'r0_km' in rr:keys=['r0_km','theta0_deg','v_mps','psi_deg']
  else:keys=fields
  expected=[float(rr[k]) if k in rr and rr[k]!='' else None for k in keys]
 else:
  if str(p) not in jsoncache:jsoncache[str(p)]=read(p)
  obj=jsoncache[str(p)]
  if 'states' in obj:
   st=next(a['state'] for a in obj['states'] if a['id']==idx);expected=[st[0]/1000,st[1]*180/math.pi,st[3],st[4]*180/math.pi]
  else:expected=next(a['h'] for a in obj['generating_cases'] if a['id']==idx)
 actual=[float(x[k]) if x[k] else None for k in fields]
 ck('original_source_numeric_identity:'+str(i),all((a is None and b is None) or (a is not None and b is not None and abs(a-b)<1e-11) for a,b in zip(actual,expected)))
aux='results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/evaluation_only/calibration_private.npz'
shape,flat=decode(R/aux,'h_truth');ck('AUX_layout',(shape==(512,4) and len(flat)==2048))
for i in range(512):
 row=v2[len(old)+i];ck('AUX_row:'+str(i),row['source']==aux and row['record']==str(i) and tuple(float(row[k]) for k in fields)==flat[4*i:4*i+4])
ck('AUX_no_duplicate_archived_truths',len({flat[i:i+4] for i in range(0,len(flat),4)})==512)
scores=rows(R/'results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/CALIBRATION_SCORES.csv');ck('AUX_all_indices_really_used',len(scores)==1024 and {int(x['calibration_index']) for x in scores}==set(range(512)))
t2='results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/CALIBRATION_GEOMETRIES.csv';tt=rows(R/t2)
for i,x in enumerate(tt):
 row=v2[728+i];ck('T2_static_row:'+str(i),row['source']==t2 and row['kind']=='static_position_projection' and row['v_mps']==row['psi_deg']=='' and float(row['r_km'])==float(x['r_km']) and float(row['theta_deg'])==float(x['theta_deg']))
ck('only_verified_additions',len(v2)==len(old)+512+len(tt))
# Guard implementation 1: normalized coordinate distances. No target-panel generator.
cfg=read(D/'DESIGN_FREEZE.json')['panel'];origin=cfg['grid_origin'];step=cfg['grid_step'];fraction=cfg['grid_fraction_guard'];w4=cfg['old_and_new_four_state_neighborhood_half_width'];w2=cfg['static_position_guard']
def test1(x,y):
 widths=w2 if y[2] is None else w4
 for sign in [1,-1]:
  b=[y[0],sign*y[1],y[2],None if y[3] is None else sign*y[3]]
  if all((abs((x[j]-b[j]+180)%360-180) if j in (1,3) else abs(x[j]-b[j]))<=widths[j] for j in range(len(widths))):return True
 return False
def test2(x,y):
 dimensions=[0,1] if y[2] is None else [0,1,2,3]
 limits=dict(zip(dimensions,w2 if len(dimensions)==2 else w4))
 for flip in (False,True):
  passed=True
  for j in dimensions:
   b=-y[j] if flip and j in (1,3) else y[j];delta=abs(x[j]-b)
   if j in (1,3):delta=math.degrees(math.atan2(math.sin(math.radians(delta)),math.cos(math.radians(delta))));delta=abs(delta)
   if delta>limits[j]+1e-13:passed=False
  if passed:return True
 return False
historical=[[float(x[k]) if x[k] else None for k in fields] for x in v2]
# All history against all history, no new state generation. Static missing dimensions are filled
# only in the comparison argument to test their independence; never persisted as target truth.
agree=True;selfpass=True;mirrorpass=True;gridagree=True
for row in historical:
 x=[a if a is not None else 0. for a in row];m=[x[0],-x[1],x[2],-x[3]]
 selfpass &= test1(x,row);mirrorpass &= test1(m,row)
 g1=any(abs((x[j]-origin[j])/step[j]-round((x[j]-origin[j])/step[j]))<=fraction for j in range(4))
 g2=any(abs(x[j]-(origin[j]+math.floor((x[j]-origin[j])/step[j]+.5)*step[j]))<=fraction*step[j] for j in range(4))
 gridagree &= g1==g2
 for y in historical:
  if test1(x,y)!=test2(x,y):agree=False
ck('all_744_squared_neighborhood_comparisons',agree,str(len(historical)**2));ck('joint_mirror_and_self_inclusion',selfpass and mirrorpass);ck('independent_coarse_grid_comparison',gridagree)
# Algebra-only boundary and circular controls, unrelated to future panel.
ck('circular_mirror_fixture',test1([10,179.999,8,179.999],[10,-180.001,8,-180.001]))
ck('static_ignores_motion_fixture',test1([10,2,999,99],[10,2,None,None]))
ck('all_coordinate_not_any_neighbor_fixture',not test1([10.2,2,8,5],[10,2,8,5]))
ck('joint_not_independent_mirror_fixture',not test1([10,2,8,-5],[10,2,8,5]))
# Private files are rehashed independently, CAP state rows decoded using standard library.
cap=R/'results/ENV_CAP_CAP06';pp={a['id']:a['state'] for a in read(cap/'P0_FROZEN_PANEL.json')['states']};rr={x['token']:x for x in rows(cap/'P2_RECORD_INDEX.csv')};ext=rows(O/'PRIVATE_ARCHIVE_IDENTITY_CHECKS.csv');identity=True;mapping=True;private_count=0
for x in ext:
 p=Path(x['path']);ok=p.is_file() and sha(p)==x['expected_sha256'];identity &=ok
 if 'CAP06' in str(p) and ok:
  shape,s=decode(p,'state');rec=rr[p.stem];private_count+=1
  if rec['stress']=='MAIN_TRACK':
   expect=list(pp[rec['state_id']]);sgn=1 if int(rec['repeat'])%2==0 else -1;expect[1]*=sgn;expect[4]*=sgn;mapping &=tuple(expect)==s
  else:mapping &=all(a==0 for a in s)
ck('external_seals_independent',identity,str(len(ext)));ck('CAP_private_roles_independent',mapping and private_count==4416)
# Independent source completeness judgement: do not infer azimuth for radial-only truth.
mfp=rows(R/'results/R3_RC3_method_recovery/MFP_UPPER_BOUND_summary.csv');unique=sorted(set((float(x['r_true_km']),float(x['z_true_m'])) for x in mfp))
code=ast.parse((R/'r3_0_method_recovery.py').read_text(encoding='utf-8-sig'));fn=next(x for x in ast.walk(code) if isinstance(x,ast.FunctionDef) and x.name=='pressure_spectrum');args=[a.arg for a in fn.args.args]
ck('range_only_gap_independently_confirmed',len(unique)==3 and args==['mode','freqs','r_m','z_s','z_r'] and not any('theta' in k or 'psi' in k for k in mfp[0]))
ck('V2_binding_verified',read(O/'HISTORY_REPAIR_INPUT_BINDING.json')['V2_sha256']==sha(O/'HISTORICAL_STATE_EXCLUSIONS_V2.csv'))
ck('no_new_panel_or_science_artifacts',not any((D/n).exists() for n in ['PANEL.csv','PUBLIC_OBSERVATIONS','MAIN_EXECUTION_STARTED.json','COLD_EXECUTION_STARTED.json']))
failed=[x for x in checks if not x['PASS']]
emit('INDEPENDENT_REPREFLIGHT_RESULT.json',dict(status='PASS' if not failed else 'FAIL',scope='Verified repair identities and rule invariance only; source completeness separately UNRESOLVED',checks=checks,PASS_count=sum(x['PASS'] for x in checks),FAIL_count=len(failed),source_completeness='UNRESOLVED',source_gap='range-only actual MFP generating truth has no frozen 4D/2D exclusion mapping; cannot certify remaining discovery inventory by implication',formal_scientific_admission=False,new_panel=0,new_innovations=0,optimizer_calls=0,Fisher_calls=0,MC_calls=0,formal_cold_calls=0,elapsed_s=time.perf_counter()-t0,STOP=True))
print(json.dumps(dict(PASS=sum(x['PASS'] for x in checks),FAIL=len(failed),source_completeness='UNRESOLVED',STOP=True)))
