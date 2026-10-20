"""Controlled registry repair; archived inputs only. No research-code imports/RNG."""
from pathlib import Path
import csv,json,hashlib,math,time,subprocess,collections
import numpy as np
R=Path(__file__).resolve().parents[3];O=Path(__file__).resolve().parent;D=O.parent
A='dc91982e67921eeb7a6aac4e12c482a980d5ad89';STOP='3c7291bd6d4683796deccb352d1e8c105db63d87'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def emit(n,x):
 with (O/n).open('x',encoding='utf-8',newline='') as f:
  if n.endswith('.json'):json.dump(x,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
  else:
   w=csv.DictWriter(f,fieldnames=list(x[0]),lineterminator='\n');w.writeheader();w.writerows(x)
emit('REPAIR_STARTED.json',dict(parent=STOP,original_A=A,panel_generated=False,formal_execution=False))
t0=time.perf_counter();assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=R).decode().strip()==STOP
inv={x['path']:x for x in rows(O/'SNAPSHOT_INVENTORY.csv')}
old=rows(D/'HISTORICAL_STATE_EXCLUSIONS.csv');assert len(old)==216
new=[dict(x) for x in old];cov=[];checks=[];external=[]
def ck(k,v,detail=''):
 checks.append(dict(check=k,PASS=bool(v),detail=detail))
def source(p,category,shape,units,index,reason,evidence,status='RESOLVED'):
 p=str(p).replace('\\','/');cov.append(dict(path=p,sha256=sha(R/p),purpose=category,data_shape=shape,units=units,truth_index=index,decision=status,reason=reason,evidence=evidence))
for p in dict.fromkeys(x['source'] for x in old):
 rr=[x for x in old if x['source']==p]
 source(p,'HISTORICAL_GENERATING_STATES',str(len(rr))+' registered records','km,degree,m/s,degree; static leaves v/psi empty','original registry record column','INCLUDE: preserve original rows byte-equivalent as CSV values','CZ-D2 original source identities and original panel policy')
aux='results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/evaluation_only/calibration_private.npz'
with np.load(R/aux,allow_pickle=False) as z:h=z['h_truth'].copy()
ck('AUX512_shape_unique_finite',h.shape==(512,4) and np.isfinite(h).all() and len(set(map(tuple,h)))==512)
for i,x in enumerate(h):new.append(dict(source=aux,record=str(i),kind='horizontal_four_state',**dict(zip(['r_km','theta_deg','v_mps','psi_deg'],map(repr,map(float,x))))))
source(aux,'THRESHOLD_CALIBRATION_TRUTH','h_truth (512,4)','km,degree,m/s,degree','h_truth[0:512,:]','INCLUDE: calibration states were used to generate 1024 resource scores','aux_d1_generate.py:44; THRESHOLD_FREEZE.json; CALIBRATION_SCORES.csv')
t2='results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/CALIBRATION_GEOMETRIES.csv';cal=rows(R/t2)
ck('T2_static16_unique',len(cal)==16 and len({(x['r_km'],x['theta_deg']) for x in cal})==16)
for i,x in enumerate(cal):new.append(dict(source=t2,record=str(i),kind='static_position_projection',r_km=x['r_km'],theta_deg=x['theta_deg'],v_mps='',psi_deg=''))
source(t2,'STATIC_CALIBRATION_POSITION','16 rows; r,theta,z','km,degree,m; v/psi absent','CSV data row0..15','INCLUDE: used for noise-reference B0 fields and cross-mesh CAL fields despite later receiver STOP','t2_tilt_g0_prepare.py:17-25; t2_tilt_g0_forward.py:56-59; NOISE_REFERENCE_FREEZE.json; REPLICA_FREEZE.json')
# Verify existing private generator archives without reconstructing any observations.
cap=R/'results/ENV_CAP_CAP06';panel={x['id']:x for x in read(cap/'P0_FROZEN_PANEL.json')['states']};records={x['token']:x for x in rows(cap/'P2_RECORD_INDEX.csv')}
for row in rows(cap/'P2_RAW_DATA_INDEX.csv'):
 if row['domain']!='truth_private':continue
 p=Path(row['path']);exists=p.is_file();digest=sha(p) if exists else '';ok=exists and digest==row['sha256'] and p.stat().st_size==int(row['bytes']);identity=False;purpose='UNRESOLVED'
 if ok:
  with np.load(p,allow_pickle=False) as z:s=z['state'].copy()
  record=records[p.stem]
  if record['stress']=='MAIN_TRACK':
   e=np.array(panel[record['state_id']]['state']);sign=1 if int(record['repeat'])%2==0 else -1;e[[1,4]]*=sign
   identity=bool(np.array_equal(s,e));purpose='REUSED_PANEL_AND_JOINT_MIRROR'
  else:
   identity=bool(np.array_equal(s,np.zeros(5)));purpose='NON_CV_DIRECTION_STRESS_ZERO_STATE_PLACEHOLDER_NOT_TARGET_TRUTH'
 external.append(dict(path=str(p),sha256=digest,expected_sha256=row['sha256'],bytes=p.stat().st_size if exists else 0,identity_PASS=ok,values_PASS=identity,purpose=purpose))
ck('CAP06_all4416_private_seals_and_state_roles',len(external)==4416 and all(x['identity_PASS'] and x['values_PASS'] for x in external))
source('results/ENV_CAP_CAP06/P2_RAW_DATA_INDEX.csv','PRIVATE_ARCHIVE_INDEX','4416 private +4416 public references','state=[m,rad,m,m/s,rad]','private/state; main-track panel vs stress placeholder','REUSE: main-track states map to original16 with joint mirrors; stress zero state is not a CV target; direction stress has r=50km coarse-grid protection at t0','code/generate_p2.py:42-76; P2_RECORD_INDEX.csv; P0_FROZEN_PANEL.json')
source('results/ENV_CAP_CAP06/P1_RAW_DATA_INDEX.csv','DIMENSIONLESS_MATH_TOY','960 configurations; three archive roles','dimensionless abstract coefficients, no physical horizontal state','ideal_shape_db (9,4)','EXCLUDE: abstract algebra toy has no r/theta/v/psi generating truth','code/p1_toy.py:generate/basis')
# SP1: actual static record labels must match known panel, not treated as fresh truths.
sp=R/'results/R4_SINGLE_HLA_STATIC_SPATIAL_RANGE_PILOT';seenpaths=set();spcounts={}
for phase in ['CALIBRATION','TEST']:
 manifest=read(sp/(phase+'_PRIVATE_ARCHIVE_MANIFEST.json'))['records'];pp=rows(sp/'CALIBRATION_PANEL.csv') if phase=='CALIBRATION' else rows(R/'results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/OFFGRID_TRUTH_EVALUATION_ONLY.csv')
 good=True
 for row in manifest:
  g=pp[int(row['geometry'])];good &= all(abs(float(row[a])-float(g[b]))<1e-12 for a,b in [('truth_r_km','r_km'),('truth_theta_deg','theta_deg'),('truth_z_m','z_label')])
  p=Path(row['private_innovations'])
  if str(p) in seenpaths:continue
  seenpaths.add(str(p));exists=p.is_file();digest=sha(p) if exists else '';ok=exists and digest==row['private_sha256']
  external.append(dict(path=str(p),sha256=digest,expected_sha256=row['private_sha256'],bytes=p.stat().st_size if exists else 0,identity_PASS=ok,values_PASS=True,purpose='SP1_SOURCE_GAIN_NOISE_INNOVATIONS;truth_labels_in_sealed_manifest'))
 spcounts[phase]=len(manifest);ck('SP1_'+phase+'_labels_to_known_panel',good)
 source(str((sp/(phase+'_PRIVATE_ARCHIVE_MANIFEST.json')).relative_to(R)),'STATIC_RECEIVER_PRIVATE_LABELS',str(len(manifest))+' records','r km,theta degree,z m','records[*].truth_r_km/truth_theta_deg/truth_z_m','REUSE: labels equal original calibration16 or H1 offgrid8; no new random positions','hla_sp1_generate.py; panel CSV; private innovation seals')
# T2 calibration archived field identity, no field recomputation.
freeze=read(R/'results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/REPLICA_FREEZE.json')
for name,digest in freeze['files'].items():
 if not name.startswith('CAL_'):continue
 p=Path('D:/ProjectStorage/WSL-CZ/T2_TILT_G0/replicas')/name;exists=p.is_file();actual=sha(p) if exists else '';shape=''
 if exists:
  a=np.load(p,mmap_mode='r',allow_pickle=False);shape=str(a.shape)
 external.append(dict(path=str(p),sha256=actual,expected_sha256=digest,bytes=p.stat().st_size if exists else 0,identity_PASS=exists and actual==digest,values_PASS=exists and a.shape[0]==16,purpose='T2_STATIC_CALIBRATION_FIELDS '+shape))
# Confirm H2 and AUX test reuse canonical existing truths, no RNG.
p=R/'results/R4_SINGLE_HLA_RADIAL_INCREMENT_PRECISION_PILOT/GENERATOR_TRUTH_AND_INNOVATIONS.npz'
with np.load(p,allow_pickle=False) as z:q=z['horizontal_truth']
canonical=rows(R/'results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/OFFGRID_TRUTH_EVALUATION_ONLY.csv');expected=np.array([[float(x[k]) for k in ['r_km','theta_deg','v_mps','psi_deg']] for x in canonical])
ck('H2_eight_existing_truths',np.array_equal(q,expected))
source(p.relative_to(R),'DEVELOPMENT_FEATURE_TRUTH','horizontal_truth (8,4)','km,deg,m/s,deg','horizontal_truth[0:8]','REUSE: equals H1 offgrid8','hla_h2_generate.py and H1 OFFGRID_TRUTH_EVALUATION_ONLY.csv')
p=R/'results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/evaluation_only/test_private.npz'
with np.load(p,allow_pickle=False) as z:q=z['h_truth']
canonical=rows(R/'results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_TRUTH_PANEL.csv');expected=np.array([[float(x[k]) for k in ['r0_km','theta0_deg','v_mps','psi_deg']] for x in canonical])
ck('AUX_test_existing12_tiled128',np.array_equal(q,np.tile(expected,(128,1))))
source(p.relative_to(R),'PERFORMANCE_TEST_TRUTH','h_truth (1536,4)','km,deg,m/s,deg','h_truth[0:1536]','REUSE: exact tiled canonical12; new noises do not make new targets','aux_d1_generate.py:63; A1_NEW_TRUTH_PANEL.csv')
# Source classification boundary requiring research-lead adjudication, not guessed coordinates.
mfp='results/R3_RC3_method_recovery/MFP_UPPER_BOUND_summary.csv';mr=rows(R/mfp);radial=sorted({(float(x['r_true_km']),float(x['z_true_m'])) for x in mr})
source(mfp,'ACTUAL_RANGE_DEPTH_GENERATING_TRUTH',str(len(mr))+' evaluation rows; '+str(len(radial))+' unique (r,z)','r km,z m; theta/v/psi absent','r_true_km,z_true_m','UNRESOLVED: actual observations generated at radial offgrid truth; frozen static guard requires (r,theta), no rule for range-only truth. Do not invent theta=0 or treat all orientations as used.','r3_0_method_recovery.py:240-255,281-300; MFP_UPPER_BOUND_rz.csv is search grid, not truth','UNRESOLVED')
source('results/R3_RC3_method_recovery/MFP_UPPER_BOUND_rz.csv','SEARCH_CANDIDATE_GRID','r,z,J columns','km,m; scores','NONE','EXCLUDE: replica search surface, not generating truth; actual generating labels in summary','r3_0_method_recovery.py:299-371')
source('results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/CATALOGUE_primary.npy','SEARCH_CANDIDATE_GRID','complete registered r/theta/z catalogue','km,degree,m','NONE','EXCLUDE: deterministic replica grid; not random calibration/test generating states','t2_tilt_g0_forward.py:complete catalogues; separate calibration loop')
source('results/ENV_CAP_CAP06/code/generate_p2.py','CALIBRATION_AND_STRESS_GENERATOR','4 orientation references;8 stress directions','degrees; no calibration range/CV state','make_cal and stress block','EXCLUDE from four-state registry: calibration angles contain no target position. Stress state=zeros is placeholder; do not register fabricated zero range/speed','source code, sealed private archive state role checks')
ck('all_external_checked_identities',all(x['identity_PASS'] and x['values_PASS'] for x in external))
assert len(new)==744
emit('HISTORICAL_STATE_EXCLUSIONS_V2.csv',new)
emit('SOURCE_COVERAGE_ADJUDICATION.csv',cov)
emit('PRIVATE_ARCHIVE_IDENTITY_CHECKS.csv',external)
emit('REPAIR_BUILD_CHECKS.json',dict(checks=checks,PASS=all(x['PASS'] for x in checks),external_files=len(external),elapsed_s=time.perf_counter()-t0))
emit('V1_V2_DIFF.json',dict(old_rows=216,new_rows=744,preserved_rows=216,additions=[dict(source=aux,rows=512,kind='horizontal_four_state'),dict(source=t2,rows=16,kind='static_position_projection')],old_sha256=sha(D/'HISTORICAL_STATE_EXCLUSIONS.csv'),new_sha256=sha(O/'HISTORICAL_STATE_EXCLUSIONS_V2.csv'),no_dedup_of_original_provenance_rows=True,no_new_panel=True))
emit('SOURCE_COMPLETENESS_GAPS.json',dict(status='UNRESOLVED',reason='Historical radial-only generating truths cannot be assigned an inclusion/exclusion mapping under frozen 4D/2D rule without a policy decision.',source=mfp,unique_radial_truths=radial,not_new_experiment=True,action='STOP; preserve verified V2 additions; no automatic release',scope='Discovery completed across snapshot; semantic adjudication stops at unresolved source. Remaining discovered hits are retained for coverage review, not silently marked PASS.',unclassified_discovery_hits='See full metadata and source usage inventories; these are discovery evidence, not certified source classifications.'))
emit('HISTORY_REPAIR_INPUT_BINDING.json',dict(original_A=A,stop_parent=STOP,old_registry='results/CZ_D2/HISTORICAL_STATE_EXCLUSIONS.csv',proposed_registry=str((O/'HISTORICAL_STATE_EXCLUSIONS_V2.csv').relative_to(R)).replace('\\','/'),old_registry_sha256=sha(D/'HISTORICAL_STATE_EXCLUSIONS.csv'),V2_sha256=sha(O/'HISTORICAL_STATE_EXCLUSIONS_V2.csv'),status='PROPOSED_V2_SOURCE_COMPLETENESS_UNRESOLVED',formal_execution_authorized=False,panel_materialization_authorized=False,unchanged_scientific_contract_sha256=sha(D/'DESIGN_FREEZE.json'),rules='Original 4D neighborhood, static2D projection, joint theta/psi mirror and coarse guard unchanged; no new range-only rule',new_panel=0,new_observations=0,new_MC=0))
print(json.dumps(dict(V2_rows=len(new),external_checked=len(external),checks_PASS=all(x['PASS'] for x in checks),completeness='UNRESOLVED',STOP=True)))
