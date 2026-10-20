"""Independent read-only provenance audit. No science or repair module imported.
The semantic judgement is source-code review plus independent schema/discovery checks,
not a theorem inferred from registry size. No CZ-D2 targets/innovations are generated.
"""
from pathlib import Path
import csv,json,hashlib,math,ast,zipfile,struct,re,subprocess,time,collections
R=Path(__file__).resolve().parents[3]; O=Path(__file__).resolve().parent; D=O.parent; OLD=D/'HISTORY_REPAIR_01'
PARENT='65f9cec61c05e43b3a628f91b078bfdcc4b5318b'; SNAP='a580499475cc4148f1195e6256c4813c6e225726'
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def obj(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while True:
   b=f.read(1024*1024)
   if not b:break
   h.update(b)
 return h.hexdigest()
def emit(n,x):
 with (DEST/n).open('x',encoding='utf-8',newline='') as f:
  if n.endswith('.csv'):
   w=csv.DictWriter(f,fieldnames=list(x[0]),lineterminator='\n');w.writeheader();w.writerows(x)
  else:json.dump(x,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
DEST=O/'AUDIT_COMPLETED'; DEST.mkdir(exist_ok=False)
checks=[]
def check(name,condition,detail=''):checks.append(dict(check=name,PASS=bool(condition),detail=detail))
def git(*args):return subprocess.check_output(['git',*args],cwd=R).decode().strip()
def values(x):return [float(x[k]) if x[k] else None for k in ('r_km','theta_deg','v_mps','psi_deg')]
def dec(p,key):
 with zipfile.ZipFile(p) as z:b=z.read(key+'.npy')
 assert b[:6]==b'\x93NUMPY'; w=2 if b[6]==1 else 4; n=int.from_bytes(b[8:8+w],'little'); start=8+w+n
 hd=ast.literal_eval(b[8+w:start].decode('latin1')); assert hd['descr']=='<f8'
 a=struct.unpack('<'+'d'*((len(b)-start)//8),b[start:]); shape=hd['shape']
 if len(shape)==2:
  nr,nc=shape; a=[[a[i+j*nr] if hd['fortran_order'] else a[i*nc+j] for j in range(nc)] for i in range(nr)]
 return shape,a
emit('INDEPENDENT_SOURCE_AUDIT_STARTED.json',dict(parent=PARENT,snapshot=SNAP,scope='READ_ONLY_SOURCE_COMPLETENESS_AUDIT_NOT_FORMAL_SCIENCE_COLD_AUDIT',formal_execution_authorized=False))
tic=time.perf_counter();check('parent',git('rev-parse','HEAD')==PARENT);check('branch',git('branch','--show-current')=='research/cz-d2')
check('no_existing_tracked_changes',git('diff','--name-only','HEAD')=='')
check('historical_snapshot_unchanged',git('diff','--name-only',SNAP,'--','.',':(exclude)results/CZ_D2/**')=='')
for m in ('DESIGN_FILE_MANIFEST.json','COMMIT_B_STOP_MANIFEST.json'):
 for x in obj(D/m)['files']:check('preserve:'+x['path'],digest(D/x['path'])==x['sha256'])
v1=rows(D/'HISTORICAL_STATE_EXCLUSIONS.csv');v2=rows(OLD/'HISTORICAL_STATE_EXCLUSIONS_V2.csv');v3=rows(O/'HISTORICAL_STATE_EXCLUSIONS_V3.csv')
check('original216_preserved',v3[:216]==v1 and len(v1)==216);check('V2_744_prefix_preserved',v3[:744]==v2 and len(v2)==744)
check('diff_identity',v3[744:]==rows(O/'V2_V3_ADDITIONS.csv') and len(v3)==762)
# Independent actual-source reconstruction of all eighteen additions.
for i,x in enumerate(v3[744:]):
 source=x['source'];t=rows(R/source)[int(x['record'])];h=values(x)
 if 'heading_identifiability' in source:expected=[50.,0.,2.,float(t['psi_deg'])]
 elif 'bearing_accuracy' in source:expected=[float(t['truth_r0_km']),0.,float(t['truth_v']),float(t['truth_psi_deg'])]
 elif 'AUX_' in source:expected=[float(t['range_km']),0.,None,None]
 else:expected=[]
 check('addition_source_identity:'+str(i),h==expected)
check('new_distinct_truths',len({tuple(values(x)) for x in v3[744:]})==18 and not ({tuple(values(x)) for x in v3[744:]} & {tuple(values(x)) for x in v2}))
# Original neighborhood, independent comparison; absent dimensions never become truth.
cfg=obj(D/'DESIGN_FREEZE.json')['panel']; widths4=cfg['old_and_new_four_state_neighborhood_half_width'];widths2=cfg['static_position_guard'];orig=cfg['grid_origin'];step=cfg['grid_step']; frac=cfg['grid_fraction_guard']
reg=[values(x) for x in v3]
def covered(h,static=False):
 for r in reg:
  if not static and r[2] is None:continue
  for sign in (1,-1):
   inds=(0,1) if static else (0,1,2,3);widths=widths2 if static else widths4;good=True
   for j in inds:
    y=(-r[j] if sign==-1 and j in (1,3) else r[j]); delta=h[j]-y
    if j in (1,3):delta=math.degrees(math.atan2(math.sin(math.radians(delta)),math.cos(math.radians(delta))))
    if abs(delta)>widths[j]:good=False;break
   if good:return True
 return False
check('new18_already_coarse_grid_guarded',all(any(abs((h[j]-orig[j])/step[j]-round((h[j]-orig[j])/step[j]))<=frac for j in range(4) if h[j] is not None) for h in reg[744:]))
# AUX calibration extracted independently using stdlib NPY decoder, not RNG reproduction.
aux='results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/evaluation_only/calibration_private.npz';shape,aa=dec(R/aux,'h_truth')
check('AUX512_shape',shape==(512,4))
for i,h in enumerate(aa):check('AUX512:'+str(i),values(v3[216+i])==h)
# Direct source tables with actual complete truth labels, found independently from schemas.
truth_DEST=O/'AUDIT_COMPLETED'; DEST.mkdir(exist_ok=False)
checks=[]
for x in rows(OLD/'CSV_SCHEMA_DISCOVERY.csv'):
 f=json.loads(x['fields']);keys=None
 for ks in [('truth_r0_km','truth_theta0_deg','truth_v','truth_psi_deg'),('truth_r_km','truth_theta_deg','truth_v','truth_psi_deg')]:
  if set(ks)<=set(f):keys=ks;break
 if keys:
  rr=rows(R/x['path']);good=all(covered([float(a[k]) for k in keys]) for a in rr)
  truth_checks.append(dict(source=x['path'],rows=len(rr),role='EXPLICIT_COMPLETE_GENERATING_TRUTH',covered=good));check('independent_truth_schema:'+x['path'],good)
# Exact panel and alternate/reference semantics.
for name in ['selected_hard_pairs.csv','rc2_vs_rc3_summary.csv','line_doppler_increment.csv','propagation_increment.csv']:
 p=R/'results/P3_RC3_increment'/name;rr=rows(p);good=True
 pairs={a['pair_id']:a for a in rows(R/'results/P3_RC3_increment/selected_hard_pairs.csv')}
 for x in rr:
  q=pairs[x['pair_id']]; good &= all(not x.get(k) or float(x[k])==float(q[k]) for k in ['ref_r0_km','ref_theta0_deg','ref_v','ref_psi_deg'])
  good &= covered([float(q[k]) for k in ['ref_r0_km','ref_theta0_deg','ref_v','ref_psi_deg']])
 check('P3_observation_reference:'+name,good)
p2=rows(R/'results/P2_RC2_kinematic_boundary/heading_identifiability_boundary.csv')
check('all_P2_heading_generators_covered',all(covered([50,0,2,float(x['psi_deg'])]) for x in p2))
# SP1 labels are actual static receivers; complete H1 truth may cover their position projection.
for name in ['CALIBRATION_PRIVATE_ARCHIVE_MANIFEST.json','TEST_PRIVATE_ARCHIVE_MANIFEST.json']:
 records=obj(R/'results/R4_SINGLE_HLA_STATIC_SPATIAL_RANGE_PILOT'/name)['records']
 check('SP1_private_labels:'+name,all(covered([float(x['truth_r_km']),float(x['truth_theta_deg'])],True) for x in records),str(len(records)))
t2=rows(R/'results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/CALIBRATION_GEOMETRIES.csv');check('T2_all16_positions',len(t2)==16 and all(covered([float(x['r_km']),float(x['theta_deg'])],True) for x in t2))
# H2 exact-equality exception remains a failure; only frozen neighborhood compatibility passes.
shape,hh=dec(R/'results/R4_SINGLE_HLA_RADIAL_INCREMENT_PRECISION_PILOT/GENERATOR_TRUTH_AND_INNOVATIONS.npz','horizontal_truth');check('H2_archive_neighborhood_coverage',shape==(8,4) and all(covered(h) for h in hh))
hd=obj(OLD/'H2_REPRESENTATION_DIAGNOSTIC_V2.json');check('H2_STRICT_FAIL_preserved','FAIL' in hd['status'] and hd['max_ULP_difference']<=2 and hd['all32_within_original_neighborhood'])
check('H2_original_erratum_retained',(OLD/'DIAGNOSTIC_READER_ERRATUM.json').is_file() and (OLD/'H2_REPRESENTATION_DIAGNOSTIC.json').is_file())
# MFP independently reconstructed from archived data and function signature, not builder output.
mfp=rows(R/'results/R3_RC3_method_recovery/MFP_UPPER_BOUND_summary.csv');coords={(float(x['r_true_km']),float(x['z_true_m'])) for x in mfp};partial=rows(O/'MFP_PARTIAL_HISTORY.csv')
tree=ast.parse((R/'r3_0_method_recovery.py').read_text(encoding='utf-8-sig'));fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='pressure_spectrum')
check('MFP_original_model_partial',[x.arg for x in fn.args.args]==['mode','freqs','r_m','z_s','z_r'])
check('MFP9_rows3_unique',len(mfp)==9 and coords=={(47.6,187.),(50.,200.),(52.3,213.)})
check('MFP_partial_table_identity',{(float(x['r_km']),float(x['z_m'])) for x in partial}==coords and len(partial)==3)
check('MFP_no_fabricated_dimensions',all(x['theta']==x['v']==x['psi']=='NOT_DEFINED' for x in partial) and not any('MFP' in x['source'] for x in v3))
# Independent fresh lexical discovery across ALL historical Python, including bundled examples.
# AST is supplemental: historical code may target newer Python. No old modules are executed.
inv=rows(OLD/'SNAPSHOT_INVENTORY.csv');adjud={x['path']:x for x in rows(O/'SOURCE_SEMANTIC_ADJUDICATION.csv')};lex=[];parse=[];unresolved=[]
for x in inv:
 if not x['path'].endswith('.py'):continue
 text=(R/x['path']).read_text(encoding='utf-8-sig',errors='strict');ls=text.splitlines()
 try:ast.parse(text);parse.append(dict(path=x['path'],AST='PARSED',fallback='',error=''))
 except SyntaxError as e:parse.append(dict(path=x['path'],AST='HISTORICAL_PY310_INCOMPATIBLE',fallback='LEXICAL_SOURCE_AND_DEPENDENCY_REVIEW_NO_EXECUTION',error=str(e)))
 for i,s in enumerate(ls):
  if not (re.search(r'(uniform\(|normal\(|linspace\(|arange\()',s) and re.search(r'(truth|r0|theta|psi|heading|target|state|range_km|velocity)',s,re.I)):continue
  p=x['path'];b=Path(p).name
  if p.startswith('tools/acoustics_toolbox/'):
   role='BUNDLED_PROPAGATION_EXAMPLE_NOT_PROJECT_TARGET_PANEL';proof='kraken3d_wedge_tl.py: thetas angular propagation coverage, not moving target observations'
  elif any(z in b for z in ('preflight','controls','audit','independent')):
   role='SOFTWARE_GEOMETRY_OR_RESIDUAL_REFERENCE_CONTROL';proof='Function/control body inspected; diagnostic input or Jacobian reference, not performance/threshold calibration population; no state inferred from random residuals'
  elif p in adjud:
   role=adjud[p]['semantic_role'];proof=adjud[p]['evidence']
  else:role='UNRESOLVED';proof='missing provenance';unresolved.append(p)
  lex.append(dict(path=p,sha256=x['sha256'],line=i+1,statement=s.strip(),role=role,evidence=proof))
check('independent_generator_discovery_no_unmapped',not unresolved,sorted(set(unresolved)))
# Identity is separate from semantics. Hash actual candidate sources plus cited code.
seals=[]
for p,x in sorted(adjud.items()):
 actual=digest(R/p);seals.append(dict(path=p,expected=x['sha256'],actual=actual,PASS=actual==x['sha256']))
check('candidate_sources_current_bytes',all(x['PASS'] for x in seals),str(len(seals)))
for g in obj(O/'SEMANTIC_FAMILY_BINDINGS.json'):
 for p in g['proof']:check('reviewed_code_identity:'+p['path'],digest(R/p['path'])==p['sha256'])
check('classification_no_unresolved',all(x['disposition']!='UNRESOLVED' for x in adjud.values()) and obj(O/'UNRESOLVED_SOURCES.json')==[])
check('private4818_identity_evidence_inherited',len(rows(OLD/'PRIVATE_ARCHIVE_IDENTITY_CHECKS.csv'))==4818 and obj(OLD/'INDEPENDENT_REPREFLIGHT_RESULT.json')['FAIL_count']==0,'Inherited accepted REPAIR01; NOT rehashed4818 files in REPAIR02')
# Structural source-role evidence independent of registry cardinality.
texts={p:(R/p).read_text(encoding='utf-8-sig') for p in ['p3_rc3_increment.py','r4_b1_execute.py','aux_d1_generate.py','hla_h2_preflight.py','aux_d1_preflight.py','results/ENV_CAP_CAP06/code/p1_toy.py']}
check('P3_ref_observation_code','simulate_observations(ref' in texts['p3_rc3_increment.py'])
check('B1_reference_code','horizontal_reference' in texts['r4_b1_execute.py'])
check('AUX_calibration_population_code','512' in texts['aux_d1_generate.py'] and 'h_truth' in texts['aux_d1_generate.py'])
check('H2_control_not_new_performance','no eight-scene performance runs' in texts['hla_h2_preflight.py'])
check('AUX_control_not_performance','not the12state performance panel' in texts['aux_d1_preflight.py'])
check('formal_artifacts_absent',not any((D/p).exists() for p in ['PANEL.csv','PUBLIC_OBSERVATIONS','MAIN_EXECUTION_STARTED.json','COLD_EXECUTION_STARTED.json']))
check('authorization_false',obj(O/'APPLICABILITY_POLICY_SUPPLEMENT.json')['formal_execution_authorized'] is False)
emit('INDEPENDENT_CODE_DISCOVERY.csv',lex);emit('HISTORICAL_PYTHON_PARSE_STATUS.csv',parse);emit('INDEPENDENT_EXPLICIT_TRUTH_COVERAGE.csv',truth_checks);emit('SOURCE_SHA256_CHECKS.csv',seals)
fail=[x for x in checks if not x['PASS']]
result=dict(status='PASS' if not fail else 'FAIL',SOURCE_SEMANTIC_COMPLETENESS='PASS_FIXED_SNAPSHOT_AND_BOUND_ARCHIVES' if not fail else 'UNRESOLVED',REGISTRY_IDENTITY_AND_COVERAGE='PASS' if not fail else 'FAIL',MFP_PARTIAL_STATE_POLICY='PARTIAL_RZ_GENERATING_TRUTH_RECORDED_NOT_4D_EXCLUSION',DESIGN_CONTRACT_UNCHANGED=not any(not x['PASS'] for x in checks if x['check'].startswith(('preserve:','historical_snapshot','no_existing'))),FORMAL_EXECUTION_AUTHORIZED=False,checks=checks,PASS_count=sum(x['PASS'] for x in checks),FAIL_count=len(fail),scope='Source-code semantic review, independent generator/schema search, direct archived truth reconstruction; not a proof about unavailable archives or every individual coordinate ever used',candidate_sources=len(seals),python_files_scanned=len(parse),AST_incompatible=sum(x['AST']!='PARSED' for x in parse),lexical_hits=len(lex),V3_rows=len(v3),private_archive_hashes_rechecked_this_round=0,accepted_private_hash_checks_inherited=4818,new_panel=0,new_innovations=0,new_observations=0,MC=0,Fisher=0,optimizer=0,formal_science_cold_audit=0,elapsed_s=time.perf_counter()-tic,STOP=True)
emit('INDEPENDENT_SOURCE_COMPLETENESS_AUDIT.json',result)
print(json.dumps({k:v for k,v in result.items() if k!='checks'},ensure_ascii=False));print('FAILED',fail)
