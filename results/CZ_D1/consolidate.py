"""CZ-D1 deterministic read-only consolidation. Never imports historical science code."""
from pathlib import Path
import csv,json,math,hashlib,subprocess,collections,re,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]; OUT=Path(__file__).resolve().parent
PARENT='ca86d324afc1689af6c96cd35145e03625086bae'
SOURCES={}; CHECKS=[]; MASTER=[]; QUANT=[]; PAIRED=[]
def git(*a): return subprocess.check_output(['git','-c','safe.directory='+str(ROOT),'-c','core.quotepath=false',*a],cwd=ROOT).decode('utf-8').strip()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def source(p):
 p=Path(p); p=p if p.is_absolute() else ROOT/p
 k=p.relative_to(ROOT).as_posix() if p.is_relative_to(ROOT) else str(p)
 if k not in SOURCES:SOURCES[k]={'path':k,'bytes':p.stat().st_size,'sha256':sha(p),'version_sha':PARENT,'scope':'saved artifact; not a fresh science run'}
 return p
def rows(p):return list(csv.DictReader(source(p).open(encoding='utf-8-sig',newline='')))
def obj(p):return json.loads(source(p).read_text('utf-8-sig'))
def write(name,data):
 if name.endswith('.json'):(OUT/name).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8');return
 fields=list(dict.fromkeys(k for r in data for k in r))
 with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fields);w.writeheader();w.writerows(data)
def chk(name,ok,detail=''):CHECKS.append({'check':name,'PASS':bool(ok),'detail':str(detail)})
def eq(name,a,b,tol=1e-10):chk(name,(a==b) or (math.isfinite(a) and math.isfinite(b) and abs(a-b)<=tol*max(1,abs(b))),f'{a} vs {b}')
def nr(a,q):return sorted(a)[math.ceil(q*len(a))-1] if a else None
def fmt(x):return 'INF' if x is not None and math.isinf(x) else ('' if x is None else x)
def base(exp,arch,path,**kw):
 d=dict(experiment_id=exp,phase=exp,architecture_id=arch,version_sha=PARENT,evidence_path=str(Path(path).parent).replace('\\','/'),result_file=path,
 scenario='',r_true_km='',theta_true_deg='',z_true_m='',v_true_mps='',psi_true_deg='',N_elements='',element_spacing_m='',vertical_channels='',aux_count='',aux_baseline_km='',turn_deg='',turn_start_s='',observation_s='',epoch_count='',signal_S_class='',frequencies_Hz='',sigma_bearing_deg='',bearing_bias_bound_deg='',navigation_sigma_m='',channel_gain_contract='',feature_sigma_and_unit='',extraction_source='',environmental_contract='',matched_or_mismatched='',oracle_flags='',estimator_family='',offgrid_yes_no='',n_total='',n_independent='',n_rejected='',rejection_rule='',metrics_scope='',r_P50_percent='',r_P95_percent='',theta_P50_deg='',theta_P95_deg='',z_P50_m='',z_P95_m='',v_P50_percent='',v_P95_percent='',v_P50_mps='',v_P95_mps='',psi_P50_deg='',psi_P95_deg='',joint_status='NOT_DEFINED',correct_support_retention='',cold_audit_status='SOURCE_REPORTED_NOT_RERUN',lead_audit_status='PENDING_OR_SOURCE_SCOPE_ONLY',gate='',claim_scope='',missing_inputs='',comments='',evidence_category='NOT_EVALUATED',precision_status='NOT_EVALUATED',z_status='NOT_EVALUATED')
 d.update(kw);MASTER.append(d);return d
# Inventory all historical result files, including nested stages and ancillary figures/data.
files=git('ls-tree','-r','--name-only',PARENT,'--','results').splitlines()
files=[p for p in files if not p.startswith('results/CZ_D1/')]
log=git('log',PARENT,'--format=@@%H','--name-only','--','results').splitlines(); changed={};current=''
for s in log:
 if s.startswith('@@'):current=s[2:]
 elif s.startswith('results/'):changed.setdefault(s,current)
write('CZ_D1_FILE_INVENTORY.csv',[{'path':p,'snapshot':PARENT,'last_change_sha':changed.get(p,''),'bytes':(ROOT/p).stat().st_size,'status':'PRESENT'} for p in files])
dirs={str(Path(p).parent).replace('\\','/') for p in files if any(x in Path(p).name.upper() for x in ['GPT_SYNC','DECISION','ADJUDICATION','FINAL_REPORT','WORKING_ENVELOPE'])}
dirs|={str(Path(p).parents[len(Path(p).parts)-3]).replace('\\','/') for p in files if len(Path(p).parts)>=3}
inv=[]
for d in sorted(dirs):
 direct=[p for p in files if str(Path(p).parent).replace('\\','/')==d]
 anchors=[p for p in direct if any(t in Path(p).name.upper() for t in ['GPT_SYNC','DECISION','ADJUDICATION','FINAL_REPORT','WORKING_ENVELOPE']) and Path(p).suffix.lower() in ['.json','.md']]
 excerpts=[]
 for p in anchors:
  s=source(p).read_text('utf-8-sig',errors='replace')
  if p.endswith('.json'):
   try:
    j=json.loads(s)
    if isinstance(j,dict):excerpts.extend(f'{Path(p).name}:{k}={v}' for k,v in j.items() if isinstance(v,(str,bool,int,float)) and any(t in k.lower() for t in ['decision','class','status','gate','stop','verdict','result','adjudication']))
   except Exception:pass
  else:excerpts.extend(f'{Path(p).name}:{v[:350]}' for v in s.splitlines() if any(t in v for t in ['结论','判定','NOT_EVALUATED','BLOCKED','SUPPORTED','NOT_ESTABLISHED','INCOMPLETE','STOP']))
 inv.append({'stage_path':d,'snapshot_sha':PARENT,'latest_entry_change_sha':max((changed.get(p,'') for p in anchors),default=''),'entrypoints':'|'.join(anchors),'tracked_files_recursive':sum(p.startswith(d+'/') for p in files),'recorded_status_verbatim':' | '.join(excerpts[:16]),'review_status':'INDEXED_NOT_ALL_METRICS_RECALCULATED','execution_status':'SEE_EXPLICIT_SOURCE_STATUS' if excerpts else 'STATUS_UNRESOLVED','missing':'NO_TERMINAL_ENTRYPOINT' if not anchors else ''})
write('CZ_D1_STAGE_COMMIT_STATUS.csv',inv)
# Recompute errors from saved per-record values, never run estimates.
H2='results/R4_SINGLE_HLA_RADIAL_INCREMENT_PRECISION_PILOT/'
AUX='results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/'
T3='results/R4_T3_W1_WIDE_BASELINE_DEMO/'
H3='results/R4_H3_G0_FINITE_VERTICAL_RESPONSE/'
CAP='results/ENV_CAP_CAP06/'
configs=[]
h2=rows(H2+'CONFIGURATION_RESULTS.csv'); hf=obj(H2+'DESIGN_FREEZE.json')
a=rows(AUX+'ESTIMATES.csv');t=rows(T3+'RUN_ERRORS.csv')
# Metric dictionaries use percent, degree, m/s. Rejection always INF for unconditional.
for exp,data,key,case,accepted,cols,path,arch in [
 ('H2',h2,lambda r:f"{r['method']}_{float(r['sigma']):g}",'geometry',lambda r:r['status']=='ACCEPTED',{'r':('relative_range_error',100),'theta':('absolute_bearing_error_deg',1),'v':('relative_speed_error',100),'psi':('absolute_heading_error_deg',1)},H2+'CONFIGURATION_RESULTS.csv','A2'),
 ('AUX-D1',a,lambda r:r['resource'],'case_id',lambda r:r['accepted']=='True',{'r':('r_error',100),'theta':('theta_error',1),'v':('v_error',100),'psi':('psi_error',1)},AUX+'ESTIMATES.csv','A5'),
 ('T3-W1',t,lambda r:r['config'],'case_id',lambda r:r['accepted']=='True',{'r':('range_rel',100),'theta':('bearing_deg',1),'v':('speed_rel',100),'v_abs':('speed_abs_mps',1),'psi':('heading_deg',1)},T3+'RUN_ERRORS.csv','A5')]:
 groups=collections.defaultdict(list)
 for r in data:groups[key(r),r[case]].append(r)
 byconfig=collections.defaultdict(list)
 for (cfg,scene),rr in sorted(groups.items()):
  qrow={'experiment':exp,'configuration':cfg,'scene':scene,'n':len(rr),'n_accepted':sum(accepted(r) for r in rr)}
  for scope in ['unconditional','conditional']:
   for m,(c,mul) in cols.items():
    vals=[float(r[c])*mul if accepted(r) else math.inf for r in rr] if scope=='unconditional' else [float(r[c])*mul for r in rr if accepted(r)]
    for tag,pct in [('P50',.5),('P90',.9),('P95',.95),('P99',.99)]:qrow[f'{m}_{scope}_{tag}']=fmt(nr(vals,pct))
  QUANT.append(qrow);byconfig[cfg].append(qrow)
  for scope in ['unconditional','conditional']:
   d=base(exp,('A0' if cfg in ['B_0','U0'] else arch),path,scenario=scene,configuration=cfg,n_total=len(rr),n_independent=len(rr),n_rejected=len(rr)-qrow['n_accepted'],rejection_rate=(len(rr)-qrow['n_accepted'])/len(rr),metrics_scope=scope+' per scene nearest-rank; P99 descriptive reaggregation',N_elements=8,element_spacing_m=2,vertical_channels=0,aux_count=0 if exp=='H2' or cfg=='U0' else 1,turn_deg=15,turn_start_s=600,observation_s=1200,epoch_count=121,offgrid_yes_no='YES_DEVELOPMENT_GEOMETRIES',extraction_source='SYNTHETIC_FEATURE_OBSERVATIONS_NOT_RECEIVER_EXTRACTION',oracle_flags='truth generation/evaluation only; signed bearing assumed',rejection_rule='Frozen acceptance; rejected/numerical failure INF unconditional; no trimming',evidence_category='INJECTED_FEATURE_VALUE' if exp=='H2' and cfg!='B_0' else 'CONTINUOUS_ESTIMATOR_ERROR',precision_status='RECALCULATED_SAVED_OUTPUT_ERRORS',claim_scope='IDEAL_CONDITIONAL',missing_inputs='Actual receiver extraction; side resolution; calibrated hardware; depth',lead_audit_status='LEAD_AUDITED_SCOPED_PER_USER; not engineering acceptance',pairing='Same base innovations across resource/method/quality; scene samples not additional independent geometries',gate='SCOPED_CONDITIONAL_RESULT')
   for m,col in [('r','r_P{}_percent'),('theta','theta_P{}_deg'),('v','v_P{}_percent'),('v_abs','v_P{}_mps'),('psi','psi_P{}_deg')]:
    for q in ['50','95']:
     if m in cols:d[col.format(q)]=qrow[f'{m}_{scope}_P{q}']
   if exp=='H2':d.update(sigma_bearing_deg=.1,navigation_sigma_m=0,feature_sigma_and_unit=cfg+'; sigma q m/s; 6x200s windows',environmental_contract='kinematics; no acoustic propagation',channel_gain_contract='NOT_MODELED',estimator_family='bounded nonlinear least squares; 64 sign expansions for U/UB',comments='B baseline; S known sign; U magnitude; UB unknown zero; U-mismatch ignores injected bias. Independent base total256, NOT2816.',n_independent_total=256)
   elif exp=='AUX-D1':d.update(sigma_bearing_deg=.05,navigation_sigma_m=25,bearing_bias_bound_deg='common +/-0.05; half-difference +/-0.025',aux_baseline_km=0 if cfg=='U0' else 5,channel_gain_contract='NOT_MODELED',environmental_contract='signed bearing CV; exact synchronization and association assumed',estimator_family='4D bounded soft_l1;128 starts; fixed bias not estimated',comments='1536 paired base records; output acceptance not accuracy; U0 often hits physical bounds',n_independent_total=1536)
   else:
    quality,geom=cfg.split('_');d.update(sigma_bearing_deg=.05 if quality=='Q1' else .075,navigation_sigma_m=25 if quality=='Q1' else 37.5,bearing_bias_bound_deg='common +/-0.05; half-difference +/-0.025' if quality=='Q1' else 'common +/-0.075; half-difference +/-0.0375',aux_baseline_km={'B0':5,'B1':7,'B2':math.sqrt(104),'B3':math.sqrt(234)}[geom],environmental_contract='signed bearing CV; simultaneous correct association; equal assumed quality both nodes',estimator_family='bounded 4-state+2 bias; Uniform hard support; navigation-equivalent likelihood',comments='6000 base innovations shared across8configs; B0 dual node; B3 all outputs touch bias support; no receiver demonstration',n_independent_total=6000)
  configs.append((exp,cfg,scene,rr,cols,accepted))
 for cfg,qq in byconfig.items():
  summary={'experiment':exp,'configuration':cfg,'scenes':len(qq),'n':sum(x['n'] for x in qq),'accepted':sum(x['n_accepted'] for x in qq)}
  for m in cols:
   for scope in ['unconditional','conditional']:
    v=[float(x[f'{m}_{scope}_P50']) for x in qq if x[f'{m}_{scope}_P50']!=''];w=[float(x[f'{m}_{scope}_P95']) for x in qq if x[f'{m}_{scope}_P95']!='']
    summary[f'{m}_{scope}_mean_scene_P50']=fmt(sum(v)/len(v) if v else None);summary[f'{m}_{scope}_worst_scene_P95']=fmt(max(w) if w else None)
  PAIRED.append(summary)
write('CZ_D1_RECOMPUTED_SCENE_QUANTILES.csv',QUANT);write('CZ_D1_RECOMPUTED_CONFIG_SUMMARY.csv',PAIRED)
# Compare full saved nearest-rank summaries.
for r in rows(H2+'METRICS_BY_GEOMETRY.csv'):
 cfg=f"{r['method']}_{float(r['sigma']):g}"; q=next(x for x in QUANT if x['experiment']=='H2' and x['configuration']==cfg and x['scene']==r['geometry'])
 for m,s,mul in [('r','range_relative',100),('theta','bearing_deg',1),('v','speed_relative',100),('psi','heading_deg',1)]:
  for sc in ['unconditional','conditional']:
   for qn,sn in [('P50','median'),('P90','P90'),('P95','P95')]:
    if r[f'{s}_{sc}_{sn}'] and q[f'{m}_{sc}_{qn}']!='':eq(f'H2:{cfg}:{r["geometry"]}:{m}:{sc}:{qn}',float(q[f'{m}_{sc}_{qn}']),float(r[f'{s}_{sc}_{sn}'])*mul)
for prefix,fn,ex in [(AUX,'PAIRED_RESULTS_BY_CASE.csv','AUX-D1'),(T3,'CASE_STATISTICS.csv','T3-W1')]:
 for r in rows(prefix+fn):
  if ex=='AUX-D1' and r['variant']!='PRIMARY':continue
  if ex=='T3-W1' and r['side_group']!='BOTH':continue
  cfg=r['resource'] if ex=='AUX-D1' else r['config'];q=next(x for x in QUANT if x['experiment']==ex and x['configuration']==cfg and x['scene']==r['case_id'])
  for m,s,mul in [('r','r' if ex=='AUX-D1' else 'range_rel',100),('theta','theta' if ex=='AUX-D1' else 'bearing_deg',1),('v','v' if ex=='AUX-D1' else 'speed_rel',100),('psi','psi' if ex=='AUX-D1' else 'heading_deg',1)]:
   for sc in ['unconditional','conditional']:
    ss='accepted_only' if sc=='conditional' and ex=='T3-W1' else sc
    for qn,sn in [('P50','median'),('P90','P90'),('P95','P95'),('P99','P99')]:eq(f'{ex}:{cfg}:{r["case_id"]}:{m}:{sc}:{qn}',float(q[f'{m}_{sc}_{qn}']),float(r[f'{s}_{ss}_{sn}'])*mul)
# Paired changes and P95 differences require same experiment/quality/scene.
paired=[]
for exp,data,key,idkey,acc,cols in [('AUX-D1',a,'resource','record_index',lambda r:r['accepted']=='True',{'r':'r_error','theta':'theta_error','v':'v_error','psi':'psi_error'}),('T3-W1',t,'config','draw_index',lambda r:r['accepted']=='True',{'r':'range_rel','theta':'bearing_deg','v':'speed_rel','psi':'heading_deg'})]:
 cfgs=sorted({r[key] for r in data})
 for cfg in cfgs:
  ref='U0' if exp=='AUX-D1' else cfg.split('_')[0]+'_B0'
  if cfg==ref:continue
  aa={r[idkey]:r for r in data if r[key]==ref};bb={r[idkey]:r for r in data if r[key]==cfg};chk(f'{exp}:{cfg}:paired_ids',set(aa)==set(bb))
  for m,c in cols.items():
   better=worse=tie=0
   for k in aa:
    x=float(aa[k][c]) if acc(aa[k]) else math.inf;y=float(bb[k][c]) if acc(bb[k]) else math.inf
    better+=y<x;worse+=y>x;tie+=y==x
   paired.append(dict(experiment=exp,reference=ref,candidate=cfg,parameter=m,n=len(aa),improved=better,worsened=worse,tied=tie,scope='saved same-record errors; no new solve'))
write('CZ_D1_PAIRED_CHANGES.csv',paired)
# A1 coarse offgrid diagnostic; preserve original linear quantile, not continuous recovered precision.
a1='results/R4_A1_OFFGRID_BEARING_BOUNDARY/CASE_LEVEL_RESULTS.csv';r=rows(a1); rr=[x for x in r if x['experiment_scope']=='NOMINAL_PILOT_ONLY'];chk('A1_pilot27',len(rr)==27)
d=base('R4-A1','A2',a1,scenario='9 offgrid states x3 frozen seeds',n_total=len(rr),n_independent=27,N_elements=8,element_spacing_m=2,aux_count=0,vertical_channels=0,turn_deg=15,turn_start_s=600,observation_s=1200,epoch_count=121,sigma_bearing_deg=.1,navigation_sigma_m=0,frequencies_Hz='201;235;283',feature_sigma_and_unit='perfect model relative TL; no extraction',extraction_source='MODEL_TL_AND_SYNTHETIC_BEARING',environmental_contract='E-STD matched rigid bottom',offgrid_yes_no='YES',metrics_scope='pooled27 linear percentile; selected grid point error only',gate='BLOCKED_COARSE_RC2_SELECTION_AND_OFFGRID_MODAL_MISMATCH',evidence_category='GRID_CANDIDATE_CONTRACTION',precision_status='GRID_POINT_ERROR_NOT_CONTINUOUS_CERTIFICATE',claim_scope='MECHANISM',comments='0/27 top1 bracketing-cell success; do not treat as all-algorithm impossibility')
for m,c,unit,mul in [('r','top1_rel_r','percent',100),('theta','top1_abs_theta_deg','deg',1),('v','top1_rel_v','percent',100),('psi','top1_abs_psi_deg','deg',1)]:
 for q in [50,95]:d[f'{m}_P{q}_{unit}']=float(np.percentile([float(x[c]) for x in rr],q))*mul
chk('A1_top1_success0',sum(x['top1_in_bracketing_cell']=='True' for x in rr)==0)
# R3 turn full catalogue preservation.
turn='results/R3_RC23_CLOSEDLOOP/R3_RC23_SINGLE_TURN_RANGE_ANCHOR/'
tr=rows(turn+'RC2_RC3_TURN_SET_METRICS.csv');br=rows(turn+'RC2_TURN_ONLY_SET_METRICS.csv');write('CZ_D1_R3_TURN_CATALOGUE.csv',[dict(observable='RC2+RC3',**x) for x in tr]+[dict(observable='RC2_ONLY',**x) for x in br])
for x in tr:
 base('R3-1C','A1',turn+'RC2_RC3_TURN_SET_METRICS.csv',scenario=f"{x['config']};z={x['z_true_m']};tau={x['tau_db']}",r_true_km=50,z_true_m=x['z_true_m'],v_true_mps=2,psi_true_deg=5,turn_deg=x['delta_deg'],turn_start_s=600,observation_s=1200,epoch_count=121,N_elements=8,element_spacing_m=2,aux_count=0,sigma_bearing_deg=.1,extraction_source='IDEAL_TRACKABLE_TL',environmental_contract='E-STD matched rigid bottom',metrics_scope='finite catalogue; one anchor motion; no MC precision',evidence_category='GRID_CANDIDATE_CONTRACTION',claim_scope='MECHANISM',correct_support_retention=x['true_retained'],candidate_r_width_km=x['r_width_km'],candidate_count=x['n_RC2RC3'],n_independent=1,gate='SMALL_TURN_RANGE_RIDGE_PERSISTS',comments='RC2+RC3; separate RC2-only CSV; source tau and frequency contract retained; not P95')
# H3 resource matched ratios, do not convert local scales to precision.
hi=rows(H3+'VERTICAL_INFORMATION_BY_SCENE.csv');hs=rows(H3+'RESOURCE_MATCHED_CONTROL.csv')
for x in hs:
 matches=[a for a in hi if a['scene_id']==x['scene_id'] and a['noise_model']==x['noise_model'] and a['mesh']=='160001' and a['calibration']=='C1b' and float(a['sigma'])==.01]
 b1=next(a for a in matches if a['resource']=='B1');b2=next(a for a in matches if a['resource']=='B2');ratio=float(b1['depth_information_horizontal_profiled'])/float(b2['depth_information_horizontal_profiled']);eq('H3_ratio:'+x['scene_id']+x['noise_model'],ratio,float(x['B1_to_B2']))
 base('H3-G0','A3',H3+'VERTICAL_INFORMATION_BY_SCENE.csv',scenario=x['scene_id']+';'+x['noise_model'],configuration='B1 vs B2; C1b;160001;sigma1%',N_elements=12,element_spacing_m=2,vertical_channels=4,aux_count=0,z_true_m=200,observation_s=1200,epoch_count=3,turn_deg=15,turn_start_s=600,extraction_source='IDEAL_COMPLEX_FIELD_NOT_RECEIVER',feature_sigma_and_unit='1%/5% proper complex design; main1%; absolute floor separate',channel_gain_contract='C1b element-frequency gains fixed across3snapshots; source free eachfrequency/time',environmental_contract='nominal matched E-STD; SSP/pose untested',evidence_category='FISHER_LOCAL_INFORMATION',claim_scope='MECHANISM',gate='CONDITIONAL_LOCAL_DEPTH_INFORMATION',info_gain_B1_B2=ratio,n_independent=3,metrics_scope='3 geometries x mirror andnoise; NOT independent MC',comments='Source-depth finite-difference absoluteaccuracyuncertified; real horizontal support and extraction missing')
# CAP06 archive identity and epoch statistics, no extraction re-run.
raw=rows(CAP+'RAW_ARCHIVE_INDEX.csv'); rawchecks=[]
for i,x in enumerate(raw):
 p=Path(x['path']);exists=p.is_file();ok=exists and p.stat().st_size==int(x['bytes']) and sha(p)==x['sha256']
 rawchecks.append({'path':str(p),'expected_sha256':x['sha256'],'expected_bytes':x['bytes'],'PASS':ok,'status':'VERIFIED' if ok else 'MISSING_OR_CHANGED'})
write('CZ_D1_CAP06_RAW_IDENTITY.csv',rawchecks);chk('CAP06_all_raw_SHA',all(x['PASS'] for x in rawchecks),f'{len(rawchecks)}files')
cap=rows(CAP+'P2_RECEIVER_BEARING_CASES.csv'); csummary=[]
for mode in ['CAL','UNCAL']:
 for snr in ['0','10','20']:
  rr=[x for x in cap if x['split']=='holdout' and x['stress']=='MAIN_TRACK' and x['mode']==mode and x['snr']==snr and x['present']=='True']
  vals=[float(x['axis_error_deg']) for x in rr];rec=dict(mode=mode,snr=snr,n_present=len(rr),extracted=sum(x['target_component_extracted']=='True' for x in rr),covered=sum(x['target_u_covered']=='True' for x in rr),axis_P50_deg=fmt(nr(vals,.5)),axis_P95_deg=fmt(nr(vals,.95)));csummary.append(rec)
  base('CAP06-O6','A6',CAP+'P2_RECEIVER_BEARING_CASES.csv',scenario='holdout MAIN+B0+B3 paired',configuration=mode+'_'+snr,N_elements=8,element_spacing_m=2,aux_count='0 or1; separate node extraction, not fusion',observation_s=1200,epoch_count=31,turn_deg=15,turn_start_s=600,n_total=len(rr),n_independent=256,n_rejected=len(rr)-rec['extracted'],extraction_source='ACTUAL_SYNTHETIC_8CHANNEL_SAMPLES',feature_sigma_and_unit=snr+' dB reference;1s snapshots;not iid bearing sigma',channel_gain_contract='CAL4knownangles40dB orUNCAL; +/-1dB,+/-20deg gains',environmental_contract='single-plane-wave snapshot with registered interference/defects; NOT CZ fullwave',evidence_category='RECEIVER_EXTRACTION',claim_scope='RECEIVER_CHAIN_CONDITIONAL',gate='P2_NOT_ADMITTED',metrics_scope='axis-angle receiver error, NOT target state theta; epochs correlated/paired',receiver_axis_P50_deg=rec['axis_P50_deg'],receiver_axis_P95_deg=rec['axis_P95_deg'],correct_support_retention=rec['covered']/len(rr),cold_audit_status='COLD_AUDIT_ONLY; archiveSHAand savedmetricrecalc in CZ-D1',lead_audit_status='PENDING_LEAD_AUDIT',comments='P3NOT_EVALUATED; mainCAL20coverage<95%; 320base trajectories total; fresh256same-modelholdout; uncalibrated zero extraction not allmethods impossibility')
write('CZ_D1_CAP06_RECALCULATED.csv',csummary)
# All source reads above carry precise hashes. More curated sources added in report builder.
write('CZ_D1_EVIDENCE_MASTER.csv',MASTER);write('CZ_D1_RECALC_CHECKS.csv',CHECKS)
write('CZ_D1_SOURCE_PROVENANCE.json',{'parent':PARENT,'design_sha':'6ff0d07fca848d8022fad2facd5a107a920f8136','files':list(SOURCES.values()),'cap06_raw_index':'CZ_D1_CAP06_RAW_IDENTITY.csv','history_policy':'Last-change SHA identifies artifact byte revision, not necessarily entire scientific execution; immutable snapshot cited for reproduction','scope':'Read-only recomputation from saved outputs; no science models or extractors rerun'})
write('CZ_D1_RECALC_STATUS.json',{'checks':len(CHECKS),'failed':sum(not x['PASS'] for x in CHECKS),'source_files':len(SOURCES),'master_rows':len(MASTER),'inventory_files':len(files),'stage_entries':len(inv),'new_solver':0,'new_MC':0,'new_optimizer':0})
print(json.dumps({'checks':len(CHECKS),'FAIL':sum(not x['PASS'] for x in CHECKS),'master':len(MASTER),'inventory':len(files),'stages':len(inv)}))
