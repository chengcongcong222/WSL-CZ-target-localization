"""Independent CZ-D1 saved-evidence/report audit; no historical science imports."""
import csv,json,math,hashlib,subprocess,collections
from pathlib import Path
import numpy as np
O=Path(__file__).resolve().parent;R=O.parents[1];checks=[];source_extra=[]
def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def h(p):
 s=hashlib.sha256()
 with Path(p).open('rb') as f:
  for a in iter(lambda:f.read(2**20),b''):s.update(a)
 return s.hexdigest()
def test(k,b,detail=''):checks.append({'check':k,'PASS':bool(b),'detail':str(detail)})
def equal(k,a,b):test(k,a==b or (np.isfinite(a) and np.isfinite(b) and abs(a-b)<1e-9*max(1,abs(b))),f'{a};{b}')
def write(n,rr):
 with (O/n).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,list(dict.fromkeys(k for r in rr for k in r)));w.writeheader();w.writerows(rr)
def source(p):
 p=R/p;source_extra.append({'path':p.relative_to(R).as_posix(),'bytes':p.stat().st_size,'sha256':h(p),'version_sha':'ca86d324afc1689af6c96cd35145e03625086bae','scope':'input contract and saved archive identity'});return p
P='ca86d324afc1689af6c96cd35145e03625086bae'
# Saved raw observation hashes, rather than archive summaries alone.
ids=[]
for name in ['PUBLIC_BLOCK_IDENTITIES.csv','PRIVATE_INNOVATION_IDENTITIES.csv']:
 for x in rows(source('results/R4_T3_W1_WIDE_BASELINE_DEMO/'+name)):
  p=Path(x['path']);ok=p.is_file() and p.stat().st_size==int(x['bytes']) and h(p)==x['sha256'];test('T3identity:'+str(p),ok);ids.append({'path':str(p),'sha256':x['sha256'],'PASS':ok})
a=json.loads(source('results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/PUBLIC_OBSERVATION_MANIFEST.json').read_text('utf-8-sig'))
for x in a['arrays']:
 p=R/x['path'];ok=p.is_file() and p.stat().st_size==int(x['bytes']) and h(p)==x['sha256'];test('AUXidentity:'+str(p),ok);ids.append({'path':str(p),'sha256':x['sha256'],'PASS':ok})
# Source contract identity and complete selected physical/metric contract sections retained.
contract={}
for label,p in [('H2','R4_SINGLE_HLA_RADIAL_INCREMENT_PRECISION_PILOT/DESIGN_FREEZE.json'),('AUX-D1','R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/AUX_D1_TASK_SPEC_20261012_v1.json'),('T3-W1','R4_T3_W1_WIDE_BASELINE_DEMO/T3_W1_DESIGN_FREEZE.json'),('H3','R4_H3_G0_FINITE_VERTICAL_RESPONSE/DESIGN_FREEZE.json'),('CAP06','ENV_CAP_CAP06/CAP06_FROZEN_SPEC.json'),('R3-turn','R3_RC23_CLOSEDLOOP/R3_RC23_SINGLE_TURN_RANGE_ANCHOR/CLOSEDLOOP_1C_CONFIG.json')]:
 j=json.loads(source('results/'+p).read_text('utf-8-sig'));contract[label]={'source':'results/'+p,'sections':{k:v for k,v in j.items() if k not in ['sha256','inputs','frozen_code','cache_metadata']}};test('contract_parse:'+label,True)
(O/'CZ_D1_INPUT_CONTRACTS.json').write_text(json.dumps(contract,ensure_ascii=False,indent=2),encoding='utf-8')
write('CZ_D1_ADDITIONAL_ARCHIVE_IDENTITIES.csv',ids)
# Independent numpy inverted-CDF implementation checks every recomputed scene quantile.
q=rows(O/'CZ_D1_RECOMPUTED_SCENE_QUANTILES.csv');specs={
 'H2':('R4_SINGLE_HLA_RADIAL_INCREMENT_PRECISION_PILOT/CONFIGURATION_RESULTS.csv',lambda r:f"{r['method']}_{float(r['sigma']):g}",'geometry',lambda r:r['status']=='ACCEPTED',{'r':('relative_range_error',100),'theta':('absolute_bearing_error_deg',1),'v':('relative_speed_error',100),'psi':('absolute_heading_error_deg',1)}),
 'AUX-D1':('R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/ESTIMATES.csv',lambda r:r['resource'],'case_id',lambda r:r['accepted']=='True',{'r':('r_error',100),'theta':('theta_error',1),'v':('v_error',100),'psi':('psi_error',1)}),
 'T3-W1':('R4_T3_W1_WIDE_BASELINE_DEMO/RUN_ERRORS.csv',lambda r:r['config'],'case_id',lambda r:r['accepted']=='True',{'r':('range_rel',100),'theta':('bearing_deg',1),'v':('speed_rel',100),'v_abs':('speed_abs_mps',1),'psi':('heading_deg',1)})}
for exp,(path,key,scene,accepted,cols) in specs.items():
 groups=collections.defaultdict(list)
 for r in rows(R/'results'/path):groups[key(r),r[scene]].append(r)
 for r in (x for x in q if x['experiment']==exp):
  rr=groups[r['configuration'],r['scene']];test('denominator:'+exp+r['configuration']+r['scene'],len(rr)==int(r['n']))
  for metric,(col,mul) in cols.items():
   for scope in ['conditional','unconditional']:
    v=np.array([float(x[col])*mul if accepted(x) else np.inf for x in rr if scope=='unconditional' or accepted(x)])
    for tag,pct in [('P50',.5),('P90',.9),('P95',.95),('P99',.99)]:
     if len(v):equal(':'.join([exp,r['configuration'],r['scene'],metric,scope,tag]),float(np.quantile(v,pct,method='inverted_cdf')),float(r[f'{metric}_{scope}_{tag}']))
master=rows(O/'CZ_D1_EVIDENCE_MASTER.csv');test('master_has_7_architecture_categories',all(any(a in r['architecture_id'].split(';') for r in master) for a in ['A0','A1','A2','A3','A4','A5','A6']))
for i,r in enumerate(master):
 test(f'row_source:{i}',(R/r['result_file']).is_file(),r['result_file'])
 test(f'row_joint_not_fabricated:{i}',r['joint_status']=='NOT_DEFINED')
 for k in ['r_P50_percent','r_P95_percent','theta_P50_deg','theta_P95_deg','z_P95_m','v_P50_percent','v_P95_percent','v_P95_mps','psi_P95_deg']:
  if r.get(k,'')!='':test(f'valid_number:{i}:{k}',not math.isnan(float(r[k])) and float(r[k])>=0)
 if r['experiment_id'].startswith('CAP06'):test(f'no_receiver_axis_as_state_precision:{i}',r['theta_P95_deg']=='')
 if r['experiment_id'] in ['H3-G0','T3-W1','H2','AUX-D1']:test(f'no_depth_precision:{i}',r['z_P95_m']=='')
 if r['experiment_id']=='T3-W1':test(f'T3dual:{i}',r['aux_count']=='1')
# CAP06 all-file verification record + independent evenly-spaced raw-file recheck, no signal extraction.
cap=rows(O/'CZ_D1_CAP06_RAW_IDENTITY.csv');test('CAP06_20545_identity_records',len(cap)==20545);test('CAP06_all_identity_pass',all(x['PASS']=='True' for x in cap))
for i in np.linspace(0,len(cap)-1,64,dtype=int):test('independent_CAP06_SHA:'+str(i),h(Path(cap[i]['path']))==cap[i]['expected_sha256'])
saved=rows(R/'results/ENV_CAP_CAP06/P2_EXTRACTION_SUMMARY.csv')
for x in rows(O/'CZ_D1_CAP06_RECALCULATED.csv'):
 y=next(r for r in saved if r['split']=='holdout' and r['stress']=='MAIN_TRACK' and r['snr']==x['snr'] and r['mode']==x['mode'])
 for a,b in [('n_present','present_epochs'),('extracted','extracted'),('covered','covered'),('axis_P50_deg','axis_error_P50_deg'),('axis_P95_deg','axis_error_P95_deg')]:equal('CAP06saved:'+x['mode']+x['snr']+a,float(x[a]),float(y[b]))
# Identity provenance checks include existing selected scientific artifacts.
pr=json.loads((O/'CZ_D1_SOURCE_PROVENANCE.json').read_text());ps={x['path']:x for x in pr['files']};ps.update({x['path']:x for x in source_extra});pr['files']=list(ps.values())
for x in pr['files']:test('sourceSHA:'+x['path'],h(R/x['path'])==x['sha256'])
(O/'CZ_D1_SOURCE_PROVENANCE.json').write_text(json.dumps(pr,ensure_ascii=False,indent=2),encoding='utf-8')
for f in ['CZ_D1_RECALC_CHECKS.csv','CZ_D1_SAVED_ESTIMATE_ERROR_CHECKS.csv']:test(f+':allpass',all(x['PASS']=='True' for x in rows(O/f)))
# Historical directories must remain byte-clean relative to parent.
diff=subprocess.check_output(['git','-c','safe.directory='+str(R),'diff','--name-only',P],cwd=R).decode('utf-8').splitlines();test('only_new_CZ_D1_changes',all(x.startswith('results/CZ_D1/') for x in diff),diff)
required=['CZ_D1_EVIDENCE_MASTER.csv','CZ_D1_ARCHITECTURE_METRICS.md','CZ_D1_MANEUVERS_AND_AUGMENTATIONS.md','CZ_D1_FIVE_STATE_CLAIM_MATRIX.md','CZ_D1_NEGATIVE_GATES_AND_GAPS.md','CZ_D2_SELECTION_CANDIDATES.md','CZ_D1_SOURCE_PROVENANCE.json','GPT_SYNC_CZ_D1.md']
for f in required:test('required:'+f,(O/f).is_file())
write('CZ_D1_INDEPENDENT_AUDIT_CHECKS.csv',checks)
fail=sum(not x['PASS'] for x in checks);result={'PASS':fail==0,'checks':len(checks),'FAIL':fail,'scope':'Independent code by same executing agent, not research-lead acceptance; numpy inverted-CDF from raw saved case errors, report invariants, source identities, CAP0664file spot check in addition to full main SHA check; no history science rerun','prior_source_summary_checks':5974,'saved_estimate_error_checks':215268,'additional_public_private_archive_files':len(ids),'CAP06_allraw_files':len(cap),'new_solver_calls':0,'new_MC':0,'new_optimizer':0,'new_signal_generation':0}
(O/'VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(result))
