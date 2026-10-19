"""Additional historical evidence only. Deterministic report metadata and algebra on saved estimates."""
from pathlib import Path
import csv,json,hashlib,subprocess,math,collections
import numpy as np
R=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parent
P='ca86d324afc1689af6c96cd35145e03625086bae'
def read(p):return list(csv.DictReader((R/p).open(encoding='utf-8-sig',newline='')))
def write(p,rr):
 ks=list(dict.fromkeys(k for r in rr for k in r))
 with (O/p).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,ks);w.writeheader();w.writerows(rr)
def git(*args):return subprocess.check_output(['git','-c','safe.directory='+str(R),'-c','core.quotepath=false',*args],cwd=R).decode('utf-8').strip()
prov=json.loads((O/'CZ_D1_SOURCE_PROVENANCE.json').read_text());sources={r['path']:r for r in prov['files']}
def cite(p):
 b=(R/p).read_bytes();sources[p]={'path':p,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'version_sha':P,'scope':'historical artifact read; no propagation/estimation rerun'};return p
master=read('results/CZ_D1/CZ_D1_EVIDENCE_MASTER.csv');qs=read('results/CZ_D1/CZ_D1_RECOMPUTED_SCENE_QUANTILES.csv');checks=[]
truthp=cite('results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_TRUTH_PANEL.csv');tp={r['case_id']:r for r in read(truthp)}
hp=cite('results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/OFFGRID_TRUTH_EVALUATION_ONLY.csv');ht={r['panel']:r for r in read(hp)}
for row in master:
 if row['experiment_id'] in ['H2','AUX-D1','T3-W1']:
  h=ht[row['scenario']] if row['experiment_id']=='H2' else tp[row['scenario']]
  row.update(r_true_km=h.get('r_km',h.get('r0_km')),theta_true_deg=h.get('theta_deg',h.get('theta0_deg')),v_true_mps=h['v_mps'],psi_true_deg=h['psi_deg'],z_true_m=h.get('z_label','NOT_MODELED'),lead_audit_status='PENDING_LEAD_AUDIT',user_review_status='SELECTED_NUMBERS_REVIEWED_IN_USER_TASK; no blanket upgrade of stored pending status')
  if row['experiment_id'] in ['H2','AUX-D1']:
   for q in ['50','95']:row[f'v_P{q}_mps']=float(row[f'v_P{q}_percent'])*float(h['v_mps'])/100
  if row['experiment_id']=='H2':row['feature_sigma_and_unit']=row['configuration']+'; sigma(q) in m/s; 6x200s; signed S vs folded U/UB; baseline B none'
 # Every numerical field has explicit scope, all blank state entries carry per-state status.
 for state,unit in [('r','percent'),('theta','deg'),('z','m'),('v','percent'),('psi','deg')]:row[state+'_metric_status']='RECALCULATED_SCOPED' if row.get(state+'_P95_'+unit,'')!='' else 'NOT_EVALUATED'
# Verify per-record errors from saved estimated states, without rerunning optimizer.
def circ(a,b):return abs((a-b+180)%360-180)
def check(exp,a,b):
 good=(a==b) or (math.isfinite(a) and math.isfinite(b) and abs(a-b)<1e-8*max(1,abs(b)))
 checks.append({'check':exp,'PASS':good,'absolute_difference':abs(a-b) if math.isfinite(a) and math.isfinite(b) else ''})
for exp,estpath,errpath,keys in [('H2','results/R4_SINGLE_HLA_RADIAL_INCREMENT_PRECISION_PILOT/ESTIMATES.csv','results/R4_SINGLE_HLA_RADIAL_INCREMENT_PRECISION_PILOT/CONFIGURATION_RESULTS.csv',['config_id']),('T3-W1','results/R4_T3_W1_WIDE_BASELINE_DEMO/ESTIMATOR_OUTPUTS.csv','results/R4_T3_W1_WIDE_BASELINE_DEMO/RUN_ERRORS.csv',['config','case_id','replicate'])]:
 ee={tuple(r[k] for k in keys):r for r in read(cite(errpath))}
 for e in read(cite(estpath)):
  er=ee[tuple(e[k] for k in keys)];h=ht[e['geometry']] if exp=='H2' else tp[e['case_id']]
  accepted=e['status']=='ACCEPTED' if exp=='H2' else e['accepted']=='True'
  if not accepted:continue
  if exp=='H2':cols=['r_km','theta_deg','v_mps','psi_deg'];ecs=['relative_range_error','absolute_bearing_error_deg','relative_speed_error','absolute_heading_error_deg']
  else:cols=['r_hat_km','theta_hat_deg','v_hat_mps','psi_hat_deg'];ecs=['range_rel','bearing_deg','speed_rel','heading_deg']
  r0=float(h.get('r_km',h.get('r0_km')));th=float(h.get('theta_deg',h.get('theta0_deg')));v=float(h['v_mps']);psi=float(h['psi_deg'])
  pred=[abs(float(e[cols[0]])-r0)/r0,circ(float(e[cols[1]]),th),abs(float(e[cols[2]])-v)/v,circ(float(e[cols[3]]),psi)]
  for m,a,c in zip(['r','theta','v','psi'],pred,ecs):check(exp+':'+':'.join(e[k] for k in keys)+':'+m,a,float(er[c]))
for e in read(cite('results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/ESTIMATES.csv')):
 h=tp[e['case_id']]
 for k,truth,err in [('r','r0_km','r_error'),('theta','theta0_deg','theta_error'),('v','v_mps','v_error'),('psi','psi_deg','psi_error')]:
  val=abs(float(e[k+'_hat'])-float(h[truth]))/float(h[truth]) if k in ['r','v'] else circ(float(e[k+'_hat']),float(h[truth]));check('AUX:'+e['configuration']+':'+k,val,float(e[err]))
write('CZ_D1_SAVED_ESTIMATE_ERROR_CHECKS.csv',checks)
# Older five-state heuristics are retained as negative diagnostics, not ignored because newer gates failed.
p4=cite('results/P4_5_parameter_accuracy/parameter_accuracy_trials.csv');p4g=collections.defaultdict(list)
for r in read(p4):p4g[r['scenario'],r['layer']].append(r)
p4out=[]
for (scene,layer),rr in p4g.items():
 d={'scenario':scene,'layer':layer,'n_saved_rows':len(rr),'n_distinct_truth_trials':len({r['trial'] for r in rr}),'n_distinct_noise_ids':len({r['noise'] for r in rr}),'status':'HISTORICAL_HEURISTIC_POINT_OUTPUT; coverage invalid/weak; no revised admission','rejection_status':'NO_UNIFIED_ACCEPTANCE_RULE; not unconditional certified output','scope':'CZ-D1 retrospective linear quantiles from saved points; fractional r/v normalized by each true value, not nominal state'}
 for k,field,scale in [('r','abs_err_r_km','truth_r_km'),('theta','abs_err_theta_deg',None),('z','abs_err_z',None),('v','abs_err_v','truth_v'),('psi','abs_err_psi_deg',None)]:
  vals=[float(r[field])/(float(r[scale]) if scale else 1)*(100 if scale else 1) for r in rr if r[field] and math.isfinite(float(r[field]))]
  d[k+'_n_finite']=len(vals)
  for q in [50,95]:d[k+'_P'+str(q)]=float(np.percentile(vals,q)) if vals else ''
  d[k+'_coverage']=sum(r['cov_'+('r' if k=='r' else k)]=='True' for r in rr)/len(rr)
 p4out.append(d)
write('CZ_D1_EARLY_P4_DIAGNOSTICS.csv',p4out)
# Curated route states. No unverified precision placed in numerical output columns.
curated=[
('R3_FINAL_CLOSEOUT','A2','R3_FINAL_WORKING_ENVELOPE.md','GRID_CANDIDATE_CONTRACTION','MECHANISM','FROZEN_IDEAL_DISCRETE_TWO_HEADING_CANDIDATES','Matched E-STD; three lines; one motion, three source labels; 15deg turn','Not continuous precision; depth profiled'),
('P2_RC2_kinematic_boundary','A0;A1','P2_RC2_GPT_SYNC.md','GRID_CANDIDATE_CONTRACTION','MECHANISM','RANGE_MOTION_RIDGE','Signed synthetic bearings; time and turn controls','No general HLA physical impossibility'),
('P4_5_parameter_accuracy','A2','P4_5_GPT_SYNC.md','GRID_CANDIDATE_CONTRACTION','MECHANISM','INVALID_COVERAGE_OR_WEAK','Old heuristic five-state outputs; supplemental P50/P95 retained separately','Old z RMSE/P95 are not certified depth precision'),
('P4_G4_performance_boundary','A1;A2','P4_G4_REVISION.md','GRID_CANDIDATE_CONTRACTION','MECHANISM','G4_PARTIAL_BOUNDARY_VALIDATED','Revised nested candidate filtering and retain-first display','Use revised report; original nonnested positive widths superseded'),
('R4_A1_FIX_CONTINUOUS_SEARCH','A2','GPT_SYNC.md','NOT_EVALUATED','MECHANISM','HISTORICAL_STOP','See frozen report','Search failure not physical nonidentifiability'),
('R4_A1_FIX2_ACOUSTIC_COVERAGE','A2','GPT_SYNC.md','GRID_CANDIDATE_CONTRACTION','MECHANISM','BLOCKED_UNCLOSED_ACOUSTIC_COVERAGE','Correct narrow basins exist; global capture not established','Cumulative retained successes not raw budget convergence'),
('R4_A1_SEARCH_TRACTABILITY_AUDIT','A2','GPT_SYNC.md','NOT_EVALUATED','MECHANISM','SEARCH_TRACTABILITY_NOT_ESTABLISHED','SHGO/DIRECT frozen finite budget','No extrapolation to all algorithms'),
('R4_RC2_OBSERVATION_SUPPORT','A0','GPT_SYNC.md','GRID_CANDIDATE_CONTRACTION','MECHANISM','OUTER_SUPPORT_BROAD_OR_INPUT_INCOMPLETE','Unresolved boxes retained; auxiliary observation archive missing at that stage','AUX-D1 later generated new paired observations; not historical recovery'),
('R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT','A2','GPT_SYNC.md','GRID_CANDIDATE_CONTRACTION','IDEAL_CONDITIONAL','GRID_GAIN_OFFGRID_SEARCH_LIMITED','Three trackable lines;0.25dB main; source common time variation removed','Offgrid0/32 accepted; truecost passes32/32; source-specific variation degrades'),
('R4_SINGLE_HLA_RADIAL_RECEIVER_EXTRACTION_PILOT','A2','GPT_SYNC.md','RECEIVER_EXTRACTION','MECHANISM','FEATURE_MAPPING_NOT_ESTABLISHED','Stable records only;1500m/s groupdelay fallback; source-reference control separate','Drift provider missing; REF cannot grant unknown-source precision'),
('R4_SINGLE_HLA_CARRIER_REFERENCE_CONTRAST_PILOT','A2','GPT_SYNC.md','RECEIVER_EXTRACTION','MECHANISM','RELATIVE_MOTION_MAPPING_NOT_SUPPORTED_IN_SAVED_MODEL','Saved received frequency changes; no new signal','Turning response not target-motion recovery'),
('R4_SINGLE_HLA_STATIC_SPATIAL_RANGE_PILOT','A2','GPT_SYNC.md','RECEIVER_EXTRACTION','RECEIVER_CHAIN_CONDITIONAL','RANGE_CONTRACTION_NOT_SUPPORTED;NUMERICAL_OR_GRID_LIMITED','8x2s static data;201/235/283Hz; C0/C1;20/5dB','Point peak accuracy does not mean total support shrinkage'),
('R4_H3_G1_EXTRACTION_SUPPORT_REVIEW','A3','GPT_SYNC.md','NOT_EVALUATED','MECHANISM','DEPTH_SUPPORT_AND_PROPAGATION_UNRESOLVED','Local C1b only; full horizontal support missing','No depth accuracy'),
('R4_H3_G2E_CONDITIONAL_EXTRACTION','A3','GPT_SYNC.md','RECEIVER_EXTRACTION','MECHANISM','COMPRESSED_RESPONSE_EXTRACTION_UNRELIABLE','Registered compressed response','Failure specific to representation'),
('R4_H3_P0_DERIVATIVE_REVIEW','A3','GPT_SYNC.md','NOT_EVALUATED','MECHANISM','STENCIL_CONSISTENCY_SUPPORTED;ACCURACY_UNCERTIFIED','R1/R2 internal consistency','No true derivative absolute error certificate'),
('R4_H3_P1_UNIFIED_WINDOW_SSP','A3','GPT_SYNC.md','NOT_EVALUATED','MECHANISM','MODAL_WINDOW_NOT_CERTIFIED;SECANT_LOCALITY_NOT_ESTABLISHED','CLOW agrees; CHIGH and locality fail','SSP effect on depth information not evaluated'),
('R4_T2_TILT_VS_SHORT_VERTICAL_G0','A4','GPT_SYNC.md','NOT_EVALUATED','MECHANISM','FORWARD_NUMERIC_NOT_ADMITTED','9 arrays; both modal grids; raw vector failures','No receiver tests; no physical rejection of tilted array'),
('R4_T2_PREC3_FUNC2_CONTINUOUS_MODAL_REFERENCE','A4','GPT_SYNC_PREC3_FUNC2.md','NOT_EVALUATED','MECHANISM','CLOSED_NUMERIC_ROUTE','Frozen implementation fails analytic admission','Do not reopen; not a method-family impossibility'),
('R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC','A5','GPT_SYNC.md','CONTINUOUS_ESTIMATOR_ERROR','IDEAL_CONDITIONAL','HISTORICAL_CONDITIONAL_HORIZONTAL_BASELINE','Old paired truth/innovation/output archive; no original public bearing archive','Not pairable with old single-node observations'),
('R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY','A5','GPT_SYNC.md','FISHER_LOCAL_INFORMATION','MECHANISM','CURRENT_GEOMETRY_INFORMATION_LIMIT','Known true positions/no-bias controls separated','CRLB not universal empirical P95 bound'),
('R4_AUX_GATE1_CROSSTRACK_REQUIREMENT','A5','GPT_SYNC.md','CONTINUOUS_ESTIMATOR_ERROR','IDEAL_CONDITIONAL','STATIC_BASELINE_BEARING_REQUIREMENTS','Static range geometry;5/7km; precision tradeoff','Not dynamic speed or joint5D'),
('R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET','A5','GPT_SYNC.md','CONTINUOUS_ESTIMATOR_ERROR','IDEAL_CONDITIONAL','NONIDEAL_STATIC_REQUIREMENTS','Individual error axes','Not simultaneous corner capability'),
('R4_AUX_GATE2B_JOINT_ERROR_BUDGET','A5','GPT_SYNC.md','CONTINUOUS_ESTIMATOR_ERROR','IDEAL_CONDITIONAL','JOINT_STATIC_REQUIREMENT','Static11ranges, signedcorners, doublemirror; shared30000 draws','T5/T10 named research criteria not client requirements; no interval guarantee'),
('R4_E1_G0_FREQUENCY_INFORMATION','A2;A5','GPT_SYNC.md','FISHER_LOCAL_INFORMATION','MECHANISM','FREQUENCY_INCREMENT_NOT_ESTABLISHED','Unknown source frequency/drift/node reference profiled','Do not rename same idealDoppler mechanism'),
('R4_E2_NINE_PAIR_INFORMATION_PILOT','A2','GPT_SYNC.md','NOT_EVALUATED','MECHANISM','IMPLEMENTATION_INVALID','Automatic promising classification explicitly superseded by formula audit','Do not export invalid positive Fisher evidence'),
('R4_E2_NINE_PAIR_FORMULA_REPAIR','A2','GPT_SYNC.md','NOT_EVALUATED','MECHANISM','FORMULA_REPAIR_FAILED','C1 evidence not established','No original E2 admission'),
('ENV_CAP_FULL_CAMPAIGN','CROSSCUTTING','GPT_SYNC_ENV_CAP_FINAL.md','NOT_EVALUATED','MECHANISM','NO_GO_CURRENT_CZ_METHOD_PERFORMANCE_DEMO','15environmentunits; E1coldfail; newCZcomplexprovider missing','Old environment result unchanged'),
('ENV_CAP_REC01','CROSSCUTTING','GPT_SYNC_ENV_CAP_REC01.md','NOT_EVALUATED','MECHANISM','CONDITIONAL_CZ_NUMERICS_RECOVERED_COHERENT_PROVIDER_NOT_READY','9SCS units; oneWOAfamily; ray/TL conditional','No localization precision'),
('ENV_CAP_CZROB02','CROSSCUTTING','GPT_SYNC_CZROB02.md','NOT_EVALUATED','MECHANISM','FULL_ANGLE_NUMERIC_GATE_FAILED','Interpolation diagnosis; failure preserved','Not physicalCZabsence'),
('ENV_CAP_NUM03','CROSSCUTTING','GPT_SYNC_NUM03.md','NOT_EVALUATED','MECHANISM','REFERENCE_AND_NOMINAL_FIRST_RETURN_GEOMETRY_SUPPORTED','9samefamilyunits; C-linear/PCHIP; raygeometry','No fullwave or extraction performance; numerical repair closed'),
('ENV_CAP_CZO04','A2','GPT_SYNC_CZO04.md','GRID_CANDIDATE_CONTRACTION','MECHANISM','LIMITED_CONDITIONAL_LEVEL_INFORMATION','200Hz single-depth I-SHD; sourceC0/C1/C2; depth/environment profiled discrete','General range anchor not established; free source null'),
('ENV_CAP_REQ05','CROSSCUTTING','GPT_SYNC_REQ05.md','NOT_EVALUATED','MECHANISM','RESOURCE_CONDITIONS_NO_NEW_PERFORMANCE','O1..O6 nuisance and resource matrix','No minimum equipment or accuracy certification'),
('ENV_CAP_CAP06','A3;A5;A6','GPT_SYNC_CAP06_FINAL.md','RECEIVER_EXTRACTION','RECEIVER_CHAIN_CONDITIONAL','NO_GO_FUSION_UNDER_FROZEN_RECEIVER_CONTRACT','P1toy only;P2 actual synthetic extraction;P3 not run','All requested combined configurations not tested positioning; no depth or5D')]
cur=[]
for stage,arch,name,cat,claim,gate,contract,limit in curated:
 p='results/'+stage+'/'+name
 if not (R/p).exists():
  ps=list((R/'results'/stage).glob('*GPT_SYNC*'));p=str(ps[0].relative_to(R)).replace('\\','/') if ps else p
 exists=(R/p).exists()
 if exists:
  cite(p);text=(R/p).read_text('utf-8-sig');last=git('log','-1','--format=%H',P,'--',p)
 else:text='SOURCE_MISSING';last=''
 row={k:'' for k in master[0]};row.update(experiment_id=stage,phase=stage,architecture_id=arch,version_sha=P,evidence_path='results/'+stage,result_file=p,artifact_last_change_sha=last,evidence_category=cat,claim_scope=claim,gate=gate,environmental_contract=contract,comments=limit,precision_status='NOT_RECALCULATED_SUMMARY_ONLY',joint_status='NOT_DEFINED',lead_audit_status='PENDING_LEAD_AUDIT_UNLESS_EXPLICIT_USER_SCOPED_ACCEPTANCE',cold_audit_status='SOURCE_REPORTED_NOT_RERUN',missing_inputs='See source; metric cells not populated from prose',metrics_scope='route-level evidence statement; not an accuracy row')
 for s in ['r','theta','z','v','psi']:row[s+'_metric_status']='NOT_EVALUATED_IN_THIS_ROW'
 master.append(row);cur.append({'stage':stage,'architecture':arch,'source':p,'last_change_sha':last,'category':cat,'gate':gate,'contract':contract,'limit':limit,'source_exists':exists})
# R3 final two candidates are reconstructed from saved set, not remembered narrative.
r3p=cite('results/R3_FINAL_PARAMETER_CONTRACTION/FINAL_SURVIVOR_STATE_TABLE_FIXED.csv');r3=read(r3p)
r3groups=collections.defaultdict(list)
for x in r3:r3groups[x['stage'],x['subset'],x['z_true']].append(x)
r3out=[]
for key,rr in r3groups.items():
 r3out.append({'stage':key[0],'subset':key[1],'source_depth_m':key[2],'count':len(rr),'r_labels_km':str(sorted({float(x['r']) for x in rr})),'theta_labels_deg':str(sorted({float(x['theta']) for x in rr})),'v_labels_mps':str(sorted({float(x['v']) for x in rr})),'psi_labels_deg':str(sorted({float(x['psi']) for x in rr})),'scope':'Frozen discrete surviving set; no continuous confidence or P95'})
write('CZ_D1_R3_FINAL_LABELS.csv',r3out)
# Fix inventory SHA: choose latest entrypoint in Git chronological traversal, never lexical hash order.
inv=read('results/CZ_D1/CZ_D1_STAGE_COMMIT_STATUS.csv');order=git('log',P,'--format=%H').splitlines();rank={s:i for i,s in enumerate(order)}
for x in inv:
 ps=x['entrypoints'].split('|') if x['entrypoints'] else []
 shas=[git('log','-1','--format=%H',P,'--',p) for p in ps]
 x['latest_entry_change_sha']=min(shas,key=lambda s:rank.get(s,10**9)) if shas else ''
 if not x['recorded_status_verbatim'] and ps:x['execution_status']='SOURCE_PRESENT_REQUIRES_MANUAL_STATUS';x['missing']='NO_AUTOMATIC_STATUS_EXTRACTION; not missing source'
write('CZ_D1_STAGE_COMMIT_STATUS.csv',inv)
allfiles=git('ls-tree','-r','--name-only',P).splitlines();write('CZ_D1_PROJECT_FILE_INVENTORY.csv',[{'path':p,'snapshot':P,'bytes':(R/p).stat().st_size,'kind':'RESULT' if p.startswith('results/') else 'CODE_DOC_ASSET_OR_DESIGN'} for p in allfiles])
# All design/contract docs discoverable, without labelling design as completed experiment.
designs=[p for p in allfiles if (p.startswith('docs/') or p.startswith('results/')) and any(s in Path(p).name.upper() for s in ['DESIGN','PROTOCOL','TASK','PLAN','FREEZE']) and Path(p).suffix.lower() in ['.md','.json','.txt']]
write('CZ_D1_DESIGN_AND_UNEXECUTED_INDEX.csv',[{'path':p,'snapshot':P,'execution_status':'DESIGN_ARTIFACT_ONLY; execution must be established from corresponding results'} for p in designs])
write('CZ_D1_EVIDENCE_MASTER.csv',master);write('CZ_D1_CURATED_ROUTES.csv',cur)
prov['files']=list(sources.values());prov['additional_error_checks']=len(checks);prov['additional_error_failures']=sum(not x['PASS'] for x in checks)
(O/'CZ_D1_SOURCE_PROVENANCE.json').write_text(json.dumps(prov,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'master':len(master),'curated':len(cur),'error_checks':len(checks),'error_failures':sum(not x['PASS'] for x in checks),'source_files':len(sources),'full_project_files':len(allfiles),'design_artifacts':len(designs)}))
