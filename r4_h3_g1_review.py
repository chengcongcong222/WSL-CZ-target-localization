"""H3-G1 documentation and read-only profiling of saved scores. No numerical experiment."""
from pathlib import Path
from decimal import Decimal
import json,csv,hashlib,subprocess,collections,ast
PARENT='0cfa6b6c7c6ff231d602522144c5a7b992cf03ad'
O=Path('results/R4_H3_G1_EXTRACTION_SUPPORT_REVIEW')
G=Path('results/R4_H3_G0_FINITE_VERTICAL_RESPONSE')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,data):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def sha(p,raw=False):
 b=Path(p).read_bytes();return hashlib.sha256(b if raw else b.replace(b'\r\n',b'\n')).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],text=True,encoding='utf-8').strip()
def main():
 assert git('rev-parse','HEAD')==PARENT
 assert git('ls-remote','origin','refs/heads/main').split()[0]==PARENT
 O.mkdir(parents=True,exist_ok=True)
 taskdir=Path('D:/CodexData/home/attachments/6fc1888e-c2ad-475b-944e-69b1ff8c0d60')
 task=next(x for x in taskdir.iterdir() if x.is_file())
 (O/'RESEARCH_LEAD_TASK.md').write_bytes(task.read_bytes())
 docs=read(O/'DESIGN_DOCUMENTS_PAYLOAD.json')
 for name,txt in docs.items():(O/name).write_text(txt,encoding='utf-8')
 source=G/'FINITE_HORIZONTAL_LABEL_DIAGNOSTIC.csv'
 score=rows(source);assert len(score)==4500
 scored=[]
 for n,r in enumerate(score,2):
  scored.append(dict(r,source_line=n))
 groups=collections.defaultdict(list)
 for r in scored:groups[r['scene_id'],r['resource'],r['noise_model']].append(r)
 profiles=[];summary=[];fixed=[]
 for key,rs in sorted(groups.items()):
  sid,res,noise=key
  for cal in ['C0','C1b']:
   field=cal+'_joint_chart_scale_1pct'
   nodes=collections.defaultdict(list)
   for r in rs:nodes[r['range_offset_m'],r['speed_offset_mps']].append(r)
   for node,rr in sorted(nodes.items()):
    sort=sorted(rr,key=lambda r:(Decimal(r[field]),Decimal(r['depth_label_m'])))
    fixed.append(dict(scene_id=sid,resource=res,noise_model=noise,calibration=cal,range_offset_m=node[0],speed_offset_mps=node[1],best_depth_label_m=sort[0]['depth_label_m'],second_depth_label_m=sort[1]['depth_label_m'],best_score_1pct=sort[0][field],second_score_1pct=sort[1][field],source_line=sort[0]['source_line'],scope='FIXED_NODE_CONDITIONAL_RANKING_NOT_RECOVERY'))
   for leave in [False,True]:
    subset=[r for r in rs if not leave or Decimal(r['range_offset_m'])!=0 or Decimal(r['speed_offset_mps'])!=0]
    carriers=[]
    for depth in ['180.0','190.0','200.0','210.0','220.0']:
     cand=[r for r in subset if r['depth_label_m']==depth]
     carrier=min(cand,key=lambda r:(Decimal(r[field]),Decimal(r['range_offset_m']),Decimal(r['speed_offset_mps'])))
     carriers.append(carrier)
     profiles.append(dict(scene_id=sid,resource=res,noise_model=noise,calibration=cal,center_excluded=leave,horizontal_node_count=len(cand),depth_label_m=depth,profile_score_1pct=carrier[field],profile_score_5pct_descriptive=str(Decimal(carrier[field])/5),range_offset_m=carrier['range_offset_m'],speed_offset_mps=carrier['speed_offset_mps'],source_line=carrier['source_line'],source_sha256=sha(source),scope='FINITE_SAVED_NOMINAL_TANGENT_CHART_MIN_NOT_GLOBAL_LIKELIHOOD'))
    ordered=sorted(carriers,key=lambda r:(Decimal(r[field]),Decimal(r['depth_label_m'])))
    full=sorted(subset,key=lambda r:(Decimal(r[field]),Decimal(r['depth_label_m']),Decimal(r['range_offset_m']),Decimal(r['speed_offset_mps'])))
    b,c=ordered[:2];gap=Decimal(c[field])-Decimal(b[field])
    summary.append(dict(scene_id=sid,resource=res,noise_model=noise,calibration=cal,center_excluded=leave,registered_horizontal_nodes=24 if leave else 25,best_depth_label_m=b['depth_label_m'],best_score_1pct=b[field],best_dr_m=b['range_offset_m'],best_dv_mps=b['speed_offset_mps'],second_depth_label_m=c['depth_label_m'],second_score_1pct=c[field],second_depth_dr_m=c['range_offset_m'],second_depth_dv_mps=c['speed_offset_mps'],score_gap=str(gap),gap_over_best=str(gap/Decimal(b[field])) if Decimal(b[field])!=0 else 'INF',nearest_other_candidate_depth_m=full[1]['depth_label_m'],nearest_other_candidate_score_1pct=full[1][field],nearest_other_dr_m=full[1]['range_offset_m'],nearest_other_dv_mps=full[1]['speed_offset_mps'],source_line=b['source_line'],second_source_line=c['source_line'],scope='DESCRIPTIVE_FINITE_PROFILE_NOT_RECOVERY_OR_CONFIDENCE'))
 write(O/'FINITE_DEPTH_PROFILES.csv',profiles)
 write(O/'FINITE_PROFILE_SUMMARY.csv',summary)
 write(O/'FIXED_NODE_LABEL_RANKING.csv',fixed)
 support=[
 dict(support='O',artifact=str(G/'DESIGN_FREEZE.json'),rule='original true center',coverage='ORACLE_NOT_OBSERVATION_SUPPORT',admitted=False),
 dict(support='F',artifact=str(source),rule='25 truth-offset r/v nodes; fixed theta/psi',coverage='FINITE_DIAGNOSTIC_ONLY',admitted=False),
 dict(support='U0_COARSE',artifact='results/R4_A1_OFFGRID_BEARING_BOUNDARY/R4_A1_CONFIG.json',rule='min SSE+13.3 sigma_rad^2; no-noise +1e-12',coverage='0/270 bracketing retention at 0.02 and 0.05deg; different panel',admitted=False),
 dict(support='U0_CONTINUOUS',artifact='r4_a1_fix_continuous.py',rule='32 starts plus4096 sampled r/v/psi bearing feasible points; no-noise1e-20',coverage='NOT_EXHAUSTIVE_CONTINUOUS_REGION',admitted=False),
 dict(support='U_AUX_STATIC',artifact='results/R4_AUX_GATE2B_JOINT_ERROR_BUDGET/GPT_SYNC.md',rule='static tested corners and11 range points',coverage='NOT_DYNAMIC_MULTIBRANCH_SUPPORT',admitted=False),
 dict(support='U_AUX_DYNAMIC',artifact='results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/GPT_SYNC.md',rule='observation-only point estimate; soft_l1',coverage='NAIVE_LOCAL_COVARIANCE_NOT_JOINT_COVERAGE',admitted=False)]
 write(O/'SUPPORT_PROVENANCE.csv',support)
 extras=[
 'r4_h3_g1_review.py',
 'r4_a1_fix_continuous.py','r4_a1_new_estimator.py','r4_a1_offgrid_bearing.py','r4_e1_g0.py','r4_e2_g0.py',
 'results/R4_E1_G0_FREQUENCY_INFORMATION/E1_G0_SCENE_BINDINGS.json',
 'results/R4_A1_OFFGRID_BEARING_BOUNDARY/R4_A1_CONFIG.json',
 'results/R4_A1_OFFGRID_BEARING_BOUNDARY/GPT_SYNC.md',
 'results/R4_A1_OFFGRID_BEARING_BOUNDARY/PANEL_RC2_BOUNDARY_SUMMARY.csv',
 'results/R4_A1_OFFGRID_BEARING_BOUNDARY/SURVIVOR_SET_QUALITY.csv',
 'results/P2_RC2_kinematic_boundary/P2_RC2_GPT_SYNC.md',
 'results/R4_A1_FIX2_ACOUSTIC_COVERAGE/METHOD_FREEZE.json',
 'results/R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_FINAL_DECISION.json',
 'results/R4_AUX_GATE2B_JOINT_ERROR_BUDGET/GPT_SYNC.md',
 'results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/GPT_SYNC.md',
 'results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_TRUTH_PANEL.csv',
 'results/R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY/SPEED_INFO_DECISION.json',
 'results/R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT/GPT_SYNC.md',
 'results/R4_H3_DEPTH_OBSERVABLE_REVIEW/H3_OBSERVABLE_AND_NUISANCE_MATRIX.md',
 'results/R4_H3_DEPTH_OBSERVABLE_REVIEW/PRIMARY_SOURCE_PROVENANCE.csv',
 str(task)]
 inputs=extras+[str(x) for x in G.iterdir() if x.is_file()]
 bindings={n:dict(kind='RAW' if Path(n).suffix=='.npz' or n==str(task) else 'LF_TEXT',sha256=sha(n,Path(n).suffix=='.npz' or n==str(task))) for n in inputs}
 dump(O/'READ_ONLY_INPUT_BINDINGS.json',dict(parent_SHA=PARENT,bindings=bindings,hash_only_inherited_modal_inputs=True,original_scientific_files_modified=False))
 oldfreeze=read(G/'DESIGN_FREEZE.json')
 counters={res:{noise:sum(r['resource']==res and r['noise_model']==noise and r['calibration']=='C1b' and r['center_excluded'] and Decimal(r['best_depth_label_m'])==200 for r in summary) for noise in ['RELATIVE','ABSOLUTE_FLOOR']} for res in ['B0','B1','B2']}
 decision=dict(stage='R4_H3_G1_EXTRACTION_AND_SUPPORT_READINESS_REVIEW',parent_SHA=PARENT,primary_decision='H3_DEPTH_SUPPORT_AND_PROPAGATION_UNRESOLVED',additional_classifications=['H3_EXTRACTION_SOURCE_CONDITION_INCOMPLETE','H3_EXTRACTION_CALIBRATION_CONDITIONAL'],unconditional_minimal_pilot_ready=False,nominal_extraction_design='BOUNDED_FREQUENCY_SAMPLE_DRAFT_FOR_SEPARATE_FREEZE_REVIEW',source_occupancy='NOT_ESTABLISHED_FOR_ACTUAL_UUV',source_frequency_and_coherence='NOT_ESTABLISHED_FOR_ACTUAL_UUV',calibration_stability='NOT_ESTABLISHED',physical_propagation_nuisance='PHYSICAL_PROPAGATION_NUISANCE_CONTRACT_NOT_ESTABLISHED',horizontal_support='FULL_RC2_HORIZONTAL_SUPPORT_NOT_ESTABLISHED',C1g_FREE_MODAL_GAIN='ZERO_DEPTH_CONTROL_RETAINED',C1g_constrained='NOT_ADMITTED_NO_INDEPENDENT_PHYSICAL_CONTRACT',C1e='NOT_EVALUATED',inherited_H3_G0='CONDITIONAL_LOCAL_DEPTH_INFORMATION_ACCEPTED',inherited_scope='DETERMINISTIC_MECHANISM_UPPER_BOUND_ONLY',identifiers=6,independent_main_geometries=3,read_only_profile_rows=len(profiles),read_only_profile_summaries=len(summary),fixed_node_rankings=len(fixed),center_excluded_C1b_source200_best=counters,finite_profile_scope='FINITE_SAVED_NOMINAL_TANGENT_CHART_MIN_NOT_GLOBAL_LIKELIHOOD',H3_extracted_observable='NOT_ESTABLISHED',H3_full_support_depth='NOT_ESTABLISHED',E2='PAUSED_WITH_FORMULA_REPAIR_FAILURE',original_E2_G0='FAIL_UNCHANGED',speed_route='CURRENT_1200S_GEOMETRY_INFORMATION_LIMIT_CONFIRMED_BY_CRLB_RETAINED',H3_A='RETAINED_NOT_EXECUTED',R4_percent=0,new_KRAKEN=0,new_FIELD=0,new_MC=0,new_audio=0,new_localization_runs=0,new_G0_science_runs=0,next='SEPARATE_NOMINAL_MINIMAL_EXTRACTION_DESIGN_FREEZE_REVIEW',next_execution='NOT_AUTHORIZED',stop_after_single_review_commit=True)
 dump(O/'H3_G1_READINESS_DECISION.json',decision)
 dump(O/'MINIMAL_PILOT_DRAFT.json',dict(status='DRAFT_NOT_EXECUTION_AUTHORIZATION',source_model='CONDITIONAL_GAUSSIAN_LINE_COEFFICIENTS_NOT_ACTUAL_UUV_EVIDENCE',scenes=['H01','H06','H12'],resources=['B0','B1','B2'],snapshots_s=[0,600,1200],primary_mesh=160001,secondary_mesh='80001_SAME_SAMPLES_NOT_INDEPENDENT_CONFIRMATION',frequency_packages=dict(F1=[200],F3=[150,200,250],F13=oldfreeze['frequencies_hz']),K=[1,8,32],realizations=16,noise_background_samples=256,noise_models=['RELATIVE','ABSOLUTE_FLOOR'],sigma=[.01,.05],absolute_reference=.00010824612978073539,fs_hz=1000,window='PERIODIC_HANN',segment_seconds=2,n_samples=2000,overlap=0,nominal_ENBW_hz=.75,noise_frequency_correlation='FULL_WINDOW_GRAM_NOT_DIAGONAL',calibration='C1b_UNKNOWN_FIXED_ACROSS_SNAPSHOTS',source_power='UNKNOWN_EACH_FREQ_SNAPSHOT_NOMINAL_ONE',main_statistic='FULL_JOINT_CSD_PLUS_NOISE_BACKGROUND_LIKELIHOOD',normalized_statistic_information='NOT_ESTABLISHED_WITHOUT_MATCHED_FINITE_SAMPLE_MODEL',raw_cell_count=5184,wall_cap_seconds=3600,storage_cap_GB=1,drift_phase_budget_example_rad=.01,clock_delay_example_microseconds=6.37,actual_clock_calibration='NOT_ESTABLISHED',new_propagation_calls=0,source_conditions='DESIGN_ASSUMPTIONS_NOT_REAL_TARGET_CONFIRMATION',horizontal_support='FULL_RC2_HORIZONTAL_SUPPORT_NOT_ESTABLISHED',propagation_nuisance='PHYSICAL_PROPAGATION_NUISANCE_CONTRACT_NOT_ESTABLISHED',next_execution='NOT_AUTHORIZED',continuous_depth_P95='NOT_OPENED',R4_percent=0))
 checks=[]
 def ck(name,ok,detail=''):checks.append(dict(check=name,PASS=bool(ok),detail=str(detail)))
 ck('raw_score_rows_4500',len(score)==4500)
 ck('profile_rows_720',len(profiles)==720)
 ck('profile_summary_rows_144',len(summary)==144)
 ck('fixed_node_rankings_1800',len(fixed)==1800)
 ck('all_input_bindings_unchanged',all(sha(n,r['kind']=='RAW')==r['sha256'] for n,r in bindings.items()))
 ck('all_inherited_G0_bindings_unchanged',all(sha(n,r['kind']=='RAW')==r['sha256'] for n,r in oldfreeze['bindings'].items()))
 ck('G0_formula_controls_remain1257_PASS',len(rows(G/'CONTROL_CHECKS.csv'))==1257 and all(r['PASS']=='True' for r in rows(G/'CONTROL_CHECKS.csv')))
 ck('G0_independent1272_PASS',read(G/'VALIDATION.json')['PASS']==1272 and read(G/'VALIDATION.json')['FAIL']==0)
 ck('G0_numeric180_stable',len(rows(G/'NUMERICAL_STABILITY.csv'))==180 and all(r['status']=='STABLE' for r in rows(G/'NUMERICAL_STABILITY.csv')))
 for r in profiles:
  line=int(r['source_line']);orig=score[line-2];col=r['calibration']+'_joint_chart_scale_1pct'
  cc=[x for x in score if x['scene_id']==r['scene_id'] and x['resource']==r['resource'] and x['noise_model']==r['noise_model'] and x['depth_label_m']==r['depth_label_m'] and (not r['center_excluded'] or Decimal(x['range_offset_m'])!=0 or Decimal(x['speed_offset_mps'])!=0)]
  # independent exact Decimal minimum; no new propagation/state score.
  minimum=min(Decimal(x[col]) for x in cc)
  ck('profile_source:'+str(len(checks)),r['profile_score_1pct']==orig[col] and Decimal(r['profile_score_1pct'])==minimum and len(cc)==r['horizontal_node_count'] and r['range_offset_m']==orig['range_offset_m'] and r['speed_offset_mps']==orig['speed_offset_mps'])
 for r in summary:
  rs=[x for x in profiles if all(x[k]==r[k] for k in ['scene_id','resource','noise_model','calibration','center_excluded'])]
  zz=sorted(rs,key=lambda x:(Decimal(x['profile_score_1pct']),Decimal(x['depth_label_m'])))
  ck('summary_minima:'+str(len(checks)),r['best_depth_label_m']==zz[0]['depth_label_m'] and r['second_depth_label_m']==zz[1]['depth_label_m'] and Decimal(r['score_gap'])==Decimal(zz[1]['profile_score_1pct'])-Decimal(zz[0]['profile_score_1pct']))
 for res in ['B0','B1','B2']:
  for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
   cnt=sum(r['resource']==res and r['noise_model']==noise and r['calibration']=='C1b' and Decimal(r['best_depth_label_m'])==200 for r in fixed)
   expected={'B0':{'RELATIVE':50,'ABSOLUTE_FLOOR':34},'B1':{'RELATIVE':62,'ABSOLUTE_FLOOR':64},'B2':{'RELATIVE':44,'ABSOLUTE_FLOOR':38}}[res][noise]
   ck('inherited_fixed_label_count:'+res+noise,cnt==expected)
 for case in ['H01','H06','H12']:
  m=[r for r in profiles if r['scene_id']==case+'_M']
  pp=[r for r in profiles if r['scene_id']==case+'_P']
  key=lambda r:(r['resource'],r['noise_model'],r['calibration'],r['center_excluded'],r['depth_label_m'],r['profile_score_1pct'],r['range_offset_m'],r['speed_offset_mps'])
  ck('mirror_duplicated_main_geometry:'+case,sorted(map(key,m))==sorted(map(key,pp)))
 required=['RECEIVED_SIGNAL_REQUIREMENTS.md','FINITE_SAMPLE_OBSERVABLE_CHAIN.md','PROPAGATION_NUISANCE_REQUIREMENT.md','HORIZONTAL_SUPPORT_AUDIT.md','FINITE_LABEL_PROFILE_SCOPE.md','LOW_SNR_EXTRACTION_LIMITS.md','MINIMAL_EXTRACTION_PILOT_DESIGN.md','H3_G1_READINESS_DECISION.json','GPT_SYNC.md']
 for name in required:ck('required:'+name,(O/name).is_file() and (O/name).stat().st_size>0)
 ck('no_wrong_information_or_project_admission',not decision['unconditional_minimal_pilot_ready'] and decision['R4_percent']==0 and decision['next_execution']=='NOT_AUTHORIZED')
 ck('pilot_resource_count',3*3*3*2*2*3*16==5184)
 ck('no_generated_propagation_or_recording_files',not any(x.suffix in ['.mod','.npz','.wav','.flac','.shd'] for x in O.iterdir()))
 allowed={'pathlib','decimal','json','csv','hashlib','subprocess','collections','ast'}
 tree=ast.parse(Path('r4_h3_g1_review.py').read_text(encoding='utf-8'))
 mods={n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)}|{x.name for n in ast.walk(tree) if isinstance(n,ast.Import) for x in n.names}
 ck('review_script_standard_library_only',mods<=allowed)
 write(O/'DOCUMENTATION_CHECKS.csv',checks)
 val=dict(checks=len(checks),PASS=sum(x['PASS'] for x in checks),FAIL=sum(not x['PASS'] for x in checks),scope='DOCUMENTATION_PROVENANCE_AND_EXACT_SAVED_SCORE_AGGREGATION_ONLY',new_experiments=0,read_only_existing_score_rows=4500,not_scientific_performance_validation=True,original_G0_unchanged=True)
 dump(O/'DOCUMENTATION_VALIDATION.json',val);assert val['FAIL']==0
 master=Path('results/R4_MASTER')
 marker='## H3-G1 extraction and support readiness review'
 text=(master/'R4_PLAN.md').read_text(encoding='utf-8')
 if marker not in text:
  with (master/'R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\n'+marker+'\n\nH3_DEPTH_SUPPORT_AND_PROPAGATION_UNRESOLVED; source and calibration conditional gaps retained. Accepted G0 conditional local mechanism unchanged. Read-only finite score profiles only; no new experiments. A bounded nominal finite-frequency-sample design is proposed, not released. No full RC2 support or constrained propagation-gain contract; no actual extracted observable. H3-B remains priority candidate; H3-A and speed bottleneck retained. R4=0%. One review commit/push/verify and STOP; next execution NOT_AUTHORIZED.\n')
 ledger=master/'R4_EVIDENCE_LEDGER.csv'
 if decision['stage'] not in ledger.read_text(encoding='utf-8'):
  with ledger.open('a',encoding='utf-8',newline='') as f:csv.writer(f).writerow([decision['stage'],'../R4_H3_G1_EXTRACTION_SUPPORT_REVIEW/H3_G1_READINESS_DECISION.json',decision['primary_decision'],'DESIGN_AND_READ_ONLY_SAVED_SCORE_PROFILE;NO_EXPERIMENT;NO_CREDIT','PENDING_RESEARCH_LEAD_AUDIT;STOP',0])
 dump(master/'R4_H3_G1_EXTRACTION_SUPPORT_PROGRESS.json',dict(stage=decision['stage'],parent_SHA=PARENT,decision=decision['primary_decision'],additional_classifications=decision['additional_classifications'],H3_G0='CONDITIONAL_LOCAL_DEPTH_INFORMATION_ACCEPTED',R4_percent=0,next_execution='NOT_AUTHORIZED',review_status='COMPLETED_STOP_AFTER_COMMIT'))
 manifest=[dict(file=str(x),kind='LF_TEXT',sha256=sha(x)) for x in sorted(O.iterdir()) if x.is_file() and x.name!='DELIVERABLE_MANIFEST.json']
 dump(O/'DELIVERABLE_MANIFEST.json',dict(files=manifest,parent_SHA=PARENT,scope='REVIEW_DELIVERABLES_NOT_NEW_SCIENTIFIC_ADMISSION',R4_percent=0))
 print(json.dumps(dict(decision=decision['primary_decision'],validation=val,read_only_center_exclusion=counters,deliverables=len(manifest)),indent=2))
if __name__=='__main__':main()
