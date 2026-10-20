"""Read-only semantic closeout. Fixed historical snapshot; no scientific imports."""
from pathlib import Path
import csv,json,hashlib,collections,re,math
R=Path(__file__).resolve().parents[3];O=Path(__file__).resolve().parent;OLD=O.parent/'HISTORY_REPAIR_01'
def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(n,x):
 with (O/n).open('x',encoding='utf-8',newline='') as f:
  if n.endswith('.json'):json.dump(x,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
  else:
   w=csv.DictWriter(f,fieldnames=list(x[0]),lineterminator='\n');w.writeheader();w.writerows(x)
V2=rows(OLD/'HISTORICAL_STATE_EXCLUSIONS_V2.csv');v3=[dict(x) for x in V2];axes=['r_km','theta_deg','v_mps','psi_deg'];evidence=[];origin_rows=[];added=[]
def add_origin(path,index,h,purpose,proof):
 h=[None if x is None else float(x) for x in h];kind='static_position_projection' if h[2] is None else 'horizontal_four_state';found=None
 for i,x in enumerate(v3):
  a=[float(x[k]) if x[k] else None for k in axes]
  if all((u is None and v is None) or (u is not None and v is not None and abs(u-v)<=1e-11) for u,v in zip(h,a)):found=i;break
 if found is None:
  row=dict(source=path,record=str(index),kind=kind,**{k:'' if v is None else repr(v) for k,v in zip(axes,h)});v3.append(row);added.append(row);found=len(v3)-1
 origin_rows.append(dict(source=path,source_sha256=sha(R/path),truth_index=str(index),purpose=purpose,units='km,degree,m/s,degree; static v/psi absent',**{k:'' if v is None else repr(v) for k,v in zip(axes,h)},registry_row_zero_based=found,coverage='V2' if found<len(V2) else 'V3_ADDITION',code_evidence=proof))
# Existing archived outputs, not regenerated truths.
p='results/P2_RC2_kinematic_boundary/bearing_accuracy_time_boundary.csv'
for i,x in enumerate(rows(R/p)):add_origin(p,i,[x['truth_r0_km'],0,x['truth_v'],x['truth_psi_deg']],'BEARING_OBSERVATION_GENERATING_TRUTH','p2_rc2_kinematic_boundary.py:scenario_truth; P2_RC2_CONFIG.json/scenarios all theta0=0')
p='results/P2_RC2_kinematic_boundary/heading_identifiability_boundary.csv';cfg=read(R/'results/P2_RC2_kinematic_boundary/P2_RC2_CONFIG.json')['baseline_truth']
for i,x in enumerate(rows(R/p)):add_origin(p,i,[cfg['r0_km'],cfg['theta0_deg'],cfg['v_mps'],x['psi_deg']],'HEADING_SCAN_ACTUAL_NOISY_OBSERVATION_TRUTH','p2_rc2_kinematic_boundary.py:573; truth_dict(50e3,0,2,psi_deg); scan generated observations, not replica-only grid')
p='results/P3_RC3_increment/selected_hard_pairs.csv'
for i,x in enumerate(rows(R/p)):add_origin(p,i,[x[k] for k in ['ref_r0_km','ref_theta0_deg','ref_v','ref_psi_deg']],'REFERENCE_SIDE_GENERATING_TRUTH','p3_rc3_increment.py:run_one_pair; simulate_observations(ref); alt only scored')
p='results/P4_G4_performance_boundary/P4_G4_CONFIG.json';c=read(R/p)
for key,x in [('B0',c['B0'])]+[('scenarios/'+k,v) for k,v in c['scenarios'].items()]:add_origin(p,key,[x['r0_m']/1000,x['theta0_deg'],x['v'],x['psi_deg']],'NOMINAL_AND_INTEGRATED_GENERATING_TRUTH','p4_g4_performance_boundary.py:truth_from/run_integrated; p4_g4_revision.py')
for p in ['results/R3_RC23_CLOSEDLOOP/R3_RC23_SINGLE_WINDOW_SET_CONTRACTION/CLOSEDLOOP_1A_CONFIG.json','results/R3_RC23_CLOSEDLOOP/R3_RC23_TWO_WINDOW_TEMPORAL_CONTRACTION/CLOSEDLOOP_1B_CONFIG.json']:
 x=read(R/p)['truth'];add_origin(p,'/truth',[x['r0_m']/1000,x['theta0_deg'],x['v'],x['psi_deg']],'RC23_GENERATING_TRUTH','r3_closedloop_1a_set_contraction.py; r3_closedloop_1b_two_window.py; downstream1C/1D/1E reuse TRUTH')
p='results/R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY/B1_DESIGN_FREEZE.json'
for i,x in enumerate(read(R/p)['cases']):add_origin(p,f'/cases/{i}/horizontal_reference',x['horizontal_reference'],'GENERATING_LEVEL_REFERENCE','r4_b1_execute.py:reference=case[horizontal_reference]; raw_levels(reference); B1A/B1B conditioning lattice is candidate-only')
# Static AUX target explicitly is (R,0); zeros here are defined by model, unlike MFP.
auxfiles=[('R4_ROUTE_RESET_AUXILIARY_NODE_GATE0','AUX_GEOMETRY_RESULTS.csv','r4_aux_gate0_geometry.py'),('R4_AUX_GATE1_CROSSTRACK_REQUIREMENT','AUX_CROSSTRACK_RESULTS.csv','r4_aux_gate1_crosstrack.py'),('R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET','AUX_BIAS_RESULTS.csv','r4_aux_gate2a.py'),('R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET','AUX_DEPLOYMENT_RESULTS.csv','r4_aux_gate2a.py'),('R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET','AUX_POSITION_RESULTS.csv','r4_aux_gate2a.py'),('R4_AUX_GATE2B_JOINT_ERROR_BUDGET','AUX_T5_JOINT_RESULTS.csv','r4_aux_gate2b.py'),('R4_AUX_GATE2B_JOINT_ERROR_BUDGET','AUX_T10_JOINT_RESULTS.csv','r4_aux_gate2b.py')]
for family,name,code in auxfiles:
 p=f'results/{family}/{name}';seen=set()
 for i,x in enumerate(rows(R/p)):
  if x.get('geometry_status')=='INFEASIBLE_TRIANGLE':continue
  r=float(x['range_km'])
  if r in seen:continue
  seen.add(r);add_origin(p,i,[r,0,None,None],'STATIC_BEARING_CALIBRATION_OR_PERFORMANCE_GENERATING_POSITION',code+': target=(R,0), MAIN=(0,0); auxiliary placement side is not target azimuth; no target velocity defined')
put('APPLICABLE_SOURCE_ROW_COVERAGE.csv',origin_rows);put('HISTORICAL_STATE_EXCLUSIONS_V3.csv',v3);put('V2_V3_ADDITIONS.csv',added)
put('REGISTRY_DIFF.json',dict(V1_rows=216,V2_rows=len(V2),V3_rows=len(v3),added_rows=len(added),V2_prefix_exact=v3[:len(V2)]==V2,V2_sha256=sha(OLD/'HISTORICAL_STATE_EXCLUSIONS_V2.csv'),V3_sha256=sha(O/'HISTORICAL_STATE_EXCLUSIONS_V3.csv'),policy='No new threshold; explicitly register additional actual grid-aligned truths already blocked by old coordinate-grid guard',generation_calls=0))
# MFP policy and partial-state evidence, grouped across three source signal conditions.
p='results/R3_RC3_method_recovery/MFP_UPPER_BOUND_summary.csv';rr=rows(R/p);partial=[]
for (r,z) in sorted({(float(x['r_true_km']),float(x['z_true_m'])) for x in rr}):
 ids=[i for i,x in enumerate(rr) if float(x['r_true_km'])==r and float(x['z_true_m'])==z]
 partial.append(dict(source=p,sha256=sha(R/p),rows=json.dumps(ids),r_km=r,z_m=z,theta='NOT_DEFINED',v='NOT_DEFINED',psi='NOT_DEFINED',signal_conditions=';'.join(rr[i]['source'] for i in ids),classification='PARTIAL_RZ_GENERATING_TRUTH_RECORDED_NOT_4D_EXCLUSION',reason='Actual pressure generation; no horizontal bearing/velocity/heading identity and no MFP pressure observation reused by CZ-D2',code='r3_0_method_recovery.py:run_mfp_upper_bound/pressure_spectrum'))
put('MFP_PARTIAL_HISTORY.csv',partial)
put('APPLICABILITY_POLICY_SUPPLEMENT.json',dict(authority='Research lead current REPAIR-02 instruction',original_A_unchanged=True,full_r_theta_v_psi='original4D+jointmirror',static_r_theta='original2D+azimuthmirror',partial_r_z='PARTIAL_RZ_GENERATING_TRUTH_RECORDED_NOT_4D_EXCLUSION',candidate_grid='NOT_GENERATING_TRUTH',numeric_reference_perturbations='NOT_GENERATING_TRUTH_UNLESS_ACTUALLY_USED_AS_OBSERVATION_GENERATOR',noise='NOT_TARGET_STATE',dimensionless_math_control='NOT_PHYSICAL_TARGET_STATE',equipment_direction_calibration='NO_FABRICATED_RANGE_OR_VELOCITY',uncertain_identity='UNRESOLVED_AND_STOP',no_new_range_only_guard=True,does_not_claim_individually_unseen_coordinate_values=True,formal_execution_authorized=False))
# Family dependency declarations reviewed against generator source and archived contracts.
groups=[]
def group(gid,prefixes,role,origins,proof):
 for p in proof:assert (R/p).is_file(),p
 groups.append(dict(id=gid,prefixes=prefixes,role=role,origin=origins,proof=[dict(path=p,sha256=sha(R/p)) for p in proof]))
group('P2',['P2_RC2_'],'mixed actual truth vs replica grid','V2 P2 offgrid; V3 P2 actual heading/bearing truths',['p2_rc2_kinematic_boundary.py','results/P2_RC2_kinematic_boundary/P2_RC2_CONFIG.json'])
group('P3',['P3_RC3_'],'ref generates observations; alt is candidate','V3 selected_hard_pairs ref rows only',['p3_rc3_increment.py','results/P3_RC3_increment/selected_hard_pairs.csv'])
group('P4',['P4_'],'random performance truth and fixed nominal truth','V2 P4 trials90; V3 P4_G4 nominal; candidate grids excluded',['p4_5_parameter_accuracy.py','p4_g4_performance_boundary.py','p4_g4_revision.py'])
group('R3_PARTIAL',['R3_A1_','R3_B1_','R3_C1_','R3_C2_','R3_RC3_','RANGE_UB'],'radial-depth/source-frequency/mode references, no full horizontal trajectory identity','Partial radial/depth history retained; no fabricated bearing/heading; search and numerical field grids not generating4Dtruth',['r3_0_method_recovery.py','r3_a1_fidelity.py','r3_b1_mmac_observability.py','r3_c1_depth_motion.py','r3_c2_2c_paper_repro.py','r3_c2_4b_2_paper_spectrum.py','r3_range_ub_1b_tdoa_oracle.py','r3_rc3_reanchor_2_profiled_discrim.py'])
group('R3_CLOSED',['R3_RC23_','R3_FINAL_'],'full nominal trajectory and downstream candidate sets','[50,0,2,5] from frozen truth; heading4 survivor is not a new generating truth',['r3_closedloop_1a_set_contraction.py','r3_closedloop_1b_two_window.py','r3_closedloop_1c_single_turn.py','r3_closedloop_1d_freq_subset.py','r3_closedloop_1e_a_amplitude.py','r3_closedloop_1e_b_tracked_freq.py','r3_closedloop_1e_c_ssp.py','r3_final_parameter_contraction_fix.py'])
group('A1',['R4_A1_'],'P/Q/F/H generators, fit endpoints and caches','V2 both Q implementations, P9,F8,H12; A1-new and speed reuse H12; unused starts/candidate cells excluded',['r4_a1_offgrid_bearing.py','r4_a1_fix_continuous.py','r4_a1_fix2_coverage.py','r4_a1_new_dynamic.py','r4_a1_search_tractability_development.py','r4_speed_diagnostic.py','r4_speed_information.py'])
group('AUX_STATIC',['R4_ROUTE_RESET_AUXILIARY','R4_AUX_GATE'],'actual static bearing observations at target(R,0)','V3 observed ranges50..60,theta0; normal draws and node-position perturbations are not target state',['r4_aux_gate0_geometry.py','r4_aux_gate1_crosstrack.py','r4_aux_gate2a.py','r4_aux_gate2b.py'])
group('AUX_D1',['R4_AUX_D1_'],'new historical random calibration and H12 test reuse','V2 AUX512, H12 repeated128; starts/endpoints are not truths',['aux_d1_generate.py','aux_d1_model.py','results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/PRIVATE_GENERATOR_MANIFEST.json'])
group('B1',['R4_B1'],'nominal truth plus conditioning candidates','V3 B1 horizontal_reference=[50,0,2,5]; altered r/v grid is a candidate used against fixed observation',['r4_b1_execute.py','r4_b1a_execute.py','r4_b1b_analyze.py'])
group('E_H3',['R4_E1_','R4_E2_','R4_H3_'],'H12/mirrors fixed scene inputs; finite differences and modal arrays','V2 H12 and joint mirrors. H3-G2E selected H01/H06/H12 reuse G0; perturbed candidate/Jacobian values not independent truth',['r4_e1_g0.py','r4_e2_g0.py','r4_e2_nine_pair_pilot.py','r4_h3_g0.py','r4_h3_g2e.py','r4_h3_p0.py','r4_h3_p1.py','results/R4_E1_G0_FREQUENCY_INFORMATION/E1_G0_SCENE_BINDINGS.json'])
group('RC2',['R4_RC2_'],'archived single-node P observations; dual-node input unavailable','V2 P-panel; interval boxes are candidate supports; no replacement synthetic truths',['r4_rc2_support.py','r4_rc2_interval.py'])
group('ENVELOPE',['R4_ROUTE_REDESIGN_'],'design/pre-run interval controls; development locked','No registered new scientific panel execution; arithmetic cells and model argument points not test truths',['r4_cz_envelope_development.py','r4_cz_envelope_model.py','r4_cz_envelope_closure_audit.py'])
group('H1',['R4_SINGLE_HLA_SOURCE_'],'main6 and offgrid8 actual generated observations','V2 H1 MAIN and OFFGRID CSV; full grid and multistarts are candidates',['hla_h1_execute.py','hla_h1_core.py'])
group('H2_RX',['R4_SINGLE_HLA_RADIAL_','R4_SINGLE_HLA_CARRIER_'],'H1 offgrid8 reuse and noises/source parameters','V2 H1 offgrid; accepted H2<=2ULP representation exception retained; RX1Z uses saved records',['hla_h2_generate.py','hla_rx1_execute.py','hla_rx1_generate.py','hla_rx1z_evaluate.py'])
group('SP1',['R4_SINGLE_HLA_STATIC_'],'16 random static calibration plus H1 offgrid8 test','V2 SP1 calibration16 and H1offgrid; private manifests accepted4818check chain',['hla_sp1_preflight.py','hla_sp1_generate.py','results/R4_SINGLE_HLA_STATIC_SPATIAL_RANGE_PILOT/CALIBRATION_PRIVATE_ARCHIVE_MANIFEST.json','results/R4_SINGLE_HLA_STATIC_SPATIAL_RANGE_PILOT/TEST_PRIVATE_ARCHIVE_MANIFEST.json'])
group('T2',['R4_T2_'],'static calibration/reference fields; performance stopped','V2 T2 CAL16; H12 selected scenes; N1/N2 probes and precision sentinels are numerical references, not generated performance observations',['t2_tilt_g0_prepare.py','t2_tilt_g0_forward.py','t2_tilt_g0_generate.py','t2_n1_analyze.py','t2_n2_evaluate.py','t2_n3_run.py'])
group('T3',['R4_T3_'],'H12 reused across dual-node schemes','V2 H12; new noise/config does not create unseen target state',['results/R4_T3_W1_WIDE_BASELINE_DEMO/t3_w1_generate.py','results/R4_T3_W1_WIDE_BASELINE_DEMO/T3_W1_DESIGN_FREEZE.json'])
group('CZO',['ENV_CAP_CZO04'],'36 actual generating cases, depth/environment labels separate','V2 DESIGN_FREEZE.generating_cases; per-case state arrays share these labels',['results/ENV_CAP_CZO04/scripts/analysis.py','results/ENV_CAP_CZO04/DESIGN_FREEZE.json'])
group('CAP',['ENV_CAP_CAP06'],'development/holdout16+mirrors; P1 mathtoy; directional stress','V2 P0 panel16 and mirrors; stress zeros placeholders not target; calibration angles lack range. Accepted private archive mapping retained',['results/ENV_CAP_CAP06/code/generate_p2.py','results/ENV_CAP_CAP06/code/p1_toy.py','results/ENV_CAP_CAP06/P0_FROZEN_PANEL.json','results/ENV_CAP_CAP06/P2_RAW_DATA_INDEX.csv'])
group('ENV_GEOMETRY',['ENV_CAP_FULL_','ENV_CAP_REC01','ENV_CAP_CZROB02','ENV_CAP_NUM03'],'SSP/ray/energy reference and shallow paper controls','Partial radial/depth observations or numerical references; no 4D target generator or static r/theta performance panel',['results/ENV_CAP_NUM03/GPT_SYNC_NUM03.md','results/ENV_CAP_REC01/GPT_SYNC_ENV_CAP_REC01.md','results/ENV_CAP_FULL_CAMPAIGN/GPT_SYNC_ENV_CAP_FINAL.md'])
group('DOC_REVIEW',['CZ_D1','ENV_CAP_REQ05','R4_MASTER','R4_OBSERVABILITY_','R4_AUX_CLIENT_','R4_SINGLE_HLA_CLOSEOUT_'],'historical aggregation/review, no new target observations','Follow canonical upstream origin; evidence row count is not number of new truths',['results/CZ_D1/GPT_SYNC_CZ_D1.md'])
put('SEMANTIC_FAMILY_BINDINGS.json',groups)
# Candidate source universe: all detected state/truth/calibration/holdout mentions or state-shaped schemas.
inv={x['path']:x for x in rows(OLD/'SNAPSHOT_INVENTORY.csv')};candidates=collections.defaultdict(set)
for n,key in [('GENERATOR_USAGE_DISCOVERY.csv','line'),('JSON_STATE_KEY_DISCOVERY.csv','pointer')]:
 for x in rows(OLD/n):candidates[x['path']].add(n)
for x in rows(O/'STRUCTURED_STATE_OBJECTS.csv'):candidates[x['path']].add('STRUCTURED_STATE_OBJECTS')
for x in rows(OLD/'CSV_SCHEMA_DISCOVERY.csv'):
 f=json.loads(x['fields']);text='|'.join(f)
 if x['candidate_truth_fields']=='True' or x['named_panel']=='True' or (re.search(r'(^|\|)(r|r0|r_km|r0_km)(\||$)',text) and re.search(r'theta|psi|v_mps',text)):candidates[x['path']].add('CSV_STATE_SCHEMA')
for x in rows(OLD/'ARCHIVE_ARRAY_HEADERS.csv'):
 if x['candidate_state_key']=='True' or re.search(r'(^|/)(state|states|h|target|parameters)\.npy',x['member']) or x['member']=='':candidates[x['path']].add('ARRAY_STATE_SCHEMA')
# Root code names route to reviewed dependency families; helper/audit code cannot be blanket treated as new truth.
rootmap=[('p2_','P2'),('p3_','P3'),('p4_','P4'),('r3_closedloop','R3_CLOSED'),('r3_final','R3_CLOSED'),('r3_','R3_PARTIAL'),('r4_a1_','A1'),('r4_speed_','A1'),('aux_d1_','AUX_D1'),('r4_aux_','AUX_STATIC'),('r4_b1','B1'),('r4_e','E_H3'),('r4_h3','E_H3'),('r4_rc2','RC2'),('r4_cz_','ENVELOPE'),('hla_h1','H1'),('hla_h2','H2_RX'),('hla_rx','H2_RX'),('hla_sp1','SP1'),('t2_','T2'),('prec','T2')]
classified=[];unresolved=[]
for path,triggers in sorted(candidates.items()):
 g=None;reason='';role=''
 if path.startswith('tests/'):
  role='SOFTWARE_ARITHMETIC_FIXTURE_NOT_RESEARCH_TRUTH';reason='Unit/control fixture: no registered performance/calibration observation archive; counterfactual labels only affect diagnostics.';gid='TEST_FIXTURE'
 elif path.startswith('results/'):
  family=path.split('/')[1];matches=[a for a in groups if any(family.startswith(v) for v in a['prefixes'])]
  if len(matches)==1:g=matches[0];gid=g['id']
  else:gid='UNRESOLVED'
 else:
  gid=next((b for a,b in rootmap if path.startswith(a)),'UNRESOLVED');g=next((a for a in groups if a['id']==gid),None)
 if g:role=g['role'];reason=g['origin']
 if gid=='UNRESOLVED':unresolved.append(path)
 classified.append(dict(path=path,sha256=inv[path]['sha256'],discovery_triggers=';'.join(sorted(triggers)),family=gid,semantic_role=role,disposition='UNRESOLVED' if gid=='UNRESOLVED' else 'CLASSIFIED_WITH_GENERATOR_PROVENANCE',exclusion_basis=reason,evidence='SEMANTIC_FAMILY_BINDINGS.json#'+gid))
put('SOURCE_SEMANTIC_ADJUDICATION.csv',classified);put('UNRESOLVED_SOURCES.json',unresolved)
print(json.dumps(dict(V3_rows=len(v3),new_rows=len(added),applicable_evidence_rows=len(origin_rows),candidate_sources=len(classified),unresolved=unresolved)))
