"""Full identity precheck and initial freeze, before any scientific evaluation."""
from model import *
import subprocess,time,ast

def git(*a):return subprocess.check_output(['git','-c','safe.directory='+str(R),*a],cwd=R).decode('utf-8').strip()
assert git('rev-parse','HEAD')=='ebdf94c2927bd1b5d4d459befaf96fbdc2bf454a'
assert git('branch','--show-current')=='research/env-cap-czo04'
assert git('diff','--name-only')=='' and git('diff','--cached','--name-only')==''
inputs={}
for a in load(R/'results/ENV_CAP_NUM03/DESIGN_FREEZE.json')['execution_inputs']:inputs[a['path']]=dict(a)
with (R/'results/ENV_CAP_NUM03/FILE_HASHES_SHA256.csv').open(encoding='utf-8-sig') as f:
 for a in csv.DictReader(f):
  p=R/'results/ENV_CAP_NUM03'/a['path'];inputs[str(p)]=dict(path=str(p),sha256=a['sha256'],source='NUM03_result_manifest')
extra=['results/R4_ROUTE_REDESIGN_GATE_1_R2_OBSERVABLE/GPT_SYNC.md','results/R4_ROUTE_REDESIGN_GATE_1_R2_OBSERVABLE/R2_MEASUREMENT_CHAIN_AUDIT.md','results/R4_ROUTE_REDESIGN_GATE_2_MEASUREMENT_CONTRACT/CZ_ENVELOPE_MEASUREMENT_CONTRACT.md','results/R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN/GPT_SYNC_CORRECTED.md','results/R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN/CZ_CLOSURE_BUDGET_LIMITATION.md','results/R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN/CZ_GRID_FREEZE.json','results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/DESIGN_FREEZE.json','hla_h1_core.py','results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_TRUTH_PANEL.csv']
for v in extra:
 p=R/v;inputs[str(p)]=dict(path=str(p),sha256=sha(p),source='fixed_parent_historical_source')
checks=[]
for k,a in inputs.items():checks.append(dict(path=k,sha256=a['sha256'],PASS=Path(k).is_file() and sha(k)==a['sha256']))
assert all(a['PASS'] for a in checks)
layout=load(O/'INPUT_LAYOUT_PREFLIGHT.json')
for v in layout['files']:assert any(a['path']==v['path'] and a['sha256']==v['sha256'] for a in inputs.values()),v['path']
for p in (O/'scripts').glob('*.py'):ast.parse(p.read_text(encoding='utf-8-sig'))
assert load(O/'SYNTHETIC_CONTROLS.json')['PASS'] and load(O/'IMPLEMENTATION_PREFLIGHT_CONTROLS.json')['PASS']
assert all(not p.exists() for p in [W/'ANALYSIS_LAUNCH.json',W/'COLD_LAUNCH.json',O/'DESIGN_FREEZE.json'])
dump(O/'PREFLIGHT.json',dict(PASS=True,all_identity_checks=len(checks),historical_tracked_diff='',synthetic_controls=35,independent_implementation_controls=len(load(O/'IMPLEMENTATION_PREFLIGHT_CONTROLS.json')['checks']),new_science_started=False,parent=git('rev-parse','HEAD'),branch=git('branch','--show-current')))
table(O/'INPUT_IDENTITY_PRECHECK.csv',checks)
for p in sorted(O.rglob('*')):
 if p.is_file():inputs[str(p)]=dict(path=str(p),sha256=sha(p),source='CZO04_frozen_contract_code_or_precheck')
cfg=dict(stage='ENV_CAP_CZO04',parent_SHA=git('rev-parse','HEAD'),branch=git('branch','--show-current'),freeze_unix=time.time(),execution_inputs=list(inputs.values()),science_wall_max_s=3600,process_tree_RSS_max_bytes=4*1024**3,new_D_max_bytes=2*1024**3,new_Git_max_bytes=64*1024**2,formal_analyses=1,cold_audits=1,retry=0,new_solvers=0,new_MC=0,new_receive_records=0,receiver_depth_m=200,frequency_hz=200,main_feature='121 RAW SINGLE_DEPTH DB LEVELS; NO SMOOTHING; LINEAR RANGE INTERPOLATION',generating_cases=states(),grid_shape=[16,21,11,31],grid_origin=[45,-5,1,-15],grid_step=[1,.5,.2,1],depth_labels_m=[150,200,250],environment_labels=['JW','JA','ANNUAL'],noise_sigma_db=[.1,.25,.5],rho=[0,.9],primary=dict(condition='C1',sigma_db=.25,rho=0),local_numerics=dict(J_grid=.1,effective_info_grid=.2,FD_step=.02,template_RMS_db=.0625),state_scale=SCALE.tolist(),FD_step=STEPS.tolist(),numerical_rank_relative_threshold=1e-10,weak_or_out_of_grid='UNRESOLVED_NOT_EXCLUDED',contract_sha256=sha(O/'RESEARCH_CONTRACT.md'))
dump(O/'DESIGN_FREEZE.json',cfg)
print(json.dumps(dict(PASS=True,historical_inputs=len(checks),frozen_inputs=len(inputs),scenarios=len(cfg['generating_cases']),grid_nodes=len(grid()))))
