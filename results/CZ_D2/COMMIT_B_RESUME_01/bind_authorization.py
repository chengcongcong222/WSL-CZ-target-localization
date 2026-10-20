"""New explicit authorization + V3 binding only. Does not alter original contract."""
from common import *
import subprocess
assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip()==PARENT
assert subprocess.check_output(['git','branch','--show-current'],cwd=ROOT).decode().strip()=='research/cz-d2'
assert sha(V3)==V3HASH and len(csvread(V3))==762
for man in ['DESIGN_FILE_MANIFEST.json','COMMIT_B_STOP_MANIFEST.json']:
 for x in js(D/man)['files']:assert sha(D/x['path'])==x['sha256']
inputs=[D/n for n in ['EXPERIMENT_CONTRACT.md','DESIGN_FREEZE.json','PANEL_AND_HOLDOUT_POLICY.md','GATE_AND_STATISTICS.md','AUDIT_DESIGN.md','DESIGN_FILE_MANIFEST.json']]+[REF/n for n in ['t3_w1_model.py','t3_w1_generate.py','t3_w1_information.py']]+[V3,D/'HISTORY_REPAIR_02/SOURCE_COMPLETENESS_FINAL.json',D/'HISTORY_REPAIR_02/MFP_PARTIAL_HISTORY.csv']
dump(O/'EFFECTIVE_INPUT_BINDING.json',dict(design_SHA=DESIGN,repair02_SHA=PARENT,registry=dict(path=str(V3),sha256=V3HASH,rows=762),only_replaced_input='historical exclusion registry V1 -> V3 by research lead authorization',inputs=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in inputs],original_panel_and_statistics_unchanged=True,history_scope='FIXED_SNAPSHOT_AND_BOUND_ARCHIVES',MFP='PARTIAL_RZ_RECORDED_NOT_4D_EXCLUSION',H2='STRICT_FLOAT_EQUALITY_FAIL_RETAINED',production_allowed_after_preflight=True),True)
dump(O/'RESUME_AUTHORIZATION.json',dict(authority='Research lead current explicit CZ-D2 Commit B resume instruction',basis=PARENT,design=DESIGN,history_gate='CZ_D2_HISTORY_SOURCE_GATE_ACCEPTED_SCOPED',new_lock_directory=str(O),main_calls=1,cold_calls=1,retries=0,history_unchanged=True,CZ_D3_authorized=False,R4_percent=0),True)
print('V3 binding sealed; no panel generated')
