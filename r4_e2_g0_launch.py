"""Pre-run infrastructure correction: attributes binding registered raw CRLF bytes.
No mathematical/provider/noise/array/candidate code changes. Guard rejection before
EXECUTION_STARTED means zero scientific executions. Verify exact raw bytes instead.
"""
from pathlib import Path
import json,hashlib
import r4_e2_g0 as screen
path='results/R4_E2_G0_HLA_DIFFERENCE_INFORMATION/.gitattributes'
freeze=screen.read(screen.OUT/'E2_G0_DESIGN_FREEZE.json')
assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==freeze['bindings'][path]
old_sha=screen.sha
def corrected_sha(p,raw=False):
 return old_sha(p,True) if str(p).replace('\\','/')==path else old_sha(p,raw)
screen.dump(screen.OUT/'PRE_EXECUTION_INFRASTRUCTURE_CORRECTION.json',dict(reason='attributes binding was registered raw CRLF, verifier defaulted LF-canonical; raw binding exactly matches',file=path,raw_sha256=freeze['bindings'][path],scientific_code_sha256=old_sha('r4_e2_g0.py'),first_launch='guard rejected before EXECUTION_STARTED; scientific execution count0',mathematical_changes=False,extra_scientific_execution=False,launcher_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()))
screen.sha=corrected_sha
screen.execute()