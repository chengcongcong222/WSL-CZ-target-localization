from hla_rx1z_common import *
import subprocess,ast

def main():
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==PARENT
 assert subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0]==PARENT
 assert sha(OUT/'DELEGATED_TASK.md')=='439da556ffbd3e610d1450e92fd4d455c831d0f5d27bb44e8bdb6fd02983a06c'
 assert not (OUT/'EXECUTION_STARTED.json').exists()
 saved=rows(OUT/'PUBLIC_SAVED_CARRIER_ROWS.csv');index=rows(OUT/'PUBLIC_RECORD_INDEX.csv');assert len(saved)==2304 and len(index)==128
 assert len({(r['record_id'],r['line'],r['window']) for r in saved})==2304
 for r in index:assert Path(r['path']).exists() and sha(r['path'])==r['sha256']
 for r in rows(OUT/'INPUT_PROVENANCE.csv'):assert sha(ROOT/r['path'])==r['working_sha256']
 sources=sorted(ROOT.glob('hla_rx1z_*.py'))
 for p in sources:ast.parse(p.read_text(encoding='utf-8-sig'))
 dump(OUT/'PREFLIGHT.json',{'PASS':True,'time_utc':stamp(),'parent':PARENT,'raw_records_hash_verified':128,'carrier_rows':2304,'no_prior_execution':True,'source_parse_PASS':True})
 inputs=['DELEGATED_TASK.md','START_INSTRUCTION.txt','PUBLIC_RECEIVER_METADATA.json','PUBLIC_RECORD_INDEX.csv','PUBLIC_SAVED_CARRIER_ROWS.csv','INPUT_PROVENANCE.csv','NEW_TASK_AND_PARENT.json','SOURCE_AND_PERMISSION_LEDGER.md']
 dump(OUT/'DESIGN_FREEZE.json',{'parent':PARENT,'time_utc':stamp(),'task_sha256':sha(OUT/'DELEGATED_TASK.md'),'code_sha256':{p.name:sha(p) for p in sources},'input_sha256':{n:sha(OUT/n) for n in inputs},'conditions':'128 saved STABLE records only; main 20dB, pressure 5dB; 8 geometries x 8 paired repeats','valid_carrier':'finite LO+offset and original detection_ratio >=6; old residual acceptance separate','feature':'signed adjacent difference Hz; -c_ref*d/mean six valid receive frequencies, no F0','scopes':{'TURN':[2],'SAME_LEG':[0,1,3,4],'ALL':[0,1,2,3,4]},'recompute_tolerance_hz':1e-8,'digital_reference_delta_hz':0.031,'digital_reference_records':['R0000','R0001'],'mapping_rule':{'supported':'at 20dB at least 5/8 geometries have >=2/3 lines beating zero RMSE in BOTH SAME_LEG and TURN','not_supported':'at 20dB zero geometries satisfy SAME_LEG >=2/3 lines beating zero; otherwise UNRESOLVED','meaning':'descriptive route rule, not project precision PASS or coverage guarantee'},'tracking_rule':'reproducible if full rows, same settings recompute <=1e-8 Hz and digital reference PASS; no independent modal peak identity certificate','truth_timing':'separate evaluator only after FEATURE_OUTPUT_MANIFEST hashes frozen','statistics':'all planned entries retained, missing INF; conditional and unconditional statistics; no independent line/difference or paired-SNR sample claim','figures':['first-ID20dB carrier and old residual trajectory','all geometry line errors vs true-change/zero predictor, turn separate'],'budget':{'wall_total_s':1800,'RSS_bytes':4294967296,'BLAS_threads':4,'new_repo_bytes':67108864,'compute_processes':1},'new_KRAKEN_FIELD':0,'new_signals_noise_states':0,'stopping':'once only; no peak tuning, no DRIFT, no H2 feedback; Commit B and STOP'})
 print('PREFLIGHT AND DESIGN FREEZE PASS')
if __name__=='__main__':main()
