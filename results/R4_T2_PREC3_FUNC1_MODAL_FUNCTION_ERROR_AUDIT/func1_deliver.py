"""Delivery aggregation from sealed FUNC1 results; no vector analysis rerun."""
from pathlib import Path
import json,csv,hashlib,collections,shutil,io,statistics
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3_FUNC1')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def rows(p):
 with Path(p).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def table(p,x):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(x[0]),lineterminator='\n');w.writeheader();w.writerows(x)
def main():
 for name in ['MAIN_OUTPUT_SEAL.json','SCIENCE_OUTPUT_SEAL.json']:
  for r in read(OUT/name)['files']:assert sha(r['path'])==r['sha256'],r['path']
 cfg=read(OUT/'FUNC1_DESIGN_FREEZE.json');dec=read(OUT/'FUNC1_DECISION.json');cold=read(OUT/'FUNC1_COLD_AUDIT.json');budget=read(OUT/'EXECUTION_BUDGET.json')
 assert cold['PASS'] and cold['FAIL']==0 and cold['checks']==12792
 assert dec['original_reconstructions']==4 and dec['coupling_vectors']==80 and budget['main_calls']==1 and budget['cold_calls']==1 and budget['science_retries']==0
 dense=rows(OUT/'FUNC1_DENSE_SAMPLING_COMPARISON.csv');phase=rows(OUT/'FUNC1_NORMALIZATION_AND_PHASE_AUDIT.csv');cross=rows(OUT/'FUNC1_FE_FD_CROSS_GRID.csv');bands=rows(OUT/'FUNC1_DEPTH_RESOLVED_ERROR.csv');coupling=rows(OUT/'FUNC1_LOCAL_MODAL_COUPLING.csv');case=[];stability=[]
 for mode in [1,481]:
  for N in [80000,160000]:
   d={r['grid']:float(r['e_sign']) for r in dense if r['mode']==str(mode) and r['N']==str(N)};p=next(r for r in phase if r['mode']==str(mode) and r['N']==str(N) and r['grid']=='S3');delta=abs(d['S2']-d['S3']);stability.append(delta)
   case.append({'mode':mode,'N':N,'old251_unweighted_e':d['SSP_UNWEIGHTED'],'old251_trapezoid_e':d['SSP_TRAPZ'],'dense025_trapezoid_e':d['S3'],'dense_vs_original_relative_increase':(d['S3']-d['SSP_UNWEIGHTED'])/d['SSP_UNWEIGHTED'],'dense_complex_projected_e':float(p['e_complex']),'dense_stability_abs_delta':delta,'historical0p002_failed_comparison_persists':d['SSP_UNWEIGHTED']>.002 and d['S3']>.002})
 table(OUT/'FUNC1_CASE_SUMMARY.csv',case)
 grouped=[]
 for mode in [1,481]:
  for N in [80000,160000]:
   for res in cfg['resources']:
    r=[v for v in coupling if v['mode']==str(mode) and v['N']==str(N) and v['resource']==res];raw=[float(v['raw_relative_difference']) for v in r];shape=[float(v['array_shape_difference']) for v in r]
    grouped.append({'mode':mode,'N':N,'resource':res,'probes':len(r),'raw_relative_min':min(raw),'raw_relative_median':statistics.median(raw),'raw_relative_max':max(raw),'shape_difference_max':max(shape),'weak_source_or_receiver_cases':sum(v['weak_source_or_receiver']=='True' for v in r),'weak_coupling_norm_cases':sum(v['weak_coupling_norm']=='True' for v in r),'interpretation':'DESCRIPTIVE_SINGLE_MODE_Q_ONLY'})
 table(OUT/'FUNC1_COUPLING_SUMMARY.csv',grouped)
 mode961=next(r for r in cross if r['comparison']=='SAME_METHOD_CROSS_GRID' and r['method']=='FE' and r['mode']=='961')
 final={'stage':'T2-PREC3-FUNC1','classifications':dec['classifications'],'scientific_scope':'ARCHIVED_FUNCTIONS_ONLY','original_reconstruction':'PASS4_OF4','dense_quadrature_stability':'PASS4_OF4','max_dense_stability_delta':max(stability),'cases':case,'source_vectors':10,'cold_checks':cold['checks'],'cold_FAIL':0,'mode961_FE_cross_grid_e':float(mode961['e_sign']),'mode961_FD_dense_and_coupling':'NOT_EVALUATED_FD_VECTOR_NOT_ARCHIVED','single_mode_coupling_vectors':80,'channels':800,'weak_mode1_coupling_vectors':sum(r['mode']=='1' and r['weak_coupling_norm']=='True' for r in coupling),'new_solvers':0,'R4_percent':0,'historical_PREC3_B':'FAIL_UNCHANGED','CHIGH_and_real_bottom':'NOT_CERTIFIED','T2_A_B_actual_performance':'NOT_EVALUATED','research_lead_audit':'PENDING','STOP':True,'next_stage_authorized':False}
 write(OUT/'FUNC1_FINAL_ADJUDICATION.json',final)
 write(OUT/'FUNC1_DATA_DICTIONARY.json',{'errors':'dimensionless; multiply100 for percent','sampling_bias':'20m original has unweighted norm; SSP_TRAPZ separately isolates endpoint weighting before dense sampling','metric_e_sign':'reference is second vector, FD in FE-vs-FD; raw archived values, sign selection only','metric_e_complex':'reference target norm is first vector; alpha minimizes weighted residual; diagnostic only','same_method_cross_grid_metric_names':'In generic metric fields FE_* means first/coarse vector and FD_* means second/fine reference, BOTH belong to method column. No FD961 vector is implied.','coupling':'phi(source)*phi(receiver), no radial propagation factor or full-mode summation; no per-channel fitted gain','equal_depth_arrays':'single-mode Q vectors at same receiver depth are collinear by construction; nearzero shape errors are a structural control, not full-field accuracy','weak_flags':'fixed design flags, no realistic SNR inferred; raw finite ratios retained, exactly zero reference undefined','own_method_order':'single cross-grid pair does not estimate individual-method convergence order','full_field_and_localization':'NOT_EVALUATED'})
 logs=OUT/'runtime_logs';logs.mkdir(exist_ok=True)
 for p in (LOCAL/'runtime_logs').glob('*'):
  if p.is_file():shutil.copyfile(p,logs/p.name)
 ledger=ROOT/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv';old=(LOCAL/'WORKTREE_LEDGER_BEFORE.csv').read_bytes();assert ledger.read_bytes()==old
 buf=io.StringIO();csv.writer(buf,lineterminator='\n').writerow(['T2-PREC3-FUNC1','../R4_T2_PREC3_FUNC1_MODAL_FUNCTION_ERROR_AUDIT/FUNC1_FINAL_ADJUDICATION.json','ARCHIVED_MODAL_FUNCTION_DIAGNOSTIC_COMPLETE','VECTORS10;REPLAY4;DENSE20;COUPLING80;NO_NEW_SOLVES;FD961_NOT_ARCHIVED;R4_delta=0','COLD12792_PASS_0FAIL;HISTORICAL_PREC3_B_FAIL_UNCHANGED;RESEARCH_LEAD_PENDING;STOP',0]);ledger.write_bytes(old+buf.getvalue().encode('utf-8'))
 write(OUT/'LEDGER_APPEND_PROVENANCE.json',{'old_sha256':hashlib.sha256(old).hexdigest(),'old_byte_prefix_preserved':True,'appended_rows':1,'credited_weight_percent':0})
 write(OUT/'STOP.json',{'STOP':True,'R4_percent':0,'next_stage_authorized':False})
 print(json.dumps(final,indent=2))
if __name__=='__main__':main()
