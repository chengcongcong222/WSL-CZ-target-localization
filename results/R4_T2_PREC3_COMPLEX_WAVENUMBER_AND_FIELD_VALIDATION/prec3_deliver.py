"""Delivery only: aggregate sealed outputs; no new physical computations."""
from pathlib import Path
import json,csv,hashlib,shutil,io,sys,collections
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def main():
 for name in ['MAIN_SEAL.json','FINAL_SCIENCE_SEAL.json']:
  for r in read(OUT/name)['files']:assert sha(r['path'])==r['sha256'],r['path']
 cfg=read(OUT/'PREC3_DESIGN_FREEZE.json');a=rows(OUT/'PREC3_REAL_K_INDEPENDENT_VALIDATION.csv');b=rows(OUT/'PREC3_COMPLEX_K_VALIDATION.csv');cold=read(OUT/'PREC3_COLD_AUDIT.json');cc=rows(OUT/'PREC3_COLD_CHECKS.csv');ivp=read(OUT/'PREC3_COLD_PRUEFER.json');budget=read(OUT/'EXECUTION_BUDGET.json');dec=read(OUT/'PREC3_DECISION.json')
 assert len(a)==15 and len(b)==15 and cold['checks']==487 and cold['FAIL']==4
 assert all(r['check']=='independent_complex_mode_shape' for r in cc if r['PASS']=='False')
 assert budget['main_calls']==1 and budget['cold_calls']==1 and budget['retries']==0
 for r in cfg['frozen_code']:assert sha(OUT/r['path'])==r['sha256']
 summary={'A_real_roots':{'completed':15,'preregistered':15,'PASS':15,'max_Galerkin_vs_PREC2':max(float(r['GALERKIN_vs_PREC2']) for r in a),'max_Galerkin_R4_gap':max(float(r['GALERKIN_R4_gap']) for r in a),'max_discrete_mass_dual_residual_k_bound':max(float(r['max_dual_discrete_residual_k_bound']) for r in a)},'B_complex_root_comparisons':{'completed':15,'preregistered':15,'individual_root_flags_PASS':15},'main_independent_complex_function_modes':{'completed':2,'preregistered':3,'third_cancelled_at_second_mode_gate':True},'fresh_cold_Pruefer':{'completed':3,'preregistered':3,'max_difference_vs_Galerkin':max(r['vs_Galerkin'] for r in ivp),'max_two_level_difference':max(r['level_difference'] for r in ivp)},'cold_checks':{'completed':487,'PASS':483,'FAIL':4,'false_categories':['independent_complex_mode_shape'],'source_identity454_PASS':True},'C_fixed_sentinels':{'completed':0,'preregistered':60,'status':'NOT_EVALUATED_GATE_B_FAIL'},'complete_C3':'NOT_EVALUATED_NOT_AVAILABLE'}
 for k in ['max_real_path_disagreement','max_imag_path_disagreement','max_real_shift_grid_difference','max_imag_shift_grid_difference','MOD_imag_vs_independent']:summary['B_complex_root_comparisons'][k]=max(float(r[k]) for r in b)
 status={'stage':'T2-PREC3','classifications':['PREC3_A_REAL_K_SUPPORTED_SUBSET_ONLY','PREC3_B_COMPLEX_K_NOT_CERTIFIED'],'Gate_A':'PASS_SUBSET_CONDITIONAL','Gate_B':'FAIL_FROZEN_MODAL_FUNCTION_COMPARISON','Gate_C':'NOT_EVALUATED','failure_scope':'201Hz mode481 main paired FE-vs-FD normalized modal functions; cold confirms modes481+961, not complex-wavenumber disagreement','counts_and_maxima':summary,'cold_returncode1_reason':'intentional COLD_FAIL_NO_RETRY after all4scientific failures saved; no crash preceding results; no scientific retry','public_source_vs_binary_build':'NOT_ESTABLISHED','CHIGH_completeness_real_bottom_full_lossy_field':'NOT_CERTIFIED','T2_A_B_actual_performance':'NOT_EVALUATED','historical_G0_N1_N2_N3_PREC1_PREC2':'UNCHANGED','R4_percent':0,'research_lead_audit':'PENDING','STOP':True,'next_stage_authorized':False}
 write(OUT/'PREC3_FINAL_ADJUDICATION.json',status)
 # Main raw decision is immutable; explicit supplementary interpretation avoids relabeling.
 with (OUT/'PREC3_FIELD_SHAPE_COMPARISON.csv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['status','completed','preregistered','reason','field_error']);w.writeheader();w.writerow({'status':'NOT_EVALUATED','completed':0,'preregistered':60,'reason':'Gate_B_FUNCTION_INDEPENDENCE_FAIL','field_error':''})
 write(OUT/'NOT_EVALUATED_BRANCHES.json',{'main_B_independent_FD_mode961':'CANCELLED_AFTER_MODE481; separately examined only by registered cold audit','C_all60sentinels':'NOT_EVALUATED','C3':'NOT_EVALUATED','C_full_field_accuracy':'NOT_EVALUATED','external_solver_calls':0,'new_waveform_MC_DFT_FIM_localization_calls':0})
 logs=OUT/'runtime_logs';logs.mkdir(exist_ok=True)
 for p in (LOCAL/'logs').glob('*'):
  if p.is_file():shutil.copyfile(p,logs/p.name)
 large=[]
 for p in (LOCAL/'vectors').glob('*.npz'):large.append({'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
 write(OUT/'PREC3_LARGE_DATA_INDEX.json',{'scope':'actual vectors for completed calculations, retained on D; unperformed field branch absent','files':large})
 pre=read(OUT/'PREC3_PREFLIGHT.json')
 write(OUT/'PREC3_PREFLIGHT_SUPPLEMENT.json',{'saved_preflight_checks':len(pre['checks']),'saved_FAIL':0,'additional_fixed15_check':True,'total_including_preparation_final_check':len(pre['checks'])+1,'scientific_calls_before_freeze':0,'initial_operational_status_stop_preserved':True})
 # Frozen failed Gate is not changed into a success by the positive k-only evidence.
 write(OUT/'PREC3_SCOPE_AND_FAILURES.json',{'main_decision_byte_unchanged':True,'source_and_extrapolation_and_scalar_k_cold_checks_PASS':True,'modal_function_comparison':'relative L2 over original251SSP nodes after bilinear normalization, minimum sign alignment; not array field error','function_reference':'independent freshly computed FE/FD fields, NOT originalMOD function accuracy','wave_number_vs_function':'extrapolated eigenvalues can agree while raw-grid eigenfunctions differ; no certified remaining full-field error','two_level_finer_improvement_is_not_absolute_error_bound':True,'weak_loss_scaling':'errors below1e-11 floor cannot demonstrate actual quadratic asymptotics','density':'constant1 g/cm3','Thorp':'nominal model absorption, not measured ocean loss','boundary':'pressure release top/rigid bottom; real bottom unknown','source_outputs':'no all-mode sidecar produced; subset computations only'})
 # Append only one zero-credit ledger row, retaining the historical working byte prefix.
 ledger=ROOT/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv';original=(LOCAL/'WORKTREE_LEDGER_BEFORE.csv').read_bytes();assert ledger.read_bytes()==original
 buf=io.StringIO();csv.writer(buf,lineterminator='\n').writerow(['T2-PREC3','../R4_T2_PREC3_COMPLEX_WAVENUMBER_AND_FIELD_VALIDATION/PREC3_FINAL_ADJUDICATION.json','PREC3_B_COMPLEX_K_NOT_CERTIFIED','A15_OF15;B_K15_OF15;MODAL_FUNCTION_GATE_FAIL;C0_OF60;R4_delta=0','COLD487_CHECKS_4_FUNCTION_FAIL;HISTORICAL_UNCHANGED;RESEARCH_LEAD_PENDING;STOP',0]);ledger.write_bytes(original+buf.getvalue().encode('utf-8'))
 write(OUT/'LEDGER_APPEND_PROVENANCE.json',{'old_prefix_sha256':hashlib.sha256(original).hexdigest(),'old_prefix_preserved':True,'rows_added':1,'credited_weight_percent':0})
 write(OUT/'STOP.json',{'STOP':True,'R4_percent':0,'next_stage_authorized':False})
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
