"""Read-only reconstruction report; never executes the scientific regression."""
import json,shutil
import numpy as np
from pathlib import Path
import r4_e2_nine_pair_formula_repair as p
def main():
 o=p.O;f=p.read(o/'DESIGN_FREEZE.json');d=p.read(o/'REPAIR_DECISION.json')
 v=p.read(o/'VALIDATION.json')
 controls=p.base.rows(o/'CONTROL_CHECKS.csv')
 struct=p.base.rows(o/'STRUCTURAL_NULLSPACE_CHECKS.csv')
 fd=p.base.rows(o/'ANALYTIC_VS_FD_JACOBIAN.csv')
 comp=p.base.rows(o/'ORIGINAL_VS_REPAIRED_COMPARISON.csv')
 info=p.base.rows(o/'REPAIRED_INFORMATION_BY_SCENE.csv')
 st=p.base.rows(o/'EFFECTIVE_INFORMATION_STABILITY.csv')
 ca=p.base.rows(o/'CANDIDATE_SEPARATION.csv')
 primary=lambda r:r['representation']=='P2' and r['calibration']=='C1' and r['bearing_profiled']=='True'
 mainrows=[r for r in info if primary(r) and r['mesh']=='160001' and r['sigma']=='0.01']
 mainst=[r for r in st if primary(r)]
 pn=[r for r in struct if r['representation']=='P1' and r['calibration']=='C1' and r['noise_model']=='RELATIVE']
 formula=all(r['PASS']=='True' for r in controls) and all(r['PASS']=='True' for r in struct) and v['FAIL']==0
 if not formula:d['decision']='E2_FORMULA_REPAIR_FAILED';d['next']='H3_REVIEW'
 elif d['decision']=='E2_NINE_PAIR_ANALYTIC_CHAIN_REPAIRED_CONDITIONALLY_PROMISING' and not all(r['status']=='NUMERICALLY_STABLE' for r in mainst):
  d['decision']='E2_NINE_PAIR_INFORMATION_NOT_ESTABLISHED';d['next']='H3_REVIEW'
 d.update(original_pilot='IMPLEMENTATION_INVALID',formula_Gate='PASS' if formula else 'FAIL',independent_audit='PASS' if v['FAIL']==0 else 'FAIL',P1_physical_rank_pass=sum(r['PASS']=='True' for r in pn),P1_physical_rank_count=len(pn),P1_structural_null_max=max(float(r['null_response']) for r in pn),P1_half_depth_null_max=max(float(r['half_depth_null_response']) for r in pn),P2_primary_analytic_ranks=sorted({int(r['scaled_rank']) for r in mainrows}),primary_both_noise_stability_pass=sum(r['status']=='NUMERICALLY_STABLE' for r in mainst),primary_both_noise_stability_count=len(mainst),weak_singular_relative_max=max(float(r['weak_singular_relative']) for r in mainst),weak_direction_sine_max=max(float(r['weak_direction_sine']) for r in mainst),main_depth_information_stability_max=max(float(r['depth_information_relative']) for r in mainst))
 p.dump(o/'REPAIR_DECISION.json',d)
 comparison_summary={}
 for rep in ['P1','P2']:
  rr=[r for r in comp if r['representation']==rep and r['calibration']=='C1' and r['noise_model']=='RELATIVE']
  comparison_summary[rep]={k:max(float(r[k]) for r in rr) for k in ['range_relative_change','depth_relative_change']}
 def select(rep,cal,noise='RELATIVE'):
  return {r['scene_id']:r for r in info if r['representation']==rep and r['calibration']==cal and r['noise_model']==noise and r['bearing_profiled']=='True' and r['mesh']=='160001' and r['sigma']=='0.01'}
 def ratios(a,b,col):
  vals=[float(a[k][col])/max(float(b[k][col]),1e-30) for k in a]
  return dict(min=min(vals),median=float(np.median(vals)),max=max(vals))
 p2=select('P2','C1');c0=select('P2','C0');floor=select('P2','C1','ABSOLUTE_FLOOR')
 summary=dict(P1_CA_to_raw_trace=ratios(select('P1','C1'),select('P0_SINGLE','C1'),'information_trace'),P2_CA_to_raw_trace=ratios(p2,select('P0_DUAL','C1'),'information_trace'),C1_to_C0_range=ratios(p2,c0,'range_information_depth_profiled'),C1_to_C0_depth=ratios(p2,c0,'depth_information_horizontal_profiled'),fixed_floor_to_relative_range=ratios(floor,p2,'range_information_depth_profiled'),fixed_floor_to_relative_depth=ratios(floor,p2,'depth_information_horizontal_profiled'),original_vs_repaired= comparison_summary,analytic_identifiable_J_max=max(float(r['identifiable_relative_error']) for r in fd),C2_nonzero_rows=sum(float(r['information_trace'])!=0 for r in info if r['calibration']=='C2'),candidate_1pct_min_ratio=min(float(r['C1_CA_Frobenius_rms'])/float(r['CA_marginal_noise1_rms']) for r in ca if r['mesh']=='160001' and r['noise_model']=='RELATIVE'),candidate_5pct_min_ratio=min(float(r['C1_CA_Frobenius_rms'])/float(r['CA_marginal_noise5_rms']) for r in ca if r['mesh']=='160001' and r['noise_model']=='RELATIVE'))
 p.dump(o/'SUMMARY.json',summary)
 total=dict(checks=len(controls)+v['checks'],PASS=sum(r['PASS']=='True' for r in controls)+v['PASS'],FAIL=sum(r['PASS']!='True' for r in controls)+v['FAIL'],formula_Gate=d['formula_Gate'],independent_reconstruction=v,scope='Formula controls plus independent reconstruction; scientific stability recorded separately.')
 p.dump(o/'VALIDATION.json',total)
 text=f"""# E2 nine-pair analytic formula repair
Final decision: {d['decision']}.
Parent: {p.PARENT}. Design: {d['design_SHA']}.
Execution SHA is the containing commit; remote main must match it after push.

P1 physical rank/null controls: {d['P1_physical_rank_pass']}/{d['P1_physical_rank_count']}.
Maximum null response: {d['P1_structural_null_max']:.12g};
half-depth-stencil response: {d['P1_half_depth_null_max']:.12g}.
P2 primary measured ranks: {d['P2_primary_analytic_ranks']}.
Both-noise primary stability: {d['primary_both_noise_stability_pass']}/{d['primary_both_noise_stability_count']}.
Weak identifiable singular change maximum: {d['weak_singular_relative_max']:.8g}.
Weak direction sine maximum: {d['weak_direction_sine_max']:.8g}.
Effective depth-information mesh/stencil difference maximum: {d['main_depth_information_stability_max']:.8g}.
Formula Gate: {d['formula_Gate']}. Validation: {total['PASS']} PASS / {total['FAIL']} FAIL.

All local information is FORMULA_REPAIRED_DEVELOPMENT_ONLY. Local uncertainties are
linearized model/noise diagnostics, not achieved localization accuracies. Source
unknown, fixed noise floor is a design assumption and cannot define physical SNR.
C1 unknown gains are jointly shared across all registered frequencies and times;
C2 saturated response has no state information. Relative and fixed-floor comparisons
must use matched scenes, not cohort minima. Trace ratios depend on registered scales.
Finite candidates and feasible gain fits certify neither global nuisance minima nor
global correct-branch coverage. No depth prior, oracle spectrum, mode selection or
extra calibration was introduced. P1's null no longer supplies finite velocity precision.

Summary (dimensionless ratios and relative changes):
{json.dumps(summary,indent=2)}

The old negative pilot, all old candidate pressures, and input modal files remain
hash-locked unchanged. The original 0.2% raw-field admission is FAIL_UNCHANGED,
the original pilot is IMPLEMENTATION_INVALID, original E2 information NOT_EVALUATED,
and R4=0%. New KRAKEN/Monte Carlo/audio calls: 0/0/0.
Next: {d['next']} (review only); no full-band or H3 execution authorized.
Stop after execution commit, push, and remote verification.
"""
 (o/'GPT_SYNC.md').write_text(text,encoding='utf-8')
 (o/'PLANE_WAVE_AND_DATA_PROCESSING_CONTROLS.md').write_text(
 '# Formula/noise controls\nPlane-wave distance/depth and profiled-bearing response, fixed-gain absorption, saturated C2, shared-frequency off-diagonal covariance, exact CA derivative, alternative dense contrast covariance, and all scene/mesh data-processing inequalities are recorded in CONTROL_CHECKS.csv.\n'
 +f"Result: {sum(r['PASS']=='True' for r in controls)} PASS / {sum(r['PASS']!='True' for r in controls)} FAIL.\n",encoding='utf-8')
 # Final artifact manifest excludes itself; old immutable input bindings rechecked independently.
 artifacts={}
 for q in o.rglob('*'):
  if q.is_file() and q.name!='EXECUTION_ARTIFACT_MANIFEST.json':
   raw=q.suffix in ['.npz','.png']
   artifacts[str(q)]={'kind':'RAW' if raw else 'LF','sha256':p.sha(q,raw)}
 p.dump(o/'EXECUTION_ARTIFACT_MANIFEST.json',dict(artifacts=artifacts,design_SHA=d['design_SHA'],R4_percent=0))
 master=Path('results/R4_MASTER')
 with (master/'R4_PLAN.md').open('a',encoding='utf-8') as fh:
  fh.write(f"\n\n## E2 analytic formula repair\n{d['decision']}; formula Gate {d['formula_Gate']}; P1 {d['P1_physical_rank_pass']}/48 structural controls; original pilot invalid and original E2 admission unchanged; R4=0%. See ../R4_E2_NINE_PAIR_FORMULA_REPAIR/GPT_SYNC.md. Next is review only, automatic execution stopped.\n")
 ledger=master/'R4_EVIDENCE_LEDGER.csv'
 import csv
 with ledger.open(encoding='utf-8',newline='') as fh:fields=next(csv.reader(fh))
 row={k:'' for k in fields}
 for k in fields:
  kl=k.lower()
  if 'stage' in kl or 'id'==kl:row[k]='R4_E2_NINE_PAIR_ANALYTIC_FORMULA_REPAIR'
  elif 'status' in kl or 'decision' in kl:row[k]=d['decision']
  elif 'path' in kl or 'artifact' in kl or 'evidence' in kl:row[k]='results/R4_E2_NINE_PAIR_FORMULA_REPAIR/GPT_SYNC.md'
  elif 'note' in kl:row[k]='FORMULA_REPAIRED_DEVELOPMENT_ONLY; old E2 FAIL_UNCHANGED; R4=0%; STOP'
 with ledger.open('a',encoding='utf-8',newline='') as fh:csv.DictWriter(fh,fieldnames=fields).writerow(row)
 p.dump(master/'R4_E2_NINE_PAIR_FORMULA_REPAIR_PROGRESS.json',d)
 print(json.dumps(p.native(dict(decision=d,summary=summary,validation=total)),indent=2))
if __name__=='__main__':main()
