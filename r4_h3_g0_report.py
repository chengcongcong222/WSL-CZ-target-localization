"""Package saved mechanism and read-only audit evidence; no science rerun."""
import json,csv
from pathlib import Path
import r4_h3_g0 as p
def main():
 o=p.O;d=p.read(o/'H3_G0_DECISION.json');v=p.read(o/'VALIDATION.json')
 info=p.rows(o/'VERTICAL_INFORMATION_BY_SCENE.csv') if (o/'VERTICAL_INFORMATION_BY_SCENE.csv').exists() else []
 st=p.rows(o/'NUMERICAL_STABILITY.csv') if (o/'NUMERICAL_STABILITY.csv').exists() else []
 rr=p.rows(o/'RESOURCE_MATCHED_CONTROL.csv') if (o/'RESOURCE_MATCHED_CONTROL.csv').exists() else []
 mainrows=[r for r in info if r['resource']=='B1' and r['calibration']=='C1b' and r['mesh']=='160001' and float(r['sigma'])==.01]
 summary=[]
 for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
  for res in ['B0','B1','B2']:
   for cal in ['C0','C1a','C1b','C2','B3_KNOWN_SOURCE_AND_CALIBRATION']:
    a=[r for r in info if r['resource']==res and r['calibration']==cal and r['noise_model']==noise and r['mesh']=='160001' and float(r['sigma'])==.01]
    if a:summary.append(dict(noise_model=noise,resource=res,calibration=cal,depth_information_min=min(float(r['depth_information_horizontal_profiled']) for r in a),depth_information_max=max(float(r['depth_information_horizontal_profiled']) for r in a),oracle_horizontal_depth_min=min(float(r['depth_information_oracle_horizontal']) for r in a),rank_set=','.join(sorted(set(r['rank'] for r in a))),worst_local_scale_m=max(float(r['depth_local_scale_m']) if r['depth_local_scale_m']!='INF' else float('inf') for r in a)))
 if summary:p.write(o/'CALIBRATION_SUMMARY.csv',summary)
 if info:
  contrast=[]
  lookup={(r['scene_id'],r['mesh'],r['resource'],r['noise_model'],r['sigma'],r['calibration']):r for r in info}
  for r in mainrows:
   oracle=lookup[r['scene_id'],r['mesh'],r['resource'],r['noise_model'],r['sigma'],'B3_KNOWN_SOURCE_AND_CALIBRATION']
   c0=lookup[r['scene_id'],r['mesh'],r['resource'],r['noise_model'],r['sigma'],'C0']
   contrast.append(dict(scene_id=r['scene_id'],noise_model=r['noise_model'],source_loss_depth_ratio=float(c0['depth_information_horizontal_profiled'])/max(float(oracle['depth_information_horizontal_profiled']),1e-30),C1b_to_C0_depth_ratio=float(r['depth_information_horizontal_profiled'])/max(float(c0['depth_information_horizontal_profiled']),1e-30),C1b_to_C0_trace_ratio=float(r['information_trace'])/max(float(c0['information_trace']),1e-30),trace_scope='SCALED_STATE_DEPENDENT_NOT_PERFORMANCE'))
  p.write(o/'UNKNOWN_GAIN_LOSS.csv',contrast)
 low=[]
 import numpy as np
 for sc in p.read(o/'DESIGN_FREEZE.json')['scenes']:
  for res in ['B0','B1','B2']:
   a=np.load(o/f"INPUT_{sc['scene_id']}_{res}_n160001.npz")['pressure']
   low.append(dict(scene_id=sc['scene_id'],resource=res,min_amplitude=float(abs(a).min()),median_amplitude=float(np.median(abs(a))),max_floor1_sigma_over_nominal=float(.01*p.REF/abs(a).min()),absolute_noise_reference=p.REF,actual_SNR='NOT_ESTABLISHED',delta_method='LOCAL_HIGH_SNR_DESIGN_APPROXIMATION_ONLY'))
 p.write(o/'LOW_ENERGY_RECEIVING_DIAGNOSTIC.csv',low)
 d.update(high_level_classification='H3_G0_NUMERICAL_OR_FORMULA_INCOMPLETE' if d['decision'] in ['H3_G0_IMPLEMENTATION_INVALID','H3_G0_VERTICAL_NUMERICAL_FIDELITY_INCOMPLETE'] else d['decision'],original_E2_G0='FAIL_UNCHANGED',original_E2_information='NOT_EVALUATED',H3_EXTRACTED_OBSERVABLE='NOT_OPENED',R4_percent=0,independent_main_geometries=3,retained_mirror_identifiers=6,finite_support_not_full_RC2=True)
 p.dump(o/'H3_G0_DECISION.json',d)
 lines=['# H3-G0 finite vertical mechanism result','',f"Final scoped decision: {d['decision']}.",f"High-level class: {d['high_level_classification']}.",f"Design SHA: {d.get('design_SHA',p.read(o/'EXECUTION_STARTED.json')['design_SHA'])}. Execution SHA is commit B; verified remote SHA is reported after commit.",'',f"Independent saved-evidence reconstruction: {v['PASS']} PASS / {v['FAIL']} FAIL. Maximum cold pressure relative difference {v['cold_pressure_relative_max']:.9g}; Jacobian {v['independent_J_relative_max']:.9g}; effective depth QR/SVD {v['independent_depth_information_relative_max']:.9g}.",'',f"Registered input/gradient arrays: 36. Information rows: {len(info)}. Finite F label diagnostics: "+str(len(p.rows(o/'FINITE_HORIZONTAL_LABEL_DIAGNOSTIC.csv')) if (o/'FINITE_HORIZONTAL_LABEL_DIAGNOSTIC.csv').exists() else 0)+'.','',
 '## Calibration and horizontal profiling','',
 '|Noise|Resource|Calibration|Profiled Iz min/max (m^-2)|Oracle-horizontal Iz min|Ranks|Worst local scale (m)|',
 '|---|---|---|---|---|---|---|']
 for r in summary:lines.append(f"|{r['noise_model']}|{r['resource']}|{r['calibration']}|{r['depth_information_min']:.7g} / {r['depth_information_max']:.7g}|{r['oracle_horizontal_depth_min']:.7g}|{r['rank_set']}|{r['worst_local_scale_m']:.7g}|")
 lines+=['','Values above are deterministic nominal local information diagnostics; if the numerical Gate fails, they must not be interpreted as established depth information or estimator precision. Any INF depth scale denotes a null/unidentifiable direction; it is never reported as zero variance. Full singular values, weakest direction and nullspace loading are retained per row.','',
 'C1b unknown per-element/per-frequency gains are fixed across all three snapshots, and source is arbitrary at each frequency/snapshot. C1a fixes per-element gain across every frequency. C2 spans the full response and gives zero information. Arbitrary free modal gains analytically absorb depth: all modal source-200 coefficients are nonzero in the registered caches. Constrained group gain is NOT_ADMITTED: no physical constraint or group count has been invented. C1e environment/pose uncertainty is NOT_EVALUATED.','',
 '## Vertical gain and numerical closure','',
 '|Scene|Noise|B1/B2 profiled depth information|B1 10m signal diagnostic|',
 '|---|---|---|---|']
 for r in rr:lines.append(f"|{r['scene_id']}|{r['noise_model']}|{float(r['B1_to_B2']):.7g}|{float(r['B1_depth10_signal']):.7g}|")
 if st:
  m=[r for r in st if r['resource']=='B1' and r['calibration']=='C1b']
  lines+=['',f"Primary numerical rows STABLE: {sum(r['status']=='STABLE' for r in m)}/{len(m)}.",
   '|Diagnostic|Worst primary value|Frozen limit|','|---|---|---|']
  for key,limit in [('response_to_noise1_ratio',.5),('Jacobian_mesh_relative',.02),('depth_step_relative',.02),('depth_information_relative',.1)]:
   lines.append(f"|{key}|{max(float(r[key]) for r in m):.9g}|{limit}|")
  lines+=['',f"Other calibration/resource rows numerically stable: {sum(r['status']=='STABLE' for r in st)}/{len(st)}. Raw two-grid field differences are retained separately; new H3 relative-response tests do not change old E2's failed 0.2% raw-field Gate. All finite labels use exact cached source/receiver entries; local depth derivatives at 200m are step comparisons, not continuous-depth validation."]
 lines+=['','## Evidence and limits','',
 'Shared-reference covariance and direct source-nuisance information are cross-checked with dense covariance; data-processing inequality is checked against original known-source raw fields. All controls are in CONTROL_CHECKS.csv and read-only reconstruction in INDEPENDENT_AUDIT_CHECKS.csv.','',
 'F retains every 25 horizontal error nodes × five depth labels at 160001. These are nominal tangent nuisance projections of finite log/phase contrasts with source and C1b gain jointly removed across time. They are not nonlinear nuisance optimization, estimator recovery, full RC2 support, branch certification or a continuous-depth guarantee. U/U0: H3_FULL_HORIZONTAL_SUPPORT_NOT_ESTABLISHED. Six labels represent three identical-acoustics mirror pairs; no evidence multiplier is applied.','',
 'Noise is a proper-complex ideal-field delta-method design at 1%/5%, with separate frozen absolute floor. LOW_ENERGY_RECEIVING_DIAGNOSTIC.csv records weak field regions and checks whether local high-SNR approximation may be stressed. Unknown source level forbids interpreting these assumptions as actual SNR. No robustness conclusion to SSP/pose or signal extraction is drawn.','',
 f"Proposed next REVIEW only: {d.get('next','STOP')}. No next execution authorized. No new KRAKEN/FIELD/MC/audio; no full-band E2. H3_EXTRACTED_OBSERVABLE NOT_OPENED. R4=0%. Commit B/push/verify and STOP."]
 (o/'GPT_SYNC.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 for name in ['VERTICAL_INFORMATION_BY_SCENE.csv','RESOURCE_MATCHED_CONTROL.csv','NOISE_AND_CALIBRATION_SENSITIVITY.csv','NUMERICAL_STABILITY.csv','FINITE_HORIZONTAL_LABEL_DIAGNOSTIC.csv']:
  if not (o/name).exists():p.write(o/name,[dict(status='NOT_EXECUTED_AFTER_REGISTERED_PRE_FORMULA_STOP')])
 master=Path('results/R4_MASTER')
 with (master/'R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\n## H3-G0 finite vertical response closed\n\n'+d['decision']+'. Six retained labels / three distinct main geometries; exact cached depths; C1b primary; free modal-gain depth null control; no justified constrained group model. Numerical and source/calibration scope retained. Full RC2 support not established; extracted observable not opened. R4=0%. Next '+d.get('next','STOP')+' is a review proposal only. STOP after result commit.\n')
 with (master/'R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8',newline='') as f:
  csv.writer(f).writerow(['R4_H3_G0_FINITE_VERTICAL_RESPONSE','../R4_H3_G0_FINITE_VERTICAL_RESPONSE/H3_G0_DECISION.json',d['decision'],'FINITE_DETERMINISTIC_FIXED_ENVIRONMENT_UPPER_BOUND;NO_FULL_RC2_SUPPORT','PENDING_RESEARCH_LEAD_AUDIT;STOP',0])
 p.dump(master/'R4_H3_G0_FINITE_VERTICAL_RESPONSE_PROGRESS.json',dict(stage='R4_H3_G0_FINITE_VERTICAL_RESPONSE_MECHANISM',decision=d['decision'],design_SHA=d.get('design_SHA'),R4_percent=0,execution_status='STOPPED_AFTER_SINGLE_REGISTERED_EXECUTION',next_review=d.get('next','STOP'),next_execution='NOT_AUTHORIZED'))
 files=[dict(file=str(x),sha256=p.sha(x,x.suffix=='.npz'),hash_kind='RAW' if x.suffix=='.npz' else 'LF_TEXT') for x in sorted(o.iterdir()) if x.is_file() and x.name!='ARTIFACT_MANIFEST.json']
 p.dump(o/'ARTIFACT_MANIFEST.json',dict(files=files,scope='FINAL_ARTIFACTS_BEFORE_COMMIT_B',R4_percent=0))
 print(json.dumps(dict(decision=d['decision'],artifacts=len(files),validation=v,summary=summary),indent=2))
if __name__=='__main__':main()
