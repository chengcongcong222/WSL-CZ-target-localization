"""Serialization-only completion from saved tables; no field, Jacobian or trial evaluation."""
from pathlib import Path
import time
import r4_e2_nine_pair_pilot as p
def main():
 O=p.O
 if (O/'PILOT_DECISION.json').exists():raise RuntimeError('do not overwrite existing decision')
 info=p.base.rows(O/'LOCAL_INFORMATION_BY_SCENE.csv');st=p.base.rows(O/'NUMERICAL_INFORMATION_STABILITY.csv');cs=p.base.rows(O/'CANDIDATE_SEPARATION.csv');ct=p.base.rows(O/'CONTROL_CHECKS.csv')
 assert len(info)==4608 and len(st)==1152 and len(cs)==768
 def selected(r,cal):
  return r['representation']=='P2' and r['noise_model']=='RELATIVE' and r['calibration']==cal and r['bearing_profiled']=='True'
 primary=[r for r in info if selected(r,'C1') and int(r['mesh'])==160001 and float(r['sigma'])==.01]
 upper=[r for r in info if selected(r,'C0') and int(r['mesh'])==160001 and float(r['sigma'])==.01]
 ss=[r for r in st if selected(r,'C1')];su=[r for r in st if selected(r,'C0')]
 strong=lambda rs:all(float(r['range250_signal_sqrt_information'])>=1 and float(r['depth10_signal_sqrt_information'])>=1 for r in rs)
 c=[r for r in cs if int(r['mesh'])==160001 and r['noise_model']=='RELATIVE']
 candidateok=all(r['numerical_interval_positive']=='True' and float(r['C1_CA_Frobenius_rms'])>float(r['CA_marginal_noise1_rms']) and float(r['C1_joint_chart_scale_1pct'])>1 for r in c)
 localok=all(r['status']=='NUMERICALLY_STABLE' for r in ss)
 if any(r['PASS']!='True' for r in ct):decision='IMPLEMENTATION_INVALID'
 elif localok and strong(primary) and candidateok:decision='NINE_PAIR_NONBEARING_MECHANISM_PROMISING'
 elif localok and all(r['status']=='NUMERICALLY_STABLE' for r in su) and strong(upper) and not strong(primary):decision='NINE_PAIR_CALIBRATION_CONDITIONAL'
 else:decision='NINE_PAIR_INCREMENT_WEAK_OR_NOT_RESOLVED'
 start=p.read(O/'EXECUTION_STARTED.json')
 elapsed=max((O/n).stat().st_mtime for n in ['LOCAL_INFORMATION_BY_SCENE.csv','NUMERICAL_INFORMATION_STABILITY.csv','CONTROL_CHECKS.csv','CANDIDATE_SEPARATION.csv'])-start['epoch']
 p.dump(O/'OUTPUT_SERIALIZATION_RECOVERY.json',dict(reason='Native NumPy int32 in strong_C1_scenes count rejected by Python JSON encoder; reproduced from frozen metrics and stored PROFILED_H01_M, with no new field evaluation.',repair='Read saved tables; native Python float/bool/int at JSON boundary; identical frozen outcome predicates. Original scientific script/hash unchanged.',scientific_executions=1,scientific_reruns=0,new_forward_states=0,new_KRAKEN=0,frozen_scientific_code_preserved=True,original_process_exit_code=1))
 result=dict(decision=decision,parent_SHA=p.PARENT,design_SHA=start['design_SHA'],local_primary_stable_scenes=sum(r['status']=='NUMERICALLY_STABLE' for r in ss),strong_C1_scenes=sum(float(r['range250_signal_sqrt_information'])>=1 and float(r['depth10_signal_sqrt_information'])>=1 for r in primary),C1_range_information_min=min(float(r['range_information_depth_profiled']) for r in primary),C1_depth_information_min=min(float(r['depth_information_horizontal_profiled']) for r in primary),candidate_primary_resolved=sum(r['numerical_interval_positive']=='True' and float(r['C1_CA_Frobenius_rms'])>float(r['CA_marginal_noise1_rms']) for r in c),candidate_primary_count=len(c),control_checks=len(ct),control_failures=sum(r['PASS']!='True' for r in ct),original_E2_admission='FAIL_UNCHANGED',E2_G0_information='NOT_EVALUATED',pilot_status='EXPLORATORY_PILOT_ONLY',R4_percent=0,new_KRAKEN=0,new_MC=0,new_audio=0,next='FULL_BAND_AND_EXTRACTION_RESEARCH_REVIEW' if decision=='NINE_PAIR_NONBEARING_MECHANISM_PROMISING' else 'H3_REVIEW' if decision!='IMPLEMENTATION_INVALID' else 'STOP',next_stage_execution='NOT_AUTHORIZED',execution_elapsed_s_from_saved_artifact_timestamps=elapsed,independent_audit='PENDING',output_serialization_recovered=True)
 p.dump(O/'PILOT_DECISION.json',result);print(result)
if __name__=='__main__':main()
