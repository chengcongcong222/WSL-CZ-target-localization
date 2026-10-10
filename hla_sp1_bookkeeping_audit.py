"""Supplementary independent interval/denominator/statistic bookkeeping."""
from hla_sp1_core import *
import math

def main():
 rec=rows(OUT/'RANGE_RESULTS_BY_RECORD.csv');checks=[]
 def check(n,x):checks.append(dict(check=n,PASS=bool(x)))
 for r in rec:
  step=float(r['grid_step_km']);labels=json.loads(r['range_labels']);actual=json.loads(r['range_intervals_km']);chunks=[]
  for x in labels:
   lo=max(45,x-step/2);hi=min(60,x+step/2)
   if chunks and x-previous<step*1.01:chunks[-1][1]=hi
   else:chunks.append([lo,hi])
   previous=x
  check('interval_cells:'+str((r['record_id'],r['method'],step,r['modal_mesh'])),len(chunks)==len(actual) and all(abs(x-y)<1e-10 for a,b in zip(chunks,actual) for x,y in zip(a,b)))
  valid=r['effective_output']=='True';truth=float(r['truth_r_km']);expected=abs(float(r['best_r_km'])-truth) if valid else math.inf;worst=max(abs(x-truth) for pair in chunks for x in pair) if valid else math.inf
  check('unconditional_errors:'+str((r['record_id'],r['method'],step,r['modal_mesh'])),expected==float(r['top1_abs_distance_error_km']) and (abs(worst-float(r['full_cell_worst_error_km']))<1e-10 if valid else math.isinf(float(r['full_cell_worst_error_km']))))
 for row in rows(OUT/'RANGE_RESULTS_BY_GEOMETRY.csv'):
  group=[r for r in rec if all(r[k]==row[k] for k in ['method','C','reference_SNR_dB','grid_step_km','modal_mesh','geometry'])];widths=sorted(float(r['envelope_width_km']) if r['effective_output']=='True' else math.inf for r in group);err=sorted(float(r['top1_abs_distance_error_km']) for r in group);med=(widths[7]+widths[8])/2;em=(err[7]+err[8])/2
  check('geometry_stats:'+str(row),len(group)==16 and med==float(row['median_envelope_km']) and em==float(row['median_point_error_km']) and sum(r['effective_output']=='True' for r in group)==int(row['effective_outputs']) and sum(r['horizontal_label_retained']=='True' for r in group)/16==float(row['horizontal_label_retention']))
 table(OUT/'BOOKKEEPING_AUDIT_CHECKS.csv',checks);v=read(OUT/'VALIDATION.json');v['supplemental_bookkeeping_PASS_count']=sum(r['PASS'] for r in checks);v['supplemental_bookkeeping_FAIL_count']=sum(not r['PASS'] for r in checks);v['PASS']=v['PASS'] and all(r['PASS'] for r in checks);v['PASS_count']+=sum(r['PASS'] for r in checks);v['FAIL_count']+=sum(not r['PASS'] for r in checks);dump(OUT/'VALIDATION.json',v);dump(OUT/'REPORT_ONLY_CORRECTION.json',{'initial_report_error':'NameError: design undefined at task SHA rendering','correction':'post-A report script loads DESIGN_FREEZE; no frozen code/threshold/inference edits or experiment rerun','report_figures_rewritten_only':True,'scientific_results_unchanged':True});print('BOOKKEEPING',len(checks),'FAIL',sum(not r['PASS'] for r in checks))
if __name__=='__main__':main()
