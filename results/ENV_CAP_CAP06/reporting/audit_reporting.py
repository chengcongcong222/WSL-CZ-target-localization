from pathlib import Path
import csv,json,math,hashlib,numpy as np
O=Path(__file__).resolve().parents[1]
def rows(n):return list(csv.DictReader((O/n).open(encoding='utf-8-sig')))
def main():
 checks=[]
 def ck(n,b):checks.append(dict(check_id=n,PASS=bool(b)))
 freeze=json.loads((O/'CODE_FREEZE.json').read_text())
 for p,h in freeze.items():ck('frozen_code:'+p,hashlib.sha256((O/p).read_bytes()).hexdigest()==h)
 for phase in ['P1','P2','P3']:
  a=json.loads((O/(phase+'_INDEPENDENT_COLD_AUDIT.json')).read_text());ck(phase+':audit_internal_count',a['checks']==len(a['records']) and a['FAIL']==sum(not x['PASS'] for x in a['records']) and a['PASS']==all(x['PASS'] for x in a['records']))
 resources=rows('CAP06_RESOURCE_PARETO.csv');support=rows('P3_GLOBAL_SUPPORT_STATUS.csv');available='resource' in support[0] and json.loads((O/'P3_INDEPENDENT_COLD_AUDIT.json').read_text())['PASS'] and json.loads((O/'P2_INDEPENDENT_COLD_AUDIT.json').read_text())['PASS']
 for r in resources:
  if available:
   z=[x for x in support if x['split']=='holdout' and x['resource']==r['resource']];a=[float(x['range_width_m']) for x in z if x['range_width_m']];b=[float(x['speed_width_mps']) for x in z if x['speed_width_mps']];ck(r['resource']+':resource_reconstruction',len(z)==int(r['holdout_base_trajectories']) and sum(x['truth_in_retained_outer']=='True' for x in z)==int(r['outer_truth_retained']) and abs(float(np.median(a))-float(r['median_range_outer_width_m']))<1e-9 and abs(float(np.median(b))-float(r['median_speed_outer_width_mps']))<1e-9)
   dom=False
   keys=['channels','baseline_km','calibration_transmissions','median_range_outer_width_m','median_speed_outer_width_mps']
   for other in resources:
    aa=[float(other[k]) for k in keys]+[-int(other['outer_truth_retained'])];bb=[float(r[k]) for k in keys]+[-int(r['outer_truth_retained'])];dom |= all(x<=y for x,y in zip(aa,bb)) and any(x<y for x,y in zip(aa,bb))
   ck(r['resource']+':finite_frontier',(r['finite_support_nondominated']=='True')==(not dom))
  else:ck(r['resource']+':not_evaluated',r['finite_support_nondominated']=='NOT_EVALUATED')
 for x in rows('EVIDENCE_LEVEL_MATRIX.csv'):ck('claim_source:'+x['claim'],hashlib.sha256((O/x['source']).read_bytes()).hexdigest()==x['source_SHA'])
 rec=rows('P2_RECORD_INDEX.csv');ck('raw_denominator',len(rec)==4416);ck('main_records',sum(x['stress']=='MAIN_TRACK' for x in rec)==2880);ck('stress_records',sum(x['stress']!='MAIN_TRACK' for x in rec)==1536);ck('development_and_holdout_separate',set(x['state_id'] for x in rec if x['split']=='development').isdisjoint(set(x['state_id'] for x in rec if x['split']=='holdout')))
 # Each resource uses the same main token, never a re-generated main observation.
 if available:
  by={}
  for x in support:by.setdefault(x['source_main'],set()).add(x['resource'])
  ck('paired_main_identity',len(by)==320 and all(v=={'S0','S2_B0','S2_B3'} for v in by.values()))
 usage=json.loads((O/'CAP06_RESOURCE_USAGE.json').read_text());ck('budgets',usage['scientific_wall_s']<=10800 and usage['peak_process_tree_RSS']<=8*2**30 and usage['D_campaign_bytes']<=12*2**30 and usage['git_result_bytes']<=120*2**20 and usage['external_solver_calls']==0)
 d=json.loads((O/'CAP06_FINAL_DECISION.json').read_text());ck('scope',d['R4_percent']==0 and d['first_CZ_5D_performance']=='NOT_EVALUATED' and d['hardware']=='UNKNOWN')
 out=dict(PASS=all(x['PASS'] for x in checks),checks=len(checks),FAIL=sum(not x['PASS'] for x in checks),records=checks,scope='Independent report arithmetic and provenance; no repeated scientific execution')
 (O/'P4_INDEPENDENT_COLD_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print(out['PASS'],out['checks'],out['FAIL'])
if __name__=='__main__':main()
