"""Final audit verdict; preserves frozen automatic outcome and every saved scientific table."""
from pathlib import Path
import r4_e2_nine_pair_pilot as p
O=p.O
raw=p.read(O/'PILOT_DECISION.json')
assert raw['decision']=='NINE_PAIR_NONBEARING_MECHANISM_PROMISING'
p.dump(O/'FROZEN_RULE_DECISION.json',raw)
chain=p.base.rows(O/'ANALYTIC_CHAIN_RULE_AND_RANK_AUDIT.csv');extra=[]
for r in chain:
 key=r['scene_id']+'_'+r['mesh']+'_'+r['representation']
 rank=int(r['analytic_horizontal_chain_rank']);fd=int(r['frozen_fd_rank'])
 extra.append(dict(check='physical_chain_rank:'+key,PASS=rank==fd,detail='analytic='+str(rank)+'; frozenFD='+str(fd)))
 er=float(r['range_information_relative_difference']);ez=float(r['depth_information_relative_difference'])
 extra.append(dict(check='physical_chain_effective_information:'+key,PASS=er<.1 and ez<.1,detail='range='+str(er)+'; depth='+str(ez)))
p.write(O/'POST_EXECUTION_FORMULA_CONTROLS.csv',extra)
val=p.read(O/'VALIDATION.json');old=dict(val);p.dump(O/'RECONSTRUCTION_VALIDATION.json',old)
checks=p.base.rows(O/'INDEPENDENT_AUDIT_CHECKS.csv')
checks=[dict(check=r['check'],PASS=r['PASS']=='True',detail=r['detail']) for r in checks]
checks.extend(extra);p.write(O/'INDEPENDENT_AUDIT_CHECKS.csv',checks)
val.update(checks=len(checks),PASS=sum(r['PASS'] for r in checks),FAIL=sum(not r['PASS'] for r in checks),reconstruction_checks=old['checks'],reconstruction_FAIL=old['FAIL'],post_execution_formula_checks=len(extra),post_execution_formula_FAIL=sum(not r['PASS'] for r in extra),physical_chain_check='FAILED_SINGLE_HLA_BEARING_PROFILED_RANK_AND_EFFECTIVE_DEPTH',scope=old['scope']+'; exact analytic horizontal pressure derivative and geometric chain rule; single-HLA nonbearing physical rank ceiling4.')
p.dump(O/'VALIDATION.json',val)
raw.update(decision='IMPLEMENTATION_INVALID',automatic_frozen_rule_outcome='NINE_PAIR_NONBEARING_MECHANISM_PROMISING',automatic_rule_outcome_ACCEPTED=False,reason='Post-execution independent physical chain rule: P1 bearing-profiled true rank4 but frozen finite-difference rank5 in all48 scene-mesh cases; single-array effective depth difference exceeds10pct in some cases. Aggregate2pct Jacobian check does not protect weak/null directions.',independent_audit='FAIL_PHYSICAL_FORMULA_CHECK',post_execution_formula_failures=val['post_execution_formula_FAIL'],next='STOP_FOR_FORMULA_REPAIR_REVIEW',next_stage_execution='NOT_AUTHORIZED',full_band_investment='NOT_ADMITTED_BY_THIS_PILOT',H3_execution='NOT_AUTHORIZED',C1_primary_dual_chain_rank_valid=True,scope='Dual-array primary numbers are retained exploratory values; entire pilot not admitted.')
p.dump(O/'PILOT_DECISION.json',raw)
print('FINAL',raw['decision'],'VALIDATION',val['PASS'],val['FAIL'])
