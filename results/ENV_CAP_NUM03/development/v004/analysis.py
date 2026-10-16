"""NUM03 once-only scientific confirmation after synthetic method-development admission."""
from integrator import *
from scipy.optimize import minimize_scalar
import time

def folds(rows,stride=1):
 a=rows[::stride];out=[]
 for i in range(1,len(a)-1):
  v=a[i-1:i+2]
  if all(x['status']=='FINITE' and x['zero_boundary_deep_return'] for x in v) and v[1]['range_m']<v[0]['range_m'] and v[1]['range_m']<v[2]['range_m'] and 10000<=v[1]['range_m']<=120000:out.append((i*stride,v[1]))
 return out

def refine(p,zs,rows,i):
 low=rows[i-1]['angle_deg'];high=rows[i+1]['angle_deg'];opt=minimize_scalar(lambda a:p.calc(zs,a)['range_m'],bounds=(low,high),method='bounded',options={'xatol':1e-8,'maxiter':100});candidates=[p.calc(zs,opt.x),p.calc(zs,rows[i]['angle_deg'])]
 if p.kind=='C_LINEAR':
  cs=float(p.f(zs))
  for z,c in zip(p.z,p.c):
   if z>zs and c>cs:
    a=math.degrees(math.acos(cs/c))
    if low<a<high:
     try:candidates.append(p.calc(zs,a))
     except ArithmeticError:pass
 best=min(candidates,key=lambda a:a['range_m']);best.update(bracket_low_deg=low,bracket_high_deg=high,optimizer_success=bool(opt.success),bracket_left_range_m=rows[i-1]['range_m'],bracket_right_range_m=rows[i+1]['range_m'],strict_interior_minimum=bool(best['range_m']<min(rows[i-1]['range_m'],rows[i+1]['range_m']) and low<best['angle_deg']<high))
 return best

def main():
 if (W/'FORMAL_STARTED.json').exists():raise RuntimeError('No formal retry')
 cfg=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))
 for a in cfg['execution_inputs']:assert sha(a['path'])==a['sha256'],a['path']
 assert json.loads((O/'DEVELOPMENT_GATE.json').read_text(encoding='utf-8'))['PASS']
 js(W/'FORMAL_STARTED.json',{'time_unix':time.time(),'freeze_sha':sha(O/'DESIGN_FREEZE.json'),'new_solver_calls':0,'new_MC':0})
 old=json.loads((OLD/'E0_FREEZE.json').read_text(encoding='utf-8'));cases=[a for a in old['E1_cases'] if a['env_id'] in ['SCS17_JW','SCS17_JA','SCS17_ANNUAL']];allrows=[];caseout=[];catalog=[];odeout=[];coeff=[];mpplan=[]
 hist={'SCS17_JW_z200':.40845,'SCS17_JW_z250':.7268499999999999,'SCS17_ANNUAL_z250':2.3387499999999997}
 for case in cases:
  name=case['case_id'];zs=case['source_z_m'];data=list(csv.DictReader((OLD/case['profile']).open(encoding='utf-8-sig')));z=np.array([float(a['depth_m']) for a in data]);c=np.array([float(a['sound_speed_mps']) for a in data]);out={'case_id':name,'source_z_m':zs,'source_family':'WOA23_17N118E','methods':{}};print('START',name,flush=True)
  for kind in ['C_LINEAR','PCHIP']:
   p=Profile(z,c,kind);rows=[];errors=[]
   for j,a in enumerate(np.linspace(.1,20,4001)):
    try:v=p.calc(zs,float(a));v.update(case_id=name,SSP_method=kind,angle_index=j)
    except Exception as e:v={'case_id':name,'SSP_method':kind,'angle_index':j,'angle_deg':float(a),'status':'NUMERIC_UNRESOLVED','reason':repr(e),'zero_boundary_deep_return':False};errors.append(j)
    rows.append(v)
   allrows.extend(rows);fine=folds(rows);coarse=folds(rows,2);first=min(fine,key=lambda a:a[1]['range_m']) if fine else None;cf=min(coarse,key=lambda a:a[1]['range_m']) if coarse else None
   mr={'all_angles':len(rows),'numerical_failure_indices':errors,'max_quadrature_error_estimate_m':max((r['error_estimate_m'] for r in rows if r['status']=='FINITE'),default=None),'first_grid_candidate':first[1] if first else None,'coarse_grid_candidate':cf[1] if cf else None,'finite_grid_fold_count':len(fine),'refined_candidates':[]}
   for i,a in fine:
    try:rr=refine(p,zs,rows,i);rr.update(case_id=name,SSP_method=kind,grid_angle_index=i);catalog.append(rr);mr['refined_candidates'].append(rr)
    except Exception as e:mr['refined_candidates'].append({'grid_angle_index':i,'status':'REFINEMENT_UNRESOLVED','reason':repr(e)})
   good=[a for a in mr['refined_candidates'] if a.get('strict_interior_minimum')];selected=min(good,key=lambda a:a['range_m']) if good else None;mr['first_refined_candidate']=selected
   mr['primary_numeric_PASS']=bool(not errors and mr['max_quadrature_error_estimate_m']<=cfg['quadrature_estimate_limit_m'])
   mr['refinement_all_PASS']=all(a.get('optimizer_success',False) and a.get('strict_interior_minimum',False) for a in mr['refined_candidates'])
   mr['angle_resolution_PASS']=bool(first and cf and abs(first[1]['range_m']-cf[1]['range_m'])<=cfg['coarse_fine_candidate_limit_m'])
   for jj in range(len(z)-1):
    a=p.f.c[:,jj].tolist();a=[0.]*(4-len(a))+a;coeff.append({'case_id':name,'SSP_method':kind,'interval_index':jj,'z_low_m':float(z[jj]),'z_high_m':float(z[jj+1]),'cubic':a[0],'quadratic':a[1],'linear':a[2],'constant':a[3]})
   if selected:
    mpplan.append({'case_id':name,'SSP_method':kind,'angle_deg':selected['angle_deg'],'reason':'FIRST_REFINED_CANDIDATE_FIXED_SELECTION_RULE'})
   if name in hist:mpplan.append({'case_id':name,'SSP_method':kind,'angle_deg':hist[name],'reason':'PREREGISTERED_CZROB02_FAILED_ANGLE'})
   if kind=='PCHIP' and name in cfg['MP_fixed_cases']:
    for a in [1.,5.]:mpplan.append({'case_id':name,'SSP_method':kind,'angle_deg':a,'reason':'FIXED_NONSELECTED_ANGLE_CONTROL'})
   tests=[]
   if name in cfg['ODE_representative_cases'] and mr['primary_numeric_PASS'] and selected:
    angles=[selected['angle_deg']+d for d in [-.1,0.,.1]]
    if name in hist and all(abs(a-hist[name])>1e-8 for a in angles):angles.append(hist[name])
    for a in angles:
     ref=p.calc(zs,a);pair=[]
     for step in cfg['ODE_max_steps_m']:
      v,rx,yy=p.ode(zs,a,step);tag=f'{name}_{kind}_{a:.9f}_{step:g}';path=W/(tag+'.npz');np.savez_compressed(path,range_m=rx,depth_m=yy[0],angle_rad=yy[1]);v.update(case_id=name,SSP_method=kind,angle_deg=a,Snell_range_m=ref['range_m'],path_file=str(path));v['range_difference_m']=abs(v['range_m']-ref['range_m']) if v['range_m'] is not None else None;v['PASS']=bool(v['success'] and v['returned'] and v['range_difference_m']<=cfg['ODE_range_limit_m'] and v['invariant_error']<=cfg['ODE_invariant_limit']);pair.append(v);odeout.append(v);tests.append(v)
     delta=abs(pair[0]['range_m']-pair[1]['range_m']) if all(t['range_m'] is not None for t in pair) else None
     for t in pair:t['step_difference_m']=delta;t['step_PASS']=delta is not None and delta<=cfg['ODE_step_limit_m']
   mr['ODE_count']=len(tests);mr['ODE_PASS']=all(a['PASS'] and a['step_PASS'] for a in tests) if tests else None;out['methods'][kind]=mr;print(name,kind,'primary_numeric',mr['primary_numeric_PASS'],'folds',len(fine),'ODE',mr['ODE_PASS'],flush=True)
   js(W/(name+'_'+kind+'_EXIT.json'),mr)
  a=out['methods']['C_LINEAR']['first_refined_candidate'];b=out['methods']['PCHIP']['first_refined_candidate'];out['representation_range_shift_m']=abs(a['range_m']-b['range_m']) if a and b else None;out['representation_turn_shift_m']=abs(a['turn_z_m']-b['turn_z_m']) if a and b else None;caseout.append(out);js(W/(name+'_EXIT.json'),out)
 tab(O/'NUM03_RETURN_CURVES.csv',allrows);tab(O/'NUM03_CANDIDATE_CATALOGUE.csv',catalog);tab(O/'NUM03_ODE_COMPARISON.csv',odeout);tab(O/'NUM03_COEFFICIENTS.csv',coeff);tab(O/'NUM03_MP_REFERENCE_PLAN.csv',mpplan)
 js(O/'NUM03_PRIMARY_RESULTS.json',{'status':'PRE_COLD','cases':caseout,'all9_retained':len(caseout)==9,'new_solver_calls':0,'new_MC':0,'energy_attribution':'UNRESOLVED_NOT_A_GEOMETRY_VETO','Yang_MMAR_performance':'NOT_EVALUATED'})
 js(O/'ANALYSIS_EXIT.json',{'completed':True,'all9_retained':len(caseout)==9,'formal_runs':1,'new_solver_calls':0,'new_MC':0})
if __name__=='__main__':main()
