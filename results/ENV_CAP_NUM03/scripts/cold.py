"""NUM03 cold audit: no import of integrator or analysis. Full-domain independent references."""
from pathlib import Path
import csv,json,hashlib,math,time,warnings
from decimal import Decimal,localcontext
import numpy as np
from independent_reference import Reference,mp_reference
O=Path(__file__).resolve().parents[1];R=O.parents[1];OLD=R/'results/ENV_CAP_FULL_CAMPAIGN';W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_NUM03')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,a):p.write_text(json.dumps(a,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def table(p,a):
 with p.open('w',encoding='utf-8',newline='') as f:
  d=csv.DictWriter(f,list(dict.fromkeys(k for r in a for k in r)));d.writeheader();d.writerows(a)
def findfolds(a,stride):
 rows=a[::stride];result=[]
 for i in range(1,len(rows)-1):
  ss=rows[i-1:i+2]
  if all(v.get('status')=='FINITE' and v['zero_boundary'] for v in ss) and ss[1]['range_m']<ss[0]['range_m'] and ss[1]['range_m']<ss[2]['range_m'] and 10000<=ss[1]['range_m']<=120000:result.append((i*stride,ss[1]))
 return result

def golden(ref,zs,a,b):
 g=(math.sqrt(5)-1)/2;c=b-g*(b-a);d=a+g*(b-a);fc=ref.calc(zs,c);fd=ref.calc(zs,d);best=min([fc,fd],key=lambda v:v['range_m']);best_angle=c if best is fc else d
 for i in range(70):
  if fc['range_m']<fd['range_m']:b,d,fd=d,c,fc;c=b-g*(b-a);fc=ref.calc(zs,c);angle=c;v=fc
  else:a,c,fc=c,d,fd;d=a+g*(b-a);fd=ref.calc(zs,d);angle=d;v=fd
  if v['range_m']<best['range_m']:best=v;best_angle=angle
 return dict(best,angle_deg=best_angle)

def main():
 if (W/'COLD_STARTED.json').exists():raise RuntimeError('No audit retry')
 cfg=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))
 for a in cfg['execution_inputs']:assert sha(a['path'])==a['sha256'],a['path']
 for a in json.loads((O/'PRIMARY_OUTPUT_FREEZE.json').read_text(encoding='utf-8'))['files']:assert sha(a['path'])==a['sha256'],a['path']
 js(W/'COLD_STARTED.json',{'time_unix':time.time(),'code_sha':sha(__file__),'freeze_sha':sha(O/'DESIGN_FREEZE.json')})
 checks=[];reference_rows=[];mp_rows=[];candidate_rows=[];angles_expected=np.linspace(.1,20,4001);old=json.loads((OLD/'E0_FREEZE.json').read_text(encoding='utf-8'));primary=json.loads((O/'NUM03_PRIMARY_RESULTS.json').read_text(encoding='utf-8'));full=list(csv.DictReader((O/'NUM03_RETURN_CURVES.csv').open(encoding='utf-8-sig')));co=list(csv.DictReader((O/'NUM03_COEFFICIENTS.csv').open(encoding='utf-8-sig')));ode=list(csv.DictReader((O/'NUM03_ODE_COMPARISON.csv').open(encoding='utf-8-sig')));mpplan=list(csv.DictReader((O/'NUM03_MP_REFERENCE_PLAN.csv').open(encoding='utf-8-sig')));group={}
 def ck(name,ok,detail=None):checks.append({'check':name,'PASS':bool(ok),'detail':detail})
 for a in full:group.setdefault((a['case_id'],a['SSP_method']),[]).append(a)
 for case in primary['cases']:
  name=case['case_id'];conf=next(a for a in old['E1_cases'] if a['case_id']==name);source=list(csv.DictReader((OLD/conf['profile']).open(encoding='utf-8-sig')));z=[float(a['depth_m']) for a in source];c=[float(a['sound_speed_mps']) for a in source];zs=conf['source_z_m'];H=conf['H_m'];print('AUDIT',name,flush=True)
  for kind in ['C_LINEAR','PCHIP']:
   ref=Reference(z,c,kind);rows=group[(name,kind)];ind=[];mr=case['methods'][kind]
   for j in range(len(z)-1):
    rr=next(a for a in co if a['case_id']==name and a['SSP_method']==kind and int(a['interval_index'])==j);ck(name+kind+f'_coeff{j}',max(abs(ref.coeff[j][i]-float(rr[k])) for i,k in enumerate(['cubic','quadratic','linear','constant']))<=1e-10)
   ck(name+kind+'_all4001_input_angles',len(rows)==4001 and all(int(a['angle_index'])==j and abs(float(a['angle_deg'])-angles_expected[j])<=1e-14 for j,a in enumerate(rows)))
   for j,a in enumerate(rows):
    try:
     with warnings.catch_warnings(record=True) as warn:
      warnings.simplefilter('always');v=ref.calc(zs,float(a['angle_deg']))
     dr=abs(v['range_m']-float(a['range_m'])) if a['status']=='FINITE' else None;dz=abs(v['turn_z_m']-float(a['turn_z_m'])) if a['status']=='FINITE' else None;ok=dr is not None and dr<=cfg['independent_range_tolerance_m'] and dz<=cfg['independent_turn_tolerance_m'] and v['bottom_encountered']==(a['bottom_encountered']=='True') and v['error_estimate_m']<=cfg['quadrature_estimate_limit_m'];ck(name+kind+f'_full_angle{j}',ok,None if ok else {'range_delta_m':dr,'turn_delta_m':dz})
     vv={'case_id':name,'SSP_method':kind,'angle_index':j,'angle_deg':float(a['angle_deg']),'range_m':v['range_m'],'turn_z_m':v['turn_z_m'],'bottom_encountered':v['bottom_encountered'],'zero_boundary':bool(not v['bottom_encountered'] and v['turn_z_m']<H-1),'main_reference_delta_m':dr,'turn_delta_m':dz,'adaptive_error_estimate_m':v['error_estimate_m'],'warnings':len(warn),'status':'FINITE'}
    except Exception as e:ck(name+kind+f'_full_angle{j}',False,repr(e));vv={'case_id':name,'SSP_method':kind,'angle_index':j,'angle_deg':float(a['angle_deg']),'status':'UNRESOLVED','reason':repr(e)}
    reference_rows.append(vv);ind.append(vv)
   ff=findfolds(ind,1);cf=findfolds(ind,2);ck(name+kind+'_all_grid_fold_count',len(ff)==mr['finite_grid_fold_count'],{'reference':len(ff),'primary':mr['finite_grid_fold_count']})
   if mr['first_refined_candidate'] and ff:
    i,v=min(ff,key=lambda a:a[1]['range_m']);main=mr['first_refined_candidate'];fit=golden(ref,zs,ind[i-1]['angle_deg'],ind[i+1]['angle_deg']);direct=ref.calc(zs,main['angle_deg']);dr=abs(fit['range_m']-main['range_m']);dz=abs(fit['turn_z_m']-main['turn_z_m']);okay=dr<=cfg['candidate_refined_range_tolerance_m'] and dz<=cfg['candidate_refined_turn_tolerance_m'] and fit['range_m']<min(ind[i-1]['range_m'],ind[i+1]['range_m']);ck(name+kind+'_first_refined_candidate',okay,{'range_delta_m':dr,'turn_delta_m':dz});candidate_rows.append({'case_id':name,'SSP_method':kind,'main_angle_deg':main['angle_deg'],'independent_angle_deg':fit['angle_deg'],'main_range_m':main['range_m'],'independent_range_m':fit['range_m'],'range_delta_m':dr,'turn_delta_m':dz,'direct_same_angle_delta_m':abs(direct['range_m']-main['range_m']),'PASS':okay})
   else:ck(name+kind+'_first_candidate_exists',False)
   for k,a in enumerate(v for v in ode if v['case_id']==name and v['SSP_method']==kind):
    arr=np.load(a['path_file']);positions=arr['depth_m'];theta=arr['angle_rad'];vals=np.array([ref.value(float(v)) for v in positions]);iv=np.cos(theta)/vals;p0=math.cos(math.radians(float(a['angle_deg'])))/ref.value(zs);err=float(np.max(abs(iv/p0-1)));q=ref.calc(zs,float(a['angle_deg']));delta=abs(float(a['range_m'])-q['range_m']) if a['range_m'] else None;good=a['success']=='True' and a['returned']=='True' and delta is not None and delta<=cfg['ODE_range_limit_m'] and err<=cfg['ODE_invariant_limit'];ck(name+kind+f'_ODE{k}',good and good==(a['PASS']=='True') and abs(err-float(a['invariant_error']))<1e-10,delta)
   orows=[a for a in ode if a['case_id']==name and a['SSP_method']==kind]
   for a in orows:
    b=next(v for v in orows if v['angle_deg']==a['angle_deg'] and v['step_max_m']!=a['step_max_m']);delta=abs(float(a['range_m'])-float(b['range_m'])) if a['range_m'] and b['range_m'] else None;ck(name+kind+'_ODE_step_'+a['angle_deg']+'_'+a['step_max_m'],delta is not None and delta<=cfg['ODE_step_limit_m'] and a['step_PASS']=='True')
   for a in (v for v in mpplan if v['case_id']==name and v['SSP_method']==kind):
    angle=float(a['angle_deg']);r60=mp_reference(z,c,zs,angle,kind,60);r80=mp_reference(z,c,zs,angle,kind,80);v=ref.calc(zs,angle);delta=abs(v['range_m']-r80['range_m']);
    with localcontext() as dc:
     dc.prec=100;conv=float(abs(Decimal(r60['range_decimal'])-Decimal(r80['range_decimal'])))
    ok=delta<=cfg['MP_range_tolerance_m'] and conv<=cfg['MP_precision_agreement_m'] and r80['quad_error_estimate_m']<=cfg['MP_quad_error_limit_m'] and abs(v['turn_z_m']-r80['turn_z_m'])<=cfg['MP_turn_tolerance_m'];ck(name+kind+'_MP_'+a['reason']+'_'+a['angle_deg'],ok,{'double_MP80_delta_m':delta,'MP60_MP80_delta_m':conv});mp_rows.append(dict(a,MP60_range_m=r60['range_m'],MP80_range_m=r80['range_m'],MP80_decimal=r80['range_decimal'],MP60_MP80_delta_m=conv,double_MP80_delta_m=delta,turn_delta_m=abs(v['turn_z_m']-r80['turn_z_m']),quad_error_estimate_m=r80['quad_error_estimate_m'],PASS=ok))
   print(name,kind,'full_domain_done','cumulative_FAIL',sum(not a['PASS'] for a in checks),flush=True)
 ck('9_cases_72018_full_range_rows',len(primary['cases'])==9 and len(full)==72018 and len(reference_rows)==72018)
 ck('new_external_solver_MC_zero',primary['new_solver_calls']==0 and primary['new_MC']==0)
 table(O/'NUM03_FULL_DOMAIN_INDEPENDENT_REFERENCE.csv',reference_rows);table(O/'NUM03_FIRST_CANDIDATE_INDEPENDENT.csv',candidate_rows);table(O/'NUM03_MP_REFERENCE.csv',mp_rows)
 result={'PASS':all(a['PASS'] for a in checks),'checks':len(checks),'FAIL_count':sum(not a['PASS'] for a in checks),'scope':'All registered 72018 angular entries independently reconstructed; own Hermite/roots/local transformed adaptive quadrature; independent golden candidate refinement; MP60/80 selected controls; saved differential-ray checks. Finite-domain evidence, not rigorous interval error bounds or ocean/energy/performance certification.','checks_detail':checks};js(O/'NUM03_INDEPENDENT_COLD_AUDIT.json',result);print(json.dumps({k:v for k,v in result.items() if k!='checks_detail'}));return 0 if result['PASS'] else 2
if __name__=='__main__':raise SystemExit(main())
