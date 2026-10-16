"""Once-only E1 registered environment units; all inputs produced before Commit-0."""
from env_provider import *
def main():
 config=json.loads((O/'E0_FREEZE.json').read_text(encoding='utf-8'));results=[];allcurves=[];raychecks=[];compares=[]
 assert json.loads((O/'E0_PROVIDER_CONTROL_VALIDATION.json').read_text())['PASS']
 for case in config['E1_cases']:
  name=case['case_id'];D=W/'E1'/name;marker=D/'SCIENCE_STARTED.json'
  if marker.exists():raise RuntimeError('Refusing repeated science unit '+name)
  js(marker,{'unit':name,'config_sha':sha(O/'E0_FREEZE.json'),'code_sha':sha(__file__),'provider_sha':sha(Path(__file__).with_name('env_provider.py')),'input_sha':case['input_sha'],'budget_s':1200,'registered_calls':5,'seed':'DETERMINISTIC_NO_MC'})
  z,c=loadprofile(O/case['profile']);zs=case['source_z_m'];H=case['H_m']; row={'case_id':name,'env_id':case['env_id'],'source_z_m':zs,'receiver_band_m':f'{zs-25}:{zs+25}','bottom_m':H,'source_kind':case['source_kind'],'independent_check':'Snell exact C-linear integrals vs Bellhop; TL same-family refinement','evidence_level':'L0_ANALYTIC'}
  start=time.perf_counter()
  try:
   analytic=return_curve(z,c,zs,np.linspace(.1,20,2001));table(D/'SNELL_RETURN_CURVE.csv',analytic,['angle_deg','first_return_km','turn_z_m','c_phase_mps','surface_bounces','bottom_bounces'])
   in_domain=[r for r in analytic if 10<=r['first_return_km']<=120]
   if not in_domain:
    row.update(CZ_class='CZ_NOT_PRESENT_IN_TESTED_DOMAIN',reason='No bottom-free first deep return in registered 10..120km / downward0.1..20deg family; other paths not claimed absent',allow_methods=False)
   else:
    aa=np.array([r['angle_deg'] for r in analytic]);ar=np.array([r['first_return_km'] for r in analytic]);interior,_=find_peaks(-ar)
    eligible=[i for i in interior if 10<=ar[i]<=120]
    if eligible:
     ix=min(eligible,key=lambda i:ar[i]);expected=float(ar[ix]);fold=True
    else:
     ix=int(np.argmin(ar));expected=float(ar[ix]);fold=False
    row.update(ray_first_caustic_km=expected,fold_interior=fold,ray_return_min_km=float(ar.min()),ray_return_max_km=float(ar.max()),turn_z_m=float(analytic[ix]['turn_z_m']))
    bs=[]
    for lev in ['R0','R1']:
     p=D/(lev+'.env');execute(p,name);rays=readrays(p.with_suffix('.ray'));b=ray_first_returns(rays,zs,H);table(D/(lev+'_FIRST_RETURNS.csv'),b,['angle_deg','first_return_km','turn_z_m','surface_bounces','bottom_bounces'])
     errs=[]
     for r in b:
      if aa[0]<=r['angle_deg']<=aa[-1] and 10<=r['first_return_km']<=120:
       e=abs(r['first_return_km']-float(np.interp(r['angle_deg'],aa,ar)));errs.append(e)
     bs.append(errs);raychecks.append({'case_id':name,'level':lev,'n_compared':len(errs),'max_return_range_difference_km':max(errs) if errs else 'NOT_EVALUATED'})
    rayok=all(len(v)>=10 and max(v)<=.1 for v in bs)
    field=[];rr=None
    for lev in ['I0','I1','I2']:
     p=D/(lev+'.env');execute(p,name);rr,zz,pp=readshd(p.with_suffix('.shd'));tl=energycurve(rr,zz,pp,zs);f=focus(rr,tl,expected);field.append((tl,f))
     np.savez_compressed(D/(lev+'_ENERGY_CURVE.npz'),r_km=rr,TL_dB=tl,depth_m=zz)
     allcurves.extend({'case_id':name,'level':lev,'range_km':float(r),'band_energy_TL_dB':float(t)} for r,t in zip(rr,tl))
    mask=np.abs(rr-expected)<=5;metrics=[]
    for i,j in [(0,1),(1,2)]:
     delta=np.abs(field[i][0][mask]-field[j][0][mask]);fi,fj=field[i][1],field[j][1]
     met={'case_id':name,'comparison':f'I{i}_I{j}','median_TL_delta_dB':float(np.median(delta)),'P95_TL_delta_dB':float(np.quantile(delta,.95)),'peak_delta_km':abs(fi['peak_km']-fj['peak_km']),'width_delta_km':abs(fi['width_km']-fj['width_km'])}
     met['PASS']=met['median_TL_delta_dB']<=1 and met['P95_TL_delta_dB']<=2 and met['peak_delta_km']<=1 and met['width_delta_km']<=2;metrics.append(met);compares.append(met)
    f=field[-1][1];num=max(x['P95_TL_delta_dB'] for x in metrics);enhance=f['prominence_dB']>max(.5,2*num) and not f['boundary_selected']
    row.update(**f,ray_integral_gate=rayok,TL_refinement_gate=all(m['PASS'] for m in metrics),enhancement_above_numeric_uncertainty=enhance,background_flanks_km='5..10 from peak',numerical_P95_dB=num,competing_paths='Full-field surface/bottom energy not independently decomposed; verified zero-bounce deep-return rays coexist')
    passed=fold and rayok and all(m['PASS'] for m in metrics) and enhance
    row.update(CZ_class='CZ2_RAY_TL_STABLE' if passed else 'CZ_NUMERIC_UNCERTIFIED',allow_methods=passed,evidence_level='L2_CZ_RAY_NUMERIC' if passed else 'L1_ENV_SOURCE_CONDITIONAL_RAY_DIAGNOSTIC',reason='Conditional nominal range-independent geometry; no coherent pressure admission' if passed else 'One or more preregistered caustic/topology/TL numerical gates not closed')
   row['execution_status']='COMPLETED'
  except Exception as e:
   row.update(CZ_class='CZ_NUMERIC_UNCERTIFIED',allow_methods=False,execution_status='FAILED_NO_RETRY',reason=repr(e));js(D/'FAILURE.json',{'error':repr(e),'cancel_remaining_calls_for_this_unit':True})
  row['elapsed_s']=time.perf_counter()-start;results.append(row);js(D/'UNIT_EXIT.json',row);table(O/'E1_CZ_ENVIRONMENTS.csv',results);print(name,row['CZ_class'],row.get('reason'),flush=True)
 table(O/'E1_CZ_RANGE_SHAPE_DIAGNOSTICS.csv',allcurves,['case_id','level','range_km','band_energy_TL_dB']);table(O/'E1_RAY_INDEPENDENT_COMPARISON.csv',raychecks,['case_id','level','n_compared','max_return_range_difference_km']);table(O/'E1_NUMERICAL_STABILITY.csv',compares,['case_id','comparison','median_TL_delta_dB','P95_TL_delta_dB','peak_delta_km','width_delta_km','PASS'])
 unique=sorted(set(x['env_id'] for x in results if x['CZ_class']=='CZ2_RAY_TL_STABLE'));families=sorted(set(x['source_family'] for x in config['E1_cases'] if x['env_id'] in unique))
 js(O/'E1_CZ_GATE.json',{'CZ2_subcases':sum(x['CZ_class']=='CZ2_RAY_TL_STABLE' for x in results),'CZ2_env_ids':unique,'distinct_source_families':families,'cross_environment_claim_allowed':len(families)>=2,'main_CZ_admitted':bool(unique),'coherent_complex_modal_provider':'NOT_ADMITTED','not_evaluated_sources':config['not_evaluated_sources'],'no_reruns':True})
if __name__=='__main__':main()
