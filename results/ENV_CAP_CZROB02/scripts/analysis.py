"""Single CZROB02 registered analysis; deterministic geometry, no MC or acoustic executable."""
from physics import *
import time

def main():
 cfg=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))
 if (W/'FORMAL_STARTED.json').exists():raise RuntimeError('No scientific retry')
 for a in cfg['execution_inputs']:assert sha(a['path'])==a['sha256'],a['path']
 js(W/'FORMAL_STARTED.json',{'time_unix':time.time(),'freeze_sha':sha(O/'DESIGN_FREEZE.json'),'new_external_solver_calls':0,'new_MC':0})
 old=json.loads((OLD/'E0_FREEZE.json').read_text(encoding='utf-8'));rec=json.loads((REC/'REC01_CZ_PHYSICS_ADJUDICATION.json').read_text(encoding='utf-8'));cases=[a for a in old['E1_cases'] if a['env_id'] in ('SCS17_JW','SCS17_JA','SCS17_ANNUAL')]
 full=[];topology=[];foldcatalog=[];profilediag=[];nodal=[];odes=[];crossings=[];paths=[];percase=[];coefficients=[]
 for case in cases:
  name=case['case_id'];print('START',name,flush=True);zs=case['source_z_m'];H=case['H_m'];oldcase=next(a for a in rec['cases'] if a['case_id']==name)
  data=list(csv.DictReader((OLD/case['profile']).open(encoding='utf-8-sig')));z=np.array([float(a['depth_m']) for a in data]);c=np.array([float(a['sound_speed_mps']) for a in data]);assert z[0]==0 and z[-1]==H
  result={'case_id':name,'source_z_m':zs,'source_family':'WOA23_17N118E','methods':{},'old_energy_peak_km':oldcase['focus_metrics'][-1]['peak_km'],'old_primary_node_distance_m':oldcase['fold']['turn_nearest_SSP_node_m'],'old_refined_node_distance_m':oldcase['refined_angle_fold']['turn_nearest_SSP_node_m']}
  for method in ['C_LINEAR','PCHIP']:
   p=Profile(z,c,method);rows=[];delta=[];errors=[]
   for jj in range(len(z)-1):
    coeff=p.f.c[:,jj].tolist();coeff=[0.]*(4-len(coeff))+coeff;coefficients.append(dict(case_id=name,SSP_method=method,interval_index=jj,z_low_m=float(z[jj]),z_high_m=float(z[jj+1]),cubic=coeff[0],quadratic=coeff[1],linear=coeff[2],constant=coeff[3]))
   for j,a in enumerate(np.linspace(.1,20,cfg['angle_samples'])):
    try:
     coarse=p.evaluate(zs,float(a),32);fine=p.evaluate(zs,float(a),64);d=max(abs(coarse['deep_return_km']-fine['deep_return_km']),abs(coarse['upper_return_km']-fine['upper_return_km']));fine.update(case_id=name,SSP_method=method,angle_index=j,deep_return_Q32_km=coarse['deep_return_km'],upper_return_Q32_km=coarse['upper_return_km'],max_order_difference_km=d);delta.append(d)
    except Exception as e:
     fine={'case_id':name,'SSP_method':method,'angle_index':j,'angle_deg':float(a),'status':'NUMERIC_UNRESOLVED','reason':repr(e),'deep_bottom_bounce':None,'deep_return_km':None};errors.append(j)
    rows.append(fine)
   full.extend(rows);folds=extrema(rows);coarsefolds=extrema(rows,2);first=min(folds,key=lambda a:a['deep_return_km']) if folds else None;cf=min(coarsefolds,key=lambda a:a['deep_return_km']) if coarsefolds else None
   angleok=bool(first and cf and abs(first['deep_return_km']-cf['deep_return_km'])<=cfg['angle_fold_tolerance_km']);quadok=not errors and max(delta)<=cfg['quadrature_tolerance_km'];mr={'method':method,'failure_angles':errors,'finite_angles':len(rows)-len(errors),'quadrature_max_delta_km':max(delta) if delta else None,'quadrature_PASS':quadok,'folds_count':len(folds),'first_fold':first,'coarse_first_fold':cf,'angle_refinement_PASS':angleok};result['methods'][method]=mr
   for i,f in enumerate(folds):foldcatalog.append(dict(case_id=name,SSP_method=method,fold_index=i,**{k:v for k,v in f.items() if k not in ('case_id','SSP_method')}))
   # Node sampling and derivative jumps: all original intervals, not a chosen focus.
   for j,(zz,cc) in enumerate(zip(z,c)):
    profilediag.append({'case_id':name,'SSP_method':method,'node_index':j,'depth_m':float(zz),'sample_c_mps':float(cc),'interpolated_c_mps':float(p.v(zz)),'left_derivative_s_inv':float(p.df(max(0,zz-1e-6))),'right_derivative_s_inv':float(p.df(min(H,zz+1e-6)))})
    if zz>zs and cc>float(p.v(zs)):
     a=math.degrees(math.acos(float(p.v(zs))/cc))
     if .1001<=a<=19.9999:
      for side in [-1,0,1]:
       aa=a+side*cfg['node_probe_delta_deg']
       try:pp=p.evaluate(zs,aa,64);row={k:v for k,v in pp.items() if k not in ('status',)};row['status']='FINITE'
       except Exception as e:row={'status':'UNRESOLVED_OR_GRAZING','reason':repr(e)}
       nodal.append(dict(case_id=name,SSP_method=method,node_index=j,node_depth_m=float(zz),side=side,probe_angle_deg=aa,**row))
   # Representative ray ODE cases registered before the science run.
   checkrows=[]
   if name in cfg['representative_cases'] and quadok:
    angles=[oldcase['fold']['angle_deg']+da for da in cfg['ODE_angle_offsets_deg']]
    if first and all(abs(first['angle_deg']-a)>1e-8 for a in angles):angles.append(first['angle_deg'])
    for angle in angles:
     reference=p.evaluate(zs,angle,64);two=[]
     for step in cfg['ODE_max_steps_m']:
      met,rr,yy=p.ode(zs,angle,step);tag=f'{name}_{method}_a{angle:.8f}_s{step:g}';np.savez_compressed(W/(tag+'.npz'),range_m=rr,depth_m=yy[0],angle_rad=yy[1]);met.update(case_id=name,SSP_method=method,angle_deg=angle,Snell_reference_km=reference['deep_return_km'],path_file=str(W/(tag+'.npz')))
      met['Snell_difference_km']=abs(met['range_km']-reference['deep_return_km']) if met['range_km'] is not None else None
      met['PASS']=bool(met['success'] and met['returned'] and met['Snell_difference_km']<=cfg['ODE_Snell_tolerance_km'] and met['max_relative_Snell_invariant_error']<=cfg['ODE_invariant_relative_tolerance'])
      odes.append(met);two.append(met);checkrows.append(met)
     dd=abs(two[0]['range_km']-two[1]['range_km']) if all(t['range_km'] is not None for t in two) else None
     for t in two:t['step_refinement_difference_km']=dd;t['step_refinement_PASS']=dd is not None and dd<=cfg['ODE_step_tolerance_km']
   mr['ODE_cases']=len(checkrows);mr['ODE_PASS']=all(a['PASS'] and a['step_refinement_PASS'] for a in checkrows) if checkrows else None
   mr['numeric_PASS']=bool(quadok and angleok and (mr['ODE_PASS'] is not False))
   js(W/(name+'_'+method+'_EXIT.json'),mr);print(name,method,'fold',first['deep_return_km'] if first else None,'quad',quadok,'ODE',mr['ODE_PASS'],flush=True)
  cl=result['methods']['C_LINEAR'];pc=result['methods']['PCHIP'];fcl=cl['first_fold'];fpc=pc['first_fold'];both=bool(fcl and fpc);dR=abs(fcl['deep_return_km']-fpc['deep_return_km']) if both else None;dZ=abs(fcl['deep_turn_z_m']-fpc['deep_turn_z_m']) if both else None
  numerical=cl['numeric_PASS'] and pc['numeric_PASS'];stable=bool(numerical and both and dR<=cfg['representation_range_tolerance_km'] and dZ<=cfg['representation_turn_tolerance_m']);classification='GEOMETRY_NUMERIC_UNCERTIFIED' if not numerical else ('SSP_GEOMETRY_ROBUST_FOR_TESTED_REPRESENTATIONS' if stable else 'SSP_INTERPOLATION_SENSITIVE')
  result.update(numerical_PASS=bool(numerical),geometry_stable=stable,classification=classification,range_shift_km=dR,turn_depth_shift_m=dZ,energy_attribution='EXCLUSIVE_ENERGY_ATTRIBUTION_UNRESOLVED',coherent_provider_investment='DEFER')
  topology.append({'case_id':name,'old_primary_node_distance_m':result['old_primary_node_distance_m'],'old_refined_node_distance_m':result['old_refined_node_distance_m'],'C_linear_first_fold_km':fcl['deep_return_km'] if fcl else None,'PCHIP_first_fold_km':fpc['deep_return_km'] if fpc else None,'C_linear_turn_z_m':fcl['deep_turn_z_m'] if fcl else None,'PCHIP_turn_z_m':fpc['deep_turn_z_m'] if fpc else None,'range_shift_km':dR,'turn_depth_shift_m':dZ,'C_linear_folds':cl['folds_count'],'PCHIP_folds':pc['folds_count'],'C_linear_numeric_PASS':cl['numeric_PASS'],'PCHIP_numeric_PASS':pc['numeric_PASS'],'geometry_stable':stable,'classification':classification})
  original=oldcase['fold']['deep_return_km'];lo=original-5;hi=original+5
  for lev in ['R0','R1']:
   cr=path_archive(RAW/name/(lev+'.ray'),zs,H,lo,hi);crossings.extend(dict(case_id=name,ray_level=lev,**a) for a in cr);inside=[a for a in cr if a['in_original_focus_band']];zero=sum(a['identity']=='ZERO_BOUNDARY' for a in inside);bd=sum(a['identity']!='ZERO_BOUNDARY' for a in inside)
   paths.append({'case_id':name,'data':'ARCHIVED_C_LINEAR_'+lev,'focus_low_km':lo,'focus_high_km':hi,'registered_downward_rays_only':True,'all_10_120km_crossings':len(cr),'zero_boundary_focus_crossings':zero,'boundary_focus_crossings':bd,'surface_focus_crossings':sum(a['prefix_surface_contacts']>0 for a in inside),'bottom_focus_crossings':sum(a['prefix_bottom_contacts']>0 for a in inside),'REFRACTED_PATH_IDENTITY_SUPPORTED':zero>0,'boundary_path_overlap_geometric':bd>0,'BOUNDARY_PATH_CONTAMINATION_PRESENT':'NOT_ESTABLISHED_AS_ENERGY','exclusive_energy_attribution':'EXCLUSIVE_ENERGY_ATTRIBUTION_UNRESOLVED','energy_fraction':'NOT_EVALUATED','no_initial_upward_fan':'true','new_PCHIP_energy':'NOT_EVALUATED'})
  percase.append(result)
  js(W/(name+'_EXIT.json'),result)
 tab(O/'CZROB02_INTERPOLATION_COEFFICIENTS.csv',coefficients);tab(O/'CZROB02_SSP_TOPOLOGY.csv',topology);tab(O/'CZROB02_RETURN_CURVES.csv',full);tab(O/'CZROB02_ALL_FOLDS.csv',foldcatalog);tab(O/'CZROB02_SSP_NODE_DIAGNOSTICS.csv',profilediag);tab(O/'CZROB02_NODE_PROBES.csv',nodal);tab(O/'CZROB02_ODE_COMPARISON.csv',odes);tab(O/'CZROB02_RAY_PREFIX_CROSSINGS.csv',crossings);tab(O/'CZROB02_PATH_ATTRIBUTION.csv',paths)
 js(O/'CZROB02_PRIMARY_RESULTS.json',{'status':'PRE_COLD','cases':percase,'new_external_solver_calls':0,'new_MC':0,'PCHIP_energy':'NOT_EVALUATED','exclusive_path_energy':'UNRESOLVED','decision_rule':'No coherent-provider investment recommendation until both physical structure and energetic attribution prerequisites close; no new energy decomposition is available here.'})
 js(O/'ANALYSIS_EXIT.json',{'completed':True,'cases':len(percase),'all_9_retained':len(percase)==9,'geometrically_stable_cases':sum(a['geometry_stable'] for a in percase),'numeric_failure_cases':sum(not a['numerical_PASS'] for a in percase),'all_energy_attributions_unresolved':True,'new_external_solver_calls':0,'new_MC':0})
if __name__=='__main__':main()
