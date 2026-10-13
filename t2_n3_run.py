"""One sequential N3 run; first analytic hard failure cancels every later call."""
from t2_n3_common import *
import subprocess,shutil,traceback,psutil

def main():
 if (OUT/'SCIENCE_STARTED.json').exists():raise RuntimeError('SCIENCE_ONCE_ONLY')
 verify_freeze();design_sha=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0];assert design_sha==remote and design_sha!=PARENT
 start=time.monotonic();dump(OUT/'SCIENCE_STARTED.json',{'utc':utc(),'design_SHA':design_sha,'once_only':True});calls=[];cross=[];modes=[];meshrows=[];extremes=[];globalrows=[];failures=[];g=arrays();fields={};analytical_PASS=0;peak=0
 plan=read(OUT/'SOLVER_PLAN.json')
 for item in plan:
  if item['kind']=='E_STD_V' and analytical_PASS!=16:break
  cwd=LOCAL/'runs'/item['stem'];cwd.mkdir(parents=True);env=Path(item['ENV']);assert sha(env)==item['sha256'];shutil.copyfile(env,cwd/env.name);binary=AT/'kraken.exe';assert sha(binary)==next(r['sha256'] for r in read(OUT/'OLD_N2_SOURCE_IDENTITIES.json')['records'] if r['role']=='frozen_solver_binary')
  row={'call_order':item['call_order'],'kind':item['kind'],'boundary':item['boundary'],'frequency_hz':item['f'],'CHIGH':item['high'],'mesh':item['mesh'],'ENV_sha256':sha(env),'binary_sha256':sha(binary),'stem':item['stem'],'cwd':str(cwd),'timeout_s':180,'status':'STARTED','budget_before':json.dumps(guard(start))};calls.append(row);table(OUT/'SOLVER_CALLS.csv',calls)
  stdout=cwd/(item['stem']+'.stdout.txt');stderr=cwd/(item['stem']+'.stderr.txt');began=time.monotonic();rc=None
  try:
   with stdout.open('wb') as so,stderr.open('wb') as se:
    proc=subprocess.Popen([str(binary),item['stem']],cwd=cwd,stdout=so,stderr=se)
    try:
     while proc.poll() is None:
      rss=psutil.Process().memory_info().rss
      try:rss+=psutil.Process(proc.pid).memory_info().rss
      except psutil.NoSuchProcess:pass
      peak=max(peak,rss)
      if rss>8*1024**3 or time.monotonic()-began>180 or time.monotonic()-start>7200:proc.kill();proc.wait();raise TimeoutError('SOLVER_OR_HARD_BUDGET_TIMEOUT')
      guard(start);time.sleep(.15)
     rc=proc.returncode
    except BaseException:
     if proc.poll() is None:proc.kill();proc.wait()
     raise
   row.update(returncode=rc,status='SUCCESS' if rc==0 else 'NONZERO',wall_s=time.monotonic()-began,stdout_sha256=sha(stdout),stderr_sha256=sha(stderr))
   if rc!=0:raise RuntimeError('SOLVER_NONZERO')
   modp=cwd/(item['stem']+'.mod');prtp=modp.with_suffix('.prt');mod=parse_mod(modp);depths,phi,k,f=mod;assert f==item['f'];d=env_data(cwd/env.name);prt=prtp.read_text(encoding='utf-8-sig');row.update(MOD_sha256=sha(modp),PRT_sha256=sha(prtp),modes=len(k))
   all_depths=np.unique(np.r_[d['source'],d['receiver']]);depterr=max(float(np.min(abs(depths-z))) for z in all_depths)
   for z in all_depths:zi(depths,z)
   if item['kind']=='ANALYTIC':
    contract=(d['option']=='CVW' and 'THORP' not in prt.upper() and np.all(d['ssp'][:,4:]==0) and float(np.max(abs(k.imag)))<=1e-12 and ('Perfectly RIGID' in prt if item['boundary']=='PR_R' else prt.count('VACUUM')>=2))
    if not contract:raise RuntimeError('ANALYTIC_ENV_NOT_MATCHED')
    mu,ka,cp=spectrum(f,item['boundary'],item['high']);n=min(len(ka),len(k));exactcount=len(ka)==len(k);krerr=abs(k.real[:n]-ka[:n])/ka[:n];krmax=float(krerr.max())
    theoretical=.02*np.sin(depths[:,None]*mu[None,:]);phase=[]
    for j in range(max(len(k),len(ka))):
     if j>=n:modes.append({'stem':item['stem'],'mode':j+1,'matched':False,'status':'UNPAIRED_NOT_REMOVED'});continue
     a=theoretical[:,j];b=phi[:,j];factor=np.vdot(a,b)/np.vdot(a,a);dif=abs(b-factor*a);small=abs(a)<1e-6*np.max(abs(a))
     modes.append({'stem':item['stem'],'mode':j+1,'matched':True,'kr_analytic':float(ka[j]),'kr_returned_Re':float(k[j].real),'kr_returned_Im':float(k[j].imag),'kr_relative_error':float(krerr[j]),'cp_analytic':float(cp[j]),'cp_returned':float(2*np.pi*f/k[j].real),'CHIGH_boundary_margin_analytic_mps':float(item['high']-cp[j]),'CLOW_boundary_margin_analytic_mps':float(cp[j]-1500),'depth_samples':len(depths),'near_node_samples':int(small.sum()),'near_node_max_abs_error':float(dif[small].max()) if np.any(small) else 'NO_NEAR_NODE_SAMPLES','modal_function_max_abs_error':float(dif.max()),'modal_function_shape_diagnostic':shape(a,b),'single_mode_factor_Re':float(factor.real),'single_mode_factor_Im':float(factor.imag),'factor_scope':'FUNCTION_DIAGNOSTIC_ONLY_NOT_USED_IN_FIELD'})
    kvals=[];avals=[];passed=True;callfields={}
    for probe,r,t,z in PROBES:
     for res,coords in g.items():
      pos=(r,t,z);a=analytic_field(f,item['boundary'],item['high'],pos,coords);afar=analytic_field(f,item['boundary'],item['high'],pos,coords,False);b=modal_field(mod,pos,coords);err=shape(a,b);radial=shape(a,afar);accept=exactcount and krmax<=2e-5 and err<=.002 and radial<=1e-4 and depterr<=1e-5;passed&=accept
      name=probe+'_'+res;callfields[name]=b;kvals.append(b);avals.append(a)
      cross.append({'stem':item['stem'],'boundary':item['boundary'],'frequency_hz':f,'CHIGH':item['high'],'mesh':item['mesh'],'probe':probe,'resource':res,'r_km':r,'abs_theta_deg':t,'source_z_m':z,'modes_analytic':len(ka),'modes_returned':len(k),'mode_count_PASS':exactcount,'max_kr_relative_error':krmax,'max_depth_error_m':depterr,'max_abs_Im_kr':float(np.max(abs(k.imag))),'spatial_shape_error':err,'exact_analytic_vs_farfield_shape':radial,'raw_relative_complex_error':float(np.linalg.norm(b-a)/np.linalg.norm(a)),'analytic_field_norm':float(np.linalg.norm(a)),'kraken_field_norm':float(np.linalg.norm(b)),'PASS':accept})
    aa=np.concatenate(avals);bb=np.concatenate(kvals);factor=np.vdot(aa,bb)/np.vdot(aa,aa);globalrows.append({'stem':item['stem'],'single_global_factor_Re':float(factor.real),'single_global_factor_Im':float(factor.imag),'profiled_raw_relative':float(np.linalg.norm(bb-factor*aa)/np.linalg.norm(bb)),'unprofiled_raw_relative':float(np.linalg.norm(bb-aa)/np.linalg.norm(bb)),'complex_samples':len(aa),'factor_scope':'ONE_FOR_ALL_FIVE_PROBES_FOUR_RESOURCES_NOT_PER_MODE_OR_CHANNEL'})
    np.savez(cwd/'FIELDS.npz',**callfields);fields[item['boundary'],f,item['high'],item['mesh']]=callfields
    if item['mesh']==80001:
     other=fields[item['boundary'],f,item['high'],160001]
     for name,b in callfields.items():meshrows.append({'boundary':item['boundary'],'frequency_hz':f,'CHIGH':item['high'],'probe_resource':name,'shape_difference':shape(other[name],b),'raw_relative_difference':float(np.linalg.norm(b-other[name])/np.linalg.norm(other[name])),'norm80001':float(np.linalg.norm(b)),'norm160001':float(np.linalg.norm(other[name])),'diagnostic_only_no_posthoc_gate':True})
    table(OUT/'ANALYTIC_KRAKEN_MODE_AND_FIELD_CROSSCHECK.csv',cross);table(OUT/'ANALYTIC_MODAL_ERRORS.csv',modes);table(OUT/'ANALYTIC_RAW_GLOBAL_FACTOR.csv',globalrows);table(OUT/'ANALYTIC_MESH_COMPARISON.csv',meshrows,fields=None if meshrows else ['boundary','frequency_hz','probe_resource','status'])
    row['analytic_gate_PASS']=passed;table(OUT/'SOLVER_CALLS.csv',calls)
    print('CALL',item['call_order'],item['stem'],'modes',len(k),'/',len(ka),'krmax',krmax,'shape_max',max(r['spatial_shape_error'] for r in cross if r['stem']==item['stem']),'PASS',passed,flush=True)
    if not passed:raise RuntimeError('ANALYTIC_HARD_GATE_FAIL')
    analytical_PASS+=1
   else:
    assert d['option']=='CVWT' and d['bottom']=='V' and 'THORP attenuation added' in prt and prt.count('VACUUM')>=2
    old=parse_mod(mod_path('K2',f));extfields={}
    for sn in sentinels():
     pos=tuple(float(sn[x]) for x in ('r_km','abs_theta_deg','source_z_m'))
     for res,coords in g.items():
      rfield=modal_field(old,pos,coords);vfield=modal_field(mod,pos,coords);err=shape(rfield,vfield);extfields[sn['sentinel_id']+'_'+res+'_R']=rfield;extfields[sn['sentinel_id']+'_'+res+'_V']=vfield
      extremes.append({'frequency_hz':f,'sentinel_id':sn['sentinel_id'],'resource':res,'r_km':pos[0],'abs_theta_deg':pos[1],'source_z_m':pos[2],'R_modes':len(old[2]),'V_modes':len(k),'R_cp_max':float((2*np.pi*f/old[2].real).max()),'V_cp_max':float((2*np.pi*f/k.real).max()),'spatial_shape_difference':err,'raw_relative_complex_difference':float(np.linalg.norm(vfield-rfield)/np.linalg.norm(rfield)),'R_field_norm':float(np.linalg.norm(rfield)),'V_field_norm':float(np.linalg.norm(vfield)),'scope':'R_V_MATHEMATICAL_BOUNDARY_SENSITIVITY_NOT_REAL_SEABED'})
    np.savez(cwd/'FIELDS.npz',**extfields);table(OUT/'E_STD_BOTTOM_EXTREMES.csv',extremes)
    print('CALL',item['call_order'],item['stem'],'R_V maxshape',max(r['spatial_shape_difference'] for r in extremes if r['frequency_hz']==f),flush=True)
   row['budget_after']=json.dumps(guard(start));table(OUT/'SOLVER_CALLS.csv',calls)
  except Exception as e:
   if rc is None:row.update(returncode=-999,status='TIMEOUT_OR_EXCEPTION',wall_s=time.monotonic()-began)
   row['failure_reason']=str(e);failures.append({'call_order':item['call_order'],'stem':item['stem'],'reason':str(e),'future_calls':'CANCELLED_NO_RETRY'});(cwd/'exception.txt').write_text(traceback.format_exc(),encoding='utf-8');table(OUT/'SOLVER_CALLS.csv',calls);print('STOP',str(e),flush=True);break
 for name,rr in [('ANALYTIC_KRAKEN_MODE_AND_FIELD_CROSSCHECK.csv',cross),('ANALYTIC_MODAL_ERRORS.csv',modes),('ANALYTIC_MESH_COMPARISON.csv',meshrows),('E_STD_BOTTOM_EXTREMES.csv',extremes),('ANALYTIC_RAW_GLOBAL_FACTOR.csv',globalrows)]:
  if not (OUT/name).exists():table(OUT/name,rr,fields=['status'])
 dump(OUT/'RUN_FAILURES.json',failures)
 intended=[]
 for item in plan:intended.append({'call_order':item['call_order'],'stem':item['stem'],'kind':item['kind'],'status':next((r['status'] for r in calls if r['stem']==item['stem']),'NOT_EVALUATED_STOPPED_EARLIER')})
 table(OUT/'ALL_PLANNED_CALL_STATUS.csv',intended)
 raw={str(p):sha(p) for p in (LOCAL/'runs').rglob('*') if p.is_file()}
 dump(OUT/'SCIENCE_OUTPUT_FREEZE.json',{'design_SHA':design_sha,'raw_files':raw,'scientific_tables':{str(OUT/n):sha(OUT/n) for n in ['SOLVER_CALLS.csv','ANALYTIC_KRAKEN_MODE_AND_FIELD_CROSSCHECK.csv','ANALYTIC_MODAL_ERRORS.csv','ANALYTIC_MESH_COMPARISON.csv','E_STD_BOTTOM_EXTREMES.csv','ANALYTIC_RAW_GLOBAL_FACTOR.csv','RUN_FAILURES.json','ALL_PLANNED_CALL_STATUS.csv']},'analytic_configs_PASS':analytical_PASS,'calls_attempted':len(calls),'E_STD_calls':sum(r['kind']=='E_STD_V' for r in calls),'budget':guard(start),'peak_process_tree_RSS_bytes':peak,'new_SCOOTER_FIELDS_FIELD_waveforms_DFT_MC_FIM_optimizer':0})
 print('SCIENCE SEALED',len(calls),'calls',analytical_PASS,'analytic passes',flush=True)
if __name__=='__main__':main()