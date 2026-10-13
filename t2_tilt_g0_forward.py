"""Once-only forward solving and full-registered-catalogue cross-mesh admission."""
from t2_tilt_g0_core import *
import subprocess,shutil,sys

def field_shd(p):
 b=Path(p).read_bytes();rec=4*int(np.frombuffer(b[:4],'<i4')[0]);hdr=np.frombuffer(b[2*rec:3*rec],'<i4');ns,nrd,nr=map(int,hdr[4:7]);assert ns==1
 depths=np.frombuffer(b[8*rec:9*rec],'<f4')[:nrd].astype(float);rr=np.frombuffer(b[9*rec:10*rec],'<f4')[:nr].astype(float)
 value=np.stack([np.frombuffer(b[(10+j)*rec:(10+j)*rec+8*nr],'<c8').astype(complex) for j in range(nrd)])
 return depths,rr,value

def main():
 start=time.monotonic();g=geometries();prov=read(OUT/'PROVIDER_PROVENANCE.json');expected=prov['expected_output_depths_m'];calls=[];modes=[];allcontrols=[]
 for item in prov['plans']:
  guard(start);stem=item['stem'];env=OUT/'environments'/f'{stem}.env';assert sha(env)==item['env_SHA256'];shutil.copyfile(env,LOCAL/'modal'/env.name)
  began=time.monotonic();p=subprocess.run([str(AT/'kraken.exe'),stem],cwd=LOCAL/'modal',capture_output=True,timeout=90);(LOCAL/'modal'/f'{stem}.stdout.txt').write_bytes(p.stdout+p.stderr)
  calls.append({'solver':'KRAKEN','stem':stem,'mesh':item['mesh'],'frequency_hz':item['frequency_hz'],'returncode':p.returncode,'wall_s':time.monotonic()-began});table(OUT/'SOLVER_CALLS.csv',calls)
  if p.returncode!=0:raise RuntimeError('KRAKEN nonzero; no rerun')
  z,phi,kk=parse_mod(LOCAL/'modal'/f'{stem}.mod',expected);cp=2*np.pi*item['frequency_hz']/kk.real
  modes.append({'stem':stem,'mesh':item['mesh'],'frequency_hz':item['frequency_hz'],'modes':len(kk),'min_phase_speed_mps':float(cp.min()),'max_phase_speed_mps':float(cp.max()),'CLOW_margin_mps':float(cp.min()-1500),'CHIGH_margin_mps':float(1800-cp.max()),'CHIGH_certification':'UNCERTIFIED','max_output_depth_match_m':max(float(np.min(abs(z-d))) for d in expected),'MOD_SHA256':sha(LOCAL/'modal'/f'{stem}.mod'),'ENV_SHA256':sha(env),'all_modes_used':True})
  print('KRAKEN',len(calls),'/6 modes',len(kk),flush=True)
 table(OUT/'PROVIDER_MODE_WINDOW_AND_MESH.csv',modes)
 high=Provider(160001,expected=expected);low=Provider(80001,expected=expected)
 fields=[]
 for fi,f in enumerate(FREQ):
  guard(start);stem=f't2_field_f{int(f)}';base=LOCAL/'modal'/f't2_f{int(f)}_n160001';dest=LOCAL/'modal'/stem
  shutil.copyfile(base.with_suffix('.mod'),dest.with_suffix('.mod'));shutil.copyfile(base.with_suffix('.env'),dest.with_suffix('.env'));shutil.copyfile(OUT/'environments'/f'{stem}.flp',dest.with_suffix('.flp'))
  began=time.monotonic();p=subprocess.run([str(AT/'field.exe'),stem],cwd=LOCAL/'modal',capture_output=True,timeout=60);dest.with_suffix('.stdout.txt').write_bytes(p.stdout+p.stderr);calls.append({'solver':'FIELD','stem':stem,'mesh':160001,'frequency_hz':float(f),'returncode':p.returncode,'wall_s':time.monotonic()-began});table(OUT/'SOLVER_CALLS.csv',calls)
  if p.returncode:raise RuntimeError('FIELD_NOT_READY; no rerun')
  zs,rr,val=field_shd(dest.with_suffix('.shd'));pred=[];obs=[];pointmeta=[]
  for res in prov['FIELD_plan']['resources']:
   a=g[res];ss=prov['FIELD_plan']['states'];hh=high.field(ss,a)[:,fi]
   for i,(r,t,z) in enumerate(ss):
    ranges=np.sqrt((r*1000*np.cos(np.radians(t))-a[:,0])**2+(r*1000*np.sin(np.radians(t))-a[:,1])**2)
    for j in range(len(a)):
     zi=int(np.argmin(abs(zs-a[j,2])));ri=int(np.argmin(abs(rr-ranges[j])))
     if abs(zs[zi]-a[j,2])>1e-5 or abs(rr[ri]-ranges[j])>.01:raise ValueError('FIELD coordinate mismatch')
     pred.append(hh[i,j]);obs.append(val[zi,ri]);pointmeta.append([res,i,j,float(zs[zi]),float(rr[ri])])
  pred=np.array(pred);obs=np.array(obs);c=np.vdot(pred,obs)/np.vdot(pred,pred);e=float(np.linalg.norm(obs-c*pred)/np.linalg.norm(obs))
  np.savez(OUT/f'FIELD_COMPARISON_f{int(f)}.npz',pred=pred,FIELD=obs,common_complex_factor=c,returned_depths=zs,returned_ranges=rr)
  table(OUT/f'FIELD_POINTS_f{int(f)}.csv',[{'resource':p[0],'state':p[1],'channel':p[2],'returned_z_m':p[3],'returned_range_m':p[4]} for p in pointmeta])
  fields.append({'control':'FIELD_shared_scalar_all3resources_all3ranges','frequency_hz':float(f),'error':e,'limit':.002,'PASS':e<=.002,'n_points':len(pred),'source_z_m':200,'one_factor_per_frequency':True})
  print('FIELD',int(f),'shape',e,flush=True)
 table(OUT/'FIELD_CONTROLS.csv',fields)
 # Both complete catalogues, no truth-local selection and no interpolation.
 fielddir=LOCAL/'replicas';fielddir.mkdir();num=[]
 for fine in [False,True]:
  name='fine' if fine else 'primary';s=catalogue(fine);np.save(OUT/f'CATALOGUE_{name}.npy',s)
  for res in RES:
   guard(start);hh=high.field(s,g[res]);ll=low.field(s,g[res]);np.save(fielddir/f'{res}_{name}_n160001.npy',hh);np.save(fielddir/f'{res}_{name}_n80001.npy',ll);ee=shape_error(ll,hh)
   np.save(LOCAL/f'SHAPE_ERROR_{res}_{name}.npy',ee)
   for fi,f in enumerate(FREQ):
    norm=np.linalg.norm(hh[:,fi],axis=-1);j=int(np.argmax(ee[:,fi]));num.append({'resource':res,'catalogue':name,'frequency_hz':float(f),'candidates':len(s),'max_profiled_shape_error':float(ee[:,fi].max()),'median_profiled_shape_error':float(np.median(ee[:,fi])),'max_error_candidate_id':j,'max_error_r_km':float(s[j,0]),'max_error_abs_theta_deg':float(s[j,1]),'max_error_z_m':float(s[j,2]),'min_vector_norm':float(norm.min()),'max_vector_norm':float(norm.max()),'limit':.002,'PASS':bool(np.max(ee[:,fi])<=.002)})
   print('MESH',name,res,'max',float(ee.max()),flush=True)
 table(OUT/'FORWARD_CATALOGUE_NUMERICAL_STABILITY.csv',num)
 # Separate calibration source positions, exact source depths, no thresholds/test data yet.
 cp=rows(OUT/'CALIBRATION_GEOMETRIES.csv');s=np.array([[float(p[k]) for k in ['r_km','theta_deg','z_m']] for p in cp]);cal=[]
 for res in RES:
  a=high.field(s,g[res]);b=low.field(s,g[res]);np.save(fielddir/f'CAL_{res}_n160001.npy',a);np.save(fielddir/f'CAL_{res}_n80001.npy',b);e=shape_error(b,a);cal.append({'control':'independent_calibration_geometry_crossmesh','resource':res,'error':float(e.max()),'limit':.002,'PASS':bool(e.max()<=.002)})
 s=np.array([[50.123,2.17,185],[55.321,4.37,200],[59.127,.53,215]]);sym=[]
 for p in [low,high]:
  for res in RES:
   a=p.field(s,g[res]);b=p.field(s*np.array([1,-1,1]),g[res]);e=float(np.max(abs(a-b)));sym.append({'control':'theta_sign_mirror_all_resources','resource':res,'mesh':p.mesh,'error':e,'limit':1e-12,'PASS':e<=1e-12})
  aa=p.field(s,g['A0_0']);bb=p.field(s,g['B1_V4']);dd=p.field(s,g['B2_DUP']);e=max(float(np.max(abs(aa-bb[:,:,:8]))),float(np.max(abs(aa-dd[:,:,:8]))),float(np.max(abs(dd[:,:,8:]-aa[:,:,3,None]))));sym.append({'control':'shared_original8_and_DUP_pressure','resource':'B1_B2_A0','mesh':p.mesh,'error':e,'limit':1e-12,'PASS':e<=1e-12})
 allcontrols=fields+cal+sym;table(OUT/'FORWARD_CONTROLS.csv',allcontrols)
 passed=all(x['PASS'] for x in allcontrols+num)
 dump(OUT/'FORWARD_DECISION.json',{'utc':utc(),'admitted':passed,'status':'FORWARD_NUMERIC_ADMITTED_CONDITIONALLY' if passed else 'FORWARD_NUMERIC_NOT_ADMITTED','CHIGH':'CONDITIONAL_MODAL_WINDOW_LIMIT','new_KRAKEN':6,'new_FIELD':3,'max_cross_mesh_profiled_error':max(x['max_profiled_shape_error'] for x in num),'max_FIELD_profiled_error':max(x['error'] for x in fields),'failed_catalogue_frequency_controls':sum(not x['PASS'] for x in num),'failed_other_controls':sum(not x['PASS'] for x in allcontrols),'full_catalogue_states_checked':sum(x['candidates'] for x in num),'interpolation':'NONE_DIRECT_ALL_MODE_SUM','calibration_test_MC_allowed':passed,'science_wall_s':time.monotonic()-start,'budget':guard(start)})
 dump(OUT/'REPLICA_FREEZE.json',{'utc':utc(),'files':{p.name:sha(p) for p in fielddir.glob('*.npy')},'complete_raw_shape_differences':{p.name:sha(p) for p in LOCAL.glob('SHAPE_ERROR*.npy')},'main_data_not_generated_yet':True})
 print('FORWARD ADMISSION',passed,flush=True)
if __name__=='__main__':main()