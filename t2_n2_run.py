"""One bounded N2 science execution, KRAKEN first, independent SCOOTER+FIELDS second."""
from t2_n2_common import *
import subprocess,shutil,traceback

def main():
 if (OUT/'SCIENCE_STARTED.json').exists():raise RuntimeError('SCIENCE_ONCE_ONLY')
 verify_sources();start=time.monotonic();dump(OUT/'SCIENCE_STARTED.json',{'utc':utc(),'once_only':True,'new_signal_MC_FIM':0});plan=read(OUT/'SOLVER_PLAN.json');calls=[];modal=[];g=geometries();kruns={};ffpruns={};failure=[]
 def call(tool,stem,cwd,timeout,extra):
  storage_guard(start);binary=AT/(tool+'.exe');expected=next(r['sha256'] for r in read(OUT/'SOURCE_IDENTITY.json')['records'] if r.get('name')==tool);assert sha(binary)==expected
  began=time.monotonic();stdout=cwd/(stem+'.stdout.txt');stderr=cwd/(stem+'.stderr.txt');rc=None;status='EXECUTED'
  try:
   with stdout.open('wb') as so,stderr.open('wb') as se:proc=subprocess.run([str(binary),stem],cwd=cwd,stdout=so,stderr=se,timeout=timeout);rc=proc.returncode
  except subprocess.TimeoutExpired:status='TIMEOUT';rc=-999
  row={'call_id':len(calls)+1,'tool':tool,'stem':stem,'cwd':str(cwd),'returncode':rc,'status':status,'timeout_s':timeout,'wall_s':time.monotonic()-began,'binary_sha256':expected,'stdout_sha256':sha(stdout),'stderr_sha256':sha(stderr)}|extra;calls.append(row);table(OUT/'PROVIDER_CALLS.csv',calls);print('CALL',row['call_id'],tool,stem,'rc',rc,'s',row['wall_s'],flush=True)
  if rc!=0:raise RuntimeError(tool+'_NONZERO_OR_TIMEOUT')
  return row
 def mode_record(v,f):
  z,phi,k=parse_mod(mod_path(v,f));cp=2*np.pi*f/k.real;mesh,high=KVAR[v];old=parse_mod(mod_path('K0',f));oldcp=2*np.pi*f/old[2].real;near=np.argmin(abs(cp[:,None]-oldcp[None,:]),axis=0);diff=abs(cp[near]-oldcp);depthmatch=np.array([zi(z,d) for d in old[0]]);newphi=phi[depthmatch][:,near].T;oldphi=old[1].T;valid=(np.linalg.norm(newphi,axis=-1)>0)&(np.linalg.norm(oldphi,axis=-1)>0);aligned=shape(newphi[valid],oldphi[valid]) if np.any(valid) else np.array([])
  modal.append({'variant':v,'frequency_hz':f,'mesh':mesh,'CLOW':1500,'CHIGH_requested':high,'modes':len(k),'mode_delta_from_K0':len(k)-len(old[2]),'cp_min':float(cp.min()),'cp_max':float(cp.max()),'lower_boundary_margin':float(cp.min()-1500),'upper_boundary_margin':float(high-cp.max()),'old_modes_nearest_cp_max_difference':float(diff.max()),'old_modes_nearest_phi_max_shape_difference':float(aligned.max()) if len(aligned) else 'NOT_EVALUATED_ZERO_MODE_VECTORS','matching_nonzero_mode_vectors':int(valid.sum()),'matching_zero_stored_mode_vectors':int((~valid).sum()),'nearest_matching_unique_new_modes':len(np.unique(near)),'all_returned_modes_used_in_main_field':True,'MOD_sha256':sha(mod_path(v,f)),'max_registered_depth_error_m':max(float(np.min(abs(z-d))) for d in np.unique(np.concatenate([np.arange(180,221,10)]+[a[:,2] for a in g.values()])))});table(OUT/'MODE_WINDOW_AND_MATCHING.csv',modal)
 for f in FREQ:mode_record('K0',f)
 # K0 never solved again; all newly registered K configurations performed irrespective of numerical gates.
 for item in [i for i in plan['primary'] if i['family']=='KRAKEN']:
  variant=item['variant'];f=item['frequency_hz'];cwd=LOCAL/'K'/variant;cwd.mkdir(parents=True,exist_ok=True);env=Path(item['env_path']);assert sha(env)==item['env_sha256'];shutil.copyfile(env,cwd/env.name)
  try:
   call('kraken',item['stem'],cwd,item['timeout_s'],{'family':'KRAKEN','variant':variant,'frequency_hz':f,'input_sha256':item['env_sha256']});mode_record(variant,f);kruns[variant,f]='READY'
  except Exception as e:
   failure.append({'path':'KRAKEN','variant':variant,'frequency_hz':f,'reason':str(e),'action':'STOP_REMAINING_K_NO_RETRY'});dump(OUT/'RUN_FAILURES.json',failure);(cwd/(item['stem']+'.exception.txt')).write_text(traceback.format_exc(),encoding='utf-8');break
 fieldsdir=LOCAL/'K_FIELDS';fieldsdir.mkdir();bands=[];fullsummary=[];sentsummary=[];kready=['K0'];s=sstate=sstates()
 for v in ['K1','K2','K3']:
  if all(kruns.get((v,f))=='READY' for f in FREQ):kready.append(v)
 for v in kready:
  if v=='K0':continue
  provider=Provider(v)
  for cat in CATS:
   ss=states(cat)
   for res in RES:
    storage_guard(start);h,bh=provider.field(ss,g[res],1800 if v=='K1' else 2000 if v=='K2' else None);np.save(fieldsdir/f'{v}_{res}_{cat}.npy',h)
    if v in ['K1','K2']:
     np.save(fieldsdir/f'BAND_{v}_{res}_{cat}.npy',bh)
     for fi,f in enumerate(FREQ):
      ratio=np.linalg.norm(bh[:,fi],axis=-1)/np.linalg.norm(h[:,fi],axis=-1);bands.append({'variant':v,'resource':res,'catalogue':cat,'frequency_hz':f,'cp_band_lower_mps':1800 if v=='K1' else 2000,'full_sum_used_for_gate':True}|quant(ratio))
    print('ALL-MODE FIELD',v,res,cat,len(ss),'states',flush=True)
 table(OUT/'WINDOW_MODAL_CONTRIBUTION.csv',bands,fields=list(bands[0]) if bands else ['variant','status'])
 for left,right,label in [('K0','K1','1800_TO_2000'),('K1','K2','2000_TO_2400'),('K2','K3','MESH160001_TO240001')]:
  if left not in kready or right not in kready:continue
  for cat in CATS:
   ss=states(cat)
   for res in RES:
    a=np.load(OLDLOCAL/'replicas'/f'{res}_{cat}_n160001.npy',mmap_mode='r') if left=='K0' else np.load(fieldsdir/f'{left}_{res}_{cat}.npy',mmap_mode='r');b=np.load(fieldsdir/f'{right}_{res}_{cat}.npy',mmap_mode='r');ee=shape(a,b);raw=np.linalg.norm(b-a,axis=-1)/np.linalg.norm(b,axis=-1);np.savez(LOCAL/f'K_DIFFERENCE_{left}_{right}_{res}_{cat}.npz',shape=ee,raw=raw,norm_left=np.linalg.norm(a,axis=-1),norm_right=np.linalg.norm(b,axis=-1))
    for fi,f in enumerate(FREQ):
     cid=int(np.argmax(ee[:,fi]));fullsummary.append({'comparison':label,'left':left,'right':right,'resource':res,'catalogue':cat,'frequency_hz':f}|quant(ee[:,fi])|{'above_0p002':int(np.sum(ee[:,fi]>.002)),'raw_relative_max':float(raw[:,fi].max()),'raw_relative_median':float(np.median(raw[:,fi])),'worst_candidate_id':cid,'worst_r_km':float(ss[cid,0]),'worst_abs_theta_deg':float(ss[cid,1]),'worst_z_m':float(ss[cid,2]),'min_right_field_norm':float(np.linalg.norm(b[:,fi],axis=-1).min()),'PASS':bool(np.max(ee[:,fi])<=.002)})
    if cat=='fine':
     for sn in sentinels():
      coord=[float(sn[k]) for k in ['r_km','abs_theta_deg','source_z_m']];cid=int(np.flatnonzero(np.all(ss==coord,axis=1))[0])
      for fi,f in enumerate(FREQ):sentsummary.append({'comparison':label,'resource':res,'sentinel_id':sn['sentinel_id'],'frequency_hz':f,'candidate_id_fine':cid,'shape_error':float(ee[cid,fi]),'raw_relative':float(raw[cid,fi]),'r_km':coord[0],'abs_theta_deg':coord[1],'source_z_m':coord[2]})
 table(OUT/'FULL_CATALOGUE_KRAKEN_STABILITY.csv',fullsummary,fields=list(fullsummary[0]) if fullsummary else ['comparison','status']);table(OUT/'SENTINEL_KRAKEN_STABILITY.csv',sentsummary,fields=list(sentsummary[0]) if sentsummary else ['comparison','status'])
 # FFP path uses only SCOOTER GRN and FIELDS, never KRAKEN FIELD.
 ffpstates=[];ffp_arrays={};interface=read(OUT/'EXECUTION_CONTRACT.json')['static_FFP_ready'];aborted=False
 if interface:
  for item in [i for i in plan['primary'] if i['family']=='SCOOTER']:
   v=item['variant'];f=item['frequency_hz'];cwd=LOCAL/'FFP'/v/f'f{f}';cwd.mkdir(parents=True,exist_ok=True);env=Path(item['env_path']);assert sha(env)==item['env_sha256'];shutil.copyfile(env,cwd/env.name)
   try:
    call('scooter',item['stem'],cwd,item['timeout_s'],{'family':'SCOOTER','variant':v,'frequency_hz':f,'input_sha256':item['env_sha256']});grn=cwd/(item['stem']+'.grn');d=parse_shd(grn,green=True);assert abs(d['freq']-f)<=1e-6
    for depth in [180,190,200,210,220]:zi(d['source_depths'],depth)
    for a in g.values():
     for depth in a[:,2]:zi(d['receiver_depths'],depth)
    if not np.all(d['axis']>0):raise ValueError('GRN_PHASE_SPEED_AXIS_INVALID')
    grn_hash=sha(grn);ffpruns[v,f]='GRN_READY'
    for fp in [p for p in plan['FIELDS'] if p['variant']==v and p['frequency_hz']==f]:
     flp=Path(fp['flp_path']);assert sha(flp)==fp['flp_sha256'];work=cwd/fp['stem'];work.mkdir();shutil.copyfile(grn,work/(fp['stem']+'.grn'));shutil.copyfile(flp,work/flp.name)
     call('fields',fp['stem'],work,fp['timeout_s'],{'family':'FIELDS','variant':v,'frequency_hz':f,'sentinel_id':fp['sentinel_id'],'range_step_m':fp['range_step_m'],'input_sha256':fp['flp_sha256'],'GRN_source_sha256':grn_hash})
     shd=work/(fp['stem']+'.shd');dsh=parse_shd(shd);sn=next(t for t in sentinels() if t['sentinel_id']==fp['sentinel_id']);state=dict(sn,frequency_hz=f)
     for res in RES:
      h=interpolate_shd(dsh,state,g[res]);ffp_arrays[v,f,fp['sentinel_id'],fp['range_step_m'],res]=h
     np.savez(work/'INTERPOLATED_FIELD.npz',**{res:ffp_arrays[v,f,fp['sentinel_id'],fp['range_step_m'],res] for res in RES})
     requested=np.linspace(fp['range_min_km']*1000,fp['range_max_km']*1000,fp['request_count']);ar=dsh['axis'];quantization=float(np.max(np.min(abs(ar[:,None]-requested[None,:]),axis=1)))
     ffpstates.append({'variant':v,'frequency_hz':f,'sentinel_id':fp['sentinel_id'],'range_step_m':fp['range_step_m'],'SHD_path':str(shd),'SHD_sha256':sha(shd),'GRN_source_path':str(grn),'GRN_source_sha256':grn_hash,'FLP_path':str(flp),'FLP_sha256':fp['flp_sha256'],'actual_source_depths':json.dumps(dsh['source_depths'].tolist()),'actual_receiver_depths':json.dumps(dsh['receiver_depths'].tolist()),'actual_NR':len(ar),'actual_range_min_m':float(ar[0]),'actual_range_max_m':float(ar[-1]),'actual_spacing_min_m':float(np.diff(ar).min()),'actual_spacing_max_m':float(np.diff(ar).max()),'axis_distance_to_requested_grid_max_m':quantization,'depth_match_error_m':max(float(np.min(abs(dsh['receiver_depths']-depth))) for a in g.values() for depth in a[:,2]),'receiver_interpolation':'ACTUAL_SAVED_AXIS_LINEAR_NO_EXTRAPOLATION','fields_EXE_not_field_EXE':True});table(OUT/'SCOOTER_FIELDS_PROVENANCE.csv',ffpstates)
    ffpruns[v,f]='READY'
   except Exception as e:
    failure.append({'path':'SCOOTER_FIELDS','variant':v,'frequency_hz':f,'reason':str(e),'action':'STOP_FFP_PATH_NO_RETRY'});dump(OUT/'RUN_FAILURES.json',failure);(cwd/(item['stem']+'.exception.txt')).write_text(traceback.format_exc(),encoding='utf-8');aborted=True;break
 # Store full intended sentinel statuses, missing fields explicitly NOT_EVALUATED.
 ffpfolder=LOCAL/'FFP_FIELDS';ffpfolder.mkdir();statuses=[]
 for v in FVAR:
  for step in ([.05,.025] if v=='F2' else [.05]):
   for res in RES:
    arr=np.full((12,3,len(g[res])),np.nan+1j*np.nan);mask=np.zeros((12,3),bool)
    for i,sn in enumerate(sentinels()):
     for fi,f in enumerate(FREQ):
      key=(v,f,sn['sentinel_id'],step,res)
      if key in ffp_arrays:arr[i,fi]=ffp_arrays[key];mask[i,fi]=True
      statuses.append({'variant':v,'resource':res,'sentinel_id':sn['sentinel_id'],'frequency_hz':f,'range_step_m':step,'status':'COMPLETED' if mask[i,fi] else 'NOT_EVALUATED'})
    np.savez(ffpfolder/f'{v}_{res}_step{step:.3f}.npz',field=arr,completed=mask)
 table(OUT/'FFP_SENTINEL_EXECUTION_STATUS.csv',statuses)
 if not ffpstates:table(OUT/'SCOOTER_FIELDS_PROVENANCE.csv',[],fields=['variant','frequency_hz','sentinel_id','status'])
 dump(OUT/'SCIENCE_OUTPUT_FREEZE.json',{'utc':utc(),'K_complete_variants':kready,'FFP_runtime_completed':not aborted and interface and all(ffpruns.get((v,f))=='READY' for v in FVAR for f in FREQ),'FFP_static_interface_ready':interface,'calls':len(calls),'failure_count':len(failure),'raw_files':{str(p):sha(p) for p in LOCAL.rglob('*') if p.is_file() and p.suffix.lower() in ['.env','.mod','.grn','.shd','.flp','.npy','.npz','.prt','.txt']},'budget':storage_guard(start),'no_new_observations':True});print('SCIENCE ONCE COMPLETE; FFP records',len(ffpstates),'K variants',kready,flush=True)
if __name__=='__main__':main()