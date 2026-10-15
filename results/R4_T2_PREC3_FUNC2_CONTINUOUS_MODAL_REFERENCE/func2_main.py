"""One preregistered FUNC2 analysis. Gate failures stop downstream science."""
from pathlib import Path
import csv,json,hashlib,traceback,math,time
import numpy as np
import func2_dop853 as dop
import func2_gauss as gl
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3_FUNC2');GRID=np.arange(20001)*.25

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,o):Path(p).write_text(json.dumps(o,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def csvout(name,data,default=('status','reason')):
 fields=list(dict.fromkeys(k for r in data for k in r)) if data else list(default)
 with (OUT/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(data)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def metric(a,b,z=None):
 weight=np.ones(len(a)) if z is None else np.r_[np.diff(z)[0]/2,(z[2:]-z[:-2])/2,np.diff(z)[-1]/2]
 an=float(np.dot(weight,abs(a)**2));bn=float(np.dot(weight,abs(b)**2));minus=float(np.dot(weight,abs(a-b)**2));plus=float(np.dot(weight,abs(a+b)**2));sg=1 if minus<=plus else -1;en=min(minus,plus)
 alpha=np.dot(weight,b.conj()*a)/bn if bn else None
 return {'e_sign':math.sqrt(en/bn) if bn else None,'e_complex':math.sqrt(float(np.dot(weight,abs(a-alpha*b)**2))/an) if alpha is not None and an else None,'sign':sg,'alpha_real':float(alpha.real) if alpha is not None else None,'alpha_imag':float(alpha.imag) if alpha is not None else None,'reference_energy':bn,'error_energy':en,'candidate_energy':an}
def contract(cfg):
 ls=(ROOT/cfg['env_path']).read_text(encoding='utf-8-sig').splitlines();end=next(i for i,l in enumerate(ls) if l.startswith("'R'"));a=np.array([list(map(float,l.split()[:6])) for l in ls[5:end]])
 assert a.shape==(251,6) and ls[3].strip()=="'CVWT'" and float(ls[1])==201 and np.all(a[:,3]==1) and np.all(a[:,2]==0) and np.all(a[:,4:]==0)
 return a[:,0],a[:,1]
def diagnostics(k,fn,info,mode):
 both=fn(GRID);p=both[0];der=both[1]*dop.S;peak=float(max(abs(p)));norm=np.trapz(p*p,GRID);bc=max(abs(p[0])/peak,abs(der[-1])/(dop.S*peak))
 return {'k_real':float(k.real),'k_imag':float(k.imag),'bilinear_norm_real':float(norm.real),'bilinear_norm_imag':float(norm.imag),'bilinear_norm_error':float(abs(norm-1)),'boundary_residual_dimensionless':float(bc),'match_residual_dimensionless':info['match_jump'],'sturm_internal_nodes':info.get('sturm_internal_nodes'),'phi_peak':peak},both

def main():
 cfg=read(OUT/'FUNC2_DESIGN_FREEZE.json');dest=LOCAL/'references';dest.mkdir(exist_ok=True)
 decision={'stage':'T2-PREC3-FUNC2','A':'NOT_EVALUATED','B':'NOT_EVALUATED','C':'NOT_EVALUATED','D':'NOT_EVALUATED','classifications':[],'R4_percent':0,'historical_PREC3_B':'COMPLEX_K_NOT_CERTIFIED_UNCHANGED','STOP':True,'completed':{}}
 for name in ['FUNC2_ANALYTIC_WAVEGUIDE_VALIDATION.csv','FUNC2_CONTINUOUS_COMPLEX_ROOTS.csv','FUNC2_CONTINUOUS_MODE_REFERENCE.csv','FUNC2_FE_FD_VS_REFERENCE.csv','FUNC2_MODAL_ERROR_BY_DEPTH.csv','FUNC2_LOCAL_COUPLING_REFERENCE.csv','FUNC2_CHANNEL_COUPLING_REFERENCE.csv','FUNC2_REFERENCE_LOCAL_STABILITY.csv','FUNC2_LOCAL_FUNCTION_REFERENCE.csv']:
  csvout(name,[])
 save(OUT/'FUNC2_DECISION.json',decision)
 traces=[];index=[];analytic=[];roots=[];references=[];fncache={};saved={};sources=rows(ROOT/cfg['probe_path'])[:5];coords=rows(ROOT/cfg['geometry_path']);queries=np.unique(np.r_[[float(q['source_z_m']) for q in sources],[float(q['z_m']) for q in coords if q['resource'] in cfg['resources']],np.arange(180,220.25,.25)])
 def stop(label,reason):
  decision['classifications'].append(label);decision['stop_reason']=reason;save(OUT/'FUNC2_DECISION.json',decision);save(OUT/'FUNC2_ROOT_TRACES.json',traces);save(OUT/'FUNC2_RAW_DATA_INDEX.json',{'files':index,'completed_only':True});print('SCIENCE_STOP',label,reason,flush=True)
 def store(tag,k,fn,info,mode,loss):
  dg,arr=diagnostics(k,fn,info,mode);p=dest/(tag+'.npz');query=fn(queries)
  np.savez_compressed(p,z=GRID,phi=arr[0],dphi=arr[1]*dop.S,z_query=queries,phi_query=query[0],dphi_query=query[1]*dop.S)
  index.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'tag':tag,'mode':mode,'loss':loss});saved[tag]=(dg,arr,query);return dg,arr
 try:
  for method,module in [('DOP853',dop),('GAUSS',gl)]:
   for level in [0,1]:
    for mode in cfg['modes']:
     zz=np.arange(251)*20.;speed=np.full(251,1500.);mu=(mode-.5)*math.pi/5000.;guess=float(np.sqrt((dop.S)**2-mu*mu))+5e-9
     def observe(loss,k,fn,info):
      if loss not in [0,1]:return
      tag=f'analytic_{method}_L{level}_m{mode}_a{loss:g}';dg,arr=store(tag,k,fn,info,mode,loss)
      cp=1500+1j*loss*dop.potential_data([0],[1500],1)[0].imag;exact=np.sqrt((math.tau*201/cp)**2-mu*mu);truth=np.sqrt(2/5000)*np.sin(mu*GRID);m=metric(arr[0],truth,GRID);ss=metric(arr[0][::80],truth[::80]);qtruth=np.sqrt(2/5000)*np.sin(mu*queries);actual=fn(queries)[0];qerror=float(max(abs(actual-m['sign']*qtruth)))
      ok=abs(k.real-exact.real)<=cfg['k_budget'] and abs(k.imag-exact.imag)<=cfg['k_budget'] and m['e_sign']<=cfg['shape_budget'] and dg['boundary_residual_dimensionless']<=cfg['boundary_budget'] and dg['match_residual_dimensionless']<=cfg['match_budget'] and dg['bilinear_norm_error']<=cfg['norm_budget'] and (loss!=0 or info['sturm_internal_nodes']==mode-1)
      analytic.append({'method':method,'level':level,'mode':mode,'loss':loss,'exact_k_real':float(exact.real),'exact_k_imag':float(exact.imag),'real_k_error':abs(k.real-exact.real),'imag_k_error':abs(k.imag-exact.imag),'dense_shape_error':m['e_sign'],'sparse251_shape_error':ss['e_sign'],'max_source_receiver_absolute_error':qerror,**dg,'PASS':bool(ok)})
      csvout('FUNC2_ANALYTIC_WAVEGUIDE_VALIDATION.csv',analytic)
     k,fn,h,info=module.solve(zz,speed,mode,guess,level,cfg['loss_scales'],observe);traces.append({'stage':'A','method':method,'level':level,'mode':mode,'continuation':h})
     print('A',method,level,mode,'rows',len(analytic),flush=True)
     if not all(r['PASS'] for r in analytic):
      decision['A']='FAIL';decision['completed']['analytic_rows']=len(analytic);stop('FUNC2_ANALYTIC_REFERENCE_GATE_FAIL','FROZEN_ANALYTIC_TOLERANCE');return
  # Both numerical paths must themselves agree when their accuracy is tightened.
  for method in ['DOP853','GAUSS']:
   for mode in cfg['modes']:
    for loss in [0,1]:
     a=saved[f'analytic_{method}_L0_m{mode}_a{loss:g}'];b=saved[f'analytic_{method}_L1_m{mode}_a{loss:g}'];e=metric(a[1][0],b[1][0],GRID)['e_sign']
     if e>cfg['shape_budget'] or abs(a[0]['k_real']-b[0]['k_real'])>cfg['k_budget'] or abs(a[0]['k_imag']-b[0]['k_imag'])>cfg['k_budget']:
      decision['A']='FAIL';stop('FUNC2_ANALYTIC_REFERENCE_GATE_FAIL','SELF_PRECISION');return
  decision['A']='PASS';decision['completed']['analytic_rows']=len(analytic);save(OUT/'FUNC2_DECISION.json',decision)
  z,c=contract(cfg);pr=cfg['initial_roots']
  for mode in cfg['modes']:
   for method,module in [('DOP853',dop),('GAUSS',gl)]:
    for level in [0,1]:
     def observe(loss,k,fn,info):
      tag=f'estd_{method}_L{level}_m{mode}_a{loss:g}'
      if loss in [0,1]:dg,arr=store(tag,k,fn,info,mode,loss)
      else:dg={'k_real':k.real,'k_imag':k.imag,**info}
      roots.append({'method':method,'level':level,'mode':mode,'loss':loss,**dg});csvout('FUNC2_CONTINUOUS_COMPLEX_ROOTS.csv',roots)
     k,fn,h,info=module.solve(z,c,mode,pr[str(mode)],level,cfg['loss_scales'],observe);fncache[method,level,mode]=fn;traces.append({'stage':'B','method':method,'level':level,'mode':mode,'continuation':h})
     final=saved[f'estd_{method}_L{level}_m{mode}_a1'][0];lossless=saved[f'estd_{method}_L{level}_m{mode}_a0'][0]
     ok=final['boundary_residual_dimensionless']<=cfg['boundary_budget'] and final['match_residual_dimensionless']<=cfg['match_budget'] and final['bilinear_norm_error']<=cfg['norm_budget'] and k.imag<0 and lossless['sturm_internal_nodes']==mode-1 and abs(lossless['k_real']-pr[str(mode)])<=cfg['k_budget'] and abs(lossless['k_imag'])<=1e-12
     if not ok:
      decision['B']='FAIL';stop('FUNC2_HIGH_ORDER_CONTINUOUS_REFERENCE_NOT_ESTABLISHED','IDENTITY_BOUNDARY_NORMALIZATION_OR_ROOT');return
     print('B',method,level,mode,repr(k),flush=True)
   fine=saved[f'estd_DOP853_L1_m{mode}_a1'];other=saved[f'estd_GAUSS_L1_m{mode}_a1'];dis=metric(fine[1][0],other[1][0],GRID)
   changes=[]
   for method in ['DOP853','GAUSS']:
    a=saved[f'estd_{method}_L0_m{mode}_a1'];b=saved[f'estd_{method}_L1_m{mode}_a1'];m=metric(a[1][0],b[1][0],GRID);changes.append({'method':method,'shape':m['e_sign'],'real':abs(a[0]['k_real']-b[0]['k_real']),'imag':abs(a[0]['k_imag']-b[0]['k_imag'])})
   real=abs(fine[0]['k_real']-other[0]['k_real']);imag=abs(fine[0]['k_imag']-other[0]['k_imag']);ok=real<=cfg['k_budget'] and imag<=cfg['k_budget'] and dis['e_sign']<=cfg['shape_budget'] and all(x['shape']<=cfg['shape_budget'] and x['real']<=cfg['k_budget'] and x['imag']<=cfg['k_budget'] for x in changes)
   references.append({'mode':mode,'k_path_real_difference':real,'k_path_imag_difference':imag,'function_path_difference':dis['e_sign'],'DOP_self_shape':changes[0]['shape'],'GAUSS_self_shape':changes[1]['shape'],'DOP_self_real_k':changes[0]['real'],'DOP_self_imag_k':changes[0]['imag'],'GAUSS_self_real_k':changes[1]['real'],'GAUSS_self_imag_k':changes[1]['imag'],'status':'CONTINUOUS_COMPLEX_MODE_REFERENCE_SUPPORTED_SUBSET' if ok else 'NOT_CERTIFIED','PASS':bool(ok)})
   csvout('FUNC2_CONTINUOUS_MODE_REFERENCE.csv',references)
   if not ok:
    decision['B']='FAIL';stop('FUNC2_HIGH_ORDER_CONTINUOUS_REFERENCE_NOT_ESTABLISHED','INDEPENDENT_OR_SELF_REFERENCE_DISAGREEMENT');return
  decision['B']='PASS_SUBSET';decision['completed']['reference_modes']=len(references);save(OUT/'FUNC2_DECISION.json',decision)
  comparisons=[];bands=[];couplings=[];channels=[];local=[];refq=[]
  for mode in cfg['modes']:
   fn=fncache['DOP853',1,mode];ref=fn(GRID)[0];sgcache={}
   for v in [v for v in cfg['vectors'] if v['mode']==mode]:
    with np.load(v['path'],allow_pickle=False) as a:vv=np.r_[0j,a['phi']];h=float(a['h'])
    def archived(zz):return np.interp(zz,np.arange(len(vv))*h,vv.real)+1j*np.interp(zz,np.arange(len(vv))*h,vv.imag)
    old=archived(GRID);met=metric(old,ref,GRID);coarse=metric(old[::2],ref[::2],GRID[::2]);sgcache[v['method'],v['N']]=met['e_sign']
    comparisons.append({'method':v['method'],'mode':mode,'N':v['N'],**met,'step0p5_e_sign':coarse['e_sign'],'quadrature_delta':abs(met['e_sign']-coarse['e_sign']),'historical_threshold0p002_met_descriptively':met['e_sign']<=.002,'precision_scope':'FINITE_CONTINUOUS_REFERENCE_INTERNAL_AGREEMENT_NOT_RIGOROUS_ERROR_BOUND'})
    sign=met['sign'];E=met['error_energy'];cumulative=0.
    for lo in range(0,5000,200):
     ix=(GRID>=lo)&(GRID<=lo+200);energy=float(np.trapz(abs(ref[ix])**2,GRID[ix]));err=float(np.trapz(abs(old[ix]-sign*ref[ix])**2,GRID[ix]));cumulative+=err
     bands.append({'method':v['method'],'mode':mode,'N':v['N'],'z_low_m':lo,'z_high_m':lo+200,'reference_energy':energy,'squared_error_integral':err,'fraction_total_error':err/E if E else None,'cumulative_fraction':cumulative/E if E else None,'global_sign':sign,'weak_local_norm':math.sqrt(energy)<=cfg['weak_local_L2_norm'],'reference_min_abs_phi':float(min(abs(ref[ix]))),'reference_sign_crossings_real':int(np.sum(ref[ix][:-1].real*ref[ix][1:].real<0))})
    zz=np.arange(180,220.25,.25);g=fn(zz)[0];f=archived(zz);lm=metric(f,g,zz);local.append({'method':v['method'],'mode':mode,'N':v['N'],**lm,'reference_min_abs':float(min(abs(g))),'reference_max_abs':float(max(abs(g))),'weak_reference':math.sqrt(lm['reference_energy'])<=cfg['weak_local_L2_norm']})
    for probe in sources:
     zs=float(probe['source_z_m']);fs=archived([zs])[0];gs=fn([zs])[0][0]
     for resource in cfg['resources']:
      cc=sorted([q for q in coords if q['resource']==resource],key=lambda q:int(q['channel']));zr=np.array([float(q['z_m']) for q in cc]);fr=archived(zr);gr=fn(zr)[0];qa=fs*fr;qb=gs*gr;mm=metric(qa,qb);nr=float(np.linalg.norm(qb));raw=float(np.linalg.norm(qa-qb)/nr) if nr else None
      weak=bool(abs(gs)<=cfg['weak_phi_amplitude'] or np.any(abs(gr)<=cfg['weak_phi_amplitude']) or nr<=cfg['weak_coupling_norm'])
      couplings.append({'method':v['method'],'mode':mode,'N':v['N'],'probe':probe['sentinel_id'],'resource':resource,'raw_relative_error':raw,'absolute_error_norm':float(np.linalg.norm(qa-qb)),'shape_error':mm['e_complex'],'reference_norm':nr,'weak':weak,'status':'WEAK_AMPLITUDE' if weak else 'UNDEFINED' if raw is None else 'FINITE','normalization':'ORIGINAL_ARCHIVE_NOT_RENORMALIZED'})
      for i,ch in enumerate(cc):channels.append({'method':v['method'],'mode':mode,'N':v['N'],'probe':probe['sentinel_id'],'resource':resource,'channel':ch['channel'],'source_depth':zs,'receiver_depth':zr[i],'reference_Q_real':float(qb[i].real),'reference_Q_imag':float(qb[i].imag),'archive_Q_real':float(qa[i].real),'archive_Q_imag':float(qa[i].imag),'amplitude_difference':float(abs(qa[i])-abs(qb[i])),'phase_difference_rad':float(np.angle(qa[i]*qb[i].conj())) if qa[i]!=0 and qb[i]!=0 else None})
   for probe in sources:
    zs=float(probe['source_z_m'])
    for resource in cfg['resources']:
     zr=np.array([float(q['z_m']) for q in sorted([q for q in coords if q['resource']==resource],key=lambda q:int(q['channel']))]);vectors=[fncache[meth,lev,mode]([zs])[0][0]*fncache[meth,lev,mode](zr)[0] for meth,lev in [('DOP853',1),('GAUSS',1),('DOP853',0),('GAUSS',0)]];q0=vectors[0];n=float(np.linalg.norm(q0));weak=n<=cfg['weak_coupling_norm'] or abs(fn([zs])[0][0])<=cfg['weak_phi_amplitude'] or np.any(abs(fn(zr)[0])<=cfg['weak_phi_amplitude'])
     for label,q,base in [('CROSS_PATH',vectors[1],q0),('DOP_SELF',vectors[2],q0),('GAUSS_SELF',vectors[3],vectors[1])]:
      den=float(np.linalg.norm(base));raw=float(np.linalg.norm(q-base)/den) if den else None;shape=metric(q,base)['e_complex'];ok=weak or (raw is not None and raw<=cfg['coupling_budget'] and shape<=cfg['coupling_budget']);refq.append({'mode':mode,'probe':probe['sentinel_id'],'resource':resource,'comparison':label,'raw_relative':raw,'shape':shape,'absolute_difference_norm':float(np.linalg.norm(q-base)),'reference_norm':den,'weak':bool(weak),'PASS_nonweak_only':bool(ok)})
   print('C_D',mode,'comparisons',len(comparisons),'couplings',len(couplings),flush=True)
  csvout('FUNC2_FE_FD_VS_REFERENCE.csv',comparisons);csvout('FUNC2_MODAL_ERROR_BY_DEPTH.csv',bands);csvout('FUNC2_LOCAL_COUPLING_REFERENCE.csv',couplings);csvout('FUNC2_CHANNEL_COUPLING_REFERENCE.csv',channels);csvout('FUNC2_LOCAL_FUNCTION_REFERENCE.csv',local);csvout('FUNC2_REFERENCE_LOCAL_STABILITY.csv',refq)
  decision['C']='COMPLETED_SUBSET';decision['D']='SUPPORTED_SUBSET_NONWEAK_ONLY' if all(q['PASS_nonweak_only'] for q in refq) else 'LOCAL_COUPLING_NOT_ESTABLISHED';decision['completed'].update(old_functions=len(comparisons),bands=len(bands),old_couplings=len(couplings),reference_couplings=len(refq))
  if any(r['e_sign']>.002 for r in comparisons if r['mode'] in [481,961]):decision['classifications'].append('FUNC2_REFERENCE_SUPPORTED_BUT_OLD_MODAL_FUNCTION_NOT_ADMITTED')
  if decision['D']=='SUPPORTED_SUBSET_NONWEAK_ONLY':decision['classifications'].append('FUNC2_LOCAL_MODAL_COUPLING_SUPPORTED_SUBSET_ONLY')
  else:decision['classifications'].append('FUNC2_LOCAL_COUPLING_NOT_ESTABLISHED')
  decision['FD961']='NOT_EVALUATED_FD_VECTOR_NOT_ARCHIVED';save(OUT/'FUNC2_ROOT_TRACES.json',traces);save(OUT/'FUNC2_RAW_DATA_INDEX.json',{'files':index});save(OUT/'FUNC2_DECISION.json',decision);print('SCIENCE_COMPLETE',json.dumps(decision),flush=True)
 except Exception as exc:
  (OUT/'FUNC2_SCIENCE_EXCEPTION.txt').write_text(traceback.format_exc(),encoding='utf-8')
  stop('FUNC2_ANALYTIC_REFERENCE_GATE_FAIL' if decision['A']!='PASS' else 'FUNC2_HIGH_ORDER_CONTINUOUS_REFERENCE_NOT_ESTABLISHED',str(exc))
  decision['A' if decision['A']!='PASS' else 'B']='FAIL';save(OUT/'FUNC2_DECISION.json',decision)
if __name__=='__main__':main()
