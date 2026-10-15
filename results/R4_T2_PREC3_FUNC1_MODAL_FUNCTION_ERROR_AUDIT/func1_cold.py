"""Independent FUNC1 reconstruction: explicit interpolation, interval sums and QR."""
from pathlib import Path
import json,csv,hashlib,math
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):
 if not Path(p).exists():return []
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def interpolate(values,h,z):
 z=np.asarray(z,float);u=z/h;i=np.minimum(np.floor(u).astype(int),len(values)-2);fraction=u-i
 return values[i]*(1-fraction)+values[i+1]*fraction
def nodeweights(z):
 z=np.asarray(z);a=np.empty(len(z));a[0]=(z[1]-z[0])/2;a[-1]=(z[-1]-z[-2])/2
 a[1:-1]=(z[2:]-z[:-2])/2;return a
def summation(a):return math.fsum(float(v) for v in a)
def metric(f,g,z,weighted=True):
 w=nodeweights(z) if weighted else np.ones(len(z));F=summation(w*(f.real*f.real+f.imag*f.imag));G=summation(w*(g.real*g.real+g.imag*g.imag))
 minus=summation(w*abs(f-g)**2);plus=summation(w*abs(f+g)**2);sign=1 if minus<=plus else -1;E=min(minus,plus)
 if G>0:
  design=(np.sqrt(w)*g)[:,None];target=np.sqrt(w)*f;q,r=np.linalg.qr(design);alpha=(q.conj().T@target)[0]/r[0,0];residual=summation(w*abs(f-alpha*g)**2)
 else:alpha=None;residual=None
 return {'e_sign':math.sqrt(E/G) if G>0 else None,'e_complex':math.sqrt(residual/F) if residual is not None and F>0 else None,'sign':sign,'FE_L2_norm':math.sqrt(F),'FD_L2_norm':math.sqrt(G),'difference_L2_norm':math.sqrt(E),'squared_error_integral':E,'FD_energy':G,'FE_energy':F,'alpha_real':float(alpha.real) if alpha is not None else None,'alpha_imag':float(alpha.imag) if alpha is not None else None,'alpha_modulus':float(abs(alpha)) if alpha is not None else None,'alpha_phase_rad':float(np.angle(alpha)) if alpha is not None and abs(alpha)>0 else None,'fitted_residual_energy':residual,'fraction_error_removed_same_denominator':1-residual/E if residual is not None and E>0 else None}
def main():
 cfg=read(OUT/'FUNC1_DESIGN_FREEZE.json');checks=[]
 def ck(name,error,tolerance,detail=''):
  checks.append({'check':name,'error':float(error),'tolerance':float(tolerance),'PASS':bool(np.isfinite(error) and error<=tolerance),'detail':detail})
 def val(name,a,b,detail=''):
  if a is None or b is None:ck(name,0 if a is None and (b is None or b=='') else 1,0,detail);return
  target=float(b);error=abs(float(a)-target)
  if name.endswith('alpha_phase_rad'):error=abs(math.atan2(math.sin(float(a)-target),math.cos(float(a)-target)))
  tol=cfg['cold_abs_tolerance']+cfg['cold_relative_tolerance']*max(abs(float(a)),abs(target))
  if name.endswith('fraction_error_removed_same_denominator') or name=='empirical_order':tol=max(tol,cfg['cold_derived_ratio_tolerance'])
  ck(name,error,tol,detail)
 def audit_row(name,expected,actual,keys):
  for k in keys:val(name+'/'+k,expected[k],actual.get(k),str({v:actual.get(v) for v in ['mode','N','grid','probe','resource','z_low_m']}))
 data={}
 for r in cfg['inputs']:
  err=0 if sha(ROOT/r['path'])==r['sha256'] else 1;ck('Git_source_identity',err,0,r['path'])
  if err:raise RuntimeError('INPUT_IDENTITY_STOP')
 for r in cfg['vectors']:
  err=0 if sha(r['path'])==r['sha256'] else 1;ck('NPZ_identity',err,0,r['path'])
  if err:raise RuntimeError('VECTOR_IDENTITY_STOP')
  with np.load(r['path'],allow_pickle=False) as a:
   ph=a['phi'].copy();h=float(a['h']);ck('NPZ_structure',0 if set(a.files)=={'phi','h'} and ph.shape==(r['N'],) and ph.dtype==np.dtype('complex128') and h==5000/r['N'] and np.isfinite(ph).all() else 1,0,r['path']);data[r['method'],r['mode'],r['N']]=(np.r_[0j,ph],h)
   v=np.r_[0j,ph];nb=(h/3*sum(v[:-1]**2+v[:-1]*v[1:]+v[1:]**2)) if r['method']=='FE' else h*(sum(ph[:-1]**2)+ph[-1]**2/2)
   ck('archived_bilinear_normalization',abs(nb-1),1e-8,r['path'])
 def at(method,mode,N,z):
  v,h=data[method,mode,N];return interpolate(v,h,z)
 env=(ROOT/cfg['env_path']).read_text(encoding='utf-8-sig').splitlines();stop=next(i for i,l in enumerate(env) if l.startswith("'R'"));arr=np.array([list(map(float,l.split()[:6])) for l in env[5:stop]]);z20=arr[:,0]
 grid={'SSP_UNWEIGHTED':(z20,False),'SSP_TRAPZ':(z20,True),'S1':(np.arange(5001,dtype=float),True),'S2':(np.arange(10001)*.5,True),'S3':(np.arange(20001)*.25,True)}
 metric_keys=['e_sign','e_complex','sign','FE_L2_norm','FD_L2_norm','difference_L2_norm','squared_error_integral','FD_energy','FE_energy','alpha_real','alpha_imag','alpha_modulus','alpha_phase_rad','fitted_residual_energy','fraction_error_removed_same_denominator']
 old=ROOT/'results/R4_T2_PREC3_COMPLEX_WAVENUMBER_AND_FIELD_VALIDATION';hist=rows(old/'PREC3_INDEPENDENT_COMPLEX_AND_PHI.csv');oldcold=rows(old/'PREC3_COLD_CHECKS.csv');cache={}
 original=rows(OUT/'FUNC1_ORIGINAL_GATE_RECONSTRUCTION.csv')
 for rr in original:
  m,N=int(rr['mode']),int(rr['N']);x=metric(at('FE',m,N,z20),at('FD',m,N,z20),z20,False);val('original_251',x['e_sign'],rr['e_original_251'],str((m,N)))
  prior=float(next(v for v in hist if v['mode']==str(m) and v['N']==str(N))['modal_function_FE_vs_FD_relative']);ck('historical_main_1e8',abs(x['e_sign']-prior),1e-8,str((m,N)))
  prior2=float(next(v for v in oldcold if v['check']=='independent_complex_mode_shape' and v['detail']==str(m)+'/'+str(N))['error']);ck('historical_cold_1e8',abs(x['e_sign']-prior2),1e-8,str((m,N)))
 phase=rows(OUT/'FUNC1_NORMALIZATION_AND_PHASE_AUDIT.csv')
 for rr in phase:
  m,N=int(rr['mode']),int(rr['N']);z,integral=grid[rr['grid']];x=metric(at('FE',m,N,z),at('FD',m,N,z),z,integral);cache[m,N,rr['grid']]=x;audit_row('normalization_phase',x,rr,metric_keys)
 dense=rows(OUT/'FUNC1_DENSE_SAMPLING_COMPARISON.csv')
 for rr in dense:
  m,N=int(rr['mode']),int(rr['N']);x=cache[m,N,rr['grid']];orig=cache[m,N,'SSP_UNWEIGHTED']['e_sign'];expected={'e_sign':x['e_sign'],'original_e_sign':orig,'signed_difference_vs_original':x['e_sign']-orig,'absolute_difference_vs_original':abs(x['e_sign']-orig),'relative_difference_vs_original':(x['e_sign']-orig)/orig if orig>0 else None};audit_row('dense',expected,rr,list(expected))
 stable=[]
 for m,N in dict.fromkeys((int(r['mode']),int(r['N'])) for r in dense):
  if (m,N,'S3') in cache:
   delta=abs(cache[m,N,'S2']['e_sign']-cache[m,N,'S3']['e_sign']);stable.append(delta<=2e-5);ck('dense_grid_stability',delta,2e-5,str((m,N)))
 z=grid['S3'][0];cross=rows(OUT/'FUNC1_FE_FD_CROSS_GRID.csv')
 for rr in cross:
  m=int(rr['mode']);n1=int(rr['N_first']);n2=int(rr['N_reference'])
  if rr['comparison']=='SAME_METHOD_CROSS_GRID':f,g=at(rr['method'],m,n1,z),at(rr['method'],m,n2,z)
  else:f,g=at('FE',m,n1,z),at('FD',m,n2,z)
  x=metric(f,g,z);audit_row('cross_grid',x,rr,metric_keys)
 samplepath=OUT/'FUNC1_DENSE_SAMPLES.npz'
 if samplepath.exists():
  with np.load(samplepath,allow_pickle=False) as samples:
   ck('dense_sample_grid',float(np.max(abs(samples['z_m']-z))),1e-12)
   for key in samples.files:
    if key=='z_m':continue
    meth,ms,ns=key.split('_');v=at(meth,int(ms[1:]),int(ns[1:]),z);ck('dense_raw_vector',float(np.max(abs(samples[key]-v))),1e-12,key)
 for rr in rows(OUT/'FUNC1_EMPIRICAL_ORDER.csv'):
  m=int(rr['mode']);key=rr['grid'];a=cache[m,80000,key]['e_sign'];b=cache[m,160000,key]['e_sign'];val('empirical_order',math.log2(a/b),rr['empirical_order'],str((m,key)))
 bands=rows(OUT/'FUNC1_DEPTH_RESOLVED_ERROR.csv')
 for rr in bands:
  m,N=int(rr['mode']),int(rr['N']);f=at('FE',m,N,z);g=at('FD',m,N,z);sign=cache[m,N,'S3']['sign'];lo=int(rr['z_low_m']);hi=int(rr['z_high_m']);sel=(z>=lo)&(z<=hi);ww=nodeweights(z[sel]);err=summation(ww*abs(f[sel]-sign*g[sel])**2);energy=summation(ww*abs(g[sel])**2);total=cache[m,N,'S3']['squared_error_integral'];expected={'squared_error_integral':err,'FD_energy':energy,'relative_error':math.sqrt(err/energy) if energy>0 else None,'fraction_of_total_squared_error':err/total if total>0 else None};audit_row('depth_bands',expected,rr,list(expected))
 for rr in rows(OUT/'FUNC1_LOCAL_FUNCTION_COMPARISON.csv'):
  m,N=int(rr['mode']),int(rr['N']);zz=np.arange(161)*.25+180;f,g=at('FE',m,N,zz),at('FD',m,N,zz);x=metric(f,g,zz);audit_row('local_function',x,rr,metric_keys)
  extra={'min_FE_amplitude':min(abs(f)),'max_FE_amplitude':max(abs(f)),'min_FD_amplitude':min(abs(g)),'max_FD_amplitude':max(abs(g)),'FE_points_below_amp_floor':int(sum(abs(f)<=cfg['weak_phi_amplitude'])),'FD_points_below_amp_floor':int(sum(abs(g)<=cfg['weak_phi_amplitude']))};audit_row('local_amplitudes',extra,rr,list(extra))
 channels=rows(OUT/'FUNC1_COUPLING_CHANNELS.csv');coupling=rows(OUT/'FUNC1_LOCAL_MODAL_COUPLING.csv');probes=rows(ROOT/cfg['probe_path'])[:5];geom=rows(ROOT/cfg['geometry_path'])
 for rr in coupling:
  m,N=int(rr['mode']),int(rr['N']);probe=next(v for v in probes if v['sentinel_id']==rr['probe']);coords=sorted([v for v in geom if v['resource']==rr['resource']],key=lambda r:int(r['channel']));zs=float(probe['source_z_m']);zr=np.array([float(v['z_m']) for v in coords]);fs,gs=at('FE',m,N,[zs])[0],at('FD',m,N,[zs])[0];fr,gr=at('FE',m,N,zr),at('FD',m,N,zr);f,g=fs*fr,gs*gr;x=metric(f,g,np.arange(len(g)),False)
  raw=np.linalg.norm(f-g)/np.linalg.norm(g) if np.linalg.norm(g)>0 else None;expected={'raw_relative_difference':raw,'array_shape_difference':x['e_complex'],'FE_coupling_norm':x['FE_L2_norm'],'FD_coupling_norm':x['FD_L2_norm'],'alpha_real_diagnostic_only':x['alpha_real'],'alpha_imag_diagnostic_only':x['alpha_imag']};audit_row('local_coupling',expected,rr,list(expected))
  for i,coord in enumerate(coords):
   ch=next(v for v in channels if int(v['mode'])==m and int(v['N'])==N and v['probe']==rr['probe'] and v['resource']==rr['resource'] and v['channel']==coord['channel'])
   expected={'FE_phi_source_real':fs.real,'FE_phi_source_imag':fs.imag,'FD_phi_source_real':gs.real,'FD_phi_source_imag':gs.imag,'FE_phi_receiver_real':fr[i].real,'FE_phi_receiver_imag':fr[i].imag,'FD_phi_receiver_real':gr[i].real,'FD_phi_receiver_imag':gr[i].imag,'FE_Q_real':f[i].real,'FE_Q_imag':f[i].imag,'FD_Q_real':g[i].real,'FD_Q_imag':g[i].imag};audit_row('raw_channel',expected,ch,list(expected))
 # Wavelengths use the archive k table, separately distinguish evanescence.
 kr=rows(old/'PREC3_COMPLEX_K_VALIDATION.csv')
 for rr in rows(OUT/'FUNC1_VERTICAL_WAVELENGTH_BY_DEPTH.csv'):
  m=int(rr['mode']);zz=float(rr['z_m']);ci=float(arr[np.flatnonzero(arr[:,0]==zz)[0],1]);k=float(next(v for v in kr if v['freq_hz']=='201' and int(v['mode'])==m)['K_complex_real']);q=(math.tau*201/ci)**2-k*k;val('vertical_kz_squared',q,rr['kz_squared'],str((m,zz)))
  if q>0:val('vertical_wavelength',math.tau/math.sqrt(q),rr['lambda_z_m'],str((m,zz)))
  else:ck('no_real_evanescent_wavelength',0 if rr['lambda_z_m']=='' else 1,0,str((m,zz)))
 decision=read(OUT/'FUNC1_DECISION.json');complete=decision.get('dense_quadrature')=='PASS' and decision.get('historical_reconstruction')=='PASS'
 if complete:
  labels=[];biased=False
  for m in [1,481]:
   for N in [80000,160000]:
    a=cache[m,N,'SSP_UNWEIGHTED']['e_sign'];b=cache[m,N,'S3']['e_sign'];biased|=abs(b-a)>max(cfg['evaluation_bias_abs_floor'],cfg['evaluation_bias_relative']*b)
  labels.append('FUNC1_EVALUATION_GRID_BIAS_OBSERVED' if biased else 'FUNC1_EVALUATION_GRID_NOT_DOMINANT_IN_TESTED_MODES')
  if any(1.8<=math.log2(cache[m,80000,'S3']['e_sign']/cache[m,160000,'S3']['e_sign'])<=2.2 for m in [1,481]):labels.append('FUNC1_SECOND_ORDER_BEHAVIOR_OBSERVED')
  labels.append('FUNC1_LOCAL_MODAL_COUPLING_DIAGNOSTIC_ONLY');ck('independent_classification',0 if set(labels)==set(decision['classifications']) else 1,0)
  for name,actual,expected in [('original',len(original),4),('dense',len(dense),20),('phase',len(phase),20),('cross',len(cross),9),('bands',len(bands),100),('coupling',len(coupling),80),('channels',len(channels),800)]:ck('completed_count/'+name,abs(actual-expected),0)
 scope=read(OUT/'FUNC1_MODE961_DATA_SCOPE.json') if (OUT/'FUNC1_MODE961_DATA_SCOPE.json').exists() else {}
 ck('961_not_falsely_paired',0 if not any(int(r['mode'])==961 for r in coupling) else 1,0)
 if scope:ck('961_scope_marked',0 if scope['FE_FD_dense_comparison']=='NOT_EVALUATED_FD_VECTOR_NOT_ARCHIVED' else 1,0)
 with (OUT/'FUNC1_COLD_CHECKS.csv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(checks[0]),lineterminator='\n');w.writeheader();w.writerows(checks)
 fail=sum(not r['PASS'] for r in checks);save(OUT/'FUNC1_COLD_AUDIT.json',{'PASS':fail==0,'checks':len(checks),'FAIL':fail,'vectors_read':len(data),'paired_modes_dense':[1,481] if complete else [],'mode961_new_FD_solves':0,'method':'independent linear interpolation and interval weights; fsum integrals; QR coefficient; fresh archived NPZ inputs; no main imports','R4_percent':0,'STOP':True})
 print('COLD_COMPLETE',len(checks),'FAIL',fail,flush=True)
 if fail:raise RuntimeError('FUNC1_INDEPENDENT_AUDIT_FAIL_NO_RETRY')
if __name__=='__main__':main()
