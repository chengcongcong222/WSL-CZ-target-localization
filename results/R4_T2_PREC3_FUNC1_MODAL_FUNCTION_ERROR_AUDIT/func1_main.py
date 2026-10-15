"""FUNC1: read-only archived modal-vector attribution, no eigen/field solvers."""
from pathlib import Path
import json,csv,hashlib,math
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def table(p,x,fields=None):
 fields=fields or list(dict.fromkeys(k for r in x for k in r))
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(x)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_vectors(cfg):
 data={}
 for r in cfg['vectors']:
  assert sha(r['path'])==r['sha256']
  with np.load(r['path'],allow_pickle=False) as a:data[r['method'],r['mode'],r['N']]=(np.r_[0j,a['phi']],float(a['h']))
 return data
def sample(data,method,mode,N,z):
 v,h=data[method,mode,N];x=np.arange(len(v))*h
 return np.interp(z,x,v.real)+1j*np.interp(z,x,v.imag)
def weights(z):
 d=np.diff(z);w=np.zeros(len(z));w[:-1]+=d/2;w[1:]+=d/2;return w
def metrics(f,g,w):
 F=float(np.sum(w*abs(f)**2));G=float(np.sum(w*abs(g)**2));minus=float(np.sum(w*abs(f-g)**2));plus=float(np.sum(w*abs(f+g)**2));sig=1 if minus<=plus else -1;E=min(minus,plus)
 alpha=np.sum(w*g.conj()*f)/G if G>0 else None
 residual=float(np.sum(w*abs(f-alpha*g)**2)) if alpha is not None else None
 return {'e_sign':math.sqrt(E/G) if G>0 else None,'e_complex':math.sqrt(residual/F) if residual is not None and F>0 else None,'sign':sig,'FE_L2_norm':math.sqrt(F),'FD_L2_norm':math.sqrt(G),'difference_L2_norm':math.sqrt(E),'squared_error_integral':E,'FD_energy':G,'FE_energy':F,'alpha_real':float(alpha.real) if alpha is not None else None,'alpha_imag':float(alpha.imag) if alpha is not None else None,'alpha_modulus':float(abs(alpha)) if alpha is not None else None,'alpha_phase_rad':float(np.angle(alpha)) if alpha is not None and abs(alpha)>0 else None,'fitted_residual_energy':residual,'fraction_error_removed_same_denominator':1-residual/E if residual is not None and E>0 else None}
def main():
 cfg=read(OUT/'FUNC1_DESIGN_FREEZE.json');data=load_vectors(cfg);old=ROOT/'results/R4_T2_PREC3_COMPLEX_WAVENUMBER_AND_FIELD_VALIDATION';hist=rows(old/'PREC3_INDEPENDENT_COMPLEX_AND_PHI.csv');cold=rows(old/'PREC3_COLD_CHECKS.csv')
 env=(ROOT/cfg['env_path']).read_text(encoding='utf-8-sig').splitlines();j=next(i for i,l in enumerate(env) if l.startswith("'R'"));ssp=np.array([list(map(float,l.split()[:6])) for l in env[5:j]])[:,0:2];z20=ssp[:,0]
 decision={'stage':'T2-PREC3-FUNC1','archive_vectors':10,'historical_reconstruction':'NOT_EVALUATED','dense_quadrature':'NOT_EVALUATED','classifications':[],'R4_percent':0,'historical_PREC3_B':'FAIL_UNCHANGED','STOP':True}
 recon=[];dense=[];phase=[];comp={}
 for mode in [1,481]:
  for N in [80000,160000]:
   f=sample(data,'FE',mode,N,z20);g=sample(data,'FD',mode,N,z20);m=metrics(f,g,np.ones(len(z20)))
   prior=float(next(r for r in hist if r['mode']==str(mode) and r['N']==str(N))['modal_function_FE_vs_FD_relative']);pc=float(next(r for r in cold if r['check']=='independent_complex_mode_shape' and r['detail']==str(mode)+'/'+str(N))['error'])
   r={'mode':mode,'N':N,'e_original_251':m['e_sign'],'historical_main':prior,'historical_cold':pc,'main_absolute_difference':abs(m['e_sign']-prior),'cold_absolute_difference':abs(m['e_sign']-pc),'PASS':abs(m['e_sign']-prior)<=1e-8 and abs(m['e_sign']-pc)<=1e-8};recon.append(r);table(OUT/'FUNC1_ORIGINAL_GATE_RECONSTRUCTION.csv',recon)
   if not r['PASS']:
    decision.update(historical_reconstruction='FAIL',classifications=['FUNC1_ARCHIVED_FUNCTION_RECONSTRUCTION_FAIL']);write(OUT/'FUNC1_DECISION.json',decision);return
 decision['historical_reconstruction']='PASS';write(OUT/'FUNC1_DECISION.json',decision)
 griddefs=[('SSP_UNWEIGHTED',z20,False),('SSP_TRAPZ',z20,True),('S1',np.linspace(0,5000,5001),True),('S2',np.linspace(0,5000,10001),True),('S3',np.linspace(0,5000,20001),True)]
 for mode in [1,481]:
  for N in [80000,160000]:
   prior=next(r for r in recon if r['mode']==mode and r['N']==N)['e_original_251']
   for name,z,integral in griddefs:
    f=sample(data,'FE',mode,N,z);g=sample(data,'FD',mode,N,z);m=metrics(f,g,weights(z) if integral else np.ones(len(z)));comp[mode,N,name]=m
    dense.append({'mode':mode,'N':N,'grid':name,'points':len(z),'integration':'TRAPEZOID_MOD_SQUARED' if integral else 'UNWEIGHTED_DISCRETE_L2','e_sign':m['e_sign'],'original_e_sign':prior,'signed_difference_vs_original':m['e_sign']-prior,'absolute_difference_vs_original':abs(m['e_sign']-prior),'relative_difference_vs_original':(m['e_sign']-prior)/prior if prior>0 else None})
    phase.append({'mode':mode,'N':N,'grid':name,**m})
   diff=abs(comp[mode,N,'S2']['e_sign']-comp[mode,N,'S3']['e_sign']);table(OUT/'FUNC1_DENSE_SAMPLING_COMPARISON.csv',dense);table(OUT/'FUNC1_NORMALIZATION_AND_PHASE_AUDIT.csv',phase)
   if diff>2e-5:
    decision.update(dense_quadrature='FAIL',classifications=['FUNC1_DENSE_QUADRATURE_NOT_STABLE'],first_failed_mode=mode,first_failed_N=N);write(OUT/'FUNC1_DECISION.json',decision);return
 decision['dense_quadrature']='PASS';write(OUT/'FUNC1_DECISION.json',decision)
 # Sampling-risk diagnostic uses only the frozen nominal c and archived real k.
 kr=rows(old/'PREC3_COMPLEX_K_VALIDATION.csv');waves=[];risk=[]
 for mode in [1,481,961]:
  k=float(next(r for r in kr if r['freq_hz']=='201' and int(r['mode'])==mode)['K_complex_real']);q=(2*np.pi*201/ssp[:,1])**2-k*k;lam=np.full(len(q),np.nan);lam[q>0]=2*np.pi/np.sqrt(q[q>0])
  for z,c2,q2,l in zip(ssp[:,0],ssp[:,1],q,lam):waves.append({'mode':mode,'z_m':z,'c_m_s':c2,'kz_squared':q2,'region':'OSCILLATORY' if q2>0 else 'EVANESCENT' if q2<0 else 'TURNING','lambda_z_m':float(l) if q2>0 else None})
  ll=lam[q>0];risk.append({'mode':mode,'min_lambda_z_m':float(min(ll)) if len(ll) else None,'max_lambda_z_m':float(max(ll)) if len(ll) else None,'oscillatory_SSP_nodes':int(sum(q>0)),'evanescent_SSP_nodes':int(sum(q<0)),'20m_points_per_shortest_wavelength':float(min(ll)/20) if len(ll) else None,'1m_points_per_shortest_wavelength':float(min(ll)) if len(ll) else None,'0p25m_points_per_shortest_wavelength':float(min(ll)/.25) if len(ll) else None,'risk':'UNDERSAMPLED_AT20M' if len(ll) and min(ll)<40 else 'NO_NYQUIST_VIOLATION_AT_SSP_NODES','bias_direction_inferred_from_Nyquist':False})
 table(OUT/'FUNC1_VERTICAL_WAVELENGTH_BY_DEPTH.csv',waves);table(OUT/'FUNC1_SAMPLING_RISK.csv',risk)
 # Common dense samples preserve the raw functions, without renormalization.
 z=np.linspace(0,5000,20001);w=weights(z);saved={'z_m':z};cross=[];orders=[]
 for method,modes in [('FE',[1,481,961]),('FD',[1,481])]:
  for mode in modes:
   a=sample(data,method,mode,80000,z);b=sample(data,method,mode,160000,z)
   saved[f'{method}_m{mode}_N80000']=a;saved[f'{method}_m{mode}_N160000']=b
   cross.append({'comparison':'SAME_METHOD_CROSS_GRID','method':method,'mode':mode,'N_first':80000,'N_reference':160000,'empirical_order':None,'order_scope':'ONE_CROSS_GRID_PAIR_DOES_NOT_IDENTIFY_METHOD_ORDER',**metrics(a,b,w)})
 for mode in [1,481]:
  for grid,_,_ in griddefs:
   e0=comp[mode,80000,grid]['e_sign'];e1=comp[mode,160000,grid]['e_sign'];p=math.log2(e0/e1) if e0>0 and e1>0 else None;orders.append({'mode':mode,'grid':grid,'coarse_error':e0,'fine_error':e1,'error_ratio':e0/e1 if e1>0 else None,'empirical_order':p,'second_order_descriptive':p is not None and 1.8<=p<=2.2})
  for N in [80000,160000]:
   f=saved[f'FE_m{mode}_N{N}'];g=saved[f'FD_m{mode}_N{N}'];p=next(r['empirical_order'] for r in orders if r['mode']==mode and r['grid']=='S3')
   cross.append({'comparison':'PAIRED_FE_FD','method':'FE_VS_FD','mode':mode,'N_first':N,'N_reference':N,'empirical_order':p,'order_scope':'CROSS_METHOD_DISAGREEMENT_ONLY',**metrics(f,g,w)})
 table(OUT/'FUNC1_FE_FD_CROSS_GRID.csv',cross);table(OUT/'FUNC1_EMPIRICAL_ORDER.csv',orders);np.savez_compressed(OUT/'FUNC1_DENSE_SAMPLES.npz',**saved)
 bands=[];local=[]
 for mode in [1,481]:
  for N in [80000,160000]:
   f=saved[f'FE_m{mode}_N{N}'];g=saved[f'FD_m{mode}_N{N}'];m=comp[mode,N,'S3'];sig=m['sign'];E=m['squared_error_integral']
   for lo in range(0,5000,200):
    ix=(z>=lo)&(z<=lo+200);ww=weights(z[ix]);err=float(np.sum(ww*abs(f[ix]-sig*g[ix])**2));energy=float(np.sum(ww*abs(g[ix])**2))
    bands.append({'mode':mode,'N':N,'z_low_m':lo,'z_high_m':lo+200,'global_sign':sig,'squared_error_integral':err,'FD_energy':energy,'relative_error':math.sqrt(err/energy) if energy>0 else None,'fraction_of_total_squared_error':err/E if E>0 else None,'weak_reference':math.sqrt(energy)<=cfg['weak_local_L2_norm']})
   ix=(z>=180)&(z<=220);ll=metrics(f[ix],g[ix],weights(z[ix]));local.append({'mode':mode,'N':N,'z_low_m':180,'z_high_m':220,'min_FE_amplitude':float(min(abs(f[ix]))),'max_FE_amplitude':float(max(abs(f[ix]))),'min_FD_amplitude':float(min(abs(g[ix]))),'max_FD_amplitude':float(max(abs(g[ix]))),'FE_points_below_amp_floor':int(sum(abs(f[ix])<=cfg['weak_phi_amplitude'])),'FD_points_below_amp_floor':int(sum(abs(g[ix])<=cfg['weak_phi_amplitude'])),'weak_reference':ll['FD_L2_norm']<=cfg['weak_local_L2_norm'],**ll})
 table(OUT/'FUNC1_DEPTH_RESOLVED_ERROR.csv',bands);table(OUT/'FUNC1_LOCAL_FUNCTION_COMPARISON.csv',local)
 probes=rows(ROOT/cfg['probe_path'])[:5];geom=rows(ROOT/cfg['geometry_path']);coupling=[];channels=[]
 for mode in [1,481]:
  for N in [80000,160000]:
   for probe in probes:
    zs=float(probe['source_z_m']);fs=sample(data,'FE',mode,N,[zs])[0];gs=sample(data,'FD',mode,N,[zs])[0]
    for resource in cfg['resources']:
     coords=sorted([r for r in geom if r['resource']==resource],key=lambda r:int(r['channel']));zr=np.array([float(r['z_m']) for r in coords]);fr=sample(data,'FE',mode,N,zr);gr=sample(data,'FD',mode,N,zr);qf=fs*fr;qg=gs*gr;m=metrics(qf,qg,np.ones(len(qg)));raw=float(np.linalg.norm(qf-qg)/np.linalg.norm(qg)) if np.linalg.norm(qg)>0 else None
     record={'mode':mode,'N':N,'probe':probe['sentinel_id'],'resource':resource,'source_z_m':zs,'channels':len(zr),'source_FE_abs':float(abs(fs)),'source_FD_abs':float(abs(gs)),'raw_relative_difference':raw,'array_shape_difference':m['e_complex'],'FE_coupling_norm':m['FE_L2_norm'],'FD_coupling_norm':m['FD_L2_norm'],'alpha_real_diagnostic_only':m['alpha_real'],'alpha_imag_diagnostic_only':m['alpha_imag'],'weak_source_or_receiver':bool(abs(fs)<=cfg['weak_phi_amplitude'] or abs(gs)<=cfg['weak_phi_amplitude'] or np.any(abs(fr)<=cfg['weak_phi_amplitude']) or np.any(abs(gr)<=cfg['weak_phi_amplitude'])),'weak_coupling_norm':bool(m['FE_L2_norm']<=cfg['weak_coupling_norm'] or m['FD_L2_norm']<=cfg['weak_coupling_norm'])}
     coupling.append(record)
     for i,rr in enumerate(coords):channels.append({'mode':mode,'N':N,'probe':probe['sentinel_id'],'resource':resource,'channel':rr['channel'],'source_z_m':zs,'receiver_z_m':zr[i],'FE_phi_source_real':float(fs.real),'FE_phi_source_imag':float(fs.imag),'FD_phi_source_real':float(gs.real),'FD_phi_source_imag':float(gs.imag),'FE_phi_receiver_real':float(fr[i].real),'FE_phi_receiver_imag':float(fr[i].imag),'FD_phi_receiver_real':float(gr[i].real),'FD_phi_receiver_imag':float(gr[i].imag),'FE_Q_real':float(qf[i].real),'FE_Q_imag':float(qf[i].imag),'FD_Q_real':float(qg[i].real),'FD_Q_imag':float(qg[i].imag)})
 table(OUT/'FUNC1_LOCAL_MODAL_COUPLING.csv',coupling);table(OUT/'FUNC1_COUPLING_CHANNELS.csv',channels)
 write(OUT/'FUNC1_MODE961_DATA_SCOPE.json',{'FE_cross_grid':'COMPLETED_0p25M_GRID','FE_FD_dense_comparison':'NOT_EVALUATED_FD_VECTOR_NOT_ARCHIVED','single_mode_FE_FD_coupling':'NOT_EVALUATED_FD_VECTOR_NOT_ARCHIVED','historical_cold_results_scope':'HISTORICAL_COLD_RESULT_ONLY','historical_cold_rows':[r for r in cold if r['check']=='independent_complex_mode_shape' and r['detail'].startswith('961/')],'new_FD_Galerkin_Pruefer_solves':0})
 classification=[];comparison=[]
 for mode in [1,481]:
  for N in [80000,160000]:
   orig=comp[mode,N,'SSP_UNWEIGHTED']['e_sign'];val=comp[mode,N,'S3']['e_sign'];bias=abs(val-orig);observed=bias>max(cfg['evaluation_bias_abs_floor'],cfg['evaluation_bias_relative']*val)
   comparison.append({'mode':mode,'N':N,'original':orig,'dense':val,'bias_observed':observed,'sign_of_dense_minus_original':int(np.sign(val-orig)),'original_failure_persists_descriptively':orig>.002 and val>.002})
 classification.append('FUNC1_EVALUATION_GRID_BIAS_OBSERVED' if any(r['bias_observed'] for r in comparison) else 'FUNC1_EVALUATION_GRID_NOT_DOMINANT_IN_TESTED_MODES')
 if any(r['second_order_descriptive'] for r in orders if r['grid']=='S3'):classification.append('FUNC1_SECOND_ORDER_BEHAVIOR_OBSERVED')
 classification.append('FUNC1_LOCAL_MODAL_COUPLING_DIAGNOSTIC_ONLY')
 decision.update(classifications=classification,original_reconstructions=len(recon),dense_sampling_rows=len(dense),phase_rows=len(phase),cross_grid_rows=len(cross),depth_band_rows=len(bands),local_function_rows=len(local),coupling_vectors=len(coupling),coupling_channels=len(channels),original_vs_dense= comparison,FD961_status='NOT_EVALUATED_FD_VECTOR_NOT_ARCHIVED',independent_audit='PENDING',new_solves=0)
 write(OUT/'FUNC1_DECISION.json',decision);print(json.dumps(decision,indent=2),flush=True)
if __name__=='__main__':main()
