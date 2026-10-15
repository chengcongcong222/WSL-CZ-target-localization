"""Independent cold: fresh source parsing, scalar integrals, ODE residuals,
fixed-k RK45 two-sided reconstruction for mode481, and raw archive recoupling.
Does not import any FUNC2 main/reference module or read their PASS flags as evidence.
"""
from pathlib import Path
import json,csv,hashlib,math,traceback
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3_FUNC2');S=math.tau*201/1500

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def table(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def weights(z):return np.r_[(z[1]-z[0])/2,(z[2:]-z[:-2])/2,(z[-1]-z[-2])/2]
def integral(z,y):return complex(math.fsum((weights(z)*y).real),math.fsum((weights(z)*y).imag))
def e_sign(f,g,z):
 w=weights(z);a=math.fsum(w*abs(f-g)**2);b=math.fsum(w*abs(f+g)**2);den=math.fsum(w*abs(g)**2)
 return math.sqrt(min(a,b)/den) if den else None

def fresh_rk(z,c,k,query):
 a=.201**2;atten=(.0033+.11*a/(1+a)+44*a/(4100+a)+.0003*a)/8685.8896
 ccomplex=c+1j*atten*c*c/(math.tau*201)
 match=2500.
 def side(start,end,y0,record):
  knots=z[(z>min(start,end))&(z<max(start,end))];knots=knots if start<end else knots[::-1];edges=np.r_[start,knots,end];v=np.array(y0,complex);acc=0.;parts=[]
  for left,right in zip(edges[:-1],edges[1:]):
   def ode(t,y):
    cc=np.interp(t,z,ccomplex.real)+1j*np.interp(t,z,ccomplex.imag)
    return [S*y[1],(k*k-(math.tau*201/cc)**2)*y[0]/S]
   sol=solve_ivp(ode,(left,right),v,method='RK45',rtol=3e-12,atol=3e-14,max_step=.5,dense_output=record)
   if not sol.success:raise RuntimeError('COLD_RK_INTEGRATION_FAILED')
   if record:parts.append((left,right,sol.sol,acc))
   v=sol.y[:,-1];scale=np.linalg.norm(v);v=v/scale;acc+=math.log(scale)
  return v,acc,parts
 left,ll,L=side(0,match,[0,1],True);right,lr,R=side(5000,match,[1,0],True);beta=np.vdot(right,left)/np.vdot(right,right);mismatch=np.linalg.norm(left-beta*right)
 def values(zz):
  yy=np.zeros((2,len(zz)),complex)
  for pieces,log,factor,up in [(L,ll,1,True),(R,lr,beta,False)]:
   for a,b,fn,base in pieces:
    ix=(zz>=min(a,b))&(zz<=max(a,b))&((zz<=match) if up else (zz>match))
    if ix.any():yy[:,ix]=fn(zz[ix])*math.exp(base-log)*factor
  return yy
 x,w=np.polynomial.legendre.leggauss(24);q=(z[:-1,None]+10*(x+1)).ravel();pv=values(q)[0].reshape(250,24);norm=np.sqrt(10*np.sum(w*pv*pv));v=values(query)/norm
 return v,float(mismatch)

def main():
 cfg=read(OUT/'FUNC2_DESIGN_FREEZE.json');checks=[]
 def ck(name,error,tol,detail=''):
  checks.append({'check':name,'error':float(error),'tolerance':float(tol),'PASS':bool(np.isfinite(error) and error<=tol),'detail':str(detail)})
 def eq(name,a,b,tol=1e-10,detail=''):ck(name,abs(a-b),tol,detail)
 try:
  for rec in cfg['inputs']:ck('parent_source_SHA',0 if sha(ROOT/rec['path'])==rec['sha256'] else 1,0,rec['path'])
  for rec in cfg['vectors']:ck('original_vector_SHA',0 if sha(rec['path'])==rec['sha256'] else 1,0,rec['path'])
  for rec in cfg['frozen_code']:ck('frozen_code_SHA',0 if sha(OUT/rec['path'])==rec['sha256'] else 1,0,rec['path'])
  if any(not r['PASS'] for r in checks):raise RuntimeError('COLD_IDENTITY_FAIL')
  index=read(OUT/'FUNC2_RAW_DATA_INDEX.json');raw={}
  ls=(ROOT/cfg['env_path']).read_text(encoding='utf-8-sig').splitlines();jj=next(i for i,s in enumerate(ls) if s.startswith("'R'"));arr=np.array([list(map(float,v.split()[:6])) for v in ls[5:jj]]);zz,cc=arr[:,0],arr[:,1]
  analytic=table(OUT/'FUNC2_ANALYTIC_WAVEGUIDE_VALIDATION.csv');rootrows=table(OUT/'FUNC2_CONTINUOUS_COMPLEX_ROOTS.csv')
  for rec in index['files']:
   ck('saved_reference_SHA',0 if sha(rec['path'])==rec['sha256'] else 1,0,rec['tag']);a=np.load(rec['path']);data={key:a[key].copy() for key in a.files};a.close();raw[rec['tag']]=data
   z=data['z'];p=data['phi'];dp=data['dphi'];norm=integral(z,p*p);eq('bilinear_norm',norm,1,cfg['norm_budget'],rec['tag']);bc=max(abs(p[0]),abs(dp[-1])/S)/max(abs(p));ck('boundary_residual',bc,cfg['boundary_budget'],rec['tag'])
   if rec['tag'].startswith('analytic_'):
    match=[r for r in analytic if r['method'] in rec['tag'] and int(r['mode'])==rec['mode'] and int(r['level'])==int(rec['tag'].split('_L')[1][0]) and float(r['loss'])==rec['loss']]
    if len(match)!=1:raise RuntimeError('COLD_ANALYTIC_ROW_AMBIGUOUS')
    row=match[0];mode=rec['mode'];mu=(mode-.5)*math.pi/5000;loss=rec['loss'];a2=.201**2;atten=(.0033+.11*a2/(1+a2)+44*a2/(4100+a2)+.0003*a2)/8685.8896;c0=1500+1j*loss*atten*1500**2/(math.tau*201);k=np.sqrt((math.tau*201/c0)**2-mu**2);num=complex(float(row['k_real']),float(row['k_imag']));eq('analytic_real_k_error',abs(num.real-k.real),float(row['real_k_error']),1e-14,rec['tag']);eq('analytic_imag_k_error',abs(num.imag-k.imag),float(row['imag_k_error']),1e-14,rec['tag']);truth=math.sqrt(2/5000)*np.sin(mu*z);e=e_sign(p,truth,z);eq('analytic_shape_reconstructed',e,float(row['dense_shape_error']),1e-11,rec['tag']);ck('analytic_shape_gate',e,cfg['shape_budget'],rec['tag'])
   else:
    row=next(r for r in rootrows if r['method'] in rec['tag'] and int(r['mode'])==rec['mode'] and int(r['level'])==int(rec['tag'].split('_L')[1][0]) and float(r['loss'])==rec['loss']);k=complex(float(row['k_real']),float(row['k_imag']));a2=.201**2;att=(.0033+.11*a2/(1+a2)+44*a2/(4100+a2)+.0003*a2)/8685.8896;cnode=cc+1j*rec['loss']*att*cc*cc/(math.tau*201);cp=np.interp(z,zz,cnode.real)+1j*np.interp(z,zz,cnode.imag);potential=(math.tau*201/cp)**2
    # Eighth-order derivative of actual sampled function, exclude SSP derivative knots.
    coeff=np.array([-1/560,8/315,-1/5,8/5,-205/72,8/5,-1/5,8/315,-1/560]);second=sum(coeff[i]*p[i:len(p)-8+i] for i in range(9))/.25**2;rr=second+(potential[4:-4]-k*k)*p[4:-4];ix=((z[4:-4]%20)>1.1)&((z[4:-4]%20)<18.9);measure=np.linalg.norm(rr[ix])/(S*S*np.linalg.norm(p[4:-4][ix]));ck('independent_ODE_residual',measure,1e-6,rec['tag'])
  # Reconstruct all old-function errors directly from original NPZ, not main flags.
  for row in table(OUT/'FUNC2_FE_FD_VS_REFERENCE.csv'):
   m,N=int(row['mode']),int(row['N']);ref=raw[f'estd_DOP853_L1_m{m}_a1'];source=next(v for v in cfg['vectors'] if v['method']==row['method'] and v['mode']==m and v['N']==N)
   with np.load(source['path']) as a:ph=np.r_[0j,a['phi']];h=float(a['h'])
   z=ref['z'];u=z/h;i=np.minimum(np.floor(u).astype(int),len(ph)-2);f=ph[i]+(u-i)*(ph[i+1]-ph[i]);e=e_sign(f,ref['phi'],z);eq('FE_FD_vs_reference_rebuilt',e,float(row['e_sign']),1e-10,str((m,N,row['method'])))
  # One preregistered high-order mode freshly integrated with independent RK45 at fixed k.
  high=raw.get('estd_DOP853_L1_m481_a1');cold_mode=False
  if high is not None:
   r=next(r for r in rootrows if r['method']=='DOP853' and r['level']=='1' and r['mode']=='481' and float(r['loss'])==1);k=complex(float(r['k_real']),float(r['k_imag']));query=np.unique(np.r_[high['z'],high['z_query']]);v,mis=fresh_rk(zz,cc,k,query);ck('fixed_k_RK45_matching_residual',mis,1e-7,'481');on=np.searchsorted(query,high['z']);err=e_sign(v[0,on],high['phi'],high['z']);ck('independent_RK45_vs_continuous_function',err,cfg['shape_budget'],'481');sg=1 if np.linalg.norm(v[0,on]-high['phi'])<=np.linalg.norm(v[0,on]+high['phi']) else -1;qi=np.searchsorted(query,high['z_query']);ck('independent_RK45_query_values',np.linalg.norm(v[0,qi]-sg*high['phi_query'])/np.linalg.norm(high['phi_query']),cfg['coupling_budget'],'481')
   np.savez_compressed(LOCAL/'COLD_RK45_MODE481.npz',z=query,phi=v[0],dphi=v[1]*S,k=k);cold_mode=True
   for rec in table(ROOT/cfg['probe_path'])[:5]:
    for resource in cfg['resources']:
     geo=sorted([r for r in table(ROOT/cfg['geometry_path']) if r['resource']==resource],key=lambda r:int(r['channel']));zs=float(rec['source_z_m']);zr=np.array([float(r['z_m']) for r in geo]);qs=v[0,np.searchsorted(query,zs)]*v[0,np.searchsorted(query,zr)];qt=high['phi_query'][np.searchsorted(high['z_query'],zs)]*high['phi_query'][np.searchsorted(high['z_query'],zr)];ck('fresh_RK45_highmode_coupling',np.linalg.norm(qs-qt)/np.linalg.norm(qt),cfg['coupling_budget'],rec['sentinel_id']+'/'+resource)
  for row in table(OUT/'FUNC2_LOCAL_COUPLING_REFERENCE.csv'):
   mode=int(row['mode']);N=int(row['N']);record=next(v for v in cfg['vectors'] if v['method']==row['method'] and v['mode']==mode and v['N']==N);rr=raw[f'estd_DOP853_L1_m{mode}_a1'];probe=next(r for r in table(ROOT/cfg['probe_path']) if r['sentinel_id']==row['probe']);zr=np.array([float(v['z_m']) for v in sorted([v for v in table(ROOT/cfg['geometry_path']) if v['resource']==row['resource']],key=lambda v:int(v['channel']))]);zs=float(probe['source_z_m']);iz=np.searchsorted(rr['z_query'],np.r_[zs,zr]);g=rr['phi_query'][iz];qb=g[0]*g[1:]
   with np.load(record['path']) as a:ph=np.r_[0j,a['phi']];h=float(a['h'])
   u=np.r_[zs,zr]/h;ix=np.minimum(np.floor(u).astype(int),len(ph)-2);f=ph[ix]+(u-ix)*(ph[ix+1]-ph[ix]);qa=f[0]*f[1:];den=np.linalg.norm(qb);e=np.linalg.norm(qa-qb)/den if den else None
   if e is not None:eq('old_archive_coupling_rebuilt',e,float(row['raw_relative_error']),1e-9,str((mode,N,row['method'],row['probe'],row['resource'])))
  # Expected data coverage is separate from physical/sample count.
  if table(OUT/'FUNC2_CONTINUOUS_MODE_REFERENCE.csv'):
   ck('three_registered_mode_references',abs(len(table(OUT/'FUNC2_CONTINUOUS_MODE_REFERENCE.csv'))-3),0)
  result={'PASS':not any(not r['PASS'] for r in checks),'checks':len(checks),'FAIL':sum(not r['PASS'] for r in checks),'reference_arrays_rebuilt_or_residual_checked':len(raw),'new_independent_RK45_highmode_reconstruction':1 if cold_mode else 0,'scope':'source SHA; analytic controls; fresh scalar norm; sampled ODE residual; original archive interpolation; fixed-k independent RK45 mode481 and20Q vectors when available','R4_percent':0,'STOP':True}
 except Exception as e:
  result={'PASS':False,'checks':len(checks),'FAIL':sum(not r['PASS'] for r in checks)+1,'exception':str(e),'R4_percent':0,'STOP':True};(OUT/'FUNC2_COLD_EXCEPTION.txt').write_text(traceback.format_exc(),encoding='utf-8')
 with (OUT/'FUNC2_COLD_CHECKS.csv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['check','error','tolerance','PASS','detail'],lineterminator='\n');w.writeheader();w.writerows(checks)
 (OUT/'FUNC2_INDEPENDENT_COLD_AUDIT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
