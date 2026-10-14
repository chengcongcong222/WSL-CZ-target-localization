"""One cold review: distinct FD/Pruefer/QR and raw archive reconstruction."""
from pathlib import Path
import json,csv,hashlib,struct,cmath,math,time,sys
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3');sys.path.insert(0,str(ROOT))
import t2_prec2_cold_audit as continuous
LIM=1e-4/60000
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):
 if not Path(p).exists():return []
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def mod(p):
 b=Path(p).read_bytes();r=struct.unpack_from('<i',b)[0]*4;n=struct.unpack_from('<i',b,92)[0];m=struct.unpack_from('<i',b,5*r)[0]
 z=np.array(struct.unpack_from('<'+'f'*n,b,4*r))
 pp=np.array([[complex(*struct.unpack_from('<2f',b,(7+j)*r+8*i)) for j in range(m)] for i in range(n)])
 kk=np.array([complex(*struct.unpack_from('<2f',b,(7+m)*r+8*j)) for j in range(m)])
 return z,pp,kk
def fd(f,z,c,N,j):
 h=5000/N;xx=np.arange(1,N+1)*h;w=2*math.pi*f;speed=np.interp(xx,z,c);diag=(w/speed)**2-2/h**2;off=np.ones(N-1)/h**2;off[-1]=math.sqrt(2)/h**2
 lam,v=eigh_tridiagonal(diag,off,select='i',select_range=(N-j,N-j),lapack_driver='stebz');lam=lam[0];v=v[:,0]
 x=f*f/1e6;alpha=(.0033+.11*x/(1+x)+44*x/(4100+x)+.0003*x)/8685.8896
 imaginary_nodes=alpha*c*c/w;cp=speed+1j*np.interp(xx,z,imaginary_nodes);d=(w/cp)**2-2/h**2;delta=np.dot(v*v,d-diag)
 A=diags([off,d,off],[-1,0,1],format='csc')
 ev,u=eigs(A,k=1,sigma=lam+delta,which='LM',v0=v.astype(complex),ncv=10,maxiter=2000,tol=1e-12);u=u[:,0];u=u/np.sqrt(np.dot(u,u));u=u/math.sqrt(h);u[-1]*=math.sqrt(2)
 return complex(ev[0])-lam,u,h
def main():
 began=time.monotonic();checks=[]
 def ck(name,error,tol,detail=''):
  checks.append({'check':name,'error':float(error),'tolerance':tol,'PASS':bool(np.isfinite(error) and error<=tol),'detail':detail})
 cfg=read(OUT/'PREC3_DESIGN_FREEZE.json')
 for r in cfg['inputs']:ck('source_identity',0 if sha(ROOT/r['path'])==r['sha256'] else 1,0,r['path'])
 a=rows(OUT/'PREC3_REAL_K_INDEPENDENT_VALIDATION.csv');b=rows(OUT/'PREC3_COMPLEX_K_VALIDATION.csv');raw=rows(OUT/'PREC3_RAW_EIGENVALUES.csv');oldroot=ROOT/'results/R4_T2_PREC2_ESTD_EIGENVALUE_PRECISION_FEASIBILITY';arc=read(oldroot/'T2_PREC2_EXPERIMENT_SPEC_20261014_v1.json')['archive_sets']
 src=next(r for r in arc if r['freq']==201 and r['mesh']==160001);z,c=continuous.env_speed(ROOT/src['paths']['env']);ivp=[]
 for ar in a:
  selected=[r for r in raw if r['stage']=='A' and r['freq_hz']==ar['freq_hz'] and r['mode']==ar['mode']]
  hs=np.array([1/int(r['N'])**2 for r in selected]);v=np.array([float(r['lambda_real']) for r in selected]);powers=np.column_stack([np.ones(4),hs/hs[0],(hs/hs[0])**2,(hs/hs[0])**3])
  lam=np.linalg.solve(powers,v)[0];ck('independent_extrapolation',abs(math.sqrt(lam)-float(ar['K_GALERKIN_R6'])),1e-12,ar['freq_hz']+'/'+ar['mode'])
 for j in [1,481,961]:
  hit=next((r for r in a if r['freq_hz']=='201' and int(r['mode'])==j),None)
  if hit is None:continue
  guess=float(hit['K_GALERKIN_R6']);loose,n1=continuous.prufer_root(201,z,c,j,guess,1.,1e-11,600,began);tight,n2=continuous.prufer_root(201,z,c,j,guess,.5,1e-12,600,began)
  ck('cold_Pruefer_two_levels',abs(tight-loose),LIM,str(j));ck('cold_Pruefer_vs_Galerkin',abs(tight-guess),LIM,str(j));ivp.append({'mode':j,'loose':loose,'tight':tight,'vs_Galerkin':abs(tight-guess),'level_difference':abs(tight-loose)})
  print('COLD_PRUEFER',j,flush=True)
 for j in [1,481,961]:
  br=next((r for r in b if r['freq_hz']=='201' and int(r['mode'])==j),None)
  if br is None:continue
  ar=next(r for r in a if r['freq_hz']=='201' and int(r['mode'])==j);d=[]
  for N in [80000,160000]:
   dl,u,h=fd(201,z,c,N,j);d.append(dl)
   p=LOCAL/'vectors'/f'gal_f201_m{j}_N{N}.npz'
   if p.exists():
    fe=np.load(p)['phi'];grid=np.arange(N+1)*h;uu=np.interp(z,grid,np.r_[0,u.real])+1j*np.interp(z,grid,np.r_[0,u.imag]);vv=np.interp(z,grid,np.r_[0,fe.real])+1j*np.interp(z,grid,np.r_[0,fe.imag])
    ck('independent_complex_mode_shape',min(np.linalg.norm(uu-vv),np.linalg.norm(uu+vv))/np.linalg.norm(uu),.002,str(j)+'/'+str(N))
  k=np.sqrt(float(ar['K_GALERKIN_R6'])**2+(4*d[1]-d[0])/3);target=complex(float(br['K_complex_real']),float(br['K_complex_imag']));ck('independent_FD_complex_real',abs(k.real-target.real),LIM,str(j));ck('independent_FD_complex_imag',abs(k.imag-target.imag),LIM,str(j))
 # Full arrays are reconstructed independently only for actually executed C diagnostics.
 fieldpath=OUT/'PREC3_RAW_PARTIAL_FIELDS.json';shape_rows=rows(OUT/'PREC3_FIELD_SHAPE_COMPARISON.csv');fields=read(fieldpath) if fieldpath.exists() else [];geometry=rows(ROOT/'results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/ARRAY_GEOMETRY_TABLE.csv');probes=rows(ROOT/'results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/T2_N2_SENTINELS_20261013_v1.csv');stored={}
 md={}
 for ent in arc:md[ent['freq'],ent['mesh']]=mod(ROOT/ent['paths']['mod'])
 for rr in fields:
  f,mesh=rr['freq'],rr['mesh'];zz,phi,kk=md[f,mesh];phi=phi.copy();kk=kk.copy();version=rr['version']
  for br in [v for v in b if int(v['freq_hz'])==f]:
   j=int(br['mode'])-1
   if version=='C1':kk[j]=complex(float(br['K_complex_real']),float(br['K_complex_imag']))
   if version=='K_REAL_ONLY':kk[j]=float(br['K_complex_real'])+1j*kk[j].imag
   if version=='K_IMAG_ONLY':kk[j]=kk[j].real+1j*float(br['K_complex_imag'])
   if version=='C2' and f==201 and j+1 in [1,481,961]:
    v=np.load(LOCAL/'vectors'/f'gal_f201_m{j+1}_N160000.npz');xx=np.arange(len(v['phi'])+1)*float(v['h']);phi[:,j]=np.interp(zz,xx,np.r_[0,v['phi'].real])+1j*np.interp(zz,xx,np.r_[0,v['phi'].imag])
  probe=next(v for v in probes if v['sentinel_id']==rr['probe']);angle=math.radians(float(probe['abs_theta_deg']));sx=float(probe['r_km'])*1000*math.cos(angle);sy=float(probe['r_km'])*1000*math.sin(angle);si=int(np.argmin(abs(zz-float(probe['source_z_m']))));coords=sorted([v for v in geometry if v['resource']==rr['resource']],key=lambda v:int(v['channel']));vals=[]
  for coord in coords:
   rho=math.hypot(sx-float(coord['x_m']),sy-float(coord['y_m']));ri=int(np.argmin(abs(zz-float(coord['z_m']))));realparts=[];imagparts=[]
   for j,k in enumerate(kk):
    term=phi[si,j]*phi[ri,j]*cmath.sqrt(2*math.pi/(k*rho))*cmath.exp(-1j*k*rho-1j*math.pi/4);realparts.append(term.real);imagparts.append(term.imag)
   vals.append(complex(math.fsum(realparts),math.fsum(imagparts)))
  value=np.array(vals);mainvalue=np.array(rr['field_real'])+1j*np.array(rr['field_imag']);ck('field_scalar_rebuild',np.linalg.norm(mainvalue-value)/np.linalg.norm(value),1e-10,f'{f}/{rr["probe"]}/{rr["resource"]}/{mesh}/{version}');stored[f,rr['probe'],rr['resource'],mesh,version]=value
 for rr in shape_rows:
  f=int(rr['freq_hz']);key=(f,rr['probe'],rr['resource'])
  for version in ['C0','C1','C2']:
   x=stored[(*key,160001,version)];y=stored[(*key,240001,version)];q,_=np.linalg.qr(y[:,None]);error=np.linalg.norm(x-q@(q.conj().T@x))/np.linalg.norm(x);ck('independent_QR_mesh_shape',abs(error-float(rr[version+'_mesh_shape'])),1e-10,str(key)+'/'+version)
 with (OUT/'PREC3_COLD_CHECKS.csv').open('w',encoding='utf-8',newline='') as file:
  w=csv.DictWriter(file,fieldnames=list(checks[0]),lineterminator='\n');w.writeheader();w.writerows(checks)
 save(OUT/'PREC3_COLD_PRUEFER.json',ivp)
 fail=sum(not r['PASS'] for r in checks);report={'PASS':fail==0,'checks':len(checks),'FAIL':fail,'new_Pruefer_roots':len(ivp),'raw_fields_rebuilt':len(fields),'wall_s':time.monotonic()-began,'scope':'independent roots3; extrapolation15 as completed; independent FD complex3 as completed; field scalar/QR only executed branch','R4_percent':0}
 save(OUT/'PREC3_COLD_AUDIT.json',report);print(json.dumps(report),flush=True)
 if fail:raise RuntimeError('COLD_FAIL_NO_RETRY')
if __name__=='__main__':main()
