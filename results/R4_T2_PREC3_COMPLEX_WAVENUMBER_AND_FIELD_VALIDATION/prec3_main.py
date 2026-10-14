from pathlib import Path
import json,csv,hashlib,time,math,struct,sys
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh,eigs
from scipy.linalg import solve_banded,eigh_tridiagonal
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3');LIMIT=1e-4/60000
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def table(p,rr,fields=None):
 fields=fields or list(dict.fromkeys(k for r in rr for k in r))
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rr)
def tick(stage,detail):print(stage,detail,flush=True)
def env(p):
 ls=Path(p).read_text(encoding='utf-8-sig').splitlines();j=next(i for i,l in enumerate(ls) if l.startswith("'R'"));a=np.array([list(map(float,l.split()[:6])) for l in ls[5:j]])
 assert a.shape==(251,6) and a[0,0]==0 and a[-1,0]==5000 and np.all(a[:,3]==1) and np.all(a[:,2]==0) and np.all(a[:,4:]==0)
 assert ls[3].strip().strip("'")=='CVWT'
 return a[:,0],a[:,1]
def thorp(f):
 x=(f/1000)**2
 return (0.0033+0.11*x/(1+x)+44*x/(4100+x)+0.0003*x)/8685.8896
def tri_mv(d,e,v):
 a=d*v;a[:-1]+=e*v[1:];a[1:]+=e*v[:-1];return a
def gal(f,z,c,N,loss=0):
 H=float(z[-1]);h=H/N;edges=np.arange(N)*h;w=2*np.pi*f;cp=c.astype(complex)+1j*loss*thorp(f)*c*c/w
 q=np.zeros((N,3),complex);gx,gw=np.polynomial.legendre.leggauss(4)
 for xi,wi in zip((gx+1)/2,gw/2):
  zz=edges+xi*h;cz=np.interp(zz,z,cp.real)+1j*np.interp(zz,z,cp.imag);V=(w/cz)**2
  q[:,0]+=h*wi*V*(1-xi)**2;q[:,1]+=h*wi*V*xi*(1-xi);q[:,2]+=h*wi*V*xi**2
 ad=q[:,2].copy();ad[:-1]+=q[1:,0];ad-=2/h;ad[-1]+=1/h;ae=q[1:,1]+1/h
 md=np.full(N,2*h/3);md[-1]=h/3;me=np.full(N-1,h/6)
 if loss==0:ad=ad.real;ae=ae.real
 A=diags([ae,ad,ae],[-1,0,1],format='csc');M=diags([me,md,me],[-1,0,1],format='csc')
 return A,M,ad,ae,md,me,h
def greater(ad,ae,md,me,x):
 d=(ad-x*md).tolist();e=((ae-x*me)**2).tolist();p=d[0];count=int(p>0)
 for i in range(1,len(d)):
  if abs(p)<1e-280:p=-1e-280
  p=d[i]-e[i-1]/p;count+=int(p>0)
 return count
def solve_real(f,z,c,N,mode):
 A,M,ad,ae,md,me,h=gal(f,z,c,N);lo=0.;hi=(2*np.pi*f/min(c))**2
 assert greater(ad,ae,md,me,lo)>=mode and greater(ad,ae,md,me,hi)<mode
 for _ in range(44):
  mid=(lo+hi)/2
  if greater(ad,ae,md,me,mid)>=mode:lo=mid
  else:hi=mid
 ev,v=eigsh(A,k=1,M=M,sigma=(lo+hi)/2,which='LM',v0=np.ones(N),ncv=10,maxiter=2000,tol=1e-12);lam=float(ev[0]);v=v[:,0];v/=np.sqrt(v@tri_mv(md,me,v))
 residual=tri_mv(ad,ae,v)-lam*tri_mv(md,me,v)
 band=np.zeros((3,N));band[0,1:]=me;band[1]=md;band[2,:-1]=me
 dual=solve_banded((1,1),band,residual);rk=float(np.sqrt(max(0.,residual@dual))/(2*np.sqrt(lam)))
 width=max(1e-10,hi-lo)*10
 identity=greater(ad,ae,md,me,lam-width)==mode and greater(ad,ae,md,me,lam+width)==mode-1
 return lam,v,rk,identity,(A,M,ad,ae,md,me,h)
def extrap(v):
 a=np.array(v);r2=(4*a[1:]-a[:-1])/3;r4=(16*r2[1:]-r2[:-1])/15;r6=(64*r4[1]-r4[0])/63
 return r6,float(abs(np.sqrt(r4[1])-np.sqrt(r4[0])))
def mod(p):
 b=Path(p).read_bytes();rec=4*struct.unpack_from('<i',b)[0];nt=struct.unpack_from('<i',b,92)[0];m=struct.unpack_from('<i',b,5*rec)[0]
 depth=np.frombuffer(b,dtype='<f4',count=nt,offset=4*rec).astype(float);phi=np.stack([np.frombuffer(b,dtype='<c8',count=nt,offset=(7+j)*rec) for j in range(m)],axis=1).astype(complex);k=np.frombuffer(b,dtype='<c8',count=m,offset=(7+m)*rec).astype(complex)
 assert np.isfinite(phi).all() and np.isfinite(k).all() and np.all(k.real>0) and np.all(k.imag<=0)
 return depth,phi,k
def complex_solve(f,z,c,N,loss,real_data):
 lam,v,*_=real_data;A,M,ad,ae,md,me,h=gal(f,z,c,N,loss);res=tri_mv(ad,ae,v)-lam*tri_mv(md,me,v);dlam=v@res
 ev,u=eigs(A,k=1,M=M,sigma=lam+dlam,which='LM',v0=v.astype(complex),ncv=10,maxiter=2000,tol=1e-12);lc=complex(ev[0]);u=u[:,0];u/=np.sqrt(u@tri_mv(md,me,u))
 r=tri_mv(ad,ae,u)-lc*tri_mv(md,me,u);rho=float(np.linalg.norm(r)/(max(np.linalg.norm(tri_mv(md,me,u)),1e-300)*2*abs(np.sqrt(lc))))
 return dlam,lc-lam,u,rho
def sample_phi(v,h,depths):
 return np.interp(depths,np.arange(len(v)+1)*h,np.r_[0.,v.real])+1j*np.interp(depths,np.arange(len(v)+1)*h,np.r_[0.,v.imag])
def fd_complex(f,z,c,N,mode,k_guess,loss):
 h=5000/N;depth=np.arange(1,N+1)*h;w=2*np.pi*f;cz=np.interp(depth,z,c);cp=c.astype(complex)+1j*loss*thorp(f)*c*c/w
 cc=np.interp(depth,z,cp.real)+1j*np.interp(depth,z,cp.imag)
 real=(w/cz)**2-2/h**2;e=np.full(N-1,1/h**2);e[-1]*=np.sqrt(2)
 lam,v=eigh_tridiagonal(real,e,select='i',select_range=(N-mode,N-mode),eigvals_only=False,lapack_driver='stebz');lam=float(lam[0]);v=v[:,0]
 d=(w/cc)**2-2/h**2;A=diags([e,d,e],[-1,0,1],format='csc');pert=v@((d-real)*v)
 ev,u=eigs(A,k=1,sigma=lam+pert,which='LM',v0=v.astype(complex),ncv=10,maxiter=2000,tol=1e-12);u=u[:,0];u/=np.sqrt(u@u);phi=u/np.sqrt(h);phi[-1]*=np.sqrt(2)
 return complex(ev[0])-lam,phi,h
def controls():
 rr=[];f=201.;H=5000.;C=1500.;w=2*np.pi*f;mode=5;mu=(mode-.5)*np.pi/H;alpha=thorp(f)
 for loss in [0.,.25,.5,1.]:
  cp=C+1j*loss*alpha*C*C/w;exact=np.sqrt((w/cp)**2-mu**2)
  A,M,*_=gal(f,np.array([0.,H]),np.array([C,C]),2000,loss);kk,_=eigs(A.astype(complex),k=1,M=M,sigma=exact**2,which='LM',v0=np.sin(mu*np.arange(1,2001)*H/2000).astype(complex),ncv=10,maxiter=2000,tol=1e-12);numeric=np.sqrt(kk[0])
  rr.append({'control':'constant_waveguide','loss_scale':loss,'k_real':float(exact.real),'k_imag':float(exact.imag),'decay_at60km':float(abs(np.exp(-1j*exact*60000))),'numeric_k_real':float(numeric.real),'numeric_k_imag':float(numeric.imag),'numeric_vs_exact_delta':float(abs(numeric-exact)),'PASS':bool(exact.imag<=0 and abs(numeric-exact)<=LIMIT and (loss!=0 or exact.imag==0))})
 rr.append({'control':'unit_Neper_conversion','loss_scale':1.,'k_real':0.,'k_imag':0.,'decay_at60km':0.,'PASS':bool(abs(20/np.log(10)-8.6858896)<1e-7)})
 return rr

def main():
 cfg=read(OUT/'PREC3_DESIGN_FREEZE.json');base=ROOT/'results/R4_T2_PREC2_ESTD_EIGENVALUE_PRECISION_FEASIBILITY';old=rows(base/'PREC2_LOSSLESS_SUBSET_FDM.csv');arc=read(base/'T2_PREC2_EXPERIMENT_SPEC_20261014_v1.json')['archive_sets'];z,c=env(ROOT/arc[0]['paths']['env']);cache={};a=[];b=[];raw=[]
 decision={'stage':'T2-PREC3','Gate_A':'NOT_EVALUATED','Gate_B':'NOT_EVALUATED','Gate_C':'NOT_EVALUATED','A_planned':15,'B_planned':15,'C_planned':60,'R4_percent':0,'STOP':True}
 write(OUT/'PREC3_DECISION.json',decision)
 def stop(status,reason):
  decision.update(status=status,stop_reason=reason,A_completed=len(a),B_completed=len(b));write(OUT/'PREC3_DECISION.json',decision);tick('SCIENCE_STOP',status)
 for source in cfg['inputs']:
  if sha(ROOT/source['path'])!=source['sha256']:raise RuntimeError('FROZEN_INPUT_IDENTITY_FAIL '+source['path'])
 for oldrow in old:
  f=int(oldrow['freq_hz']);j=int(oldrow['mode']);vals=[];errs=[];ids=[]
  for N in cfg['A_grids']:
   data=solve_real(f,z,c,N,j);lam,v,rk,ident,ops=data;vals.append(lam);errs.append(rk);ids.append(ident)
   raw.append({'stage':'A','freq_hz':f,'mode':j,'N':N,'lambda_real':lam,'k_real':float(np.sqrt(lam)),'dual_discrete_residual_k_bound':rk,'mode_identity':ident})
   if N in cfg['B_grids']:
    cache[f,j,N]=data;vp=LOCAL/'vectors'/f'real_f{f}_m{j}_N{N}.npz';vp.parent.mkdir(exist_ok=True);np.savez_compressed(vp,phi=v,h=5000/N)
   if not ident or rk>LIMIT:
    table(OUT/'PREC3_RAW_EIGENVALUES.csv',raw);stop('PREC3_A_NOT_ESTABLISHED','A_DISCRETE_IDENTITY_OR_RESIDUAL_FAIL');return
  lr6,gap=extrap(vals);kr=float(np.sqrt(lr6));prior=float(oldrow['K_real_lossless_R6']);passed=bool(gap<=LIMIT and abs(kr-prior)<=LIMIT)
  a.append({'freq_hz':f,'mode':j,'K_MOD':oldrow['K_mod'],'K_PRT':oldrow['K_prt'],'K_PREC2_R6':prior,'K_GALERKIN_R6':kr,'GALERKIN_R4_gap':gap,'GALERKIN_vs_PREC2':abs(kr-prior),'max_dual_discrete_residual_k_bound':max(errs),'mode_identity_all_grids':all(ids),'phase_delta_at60km':60000*(kr-prior),'PASS':passed})
  table(OUT/'PREC3_REAL_K_INDEPENDENT_VALIDATION.csv',a);table(OUT/'PREC3_RAW_EIGENVALUES.csv',raw);tick('A',f'{f}/{j} PASS={passed}')
  if not passed:stop('PREC3_A_NOT_ESTABLISHED','A_CONVERGENCE_OR_INDEPENDENT_DISAGREEMENT');return
 priorivp=rows(base/'PREC2_PRUFER_COLD_ROOTS.csv')
 for r in priorivp:
  ar=next(v for v in a if v['freq_hz']==201 and v['mode']==int(r['mode']))
  if abs(ar['K_GALERKIN_R6']-float(r['K_PRUFER_tight']))>LIMIT:stop('PREC3_A_NOT_ESTABLISHED','A_ARCHIVED_PRUFER_DISAGREEMENT');return
 decision['Gate_A']='PASS_SUBSET_CONDITIONAL';write(OUT/'PREC3_DECISION.json',decision)
 con=controls();table(OUT/'PREC3_LOSS_CONTROLS.csv',con)
 if not all(r['PASS'] for r in con):decision['Gate_B']='FAIL';stop('PREC3_B_COMPLEX_K_NOT_CERTIFIED','LOSS_UNIT_OR_SIGN_CONTROL_FAIL');return
 complexcache={}
 for ar in a:
  f,j=ar['freq_hz'],ar['mode'];kr=ar['K_GALERKIN_R6'];lp=kr*kr;details=[];results={}
  for scale in cfg['loss_scales']:
   ds=[];ps=[];residual=[]
   for N in cfg['B_grids']:
    d,l,u,rho=complex_solve(f,z,c,N,scale,cache[f,j,N]);ds.append(l);ps.append(d);residual.append(rho);h=5000/N
    raw.append({'stage':'B','freq_hz':f,'mode':j,'N':N,'loss_scale':scale,'lambda_shift_complex_real':l.real,'lambda_shift_complex_imag':l.imag,'lambda_shift_pert_real':d.real,'lambda_shift_pert_imag':d.imag,'discrete_residual_k_indicator':rho})
    if scale==1:
     p=LOCAL/'vectors'/f'gal_f{f}_m{j}_N{N}.npz';p.parent.mkdir(exist_ok=True);np.savez_compressed(p,phi=u,h=h);complexcache[f,j,N]=p
    if rho>LIMIT:
     table(OUT/'PREC3_RAW_EIGENVALUES.csv',raw);decision['Gate_B']='FAIL';stop('PREC3_B_COMPLEX_K_NOT_CERTIFIED','B_DISCRETE_RESIDUAL_FAIL');return
   dre=(4*ds[1]-ds[0])/3;pre=(4*ps[1]-ps[0])/3;kcomp=complex(np.sqrt(lp+dre));kpert=complex(np.sqrt(lp+pre));kf=complex(np.sqrt(lp+ds[1]));kc=complex(np.sqrt(lp+ds[0]));delta=kcomp-kpert;conv=kf-kc
   results[scale]=(kcomp,kpert);details.append({'scale':scale,'k_complex':kcomp,'k_pert':kpert,'real_disagreement':abs(delta.real),'imag_disagreement':abs(delta.imag),'real_shift_grid_difference':abs(conv.real),'imag_shift_grid_difference':abs(conv.imag),'perturb_error_abs':abs(delta)})
  e1,ehalf,eq=[results[s][0]-results[s][1] for s in cfg['loss_scales']]
  scaling=abs(ehalf)<=.75*abs(e1)+1e-11 and abs(eq)<=.75*abs(ehalf)+1e-11
  kcomp,kpert=results[1.];stored=next(rr for rr in arc if rr['freq']==f and rr['mesh']==160001);km=mod(ROOT/stored['paths']['mod'])[2][j-1]
  passnum=all(d['real_disagreement']<=LIMIT and d['imag_disagreement']<=LIMIT and d['real_shift_grid_difference']<=LIMIT and d['imag_shift_grid_difference']<=LIMIT and d['k_complex'].imag<0 for d in details)
  row={'freq_hz':f,'mode':j,'K_complex_real':kcomp.real,'K_complex_imag':kcomp.imag,'K_pert_real':kpert.real,'K_pert_imag':kpert.imag,'original_MOD_imag':km.imag,'MOD_imag_vs_independent':abs(km.imag-kcomp.imag),'max_real_path_disagreement':max(d['real_disagreement'] for d in details),'max_imag_path_disagreement':max(d['imag_disagreement'] for d in details),'max_real_shift_grid_difference':max(d['real_shift_grid_difference'] for d in details),'max_imag_shift_grid_difference':max(d['imag_shift_grid_difference'] for d in details),'loss_scaling_with_roundoff_floor_PASS':scaling,'phase_path_difference_at60km':60000*(kpert.real-kcomp.real),'amplitude_path_ratio_at60km':float(np.exp(60000*(kpert.imag-kcomp.imag))),'PASS':passnum and scaling and abs(km.imag-kcomp.imag)<=LIMIT}
  for d in details:raw.append({'stage':'B_PATH_COMPARISON','freq_hz':f,'mode':j,'loss_scale':d['scale'],'perturb_error_abs':d['perturb_error_abs'],'real_disagreement':d['real_disagreement'],'imag_disagreement':d['imag_disagreement']})
  b.append(row);table(OUT/'PREC3_COMPLEX_K_VALIDATION.csv',b);table(OUT/'PREC3_RAW_EIGENVALUES.csv',raw);tick('B',f'{f}/{j} PASS={row["PASS"]}')
  if not row['PASS']:decision['Gate_B']='FAIL';stop('PREC3_B_COMPLEX_K_NOT_CERTIFIED','B_PATH_GRID_LOSS_SCALING_OR_ARCHIVED_IMAG_DISAGREEMENT');return
 phicontrol=[]
 for j in [1,481,961]:
  ar=next(r for r in a if r['freq_hz']==201 and r['mode']==j);br=next(r for r in b if r['freq_hz']==201 and r['mode']==j);dk=[]
  for N in cfg['B_grids']:
   dl,ph,h=fd_complex(201,z,c,N,j,ar['K_GALERKIN_R6'],1.);dk.append(dl);fep=np.load(complexcache[201,j,N]);pfe=sample_phi(fep['phi'],h,z);pfd=sample_phi(ph,h,z);err=min(np.linalg.norm(pfe-pfd),np.linalg.norm(pfe+pfd))/np.linalg.norm(pfd)
   phicontrol.append({'freq_hz':201,'mode':j,'N':N,'modal_function_FE_vs_FD_relative':float(err),'PASS':bool(err<=.002)})
   np.savez_compressed(LOCAL/'vectors'/f'fd_f201_m{j}_N{N}.npz',phi=ph,h=h)
  knew=np.sqrt(ar['K_GALERKIN_R6']**2+(4*dk[1]-dk[0])/3);difference=abs(knew-complex(br['K_complex_real'],br['K_complex_imag']));phicontrol.append({'freq_hz':201,'mode':j,'N':0,'FE_vs_FD_complex_k':float(difference),'PASS':bool(difference<=LIMIT)})
  table(OUT/'PREC3_INDEPENDENT_COMPLEX_AND_PHI.csv',phicontrol)
  if not all(r['PASS'] for r in phicontrol):decision['Gate_B']='FAIL';stop('PREC3_B_COMPLEX_K_NOT_CERTIFIED','B_INDEPENDENT_FD_OR_MODE_FUNCTION_FAIL');return
 decision['Gate_B']='PASS_SUBSET_CONDITIONAL';write(OUT/'PREC3_DECISION.json',decision)
 run_fields(arc,a,b,complexcache,decision)

def run_fields(arc,a,b,complexcache,decision):
 geometry=rows(ROOT/'results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/ARRAY_GEOMETRY_TABLE.csv');probes=rows(ROOT/'results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/T2_N2_SENTINELS_20261013_v1.csv')[:5];cc=[];storedfields=[]
 def field(zz,pp,kk,probe,coords):
  def index(v):
   ix=int(np.argmin(abs(zz-v)))
   if abs(zz[ix]-v)>1e-5:raise RuntimeError('ARCHIVED_MODAL_DEPTH_MISSING')
   return ix
  xy=float(probe['r_km'])*1000*np.array([np.cos(np.deg2rad(float(probe['abs_theta_deg']))),np.sin(np.deg2rad(float(probe['abs_theta_deg'])))])
  rho=np.linalg.norm(xy-coords[:,:2],axis=1);weights=pp[index(float(probe['source_z_m']))][:,None]*pp[[index(v) for v in coords[:,2]]].T
  kernel=np.sqrt(2*np.pi/(kk[:,None]*rho))*np.exp(-1j*kk[:,None]*rho-1j*np.pi/4)
  return np.sum(weights*kernel,axis=0)
 def shape(x,y):return float(np.linalg.norm(x-y*(np.vdot(y,x)/np.vdot(y,y)))/np.linalg.norm(x))
 for f in [201,235,283]:
  modes=[r for r in b if r['freq_hz']==f]
  for probe in probes:
   for res in ['A0_0','A_P4M','B1_V4','B2_DUP']:
    coords=np.array([[float(q[k]) for k in ['x_m','y_m','z_m']] for q in sorted([q for q in geometry if q['resource']==res],key=lambda q:int(q['channel']))]);ff={}
    for mesh in [160001,240001]:
     ent=next(v for v in arc if v['freq']==f and v['mesh']==mesh);zz,pp,kk=mod(ROOT/ent['paths']['mod']);newk=kk.copy();realonly=kk.copy();imagonly=kk.copy();newp=pp.copy()
     for br in modes:
      j=br['mode']-1;newk[j]=complex(br['K_complex_real'],br['K_complex_imag']);realonly[j]=br['K_complex_real']+1j*kk[j].imag;imagonly[j]=kk[j].real+1j*br['K_complex_imag']
      if f==201 and j+1 in [1,481,961]:
       v=np.load(complexcache[f,j+1,160000]);newp[:,j]=sample_phi(v['phi'],float(v['h']),zz)
     versions={'C0':(pp,kk),'C1':(pp,newk),'C2':(newp,kk),'K_REAL_ONLY':(pp,realonly),'K_IMAG_ONLY':(pp,imagonly)}
     for version,(pv,kv) in versions.items():
      val=field(zz,pv,kv,probe,coords);ff[mesh,version]=val;storedfields.append({'freq':f,'probe':probe['sentinel_id'],'resource':res,'mesh':mesh,'version':version,'field_real':val.real.tolist(),'field_imag':val.imag.tolist()})
    row={'freq_hz':f,'probe':probe['sentinel_id'],'resource':res,'certified_k_replaced_modes':len(modes),'independent_phi_replaced_modes':3 if f==201 else 0,'C3':'NOT_EVALUATED_NO_COMPLETE_REFERENCE','C0_mesh_shape':shape(ff[160001,'C0'],ff[240001,'C0']),'C1_mesh_shape':shape(ff[160001,'C1'],ff[240001,'C1']),'C2_mesh_shape':shape(ff[160001,'C2'],ff[240001,'C2'])}
    for version in versions:
     row[version+'_norm']=float(np.linalg.norm(ff[160001,version]));row[version+'_vs_C0_shape']=shape(ff[160001,version],ff[160001,'C0']);row[version+'_vs_C0_raw_relative']=float(np.linalg.norm(ff[160001,version]-ff[160001,'C0'])/np.linalg.norm(ff[160001,'C0']))
    row['weak_field_fixed_design_flag']=row['C0_norm']<1e-6;row['remaining_full_field_error']='UNKNOWN_NO_COMPLETE_REFERENCE';cc.append(row)
    table(OUT/'PREC3_FIELD_SHAPE_COMPARISON.csv',cc);write(OUT/'PREC3_RAW_PARTIAL_FIELDS.json',storedfields)
 decision.update(Gate_C='NOT_CERTIFIED_COMPLETE_REFERENCE_UNAVAILABLE',C_completed=len(cc),A_completed=len(a),B_completed=len(b),status='PREC3_C_FIELD_NOT_CERTIFIED',stop_reason='PARTIAL_SUBSTITUTIONS_ONLY_C3_NOT_EVALUATED');write(OUT/'PREC3_DECISION.json',decision);tick('SCIENCE_STOP',decision['status'])
if __name__=='__main__':main()
