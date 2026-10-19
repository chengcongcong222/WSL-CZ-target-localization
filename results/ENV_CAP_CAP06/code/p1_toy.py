from common import *
def basis():
 t=np.linspace(-1,1,9);c=np.array([1.,-1.,2.,-2.]);d=np.array([1.,2.,-1.,-2.])
 A=np.stack([np.cos(np.pi*t)[:,None]*c,np.sin(2*np.pi*t)[:,None]*d,np.sin(np.pi*t)[:,None]*d,t[:,None]**2*c,t[:,None]**3*d],axis=-1).reshape(36,5)
 E=(np.cos(3*np.pi*t)[:,None]*np.array([1.,1.,-1.,-1.])).reshape(36,1)
 return t,A,E

def generate():
 public=RAW/'p1/receiver_public';private=RAW/'p1/truth_private';public.mkdir(parents=True);private.mkdir(parents=True)
 t,A,E=basis();logshape=(A@np.array([.4,.2,.3,.15,.1])+E[:,0]*.1).reshape(9,4);h=10**(np.column_stack((np.zeros(9),logshape))/20)
 for rep in range(32):
  for rho in [0.,.9]:
   for snr in [-20,0,10]:
    for gi,G in enumerate(['G0','G1','G2','G3','G4']):
     rng=np.random.default_rng(2026101901+rep);source=cn(rng,(9,128))*10**((np.sin(t*np.pi)+.3*t)[:,None]/20)
     n=cn(rng,(9,5,128));bg=cn(rng,(9,5,128))
     for k in range(1,9):n[k]=rho*n[k-1]+np.sqrt(1-rho*rho)*n[k];bg[k]=rho*bg[k-1]+np.sqrt(1-rho*rho)*bg[k]
     offset=rng.uniform(-1,1,5);slope=rng.uniform(-.5,.5,5)
     gain_db=np.tile(offset,(9,1))
     if G=='G2':gain_db+=t[:,None]*slope
     if G=='G3':gain_db+=rng.normal(size=(9,5))
     gain=10**(gain_db/20);y=np.sqrt(10**(snr/10))*source[:,None,:]*h[:,:,None]*gain[:,:,None]+n
     if G=='G4':
      agc=10**((2*np.cos(t*np.pi))[:,None,None]/20);y*=agc;bg*=agc
     token=f'{rep}_{rho}_{snr}_{G}'
     np.savez_compressed(public/(token+'.npz'),samples=y.astype(np.complex64),background=bg.astype(np.complex64),gain_metadata=gain if G=='G0' else np.ones((9,5)),time=t,gain_model=np.array(G),metadata_kind=np.array('IDEAL_INDEPENDENT_CALIBRATION_CONTROL' if G=='G0' else 'NO_DIFFERENTIAL_GAIN_VALUES'))
     np.savez_compressed(private/(token+'.npz'),ideal_shape_db=logshape,gain_db=gain_db,rep=np.array(rep),rho=np.array(rho),snr=np.array(snr),model=np.array('MATH_TOY_NOT_CZ'))
 print('P1 generated 960 toy configurations',flush=True)

def extract_only():
 # Separate process dispatch; never opens private labels.
 import sys
 def guard(event,args):
  if event=='open' and isinstance(args[0],str) and 'truth_private' in args[0]:raise PermissionError('P1 receiver private read')
 sys.addaudithook(guard)
 dest=RAW/'p1/products';dest.mkdir()
 for p in sorted((RAW/'p1/receiver_public').glob('*.npz')):
  with np.load(p) as a:
   power=np.mean(abs(a['samples'].astype(complex))**2,axis=-1)-np.mean(abs(a['background'].astype(complex))**2,axis=-1)
   power/=a['gain_metadata']**2;valid=(power[:,1:]>0)&(power[:,0,None]>0)
   ratio=np.full((9,4),np.nan);np.log10(np.divide(power[:,1:],power[:,0,None],out=np.ones((9,4)),where=valid),out=ratio,where=valid);ratio*=10
   np.savez_compressed(dest/p.name,ratio=ratio,valid=valid,subtracted_power=power)

def evaluate():
 t,A,E=basis();rows=[];infos=[]
 for G in ['G0','G1','G2','G3','G4']:
  cols=[]
  if G in ['G1','G2','G4']:cols.append(np.tile(np.eye(4),(9,1)))
  if G=='G2':cols.append(np.repeat(t,4)[:,None]*np.tile(np.eye(4),(9,1)))
  if G=='G3':cols.append(np.eye(36))
  N=np.column_stack([E]+cols);u,s,_=np.linalg.svd(N,full_matrices=False);q=u[:,s>1e-10];F=A.T@(np.eye(36)-q@q.T)@A
  eff=[]
  for j in range(5):
   other=np.column_stack((N,np.delete(A,j,axis=1)));u,s,_=np.linalg.svd(other,full_matrices=False);v=A[:,j]-u[:,s>1e-10]@(u[:,s>1e-10].T@A[:,j]);eff.append(float(v@v))
  infos.append(dict(gain=G,rank=int(np.linalg.matrix_rank(F,tol=1e-9)),profiled_state_information=eff,scope='DIMENSIONLESS_TOY_MECHANISM_ONLY',all_state_and_environment_nuisance=True))
  for rep in range(32):
   for rho in [0.,.9]:
    for snr in [-20,0,10]:
     token=f'{rep}_{rho}_{snr}_{G}';a=np.load(RAW/'p1/products'/(token+'.npz'));truth=np.load(RAW/'p1/truth_private'/(token+'.npz'))
     valid=a['valid'];ratio=a['ratio'];v=valid.ravel();known=truth['ideal_shape_db'].ravel()
     if v.any():
      B=N[v];u,s,_=np.linalg.svd(B,full_matrices=False);qq=u[:,s>1e-10];error=ratio.ravel()[v]-known[v];pr=error-qq@(qq.T@error);rms=float(np.sqrt(np.mean(pr**2)))
     else:rms=None
     rows.append(dict(token=token,gain=G,snr=snr,rho=rho,repeat=rep,valid_contrasts=int(v.sum()),total_contrasts=36,nonpositive_power=int((a['subtracted_power']<=0).sum()),profiled_ratio_error_RMS_db=rms,scope='SYNTHETIC_MECHANISM_ONLY'))
 tab('P1_EXTRACTION_CASES.csv',rows);write('P1_NUISANCE_INFO.json',infos)
 tab('P1_NUISANCE_INFO.csv',[dict(gain=x['gain'],rank=x['rank'],r=x['profiled_state_information'][0],theta=x['profiled_state_information'][1],z=x['profiled_state_information'][2],v=x['profiled_state_information'][3],psi=x['profiled_state_information'][4],scope=x['scope']) for x in infos])
 write('P1_POWER_FORWARD_GATE.json',dict(G1a='MATHEMATICAL_POWER_AND_GAIN_CONTROL_SUPPORTED',G1b='O2_FIRST_CZ_FORWARD_NOT_ADMITTED',reason='Archived SCS I-class energy lacks independently certified matched CW or defined stochastic-path receiver power; no qualified new provider adopted; old T2 closed',O2_CZ_EXTRACTED_OBSERVABLE='NOT_EVALUATED',CZ_depth_performance='NOT_EVALUATED',toy_configurations=960,G3_free_gain_zero=infos[3]['rank']==0,scope='No inference from toy to CZ; no H3 complex-information transfer',external_solver_calls=0))
 tab('P1_GLOBAL_SUPPORT.csv',[dict(config=x,status='NOT_EVALUATED',reason='O2_FIRST_CZ_FORWARD_NOT_ADMITTED') for x in ['S1','S3']])
 tab('P1_SIGNAL_AND_GAIN_CONTROLS.csv',[dict(gain=x['gain'],free_gain_zero=(x['rank']==0),scope=x['scope']) for x in infos])
 idx=[]
 for p in sorted((RAW/'p1').rglob('*.npz')):idx.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
 tab('P1_RAW_DATA_INDEX.csv',idx)
if __name__=='__main__':
 import sys
 {'generate':generate,'extract':extract_only,'evaluate':evaluate}[sys.argv[1]]()
