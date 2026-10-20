"""Public-only estimator process. Original T3 estimate called without modification."""
from common import *
import sys,time,os
model=module('t3_w1_model')
# Preload scientific dependencies before restricting reads. No private imports.
access=set();denied=[]
allowed_roots=[(DATA/'public').resolve(),REF.resolve(),O.resolve()]
allowed_outputs=[(DATA/'terminals.npy').resolve(),(DATA/'worker_progress.json').resolve(),(DATA/'worker_exception.json').resolve(),(O/'ESTIMATOR_FILE_ACCESS.json').resolve()]
def gate(event,args):
 if event!='open':return
 name=args[0]
 if not isinstance(name,(str,bytes,os.PathLike)):return
 path=Path(os.fsdecode(name)).resolve();mode=args[1];flags=args[2]
 writing=(isinstance(mode,str) and any(c in mode for c in 'wax+')) or (isinstance(flags,int) and bool(flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT)))
 ok=path in allowed_outputs if writing else (path.is_relative_to(allowed_roots[0]) or path.is_relative_to(REF.resolve()))
 # Model reads no data except the supplied public arrays. Metadata access is allowlisted.
 if not writing and path in [O/'PUBLIC_ARCHIVE_SEAL.json',O/'PUBLIC_RECORD_IDENTITIES.csv',O/'PRODUCTION_CODE_FREEZE.json']:ok=True
 if not writing and path.is_relative_to(O.resolve()) and path.suffix=='.py':ok=True
 if not writing and path.suffix in ['.py','.pyc','.pyd','.dll'] and (path.is_relative_to(Path(sys.base_prefix).resolve()) or path.is_relative_to(Path('D:/ProjectStorage/WSL-CZ/T2_PREC3_FUNC2/runtime_dependencies').resolve())):ok=True
 if not ok:denied.append(str(path));raise PermissionError('ESTIMATOR_INPUT_DENIED:'+str(path))
 access.add(('WRITE:' if writing else 'READ:')+str(path))

def main():
 sys.addaudithook(gate)
 for item in js(O/'PRODUCTION_CODE_FREEZE.json')['files']:assert sha(item['path'])==item['sha256']
 seal=js(O/'PUBLIC_ARCHIVE_SEAL.json');expected={r['path']:r['sha256'] for r in seal['files']}
 out=np.lib.format.open_memmap(DATA/'terminals.npy',mode='w+',dtype=DT,shape=(48000,));out[:]=np.zeros((),dtype=DT);done=0;start=time.monotonic()
 try:
  for qi,q in enumerate(['Q1','Q2']):
   for bi,b in enumerate(['B0','B3']):
    for ci in range(12):
     path=recordpath(q,b,ci);assert sha(path)==expected[str(path)]
     with np.load(path,allow_pickle=False) as archive:t=archive['times'];P=archive['reported'];y=archive['bearing']
     for rep in range(500):
      for pi,p in enumerate(['P0','P1']):
       ix=index(qi,bi,ci,rep,pi);qual=quality(q,p);h,z=model.estimate(t,P[rep],y[rep],qual)
       row=np.zeros((),dtype=DT);row['executed']=1;row['accepted']=z['accepted'];row['h']=h;row['initial']=z['initial_h']
       for dest,key in [('solver_status','solver_status'),('nfev','nfev'),('cost','cost'),('optimality','optimality'),('fallback','initializer_fallback'),('intersections','initializer_intersections'),('rank','rank_augmented_likelihood_jac'),('ratio','smallest_largest_jac_ratio'),('rms','mean_residual_RMS_standard')]:row[dest]=z[key]
       c,d=qual['common_bound_deg'],qual['half_diff_bound_deg'];low=np.array([45,-5,1,-15,-c,-d]);span=np.array([15,10,2,30,2*c,2*d]);unit=(h-low)/span;row['near']=(unit<=1e-6)|(unit>=1-1e-6)
       original=np.array([CFG['quality'][q]['common_bound_deg'],CFG['quality'][q]['half_diff_bound_deg']]);row['active_bias_distance']=np.array([c,d])-abs(h[4:]);row['original_bias_distance']=original-abs(h[4:]);row['outside_original']=abs(h[4:])>original
       out[ix]=row;done+=1
     out.flush();dump(DATA/'worker_progress.json',dict(done=done,seconds=time.monotonic()-start,last=[q,b,ci]));print('FIT',done,'/48000',flush=True)
 except Exception as e:
  out.flush();dump(DATA/'worker_exception.json',dict(type=type(e).__name__,message=str(e),completed=done));raise
 finally:dump(O/'ESTIMATOR_FILE_ACCESS.json',dict(enforcement='Python audit-hook read/write allowlist; process receives public arrays only; not an OS sandbox claim',paths=sorted(access),denied=denied,private_reads=0 if not denied else 'CHECK_DENIED',completed=done))
if __name__=='__main__':main()
