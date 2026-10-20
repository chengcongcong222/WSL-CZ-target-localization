"""No-panel/no-optimization implementation controls then production SHA freeze."""
from common import *
import ast,sys,scipy,psutil,os
from independent_math import residual_and_jac,initial

def main():
 verify_binding();assert not DATA.exists(),'DATA_DIRECTORY_ALREADY_EXISTS';assert not (O/'PRODUCTION_CODE_FREEZE.json').exists()
 checks=[]
 def ck(n,ok,detail=''):checks.append(dict(check=n,PASS=bool(ok),detail=detail))
 ck('single_BLAS_thread_environment',all(os.environ.get(k)=='1' for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']))
 ck('numpy_version',np.__version__=='1.26.4');ck('scipy_version',scipy.__version__=='1.15.3')
 model=module('t3_w1_model');t=np.array([0.,300.,600.,900.,1200.]);P=np.zeros((5,2,2));P[:,0,0]=[0,600,1200,1779.55,2359.11];P[:,0,1]=[0,0,0,155.29,310.58];P[:,1]=P[:,0]+[3000,15000];y=np.tile([.012,-.26],(5,1));physical=np.array([51.234,.713,1.813,-3.17,.01,-.005])
 for q in ['Q1','Q2']:
  m0=model.Model(t,P,y,CFG['quality'][q]);q0=(physical-m0.low)/m0.span;f0,j0=m0.evaluate(q0)
  for p in ['P0','P1']:
   m=model.Model(t,P,y,quality(q,p));u=(physical-m.low)/m.span;f,j=m.evaluate(u);ind,ji=residual_and_jac(physical,t,P,y,quality(q,p))
   ck(q+p+'_physical_objective_invariant',np.allclose(f,f0,rtol=1e-12,atol=1e-12))
   ck(q+p+'_Jacobian_chain_rule',np.allclose(j/m.span,j0/m0.span,rtol=1e-11,atol=1e-11))
   ck(q+p+'_independent_residual',np.allclose(f,ind,rtol=1e-11,atol=1e-11));ck(q+p+'_independent_derivative',np.allclose(j,ji*m.span,rtol=1e-10,atol=1e-10))
   before=dict(CFG['quality'][q]);qual=quality(q,p);ck(q+p+'_only_bias_support_changed',all(qual[k]==before[k] for k in ['bearing_sigma_deg','nav_sigma_m']) and CFG['quality'][q]==before)
 a,fb,n=model.initialize(t,P,y);b,fb2,n2=initial(t,P,y);ck('public_initializer_independent',np.allclose(a,b[:4],rtol=1e-10,atol=1e-10) and fb==fb2 and n==n2)
 # Static access test; no file opened, and no panel or random number generated.
 import estimate_public as worker
 try:worker.gate('open',(str(DATA/'private/hidden.npz'),'r',0));deny=False
 except PermissionError:deny=True
 ck('private_path_denied',deny);worker.denied.clear()
 code=(O/'estimate_public.py').read_text(encoding='utf-8-sig');ck('single_original_estimate_call',code.count('model.estimate(')==1 and 'least_squares(' not in code)
 ck('no_private_import_or_load',not any(isinstance(n,(ast.Import,ast.ImportFrom)) and any(s in ast.unparse(n) for s in ['generate','independent','statistics','panel']) for n in ast.walk(ast.parse(code))))
 ck('no_science_yet',not DATA.exists());ck('nearest_rank_INF_fixture',math.isinf(quantile([1]*94+[math.inf]*6,.95)) and quantile([1,2,3,4],.5)==2)
 ck('storage_available',psutil.disk_usage('D:/').free>CFG['budget']['D_data_bytes']+64*1024**2)
 files=list(O.glob('*.py'))+[REF/x for x in ['t3_w1_model.py','t3_w1_generate.py','t3_w1_information.py']]
 for p in files:ast.parse(p.read_text(encoding='utf-8-sig'))
 result=dict(PASS=all(c['PASS'] for c in checks),checks=checks,new_panel=0,new_innovations=0,new_observations=0,optimizer_calls=0,Fisher_calls=0,scope='SYNTHETIC_FORMULA_AND_PERMISSION_CONTROLS_ONLY')
 dump(O/'PRODUCTION_PREFLIGHT.json',result,True)
 if not result['PASS']:raise RuntimeError('PREFLIGHT_GATE_FAIL_NO_FORMAL_EXECUTION')
 dump(O/'PRODUCTION_CODE_FREEZE.json',dict(files=[dict(path=str(p),sha256=sha(p)) for p in files],binding_sha256=sha(O/'EFFECTIVE_INPUT_BINDING.json'),before_panel=True,zero_scientific_retries=True),True)
 print('PREFLIGHT_PASS',len(checks),'production frozen; panel still absent',flush=True)
if __name__=='__main__':main()
