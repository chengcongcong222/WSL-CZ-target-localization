"""Bounded implementation controls only; not the12state performance panel."""
import io,ast,time
import numpy as np
from aux_d1_common import *
from aux_d1_model import *
from aux_d1_generate import make_record

def run():
 if (OUT/'PREFLIGHT_START.json').exists():raise RuntimeError('Controls already started; cumulative limit must be preserved')
 dump('PREFLIGHT_START.json',{'utc':utc(),'seed':2026101212,'not_test_panel':True});tic=time.perf_counter();checks=[];fits=0;gradients=0
 def ck(n,ok,detail=''):checks.append({'name':n,'PASS':bool(ok),'detail':str(detail)})
 fixed=np.load(OUT/'SHARED_FIXED_126.npy');h=np.array([51.237,.731,2.173,7.31]);t=np.arange(121)*10.;pM=np.column_stack([2*np.minimum(t,600)+2*np.maximum(t-600,0)*np.cos(np.pi/12),2*np.maximum(t-600,0)*np.sin(np.pi/12)]);P=np.stack([pM,pM+np.array([-123.,4980.])],axis=1);target=1000*h[0]*np.array([np.cos(np.deg2rad(h[1])),np.sin(np.deg2rad(h[1]))])+t[:,None]*h[2]*np.array([np.cos(np.deg2rad(h[3])),np.sin(np.deg2rad(h[3]))]);delta=target[:,None,:]-P;y=np.arctan2(delta[:,:,1],delta[:,:,0]);rng=np.random.default_rng(2026101212)
 for node_count in [1,2]:
  pp=P[:,:node_count].copy();yy=y[:,:node_count].copy();u=(h-LOW)/SPAN;model=Model(t,pp,yy);model.geometry(u);ck('independent_cartesian_'+str(node_count),np.max(abs(model.f))<1e-9)
  for trial in range(8):
   q=rng.uniform(.15,.85,4);m=Model(t,pp,yy);J=m.jac(q)
   for j in range(4):
    gradients+=1;step=np.zeros(4);step[j]=1e-6;fd=(m.residual(q+step)-m.residual(q-step))/(2e-6);diff=np.linalg.norm(fd-J[:,j])/max(np.linalg.norm(J[:,j]),1e-12);ck('gradient_'+str(node_count)+str(trial)+str(j),diff<1e-6,diff)
  ss,diag=starts(t,pp,yy,fixed)
  for q in [ss[0],(np.array([52.4,.2,2.2,6.])-LOW)/SPAN]:
   fits+=1;result=fit_one(t,pp,yy,q);ck('controlfit_'+str(fits),result[6]==1 and result[4]<1e-10,result[4]);ck('scorecost_'+str(fits),abs(result[4]-2*result[5]/(121*node_count))<1e-10)
  moved=pp.copy();moved[15,0,1]+=20.;m=Model(t,moved,yy);ck('reportednavigation_matters'+str(node_count),abs(m.residual(u)[15*node_count])>.01)
  shift=np.array([237.,-149.]);m1=Model(t,pp,yy);m2=Model(t,pp+shift,yy);ck('reference_translation'+str(node_count),np.max(abs(m1.residual(u)-m2.residual(u)))<1e-9)
 s0,_=starts(t,P[:,:1],y[:,:1],fixed);s5,_=starts(t,P,y,fixed);ck('127_sharedbytes',s0[[0,*range(2,128)]].tobytes()==s5[[0,*range(2,128)]].tobytes());ck('U0center',np.array_equal(s0[1],np.full(4,.5)))
 badP=np.stack([pM,pM],axis=1);badY=np.concatenate([y[:,:1],y[:,:1]],axis=1);_,diag=starts(t,badP,badY,fixed);ck('invalidintersection_fixedfallback',diag['resource']['fallback'])
 for side in [-1,1]:
  pub,private=make_record(h,2026101212,31,0,side);offset=private['true_node_positions'][:,1]-private['true_node_positions'][:,0];ck('fixed_globaloffset'+str(side),np.max(abs(offset-offset[0]))<1e-10)
  ck('biasfixed'+str(side),private['systematics'].shape==(3,));ck('U0U5publicbytes'+str(side),np.ascontiguousarray(pub['bearing_main']).tobytes()==np.concatenate([pub['bearing_main'],pub['bearing_aux']],axis=1)[:,:1].tobytes())
 e=np.linspace(-100,100,1001);stable=score(e);ordinary=float(np.mean(2*A*A*(np.sqrt(1+(e/A)**2)-1)));ck('stable_score_formula',math.isclose(stable,ordinary,rel_tol=1e-13,abs_tol=1e-10))
 fake=np.zeros((128,len(FIELDS)));fake[:,:4]=CENTER;fake[:,6]=1;fake[:,4]=1.;fake[5,4]=np.nextafter(1.,0);chosen,ids,reason=choose(fake,1.);ck('strictfloatselection',chosen==5)
 txt=io.StringIO();w=csv.writer(txt);w.writerow([format_float(v) for v in fake[:,4]]);back=np.array([float(v) for v in next(csv.reader(io.StringIO(txt.getvalue())))]);ck('17digitCSVroundtrip',np.array_equal(back,fake[:,4]));ck('INFnotzero',math.isinf(sorted([1.,math.inf])[1]))
 fake[:,:4]=CENTER;env=envelope(fake,np.array([0,1]));ck('finiteclustermeaning',env['cluster_count']==1 and env['scope']=='FINITE_TERMINAL_ENVELOPE_ONLY')
 tree=ast.parse((ROOT/'aux_d1_model.py').read_text(encoding='utf-8-sig'));ck('estimator_no_private_dependencies',not any(isinstance(n,(ast.Import,ast.ImportFrom)) and any(x in ast.unparse(n) for x in ['generate','evaluate','seed','truth']) for n in ast.walk(tree)))
 ck('RSSmonitor',rss()<8*1024**3);ck('controlbudgets',fits<=12 and gradients<=128)
 duration=time.perf_counter()-tic;result={'PASS':all(c['PASS'] for c in checks),'fits':fits,'gradient_checks':gradients,'checks':checks,'elapsed_s':duration,'not_performance_panel':True,'no_test_observations_generated':True};dump('PREFLIGHT.json',result);record_clock('preflight',duration,control_fits=fits,gradient_checks=gradients);print('PREFLIGHT',result['PASS'],'fits',fits,'gradients',gradients,'FAIL',[c for c in checks if not c['PASS']],flush=True)
 if not result['PASS']:raise SystemExit(1)
if __name__=='__main__':run()