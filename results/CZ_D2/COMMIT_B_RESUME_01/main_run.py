"""One registered main invocation. No scientific retries or alternative fitting."""
from common import *
import subprocess,sys

def information(states):
 g=module('t3_w1_generate');f=module('t3_w1_information');mat=[];eff=[];rows=[];resource=[]
 for qi,Q in enumerate(['Q1','Q2']):
  for bi,B in enumerate(['B0','B3']):
   for ci,h in enumerate(states):
    tr=dict(zip(['r0_km','theta0_deg','v_mps','psi_deg'],h));x=np.r_[g.cartesian(tr),0.,0.]
    for side in [-1,1]:
     t,P=g.node_geometry(CFG['geometry'][B],side);mu=(x[:2]+t[:,None]*x[2:4])[:,None,:]-P;rr=np.linalg.norm(mu,axis=-1)/1000
     resource.append(dict(quality=Q,geometry=B,case_id=f'D{ci+1:02d}',side=side,node_baseline_km=float(np.linalg.norm(P[0,1]-P[0,0])/1000),main_min_km=rr[:,0].min(),main_max_km=rr[:,0].max(),aux_min_km=rr[:,1].min(),aux_max_km=rr[:,1].max(),main_outside45_60=int(np.sum((rr[:,0]<45)|(rr[:,0]>60))),aux_outside45_60=int(np.sum((rr[:,1]<45)|(rr[:,1]>60))),not_detection_coverage=True))
     dm,dv=f.finite_difference_check(x,t,P,CFG['quality'][Q]);assert dm<=1e-7 and dv<=1e-6,('INFORMATION_FD_GATE',Q,B,ci,side,dm,dv)
     for known in [False,True]:
      row,fp,fs,e,cov=f.analyze(x,t,P,CFG['quality'][Q],known);z=np.zeros((6,6));z[:len(fs),:len(fs)]=fs;mat.append(z);eff.append(e)
      row.update(dict(quality=Q,geometry=B,case_id=f'D{ci+1:02d}',side=side,bias_known=known,matrix_index=len(rows),mean_fd_error=dm,variance_fd_error=dv,P0_P1_data_information_identical=True,rank_near_threshold='UNRESOLVED' if .5e-10<=1/row['condition_design']<=2e-10 else 'NOT_NEAR_FULL_RANK_THRESHOLD'))
      rows.append(row)
 np.savez_compressed(O/'INFORMATION_MATRICES.npz',Fscaled=np.array(mat),velocity_effective=np.array(eff));csvwrite(O/'INFORMATION.csv',rows);csvwrite(O/'RESOURCE_GEOMETRY.csv',resource)

def main():
 verify_code();assert not DATA.exists();DATA.mkdir();(DATA/'private').mkdir();(DATA/'public').mkdir()
 from generate import panel,generate
 states=panel();generate(states)
 # No access to estimator internals or truth transfer; worker receives no state arguments.
 child=subprocess.run([sys.executable,'-B',str(O/'estimate_public.py')],cwd=ROOT)
 if not (DATA/'terminals.npy').exists():raise RuntimeError('ESTIMATION_DID_NOT_START')
 from statistics_main import run
 run()
 if child.returncode:raise RuntimeError('WORKER_FAILED_NO_RETRY')
 information(states)
 dump(O/'MAIN_OUTPUT_SEAL.json',dict(files=[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(DATA.rglob('*')) if p.is_file()]+[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(O.glob('*')) if p.is_file() and p.suffix in ['.csv','.npz']]),True)
 print('MAIN_COMPLETED_ONCE',flush=True)
if __name__=='__main__':main()
