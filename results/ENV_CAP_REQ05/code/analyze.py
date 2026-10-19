"""REQ05 deterministic algebra controls; no propagation or performance simulation."""
import json, pathlib, time, numpy as np
O=pathlib.Path(__file__).resolve().parents[1]
def evaluate(spec):
 out=[]
 for c in spec['controls']:
  J=np.asarray(c['J'],float); N=np.asarray(c['N'],float).reshape(len(J),-1)
  U,s,_=np.linalg.svd(N,full_matrices=False)
  q=U[:,s>spec['svd_absolute_tolerance']]
  A=J-q@(q.T@J); F=A.T@A
  _,sv,V=np.linalg.svd(A,full_matrices=True)
  rank=int(np.sum(sv>spec['svd_absolute_tolerance']))
  out.append(dict(id=c['id'],rank=rank,nuisance_rank=q.shape[1],information=F.tolist(),
     singular_values=sv.tolist(),nullspace=V[rank:].T.tolist(),scope='DIMENSIONLESS_ALGEBRA_ONLY',
     expected_rank=c['expected_rank'],PASS=rank==c['expected_rank']))
 # D is a lossy fixed contrast, with its full shared-reference covariance.
 J=np.array([[1,0],[0,1],[2,1.]],float)
 D=np.array([[-1,1,0],[-1,0,1.]],float); C=D@D.T
 raw=J.T@J; derived=J.T@D.T@np.linalg.solve(C,D@J)
 # Removing an unknown common scalar in raw space must equal full contrast information.
 Pc=np.eye(3)-np.ones((3,3))/3
 special={'shared_reference_covariance':C.tolist(),'raw_information':raw.tolist(),
  'contrast_information':derived.tolist(),'profile_information':(J.T@Pc@J).tolist(),
  'data_processing_min_eigenvalue':float(np.linalg.eigvalsh(raw-derived).min()),
  'source_profile_agreement':float(np.max(np.abs(derived-J.T@Pc@J)))}
 # Same-frequency CSD: common source phase cancels, channel phase does not.
 h=np.array([1+2j,2-1j,1-.5j]);g=np.array([1+.2j,.8-.3j,1.1+.1j]);s=2*np.exp(.7j)
 y=s*g*h; Cyy=np.outer(y,y.conj()); expected=abs(s)**2*np.outer(g*h,(g*h).conj())
 special['csd_source_phase_cancellation_error']=float(np.max(abs(Cyy-expected)))
 # Local plane wave y_i=S exp(i k a_i sin(theta)): global source-scale removed.
 # This is a control only, not a CZ prediction.
 a=np.array([-3,-1,1,3.]);theta=.3;k=.7;eps=1e-5
 dh=1j*k*a*np.cos(theta)*np.exp(1j*k*a*np.sin(theta))
 fd=(np.exp(1j*k*a*np.sin(theta+eps))-np.exp(1j*k*a*np.sin(theta-eps)))/(2*eps)
 special['plane_wave_analytic_fd_max']=float(np.max(abs(dh-fd)))
 # Polar constant-velocity trajectory bearing derivative, including exactly zero z column.
 x=np.array([10.,.2,3.,.7,.4]);t=np.array([0.,1.,2.,3.]);pos=np.column_stack((t,.1*t*t))
 def obs(v):
  p=v[0]*np.array([np.cos(v[1]),np.sin(v[1])])+t[:,None]*v[3]*np.array([np.cos(v[4]),np.sin(v[4])])-pos
  return np.arctan2(p[:,1],p[:,0])
 p=x[0]*np.array([np.cos(x[1]),np.sin(x[1])])+t[:,None]*x[3]*np.array([np.cos(x[4]),np.sin(x[4])])-pos
 B=[]
 for ti,pi in zip(t,p):
  dp=np.array([[np.cos(x[1]),-x[0]*np.sin(x[1]),0,ti*np.cos(x[4]),-ti*x[3]*np.sin(x[4])],
   [np.sin(x[1]),x[0]*np.cos(x[1]),0,ti*np.sin(x[4]),ti*x[3]*np.cos(x[4])]])
  B.append(np.array([-pi[1],pi[0]])@dp/(pi@pi))
 B=np.array(B);fd=np.column_stack([(obs(x+eps*np.eye(5)[i])-obs(x-eps*np.eye(5)[i]))/(2*eps) for i in range(5)])
 special.update(bearing_analytic=B.tolist(),bearing_fd=fd.tolist(),bearing_fd_max=float(abs(B-fd).max()))
 return {'controls':out,'special':special}
def main():
 spec=json.loads((O/'ALGEBRA_SPEC.json').read_text());t=time.monotonic();result=evaluate(spec)
 result['elapsed_s']=time.monotonic()-t
 (O/'ALGEBRA_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
 if not all(r['PASS'] for r in result['controls']):raise SystemExit('registered algebra rank failed')
if __name__=='__main__':main()
