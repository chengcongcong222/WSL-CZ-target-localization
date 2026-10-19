import sys, json, hashlib, pathlib, csv, subprocess
R=pathlib.Path(r'D:\XiaomiMiMoProjects\WSL系统汇聚区目标定位')
O=R/'results/ENV_CAP_REQ05'
def put(name,body):
 p=O/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(body.strip()+'\n',encoding='utf-8')
put('.gitattributes','* -text')
put('code/analyze.py',r'''
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
''')
put('code/cold_audit.py',r'''
"""Independent exact-rational elimination; never imports analyze.py."""
from fractions import Fraction as Q
import json, pathlib, hashlib, csv, math, time
O=pathlib.Path(__file__).resolve().parents[1];R=O.parents[1]
def tr(A):return [list(x) for x in zip(*A)]
def mm(A,B):return [[sum(x*y for x,y in zip(a,b)) for b in tr(B)] for a in A]
def sub(A,B):return [[x-y for x,y in zip(a,b)] for a,b in zip(A,B)]
def rref(A):
 A=[row[:] for row in A];p=[];i=0
 for j in range(len(A[0])):
  k=next((k for k in range(i,len(A)) if A[k][j]),None)
  if k is None:continue
  A[i],A[k]=A[k],A[i];v=A[i][j];A[i]=[x/v for x in A[i]]
  for k in range(len(A)):
   if k!=i:
    v=A[k][j];A[k]=[x-v*y for x,y in zip(A[k],A[i])]
  p.append(j);i+=1
  if i==len(A):break
 return A,p
def inv(A):
 n=len(A);B=[a+[Q(i==j) for j in range(n)] for i,a in enumerate(A)]
 C,p=rref(B)
 assert p[:n]==list(range(n))
 return [x[n:] for x in C]
def exact(c):
 J=[[Q(x) for x in a] for a in c['J']];N=[[Q(x) for x in a] for a in c['N']]
 piv=rref(N)[1] if len(N[0]) else []
 if piv:
  B=[[a[i] for i in piv] for a in N];A=sub(J,mm(mm(B,inv(mm(tr(B),B))),mm(tr(B),J)))
 else:A=J
 return A,mm(tr(A),A),len(rref(A)[1]),len(piv)
def main():
 start=time.monotonic();spec=json.loads((O/'ALGEBRA_SPEC.json').read_text());res=json.loads((O/'ALGEBRA_RESULTS.json').read_text());checks=[]
 def ck(name,ok,detail=''):checks.append(dict(name=name,PASS=bool(ok),detail=str(detail)))
 byid={x['id']:x for x in res['controls']}
 for c in spec['controls']:
  a,F,rank,nrank=exact(c);z=byid[c['id']]
  ck(c['id']+':exact_rank',rank==z['rank']==c['expected_rank'],rank)
  ck(c['id']+':nuisance_rank',nrank==z['nuisance_rank'])
  err=max(abs(float(F[i][j])-z['information'][i][j]) for i in range(5) for j in range(5))
  ck(c['id']+':rational_information',err<spec['matrix_tolerance'],err)
  Z=z['nullspace']; err=max([abs(sum(float(a[i][j])*Z[j][k] for j in range(5))) for i in range(len(a)) for k in range(5-rank)]+[0])
  ck(c['id']+':null_residual',err<spec['matrix_tolerance'],err)
  # Independent exact null vectors for reproducible structural audit.
  rr,piv=rref(a);free=[j for j in range(5) if j not in piv];vectors=[]
  for f in free:
   v=[Q(0)]*5;v[f]=Q(1)
   for i,p in enumerate(piv):v[p]=-rr[i][f]
   vectors.append([str(e) for e in v])
  z['exact_null_vectors']=vectors
 # Exact shared reference covariance and Gaussian profile / data processing.
 J=[[Q(x) for x in a] for a in [[1,0],[0,1],[2,1]]];D=[[Q(x) for x in a] for a in [[-1,1,0],[-1,0,1]]]
 C=mm(D,tr(D));F=mm(mm(mm(tr(J),tr(D)),inv(C)),mm(D,J));s=res['special']
 ck('shared_covariance_exact',s['shared_reference_covariance']==[[2.,1.],[1.,2.]])
 ck('contrast_exact',max(abs(float(F[i][j])-s['contrast_information'][i][j]) for i in range(2) for j in range(2))<spec['matrix_tolerance'])
 ck('source_profile_equal',s['source_profile_agreement']<spec['matrix_tolerance'])
 raw=mm(tr(J),J);loss=sub(raw,F)
 ck('data_processing_exact_PSD',loss[0][0]>=0 and loss[1][1]>=0 and loss[0][0]*loss[1][1]-loss[0][1]**2>=0)
 ck('csd_phase_cancel',s['csd_source_phase_cancellation_error']<spec['matrix_tolerance'])
 ck('plane_wave_fd',s['plane_wave_analytic_fd_max']<spec['fd_tolerance'])
 # Independent finite difference of scalar atan2 trajectory, different code path and step.
 x=[10.,.2,3.,.7,.4];delta=2e-5;err=0.
 for row,t in enumerate([0.,1.,2.,3.]):
  def b(u):return math.atan2(u[0]*math.sin(u[1])+t*u[3]*math.sin(u[4])-.1*t*t,u[0]*math.cos(u[1])+t*u[3]*math.cos(u[4])-t)
  for j in range(5):
   xp=x[:];xm=x[:];xp[j]+=delta;xm[j]-=delta
   err=max(err,abs((b(xp)-b(xm))/(2*delta)-s['bearing_analytic'][row][j]))
 ck('bearing_independent_fd',err<spec['fd_tolerance'],err)
 ck('bearing_depth_zero',all(row[2]==0 for row in s['bearing_analytic']))
 # Historical inputs verified against the frozen evidence list, not historical experiment reruns.
 for row in csv.DictReader((O/'EVIDENCE_INDEX.csv').open(encoding='utf-8',newline='')):
  p=R/row['path'];ck('input:'+row['id'],hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'])
 frozen=json.loads((O/'DESIGN_FREEZE.json').read_text())
 for row in frozen['code_and_spec']:
  ck('freeze:'+row['path'],hashlib.sha256((O/row['path']).read_bytes()).hexdigest()==row['sha256'])
 # Structural nesting independent exact F(C0)-F(C1) for level: quadratic form is rank1 PSD.
 f0=exact(next(c for c in spec['controls'] if c['id']=='LEVEL_CONSTANT'))[1]
 f1=exact(next(c for c in spec['controls'] if c['id']=='LEVEL_LINEAR'))[1]
 loss=sub(f0,f1);rank=len(rref(loss)[1]);ck('level_drift_nested_PSD',rank==1 and sum(loss[i][i] for i in range(5))>0)
 out={'PASS':all(c['PASS'] for c in checks),'checks':len(checks),'failed':sum(not c['PASS'] for c in checks),'records':checks,'elapsed_s':time.monotonic()-start,'scope':'EXACT_ALGEBRA_AND_SOURCE_IDENTITY_NOT_PHYSICAL_PERFORMANCE','independent_implementation':'Fraction Gaussian elimination; no import of production evaluator','independent_researcher':False}
 (O/'EXACT_NULLSPACES.json').write_text(json.dumps({z['id']:z['exact_null_vectors'] for z in byid.values()},indent=2)+'\n')
 (O/'INDEPENDENT_COLD_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n')
 if not out['PASS']:raise SystemExit('COLD_AUDIT_FAILED')
if __name__=='__main__':main()
''')
cs=[]
def add(i,J,N,rank,meaning):cs.append(dict(id=i,J=J,N=N,expected_rank=rank,interpretation=meaning))
def empty(n):return [[] for _ in range(n)]
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
t=list(range(-2,3));J=[[i,0,i*i,1,i**3] for i in t]
add('LEVEL_CONSTANT',J,[[1] for i in t],3,'constant nuisance removes constant state direction; synthetic polynomial only')
add('LEVEL_LINEAR',J,[[1,i] for i in t],2,'range-column linear trend absorbed; no CZ numerical claim')
add('LEVEL_FREE',J,eye(5),0,'arbitrary per-epoch source plus gain absorbs every state column')
add('LEVEL_FLAT',[[0,0,1,0,0]]*3,[[1]]*3,0,'no range-dependent propagation')
J=[[1,0,1,0,0],[2,0,3,0,0],[3,0,2,0,0]]
add('RELATIVE_KNOWN_GAIN',J,empty(3),2,'two independent range/depth signatures possible, hypothetical')
add('RELATIVE_CONSTANT_GAIN',J,[[1]]*3,2,'fixed differential gain over time need not absorb all changing shape')
add('RELATIVE_FREE_GAIN',J,eye(3),0,'arbitrary per-epoch differential gain destroys contrast information')
add('RELATIVE_STATIC_GAIN',[[1,0,2,0,0]]*3,[[1]]*3,0,'static unknown gain absorbs static ratio')
add('RELATIVE_COUPLED',[[1,0,2,0,0],[2,0,4,0,0]],empty(2),1,'range-depth tradeoff remains despite nonzero sensitivity')
J=[[0,i,0,0,0] for i in [-3,-1,1,3]]
add('PHASE_PLANE_WAVE',J,[[1]]*4,1,'horizontal plane wave has bearing only after source-scale removal')
add('PHASE_FREE_CHANNEL',J,eye(4),0,'arbitrary channel phases absorb spatial phase')
J=[[i,0,i*i,0,0] for i in [1,2,3]]
add('DELAY_COMMON_TIME',J,[[1]]*3,2,'relative delay cancels emission time; conditional independent paths')
add('DELAY_FREE_PATH',J,eye(3),0,'unknown independent path delays absorb states')
J=[[1,0,0,i,0] for i in t]
add('FREQUENCY_UNKNOWN_OFFSET',J,[[1]]*5,1,'constant shift lost; changing signature may survive')
add('FREQUENCY_UNKNOWN_DRIFT',J,[[1,i] for i in t],0,'affine motion signature lost to affine reference/source drift')
add('APERTURE_FREE_SOURCE_PHASE',[[1,0,i*i,i,0] for i in t],eye(5),0,'unconstrained temporal source phase absorbs coherent phase trajectory')
J=[[1,0,0,0,0],[0,1,0,0,0],[0,0,0,1,0],[0,0,0,0,1]]
add('AUX_HORIZONTAL_UPPER_BOUND',J,empty(4),4,'abstract independent horizontal signatures, not actual bearing FIM')
add('AUX_FREE_EPOCH_BIAS',J,eye(4),0,'per-observation free angle bias absorbs every horizontal signature')
put('ALGEBRA_SPEC.json',json.dumps(dict(state=['r','theta','z','v','psi'],units='DIMENSIONLESS_SYNTHETIC_JACOBIANS_ONLY',svd_absolute_tolerance=1e-10,matrix_tolerance=1e-9,fd_tolerance=1e-7,controls=cs),indent=2))
put('code/run_once.py',r'''
"""Bounded single-execution launcher, writes only to D:, no scientific retry."""
import sys, os, subprocess, pathlib, time, json, psutil, hashlib
O=pathlib.Path(__file__).resolve().parents[1];W=pathlib.Path(r'D:\ProjectStorage\WSL-CZ\ENV_CAP_REQ05')
phase=sys.argv[1];assert phase in ['analysis','audit']
freeze=json.loads((O/'DESIGN_FREEZE.json').read_text())
for row in freeze['code_and_spec']:
 assert hashlib.sha256((O/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
marker=W/(phase+'.started');marker.open('x').write(str(time.time()))
t=time.monotonic();peak=0;previous=0
if phase=='audit':
 prior=json.loads((O/'ANALYSIS_RUN.json').read_text());assert prior['exit_code']==0
 previous=prior['wall_s']
child=subprocess.Popen([sys.executable,str(O/'code'/('analyze.py' if phase=='analysis' else 'cold_audit.py'))],env=os.environ.copy())
reason=None
while child.poll() is None:
 try:
  p=psutil.Process(child.pid);rss=sum(q.memory_info().rss for q in [p]+p.children(recursive=True) if q.is_running());peak=max(peak,rss)
 except psutil.Error:pass
 if time.monotonic()-t+previous>600 or peak>4*1024**3:
  reason='BUDGET_EXCEEDED'
  for q in psutil.Process(child.pid).children(recursive=True):q.kill()
  child.kill();break
 time.sleep(.05)
rc=child.wait();rec=dict(phase=phase,wall_s=time.monotonic()-t,peak_sampled_tree_RSS=peak,exit_code=rc,stop_reason=reason,acoustic_calls=0,MC=0)
(O/(phase.upper()+'_RUN.json')).write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec));sys.exit(rc or bool(reason))
''')
print('Created',len(cs),'registered algebra controls')
