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
