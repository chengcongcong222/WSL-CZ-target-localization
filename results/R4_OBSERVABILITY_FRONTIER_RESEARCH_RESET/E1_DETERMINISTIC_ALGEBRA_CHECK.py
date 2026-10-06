"""Only fixed dimensionless algebra checks. No project scenario input."""
from pathlib import Path
import json
import numpy as np
checks=[]
def check(name, actual, expected, passed):
    checks.append(dict(id=name,actual=actual,expected=expected,passed=bool(passed)))
def proj(N): return np.eye(N.shape[0])-N@np.linalg.pinv(N)
def rank(M): return int(np.linalg.matrix_rank(M,tol=1e-10))
F=np.array([2.,3.]); D=np.array([.9,.8]); a=1.1
err=float(np.max(np.abs(np.outer(D,F)-np.outer(D/a,F*a))))
check('A',err,'<=1e-12',err<=1e-12)
t=np.array([-1.,0.,1.])
check('B',rank(np.column_stack([t,-t])),1,rank(np.column_stack([t,-t]))==1)
Z=np.column_stack([np.ones(3),t])
N=np.zeros((6,4)); N[:3,:2]=Z; N[3:,2:]=Z
J=np.column_stack([np.r_[2+3*t,np.zeros(3)],np.r_[np.zeros(3),4-2*t]])
n=float(np.linalg.norm(proj(N)@J))
check('C',n,'<=1e-12',n<=1e-12)
n=float(np.linalg.norm(proj(N)@np.r_[1.,-2.,1.,0.,0.,0.]))
check('C2',n,'>1e-6',n>1e-6)
L=np.array([1.,2.,3.])
check('D',rank(np.column_stack([L,L])),1,rank(np.column_stack([L,L]))==1)
check('D2',rank(np.column_stack([L,np.ones(3)])),2,rank(np.column_stack([L,np.ones(3)]))==2)
j=np.array([-1.,1.]); common=np.ones((2,1))
n=float(np.linalg.norm(proj(common)@j))
n2=float(np.linalg.norm(proj(np.eye(2))@j))
check('E',dict(common_source_residual=n,free_node_residual=n2),'common>1e-6; free<=1e-12',n>1e-6 and n2<=1e-12)
acoustic=1.-1./1.; prior=1.-1./(1.+1.)
check('F',dict(acoustic=acoustic,with_external_prior=prior),dict(acoustic=0.,with_external_prior=.5),acoustic==0 and prior==.5)
n=7; k=-3; s=k*(n+1)/(n+k)
check('SUN-W',s,'printed sum -6, not 1',s==-6 and s!=1)
check('SUN-F',10*2,'printed 20 versus CV 2',10*2!=2)
p=np.array([3.,4.]); v=np.array([-1.,2.]); r=np.linalg.norm(p); e=p/r
g=-200/1500/r*(np.eye(2)-np.outer(e,e))@v
check('SUN-H',g.tolist(),'position derivative nonzero',np.linalg.norm(g)>1e-6)
result=dict(scope='FIXED_DETERMINISTIC_ALGEBRA_ONLY',scientific_runs=0,scenario_rank_runs=0,checks=checks,pass_count=sum(c['passed'] for c in checks),fail_count=sum(not c['passed'] for c in checks))
Path(__file__).with_name('E1_ALGEBRA_CHECK_RESULTS.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
assert result['fail_count']==0
print(json.dumps(result,ensure_ascii=False,indent=2))
