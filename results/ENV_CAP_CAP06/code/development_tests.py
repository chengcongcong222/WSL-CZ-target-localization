import sys,time,json
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from generate_p2 import make_cal,waveform
from frontend import extract,calibration
from support_p3 import tests,search
from common import platform,axis
checks=[]
def ck(n,b,e=0):checks.append(dict(check=n,PASS=bool(b),error=float(e)))
rng=np.random.default_rng(101);g=10**(rng.uniform(-1,1,8)/20)*np.exp(1j*rng.uniform(-.3,.3,8));beta=.0007
cal,ang=make_cal(rng,g,beta);b,gg,res=calibration(cal,ang)
ck('known calibration orientation',abs(b-beta)<1e-4,abs(b-beta));ck('known calibration response',res<.002,res)
t=np.array([0.,1200.]);pos=platform(t)
for a in [-50.,-12.,0.,12.,50.]:
 target=pos+50000*np.column_stack((np.cos(np.deg2rad(a)+axis(t)),np.sin(np.deg2rad(a)+axis(t))))
 y,bg,truth=waveform(np.random.default_rng(102),t,target,pos,beta,g,20,'POSITIVE_CONTROL')
 z=extract(y[0],bg[0],cal,ang)
 ck('signal detected '+str(a),len(z['products'])>0)
 if z['products']:
  pp=min(z['products'],key=lambda q:abs(q['f_rx_hz']-200))
  ck('u compatibility '+str(a),pp['u_lo']<=np.cos(np.deg2rad(a)-beta)<=pp['u_hi'])
  ck('mirror preserved '+str(a),len(pp['mirror_angles_deg'])==2)
noise=np.zeros((8,128),complex);z=extract(noise,np.ones_like(noise),cal,ang);ck('zero signal rejected',not z['products'])
x=np.arange(-7,8,2)
ck('far field mirror',np.max(abs(np.exp(1j*x*np.cos(.2))-np.exp(1j*x*np.cos(-.2))))<1e-14)
# Free gain and arbitrary per-epoch source null information.
A=rng.normal(size=(20,5));ck('free nuisance null',np.linalg.norm(A-np.eye(20)@A)<1e-14)
D=np.eye(5)[1:]-np.eye(5)[0];ck('common source removed',np.max(abs(D@np.ones(5)))<1e-14)
ck('shared reference covariance',np.max(abs(D@D.T-(np.eye(4)+np.ones((4,4)))))<1e-14)
# all possible angular sets retain every point of physical Cartesian rectangle
obs=[dict(t=0.,pos=np.zeros(2),heading=0.,heading_lo=0.,heading_hi=0.,beta=0.,sigma=25.,sets=[(-1.,1.)])]
k,d,w,b,s,n=search(obs,31);vol=np.sum(np.prod(k[:,1]-k[:,0],axis=1));initial=np.prod([60000-45000*np.cos(np.deg2rad(5)),120000*np.sin(np.deg2rad(5)),3-np.cos(np.deg2rad(15)),6*np.sin(np.deg2rad(15))])
ck('unresolved domain preserved',abs(vol/initial-1)<1e-12,abs(vol/initial-1))
# Uniform level/no relative motion derivatives are zero.
ck('constant level zero derivative',np.linalg.norm(np.diff(np.ones((9,5)),axis=0))==0)
ck('zero relative motion frequency change',np.linalg.norm(np.diff(np.zeros(6)))==0)
for f in Path(__file__).parent.glob('*.py'):compile(f.read_text(),str(f),'exec')
out=dict(PASS=all(c['PASS'] for c in checks),checks=checks,scope='PRE_FREEZE_SYNTHETIC_DEVELOPMENT_ONLY',repairs=['attempt1 positive control incorrectly expected detection during deliberate dropout; use all-present analytic control, no science seed changed','fixed timestamp orientation straddling 600s','preserve unresolved Cartesian boxes','independent cold raw sums double precision'])
(Path(__file__).resolve().parents[1]/'P0_DEVELOPMENT_TESTS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
