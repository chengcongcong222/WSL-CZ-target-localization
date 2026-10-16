"""Synthetic-only pre-freeze algebra controls. No archived field evaluation."""
import numpy as np
from model import *
checks=[]
def ck(n,b):checks.append({'check':n,'PASS':bool(b)})
for rho in [0,.9]:
 w=weight(rho);rng=np.random.default_rng(2026101604);j=rng.normal(size=(121,4));v=rng.normal(size=121)
 for c in ['C0','C1','C2']:
  p=project(j.T,w,c).T;f=p.T@p;raw=(w@j).T@(w@j)
  ck(f'{rho}/{c}/DPI',np.linalg.eigvalsh(raw-f).min()>-1e-9)
  ck(f'{rho}/{c}/constant',np.linalg.norm(project(np.ones(121)*7,w,c))<1e-10)
  ck(f'{rho}/{c}/flat-range',np.linalg.norm(project(np.zeros(121),w,c))==0)
  if c in ['C1','C2']:ck(f'{rho}/{c}/linear',np.linalg.norm(project(V@np.array([3.,2.]),w,c))<1e-10)
  if c=='C2':ck(f'{rho}/{c}/free',np.linalg.norm(project(v,w,c))==0)
 ck(f'{rho}/nested',np.linalg.eigvalsh(project(j.T,w,'C0')@project(j.T,w,'C0').T-project(j.T,w,'C1')@project(j.T,w,'C1').T).min()>-1e-9)
for h in [[50,0,2,5],[47,-1,1.4,-7],[59,-.5,1.8,-11]]:
 b,r=geom(h);ck(f'{h}/finite',np.isfinite(b).all() and np.isfinite(r).all());ck(f'{h}/not-use-side',np.min(b)>=0)
# Piecewise-linear feature and finite-difference scaling, independent explicit geometry.
h=np.array([52.,1.,2.2,4.]);b,r=geom(h);rr=np.linspace(40000,66000,209);yy=rr*.0003
jb,jl=jac(h,rr,yy);jb2,jl2=jac(h,rr,yy,2)
ck('synthetic-jac-step',np.max(abs(jl-jl2))<1e-7)
ck('full-grid',grid().shape==(114576,4) and np.min(grid()[:,1])==-5 and np.max(grid()[:,1])==5)
infozero=info(np.zeros((121,4)));ck('rankzero-not-zero-variance',infozero['rank']==0 and infozero['condition'] is None)
dump(O/'SYNTHETIC_CONTROLS.json',{'scope':'synthetic algebra only; no SHD scientific evaluation','checks':checks,'PASS':all(c['PASS'] for c in checks)})
print(len(checks),all(c['PASS'] for c in checks));assert all(c['PASS'] for c in checks)
