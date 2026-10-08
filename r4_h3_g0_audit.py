"""Independent cold pressure/analytic derivative and QR information reconstruction. Never reruns science."""
import json
from pathlib import Path
import numpy as np
from scipy import linalg
import r4_h3_g0 as p
import r4_e2_nine_pair_pilot_audit as binary
checks=[]
def ck(name,passed,detail=''):checks.append(dict(check=name,PASS=bool(passed),detail=str(detail)))
def qrremove(a,n):
 if not n.size:return a.copy()
 v=np.linalg.norm(n,axis=0);n=n[:,v>1e-12]/v[v>1e-12]
 if not n.size:return a.copy()
 q,r,piv=linalg.qr(n,mode='economic',pivoting=True)
 sv=linalg.svdvals(n);rank=int(np.sum(sv>sv[0]*1e-10));q=q[:,:rank]
 return a-q@(q.T@a)
def cold(mods,st,res):
 # Construct Cartesian geometry independently; no polar reconstruction.
 tt=np.array([0.,600.,1200.]);post=np.maximum(tt-600.,0)
 xy=np.array(st[:2])[None,:]+tt[:,None]*np.array(st[2:])
 main=np.column_stack([2*np.minimum(tt,600)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)])
 offsets=np.arange(-7.,8.,2.);zr=np.full(8,200.)
 if res!='B0':offsets=np.r_[offsets,[-1.,-1.,-1.,-1.]];zr=np.r_[zr,[198.,199.,201.,202.] if res=='B1' else [200.]*4]
 delta=xy[:,None,:]-main[:,None,:]-np.column_stack([offsets,np.zeros(len(offsets))])[None,:,:]
 rr=np.sqrt(delta[:,:,0]**2+delta[:,:,1]**2)
 rnorm=np.sqrt(st[0]**2+st[1]**2);vnorm=np.sqrt(st[2]**2+st[3]**2)
 geo=np.zeros((3,2,5));geo[:,:,0]=np.array(st[:2])/rnorm;geo[:,:,1]=[-st[1],st[0]]
 geo[:,:,3]=tt[:,None]*np.array(st[2:])/vnorm;geo[:,:,4]=tt[:,None]*np.array([-st[3],st[2]])
 ps=[];js=[];hs=[];roundscale=[]
 for zz,phi,k in mods:
  ix={float(z):i for i,z in enumerate(zz)};pr=phi[[ix[z] for z in zr]]
  recv=np.tile(pr,(3,1)).T
  r=rr.ravel()
  E=np.sqrt(2*np.pi/(k[:,None]*r))*np.exp(-1j*k[:,None]*r-1j*np.pi/4)
  term=phi[ix[200.],:,None]*recv*E
  pressure=np.add.reduce(term,axis=0).reshape(rr.shape);ps.append(pressure)
  rad=np.add.reduce(term*(-1j*k[:,None]-1/(2*r)),axis=0).reshape(rr.shape)
  grad=rad[:,:,None]*delta/rr[:,:,None]
  J=np.stack([np.sum(grad*geo[:,None,:,i],axis=-1)*p.SCALE[i] for i in range(5)],axis=-1)
  H=J.copy()
  for step,acc in [(1.,J),(.5,H)]:
   d=(phi[ix[200.+step]]-phi[ix[200.-step]])/(2*step)
   acc[...,2]=np.add.reduce(d[:,None]*recv*E,axis=0).reshape(rr.shape)*p.SCALE[2]
  js.append(J);hs.append(H)
 return np.stack(ps,axis=-1),np.stack(js,axis=-2),np.stack(hs,axis=-2)
def main():
 o=p.O;f=p.read(o/'DESIGN_FREEZE.json');d=p.read(o/'H3_G0_DECISION.json')
 for n,r in f['bindings'].items():ck('binding:'+n,p.sha(n,r['kind']=='RAW')==r['sha256'])
 info=p.rows(o/'VERTICAL_INFORMATION_BY_SCENE.csv') if (o/'VERTICAL_INFORMATION_BY_SCENE.csv').exists() else []
 idx={(r['scene_id'],int(r['mesh']),r['resource'],r['noise_model'],r['calibration'],float(r['sigma'])):r for r in info}
 coldmax=0.;jmax=0.;profilemax=0.;rankmax=0.
 for sc in f['scenes']:
  for grid in p.GRIDS:
   mods=[binary.mod(p.paths.path_for(ff,grid)) for ff in p.RAW]
   for res in ['B0','B1','B2']:
    a=np.load(o/f"INPUT_{sc['scene_id']}_{res}_n{grid}.npz")
    pp,jj,hh=cold(mods,sc['state'],res)
    pressure_error=float(np.max(abs(pp-a['pressure']))/np.max(abs(pp)));coldmax=max(coldmax,pressure_error)
    ck('termwise_pressure_roundoff:'+sc['scene_id']+str(grid)+res,bool(np.all(abs(pp-a['pressure'])<=a['roundoff_bound'])))
    for key,jc in [('J',jj),('J_half',hh)]:
     rel=jc/pp[...,None]
     err=float(np.linalg.norm(rel-a[key])/np.linalg.norm(rel));jmax=max(jmax,err)
     ck('independent_analytic_J:'+sc['scene_id']+str(grid)+res+key,err<1e-9,err)
    if not info:continue
    for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
     for cal in ['C0','C1a','C1b','C2','B3_KNOWN_SOURCE_AND_CALIBRATION']:
      y,n=p.system(a['pressure'],a['J'],noise,cal,cal.startswith('B3'))
      z=qrremove(y,n) if cal!='C2' else np.zeros_like(y)
      z0=p.projected(a['pressure'],a['J'],noise,cal,cal.startswith('B3'))
      e=np.linalg.norm(z.T@z-z0.T@z0)/max(np.linalg.norm(z0.T@z0),1e-30)
      ck('QR_vs_SVD_F:'+sc['scene_id']+str(grid)+res+noise+cal,e<1e-7,e)
      for sig in [.01,.05]:
       row=idx[sc['scene_id'],grid,res,noise,cal,sig];w=z/sig
       _,sv,vh=linalg.svd(w,full_matrices=False);rank=int(np.sum(sv>sv[0]*1e-10)) if sv[0]>0 else 0
       null=vh[rank:];nz=np.linalg.norm(null[:,2]) if len(null) else 0
       rem=qrremove(w[:,2:3],w[:,[0,1,3,4]])
       iv=float(np.sum(rem*rem));iz=iv/p.SCALE[2]**2 if iv>(max(sv[0],1e-30)*1e-10)**2 and nz<=1e-8 else 0.
       err=abs(iz-float(row['depth_information_horizontal_profiled']))/max(iz,float(row['depth_information_horizontal_profiled']),1e-20);profilemax=max(profilemax,err)
       ck('independent_effective_depth:'+sc['scene_id']+str(grid)+res+noise+cal+str(sig),err<1e-5 and rank==int(row['rank']),err)
  print('AUDIT',sc['scene_id'],flush=True)
 ctrl=p.rows(o/'CONTROL_CHECKS.csv')
 ck('registered_controls_pass',all(r['PASS']=='True' for r in ctrl))
 if info:
  ck('complete_information_rows',len(info)==720)
  ck('finiteF_complete',len(p.rows(o/'FINITE_HORIZONTAL_LABEL_DIAGNOSTIC.csv'))==4500)
  ck('saturated_C2_zero',all(float(r['depth_information_horizontal_profiled'])==0 and r['depth_local_scale_m']=='INF' for r in info if r['calibration']=='C2'))
 # Scientific grid stability is a separate Gate, not a reconstruction failure.
 ck('no_false_project_credit',d.get('R4_percent',0)==0)
 ck('unchanged_geometry_and_input_sources',all(r['PASS'] for r in checks if r['check'].startswith('binding:')))
 p.write(o/'INDEPENDENT_AUDIT_CHECKS.csv',checks)
 v=dict(checks=len(checks),PASS=sum(r['PASS'] for r in checks),FAIL=sum(not r['PASS'] for r in checks),cold_pressure_relative_max=coldmax,independent_J_relative_max=jmax,independent_depth_information_relative_max=profilemax,scope='COLD_READ_ONLY_FORMULA_AND_INFORMATION_RECONSTRUCTION',new_solver=0,new_MC=0)
 p.dump(o/'VALIDATION.json',v)
 if v['FAIL']:d.update(decision='H3_G0_IMPLEMENTATION_INVALID',independent_audit='FAIL',next='STOP')
 else:d['independent_audit']='PASS'
 p.dump(o/'H3_G0_DECISION.json',d)
 print(json.dumps(p.native(v),indent=2))
if __name__=='__main__':main()
