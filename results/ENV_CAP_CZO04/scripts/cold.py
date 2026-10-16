"""Independent cold reconstruction. Does not import the main model/information functions."""
from pathlib import Path
import csv,json,struct,hashlib,time
import numpy as np
from scipy.stats import chi2
O=Path(__file__).resolve().parents[1];R=O.parents[1];W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_CZO04');checks=[]
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def ck(n,b,detail=''):
 checks.append(dict(check=n,PASS=bool(b),detail=str(detail)))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tab(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

def raw(level):
 aa=[]
 for env in ['JW','JA','ANNUAL']:
  for dep in [150,200,250]:
   data=Path(f'D:/ProjectStorage/WSL-CZ/ENV_CAP/E1/SCS17_{env}_z{dep}/I{level}.shd').read_bytes();stride=4*int.from_bytes(data[:4],'little',signed=True)
   nd=struct.unpack('<7i',data[2*stride:2*stride+28]);z=struct.unpack('<81f',data[8*stride:8*stride+324]);r=np.array(struct.unpack('<881f',data[9*stride:9*stride+3524]));j=z.index(200.)
   a=struct.unpack('<1762f',data[(10+j)*stride:(10+j)*stride+7048]);en=np.array([float(a[2*k])**2+float(a[2*k+1])**2 for k in range(881)]);v=np.full(881,np.nan);v[en>0]=10*np.log10(en[en>0]);aa.append(v)
 return r,np.array(aa)
t=np.linspace(0,1200,121);nuis=np.array([np.ones(121),(t-600)/600]).T

def geometry(s):
 s=np.atleast_2d(s);tau=t.reshape(1,-1);heading=np.deg2rad(15)
 platform=np.column_stack([2*np.minimum(t,600)+2*np.maximum(t-600,0)*np.cos(heading),2*np.maximum(t-600,0)*np.sin(heading)])
 start=s[:,0,None]*1000*np.stack([np.cos(np.deg2rad(s[:,1])),np.sin(np.deg2rad(s[:,1]))],axis=1)
 velocity=s[:,2,None]*np.stack([np.cos(np.deg2rad(s[:,3])),np.sin(np.deg2rad(s[:,3]))],axis=1)
 positions=start[:,None,:]+velocity[:,None,:]*t[None,:,None]-platform[None,:,:]
 rr=np.sqrt(np.sum(positions**2,axis=2));a=abs(np.arctan2(positions[:,:,1],positions[:,:,0])-np.deg2rad(np.where(t>600,15.,0)))
 return a,rr

def lin(rr,r,y):
 ix=np.searchsorted(r,rr,side='right')-1;ix=np.clip(ix,0,len(r)-2);f=(rr-r[ix])/(r[ix+1]-r[ix]);return y[ix]*(1-f)+y[ix+1]*f

def projected(a,wh,c):
 x=(np.asarray(a).reshape(-1,121)@wh.T).T
 if c=='C2':return np.zeros_like(np.asarray(a))
 d=wh@nuis[:,:1 if c=='C0' else 2];coef=np.linalg.lstsq(d,x,rcond=None)[0];rr=(x-d@coef).T
 return rr.reshape(np.asarray(a).shape)

def efficient(j):
 # Independent least-squares profile, no information-matrix pseudoinverse.
 a=np.linalg.lstsq(j[:,1:],j[:,0],rcond=1e-10)[0];v=j[:,0]-j[:,1:]@a;return float(v@v)

def main():
 start=time.perf_counter();cfg=read(O/'DESIGN_FREEZE.json')
 for a in cfg['execution_inputs']:ck('frozen:'+a['path'],sha(a['path'])==a['sha256'])
 for a in read(O/'PRIMARY_OUTPUT_FREEZE.json')['files']:ck('primary-seal:'+a['path'],sha(a['path'])==a['sha256'])
 curves={i:raw(i) for i in [1,2]};ranges,yy=curves[2];ll=curves[1][1]
 rows=tab(O/'SINGLE_DEPTH_INPUT_LEVELS.csv')
 for v in rows:
  p=int(v['profile']);k=int(round((float(v['range_m'])-10000)/125));a=float(v['I2_level_db']) if v['I2_level_db'] else np.nan
  ck(f'raw-point:{p}/{k}',(np.isnan(a) and np.isnan(yy[p,k])) or abs(a-yy[p,k])<1e-11)
 g=np.round(np.array([[a,b,c,d] for a in range(45,61) for b in np.arange(-5,5.01,.5) for c in np.arange(1,3.001,.2) for d in range(-15,16)]),10)
 supports=tab(O/'FINITE_CANDIDATE_SUPPORT.csv');local=tab(O/'LOCAL_EFFECTIVE_INFORMATION.csv');bstats=tab(O/'BEARING_FINITE_GRID_BASELINE.csv')
 bucket={};lbucket={}
 for s in supports:bucket.setdefault(s['scene'],[]).append(s)
 for s in local:lbucket.setdefault(s['scene'],[]).append(s)
 maxcost=0;maxinfo=0;fullcases=0
 resources=tab(O/'RESOURCE_SEPARATION_SCALES.csv')
 for br in bstats:
  name=br['scene'];z=np.load(W/(name+'_SCORES.npz'));h=z['state'];b,rr=geometry(h);bg,rg=geometry(g);bc=np.sum(((bg-b)/np.deg2rad(.1))**2,axis=1);ids=np.flatnonzero(bc<=chi2.ppf(.975,121))
  ck(name+'/full-grid-identity',np.array_equal(ids,z['bearing_ids']));ck(name+'/trajectory-range',np.max(abs(rr[0]-z['truth_range_m']))<1e-8);ck(name+'/trajectory-unsigned',np.max(abs(b[0]-z['truth_unsigned_bearing_rad']))<1e-12)
  truth={i:np.array([lin(rr[0],ranges,y) for y in curves[i][1]]) for i in [1,2]}
  pred={i:np.stack([lin(rg[ids],ranges,y) for y in curves[i][1]],axis=1) for i in [1,2]}
  cache={};whs={}
  for rho in [0,.9]:
   corr=rho**np.abs(np.arange(121)[:,None]-np.arange(121)[None,:]);wh=np.linalg.solve(np.linalg.cholesky(corr),np.eye(121));whs[rho]=wh
   for c in ['C0','C1','C2']:
    for i in [1,2]:
     pp=projected(pred[i],wh,c).reshape(-1,121);tt=projected(truth[i],wh,c);dist=np.maximum(np.sum(pp**2,axis=1)[:,None]+np.sum(tt**2,axis=1)[None,:]-2*pp@tt.T,0).reshape(len(ids),9,9).min(axis=1).T
     key=f'I{i}_rho{rho}_{c}';diff=float(np.max(abs(z[key]-dist))) if dist.size else 0;maxcost=max(diff,maxcost)
     ck(name+'/'+key+'/scores',diff<1e-6,diff);cache[(i,rho,c)]=dist
  for s in bucket[name]:
   rho=float(s['rho']);sigma=float(s['sigma_db']);c=s['condition'];tp=int(s['truth_profile']);cut=0 if c=='C2' else chi2.ppf(.975,120 if c=='C0' else 119)*sigma**2
   a=cache[(2,rho,c)][tp]<=cut+1e-10;al=cache[(1,rho,c)][tp]<=cut+1e-10
   ck(name+'/'+str((tp,rho,c,sigma))+'/support',int(a.sum())==int(s['I2_retained']) and int(al.sum())==int(s['I1_retained']) and int((a!=al).sum())==int(s['disagreement']))
   ck(name+'/'+str((tp,rho,c,sigma))+'/range',not a.any() or abs(np.ptp(g[ids[a],0])-float(s['range_width_km']))<1e-10)
   if c=='C2':ck(name+'/C2=baseline',a.all() and al.all())
  for q in [v for v in resources if v['scene']==name]:
   tp=int(q['truth_profile']);rho=float(q['rho']);c=q['condition'];ds=cache[(2,rho,c)][tp];wrong=np.flatnonzero(abs(g[ids,0]-h[0])>=1);kk=wrong[np.argmin(ds[wrong])];gap=np.sqrt(ds[kk]/121)
   ck(name+'/'+str((tp,rho,c))+'/resource-gap',abs(gap-float(q['whitened_residual_RMS_db']))<1e-6)
  for s in lbucket[name]:
   tp=int(s['truth_profile']);rho=float(s['rho']);c=s['condition'];wh=whs[rho];jbs=[];jls=[]
   for k,step in enumerate([.001,.0001,.0001,.0001]):
    a=h.copy();d=h.copy();a[k]+=step;d[k]-=step;ba,ra=geometry(a);bb,rb=geometry(d)
    scale=[1,1,.2,1][k];jbs.append(((ba-bb)[0]/(2*step))*scale);jls.append((lin(ra[0],ranges,yy[tp])-lin(rb[0],ranges,yy[tp]))/(2*step)*scale)
   bj=np.array(jbs).T/np.deg2rad(.1);jl=projected(np.array(jls),wh,c).T/.25;j=np.concatenate([bj,jl]);val=efficient(j);expected=float(s['I2_joint_range_info']);err=abs(val-expected)/max(expected,1e-20);maxinfo=max(maxinfo,err)
   ck(name+'/'+str((tp,rho,c))+'/profile-info',err<1e-5,err)
   ff=np.array(json.loads(s['Fisher']));ck(name+'/'+str((tp,rho,c))+'/Fisher',np.max(abs(ff-j.T@j))<max(1e-6,np.max(abs(ff))*1e-7))
  fullcases+=1;print('cold',name,flush=True)
 ck('all36-horizontal-cases',fullcases==36)
 ck('324-generating-profiles',len({(s['scene'],s['truth_profile']) for s in supports})==324)
 # Independent strict nuisance controls.
 for rho in [0,.9]:
  wh=whs[rho];v=np.sin(t/137)+.002*t
  ck('constant-control-'+str(rho),np.linalg.norm(projected(np.ones(121)*8,wh,'C0'))<1e-10)
  ck('linear-control-'+str(rho),np.linalg.norm(projected(.003*t+7,wh,'C1'))<1e-10)
  ck('free-source-control-'+str(rho),np.linalg.norm(projected(v,wh,'C2'))==0)
 summary=read(O/'CZO04_ANALYSIS_SUMMARY.json');prim=[s for s in supports if s['rho']=='0' and s['sigma_db']=='0.25' and s['condition']=='C1'];shrink=[s for s in prim if int(s['I2_retained'])>0 and float(s['range_width_km'])<float(s['bearing_range_width_km'])]
 ck('primary-denominator',len(prim)==summary['primary_C1_sigma025_iid_cases']);ck('primary-shrink',len(shrink)==summary['primary_nonempty_range_shrink_cases'])
 result=dict(PASS=all(x['PASS'] for x in checks),checks=len(checks),FAIL=sum(not x['PASS'] for x in checks),max_raw_cost_difference=maxcost,max_effective_info_relative_difference=maxinfo,full_horizontal_cases_rebuilt=fullcases,scope='independent raw SHD struct + explicit geometry + linear interpolation + least-squares nuisance profile; all finite supports and local matrices',elapsed_s=time.perf_counter()-start)
 save(O/'CZO04_INDEPENDENT_COLD_AUDIT.json',result)
 with (O/'INDEPENDENT_AUDIT_CHECKS.csv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['check','PASS','detail']);w.writeheader();w.writerows(checks)
 print(json.dumps(result));return 0 if result['PASS'] else 2
if __name__=='__main__':raise SystemExit(main())
