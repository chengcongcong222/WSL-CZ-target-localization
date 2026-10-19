"""Conservative Cartesian boxes conditional on empirical angular product sets.
No true state, original panel or private generator output may be read.
"""
from common import *
import sys,heapq

def protected(event,args):
 if event=='open' and isinstance(args[0],str) and ('truth_private' in args[0] or 'P0_FROZEN_PANEL' in args[0]):raise PermissionError('Estimator private read')
sys.addaudithook(protected)
SCALE=np.array([16000.,11000.,2.1,1.6])
def tests(box,obs):
 # box in [x0,y0,vx,vy]; all possible positions in it and navigation/time bounds.
 lo=box[0];hi=box[1];bad=-1
 for j,d in enumerate(obs):
  t=d['t'];pad=6*d['sigma']+1.250001;lower=lo[:2]+t*lo[2:]-d['pos']-pad;upper=hi[:2]+t*hi[2:]-d['pos']+pad
  if lower[0]<=0:return True,-1 # conservative unresolved; expected domain x>0
  vals=np.arctan2(np.array([lower[1],upper[1],lower[1],upper[1]]),np.array([lower[0],lower[0],upper[0],upper[0]]))
  a=vals.min()-d['heading_hi']-d['beta']-np.deg2rad(.1)-1e-12;b=vals.max()-d['heading_lo']-d['beta']+np.deg2rad(.1)+1e-12
  cmin=np.nextafter(min(np.cos(a),np.cos(b)),-np.inf);cmax=np.nextafter(1. if a<=0<=b else max(np.cos(a),np.cos(b)),np.inf)
  if not any(v[0]<=cmax and v[1]>=cmin for v in d['sets']):bad=j;return False,bad
 return True,bad

def point_score(x,obs):
 s=0.
 for d in obs:
  p=x[:2]+d['t']*x[2:]-d['pos'];u=np.cos(np.arctan2(p[1],p[0])-(d['heading_lo']+d['heading_hi'])/2-d['beta']);s+=min(((max(v[0]-u,0,u-v[1])/(v[1]-v[0]+.01))**2 for v in d['sets']),default=0.)
 return s

def search(obs,budget=1023):
 initial=np.array([[45000*np.cos(np.deg2rad(5)),-60000*np.sin(np.deg2rad(5)),np.cos(np.deg2rad(15)),-3*np.sin(np.deg2rad(15))],[60000,60000*np.sin(np.deg2rad(5)),3,3*np.sin(np.deg2rad(15))]])
 queue=[initial];kept=[];discard=[];witness=[];count=0;best=None;score=np.inf
 while queue and count<budget:
  box=queue.pop(0);count+=1;ok,w=tests(box,obs)
  if not ok:discard.append(box);witness.append(w);continue
  mid=box.mean(axis=0);s=point_score(mid,obs)
  if s<score:score=s;best=mid.copy()
  sizes=(box[1]-box[0])/SCALE;k=int(np.argmax(sizes))
  if np.all(box[1]-box[0]<=np.array([100.,100.,.02,.02])):kept.append(box);continue
  a=box.copy();b=box.copy();m=mid[k];a[1,k]=m;b[0,k]=m;queue.extend([a,b])
 kept.extend(queue)
 return np.array(kept).reshape(-1,2,4),np.array(discard).reshape(-1,2,4),np.array(witness,int),best,score,count

def observations(products):
 obs=[]
 for z in products:
  for j,e in enumerate(z['epochs']):
   if not e['products']:continue # missing detection never excludes a state
   obs.append(dict(t=z['timestamps'][j],pos=np.array(z['reported_positions'][j]),heading=z['nominal_heading'][j],heading_lo=0. if abs(z['timestamps'][j]-600)<=.25 else z['nominal_heading'][j],heading_hi=np.deg2rad(15) if abs(z['timestamps'][j]-600)<=.25 else z['nominal_heading'][j],beta=e['beta_est_rad'],sigma=z['nav_sigma_axis_m'],sets=[(p['u_lo'],p['u_hi']) for p in e['products']]))
 return obs

def main():
 gate=read('P2_GATE.json');dest=RAW/'p3';dest.mkdir()
 if not gate['P3_eligible']:
  tab('P3_GLOBAL_SUPPORT_STATUS.csv',[dict(config='ALL',status='NOT_EVALUATED',reason='P2_G2_NOT_ADMITTED; O2_L3_NOT_ADMITTED')]);write('P3_GATE.json',dict(status='NOT_EVALUATED',eligible=[]));return
 records=list(csv.DictReader((O/'P2_RECORD_INDEX.csv').open(encoding='utf-8')));rows=[]
 mains=[r for r in records if r['node']=='MAIN' and r['snr']=='20']
 for ii,r in enumerate(mains):
  base=r['token'].removesuffix('_MAIN_20');paths=[RAW/'p2/products'/(base+'_'+node+'_20_CAL.json') for node in ['MAIN','B0','B3']];prod=[json.loads(p.read_text()) for p in paths]
  for resource,sub in [('S0',[prod[0]]),('S2_B0',[prod[0],prod[1]]),('S2_B3',[prod[0],prod[2]])]:
   obs=observations(sub);kept,discard,w,best,s,n=search(obs,read('CAP06_FROZEN_SPEC.json')['P3_box_budget'])
   token=base+'_'+resource;path=dest/(token+'.npz');np.savez_compressed(path,kept=kept,discarded=discard,witness=w,best=np.full(4,np.nan) if best is None else best,score=np.array(s))
   rows.append(dict(token=token,source_main=r['token'],split=r['split'],state_id=r['state_id'],repeat=r['repeat'],resource=resource,boxes_examined=n,retained_boxes=len(kept),excluded_boxes=len(discard),unresolved_boxes=len(kept),global_status='CONDITIONAL_OUTER_BUT_UNRESOLVED' if len(kept) else 'EMPTY_CONDITIONAL_SET',best_score=s if np.isfinite(s) else None,accepted_point=best is not None and s<=1.,raw_path=str(path),sha256=sha(path),observations=len(obs),z='NOT_EVALUATED',PHYSICAL_CZ_LINK='UNVERIFIED'))
  if ii%16==0:print('P3 trajectories',ii,flush=True)
 tab('P3_SEARCH_INDEX.csv',rows);write('P3_GATE.json',dict(status='PRODUCT_CONDITIONAL_OUTER_SEARCH_ONLY',eligible=gate['P3_eligible'],continuous_global_optimum='NOT_CERTIFIED',unresolved_retained=True,depth='NOT_EVALUATED'))
if __name__=='__main__':main()
