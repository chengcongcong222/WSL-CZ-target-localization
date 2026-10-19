from pathlib import Path
import csv,json,math,hashlib,numpy as np
O=Path(__file__).resolve().parents[1];W=Path(r'D:\ProjectStorage\WSL-CZ\ENV_CAP_CAP06/raw')
def tab(n,rows):
 with (O/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def q(v,p):
 if not v:return 'NOT_EVALUATED'
 a=sorted(v);x=a[math.ceil(len(a)*p)-1];return x if math.isfinite(x) else 'INF'
def wilson(k,n):
 if not n:return [None,None]
 z=1.959963984540054;p=k/n;return [(p+z*z/(2*n)+sg*z*math.sqrt(p*(1-p)/n+z*z/(4*n*n)))/(1+z*z/n) for sg in [-1,1]]
def main():
 groups={}
 for r in csv.DictReader((O/'P2_RECEIVER_BEARING_CASES.csv').open()):
  if r['stress']=='MAIN_TRACK':groups.setdefault((r['token'],r['mode']),[]).append(r)
 tracks=[]
 for (token,mode),rr in groups.items():
  n=sum(x['present']=='True' for x in rr);k=sum(x['target_component_extracted']=='True' for x in rr);cov=sum(x['target_u_covered']=='True' for x in rr);errs=[float(x['axis_error_deg']) for x in rr if x['present']=='True'];r=rr[0]
  tracks.append(dict(token=token,mode=mode,split=r['split'],state_id=r['state_id'],repeat=r['repeat'],side='PLUS' if int(r['repeat'])%2==0 else 'MINUS',node=r['node'],snr=r['snr'],present_epochs=n,extracted=k,covered=cov,extraction_fraction=k/n if n else None,coverage_fraction=cov/n if n else None,axis_P50_deg=q(errs,.5),axis_P95_deg=q(errs,.95),all_present_epochs_covered=cov==n,scope='one trajectory; paired across nodes and quality; evaluation-only association'))

 strata={}
 for r in csv.DictReader((O/'P2_RECEIVER_BEARING_CASES.csv').open()):
  if r['stress']=='MAIN_TRACK' and r['present']=='True':strata.setdefault((r['split'],r['snr'],r['mode'],r['false_source']),[]).append(r)
 diagnosis=[]
 for key,rr in strata.items():
  diagnosis.append(dict(split=key[0],snr=key[1],mode=key[2],distractor_present=key[3],present_epochs=len(rr),extracted=sum(x['target_component_extracted']=='True' for x in rr),covered=sum(x['target_u_covered']=='True' for x in rr),scope='POSTHOC_DESCRIPTIVE_STRATIFICATION_OF_SEALED_RESULTS; no threshold changes'))
 tab('P2_DISTRACTOR_STRATIFICATION.csv',diagnosis)
 tab('P2_TRAJECTORY_METRICS.csv',tracks)
 out=[];gg={}
 for r in tracks:gg.setdefault(tuple(r[k] for k in ['split','state_id','node','snr','mode','side']),[]).append(r)
 for key,rr in gg.items():
  k=sum(x['all_present_epochs_covered'] for x in rr);lo,hi=wilson(k,len(rr));out.append(dict(zip(['split','state_id','node','snr','mode','side'],key),n=len(rr),all_epoch_covered=k,Wilson95_lo=lo,Wilson95_hi=hi,median_track_axis_P50_deg=q([float(x['axis_P50_deg']) for x in rr],.5),P95_track_axis_P50_deg=q([float(x['axis_P50_deg']) for x in rr],.95),scope='small per-state/per-side descriptive sample; no acoustic generalization'))
 tab('P2_STATE_SIDE_SUMMARY.csv',out)
 geometry=[]
 for r in csv.DictReader((O/'P2_RECORD_INDEX.csv').open()):
  if r['stress']!='MAIN_TRACK' or r['snr']!='20':continue
  a=np.load(W/'p2/truth_private'/(r['token']+'.npz'));s=a['state'];t=a['true_times'];target=s[0]*np.array([np.cos(s[1]),np.sin(s[1])])+t[:,None]*s[3]*np.array([np.cos(s[4]),np.sin(s[4])]);d=np.linalg.norm(target-a['true_positions'],axis=1)
  geometry.append(dict(token=r['token'],split=r['split'],state_id=r['state_id'],node=r['node'],repeat=r['repeat'],min_range_km=float(d.min()/1000),max_range_km=float(d.max()/1000),epochs_outside_45_60_km=int(np.sum((d<45000)|(d>60000))),scope='geometry diagnostic only; 45-60km not universal CZ condition or detectability'))
 tab('CAP06_NODE_RANGE_RESOURCES.csv',geometry)
 print('descriptive track rows',len(tracks),'state/side rows',len(out),'node ranges',len(geometry))
if __name__=='__main__':main()
