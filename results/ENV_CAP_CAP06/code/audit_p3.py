from pathlib import Path
import json,csv,numpy as np,math,hashlib
O=Path(__file__).resolve().parents[1];W=Path(r'D:\ProjectStorage\WSL-CZ\ENV_CAP_CAP06/raw')
def main():
 c=[]
 def ck(n,b,e=0):c.append(dict(check_id=n,PASS=bool(b),error=float(e)))
 gate=json.loads((O/'P3_GATE.json').read_text())
 if not gate['eligible']:
  ck('legitimate_not_evaluated',not json.loads((O/'P2_GATE.json').read_text())['P3_eligible'] or not json.loads((O/'P2_INDEPENDENT_COLD_AUDIT.json').read_text())['PASS'])
 else:
  rows=list(csv.DictReader((O/'P3_GLOBAL_SUPPORT_STATUS.csv').open()));tol=1e-7
  root=np.array([[45000*math.cos(math.radians(5)),-60000*math.sin(math.radians(5)),math.cos(math.radians(15)),-3*math.sin(math.radians(15))],[60000,60000*math.sin(math.radians(5)),3,3*math.sin(math.radians(15))]])
  for row in rows:
   a=np.load(row['raw_path']);K=a['kept'];D=a['discarded'];leaves=np.concatenate((K,D));volume=sum(np.prod(b[1]-b[0]) for b in leaves);rootvol=np.prod(root[1]-root[0]);ck(row['token']+':partition_volume',abs(volume/rootvol-1)<1e-10,abs(volume/rootvol-1))
   # Every final leaf must form a complete recursive dyadic partition, not just equal total volume.
   keys={tuple(b.ravel()) for b in leaves};actualroot=np.array([leaves[:,0].min(axis=0),leaves[:,1].max(axis=0)]);stack=[actualroot];valid=len(keys)==len(leaves) and np.max(abs(actualroot-root))<1e-8
   while stack and valid:
    box=stack.pop();key=tuple(box.ravel())
    if key in keys:keys.remove(key);continue
    if len(stack)>5000 or np.max(box[1]-box[0])<1e-9:valid=False;break
    k=int(np.argmax((box[1]-box[0])/[16000,11000,2.1,1.6]));m=(box[0,k]+box[1,k])/2;b=box.copy();bb=box.copy();b[1,k]=m;bb[0,k]=m;stack.extend([b,bb])
   ck(row['token']+':disjoint_exhaustive_partition',valid and not keys)
   main=row['source_main'];base=main.removesuffix('_MAIN_20');nodes=['MAIN']+([] if row['resource']=='S0' else [row['resource'].split('_')[1]]);obs=[]
   for node in nodes:
    z=json.loads((W/'p2/products'/(base+'_'+node+'_20_CAL.json')).read_text())
    for j,e in enumerate(z['epochs']):
     if e['products']:obs.append((z,j,e))
   valid=True
   for b,idx in zip(D,a['witness']):
    z,j,e=obs[int(idx)];t=z['timestamps'][j];pad=6*z['nav_sigma_axis_m']+1.250001;lo=b[0,:2]+t*b[0,2:]-z['reported_positions'][j]-pad;hi=b[1,:2]+t*b[1,2:]-z['reported_positions'][j]+pad
    az=[math.atan2(y,x) for x in [lo[0],hi[0]] for y in [lo[1],hi[1]]];turn=abs(t-600)<=.25;hlo=0 if turn else z['nominal_heading'][j];hhi=math.radians(15) if turn else hlo
    low=min(az)-hhi-e['beta_est_rad']-math.radians(.1)-1e-12;high=max(az)-hlo-e['beta_est_rad']+math.radians(.1)+1e-12
    values=[math.cos(low),math.cos(high)]+([1.] if low<=0<=high else []);mn=min(values)-1e-14;mx=max(values)+1e-14
    valid &= lo[0]>0 and all(pp['u_hi']<mn or pp['u_lo']>mx for pp in e['products'])
   ck(row['token']+':all_exclusion_witnesses',valid)
   tr=np.load(W/'p2/truth_private'/(main+'.npz'))['state'];x=np.array([tr[0]*math.cos(tr[1]),tr[0]*math.sin(tr[1]),tr[3]*math.cos(tr[4]),tr[3]*math.sin(tr[4])]);cover=any(np.all(x>=b[0]-tol) and np.all(x<=b[1]+tol) for b in K)
   ck(row['token']+':coverage',cover==(row['truth_in_retained_outer']=='True'))
   h=hashlib.sha256(Path(row['raw_path']).read_bytes()).hexdigest();ck(row['token']+':raw_hash',h==row['sha256'])
   est=a['best'];values=[math.hypot(*est[:2]),math.atan2(est[1],est[0]),math.hypot(*est[2:]),math.atan2(est[3],est[2])];truth=tr[[0,1,3,4]]
   for j,key in enumerate(['range_m','bearing_deg','speed_mps','heading_deg']):
    err=abs(values[j]-truth[j]) if j%2==0 else abs((values[j]-truth[j]+math.pi)%(2*math.pi)-math.pi)*180/math.pi
    ck(row['token']+':error_'+key,(row['accepted_point']=='False' and row[key]=='INF') or (row['accepted_point']=='True' and abs(float(row[key])-err)<1e-7))
  for s in csv.DictReader((O/'P3_CONDITIONAL_ERROR_TABLE.csv').open()):
   a=[r for r in rows if r['split']==s['split'] and r['resource']==s['resource'] and (s['state_id']=='ALL' or r['state_id']==s['state_id'])];v=sorted(float(r[s['metric']]) for r in a);good=len(a)==int(s['n'])
   for quant,q in [('P50',.5),('P90',.9),('P95',.95),('P99',.99)]:
    val=v[math.ceil(q*len(v))-1];good &= (math.isinf(val) and s[quant]=='INF') or abs(val-float(s[quant]))<1e-9
   ck('quantiles:'+str((s['split'],s['resource'],s['state_id'],s['metric'])),good)
 out=dict(PASS=all(x['PASS'] for x in c),checks=len(c),FAIL=sum(not x['PASS'] for x in c),records=c,scope='Conditional interval geometry and saved distributions; not global optimum or real coverage')
 (O/'P3_INDEPENDENT_COLD_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print(out['PASS'],out['checks'])
if __name__=='__main__':main()
