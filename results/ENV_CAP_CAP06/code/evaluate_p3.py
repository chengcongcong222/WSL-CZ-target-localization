from common import *
def main():
 gate=read('P3_GATE.json')
 if not gate['eligible']:
  for name in ['P3_HOLDOUT_CASE_STATS.csv','P3_CONDITIONAL_ERROR_TABLE.csv']:tab(name,[dict(status='NOT_EVALUATED',reason='No admitted extractor')])
  return
 records=list(csv.DictReader((O/'P3_SEARCH_INDEX.csv').open(encoding='utf-8')));out=[]
 for r in records:
  a=np.load(r['raw_path']);truth=np.load(RAW/'p2/truth_private'/(r['source_main']+'.npz'));s=truth['state'];x=np.array([s[0]*np.cos(s[1]),s[0]*np.sin(s[1]),s[3]*np.cos(s[4]),s[3]*np.sin(s[4])]);boxes=a['kept']
  retained=bool(np.any(np.all(x>=boxes[:,0]-1e-7,axis=1)&np.all(x<=boxes[:,1]+1e-7,axis=1))) if len(boxes) else False
  b=a['best'];accepted=r['accepted_point']=='True';est=np.array([np.hypot(*b[:2]),np.arctan2(b[1],b[0]),np.hypot(*b[2:]),np.arctan2(b[3],b[2])]) if np.isfinite(b).all() else np.full(4,np.nan)
  tr=s[[0,1,3,4]];err=abs(est-tr);err[[1,3]]=abs((est[[1,3]]-tr[[1,3]]+np.pi)%(2*np.pi)-np.pi)*180/np.pi
  if len(boxes):
   lo=boxes[:,0].min(axis=0);hi=boxes[:,1].max(axis=0)
   rlo=np.hypot(lo[0],max(0,lo[1],-hi[1]));rhi=max(np.hypot(xx,yy) for xx in [lo[0],hi[0]] for yy in [lo[1],hi[1]])
   vlo=np.hypot(lo[2],max(0,lo[3],-hi[3]));vhi=max(np.hypot(xx,yy) for xx in [lo[2],hi[2]] for yy in [lo[3],hi[3]])
   rw=max(0,min(60000,rhi)-max(45000,rlo));vw=max(0,min(3,vhi)-max(1,vlo))
  else:rw=vw=None
  row=dict(**r,truth_in_retained_outer=retained,range_width_m=rw,speed_width_mps=vw,point_scope='BEST_VISITED_CENTER_NOT_GLOBAL_OPTIMUM',true_r_m=float(s[0]),true_v_mps=float(s[3]))
  for key,e in zip(['range_m','bearing_deg','speed_mps','heading_deg'],err):row[key]=float(e) if accepted and np.isfinite(e) else 'INF';row[key+'_conditional']=float(e) if accepted and np.isfinite(e) else None
  row['range_pct']=float(err[0]/s[0]*100) if accepted else 'INF';row['speed_pct']=float(err[2]/s[3]*100) if accepted else 'INF'
  out.append(row)
 tab('P3_GLOBAL_SUPPORT_STATUS.csv',out)
 stats=[]
 for split in ['development','holdout']:
  for res in ['S0','S2_B0','S2_B3']:
   z=[r for r in out if r['split']==split and r['resource']==res]
   for state in ['ALL']+sorted(set(r['state_id'] for r in z)):
    a=z if state=='ALL' else [r for r in z if r['state_id']==state];n=len(a);accept=sum(r['accepted_point']=='True' for r in a);cover=sum(r['truth_in_retained_outer'] for r in a)
    for metric in ['range_m','range_pct','bearing_deg','speed_mps','speed_pct','heading_deg']:
     q=qtiles([np.inf if r[metric]=='INF' else r[metric] for r in a]);cond=qtiles([r[metric] for r in a if r[metric]!='INF'])
     stats.append(dict(split=split,resource=res,state_id=state,metric=metric,n=n,accepted=accept,rejected=n-accept,truth_outer_retained=cover,coverage_Wilson95_lo=wilson(cover,n)[0],coverage_Wilson95_hi=wilson(cover,n)[1],P50=q[0],P90=q[1],P95=q[2],P99=q[3],conditional_P50=cond[0],conditional_P95=cond[2],scope='MATCHED_PLANE_WAVE_PRODUCT_CONTROL; finite search; full unresolved boxes retained; not CZ performance'))
 tab('P3_CONDITIONAL_ERROR_TABLE.csv',stats);tab('P3_HOLDOUT_CASE_STATS.csv',[r for r in stats if r['split']=='holdout'])
 pairs=[]
 for base in sorted(set(r['source_main'] for r in out)):
  d={r['resource']:r for r in out if r['source_main']==base}
  for res in ['S2_B0','S2_B3']:
   for metric in ['range_m','bearing_deg','speed_mps','heading_deg']:
    x=d['S0'][metric];y=d[res][metric];x=np.inf if x=='INF' else x;y=np.inf if y=='INF' else y
    pairs.append(dict(source_main=base,resource=res,metric=metric,improved=y<x,worsened=y>x,equal=y==x,paired=True))
 tab('P3_PAIRED_CHANGE.csv',pairs)
if __name__=='__main__':main()
