"""Post-fit evaluation only. Nearest-rank, INF and incomplete denominators explicit."""
from common import *

def run():
 raw=np.load(DATA/'terminals.npy',mmap_mode='r');states=np.array(js(DATA/'private/PANEL_AND_REJECTIONS.json')['states']);err=np.full((48000,6),np.inf)
 for q in range(2):
  for b in range(2):
   for c in range(12):
    ids=np.arange(index(q,b,c,0,0),index(q,b,c,499,1)+1);h=raw['h'][ids,:4];d=h-states[c];e=np.column_stack([abs(d[:,0]),abs(d[:,0])/states[c,0],abs((d[:,1]+180)%360-180),abs(d[:,2]),abs(d[:,2])/states[c,2],abs((d[:,3]+180)%360-180)]);ok=(raw['executed'][ids]==1)&(raw['accepted'][ids]==1);err[ids[ok]]=e[ok]
 np.savez_compressed(O/'EVALUATION_ARRAYS.npz',errors=err,executed=raw['executed'],accepted=raw['accepted'])
 np.savez_compressed(O/'TERMINALS_ARCHIVE.npz',terminals=np.array(raw))
 stats=[];diag=[];groups={};lookup={}
 for q,Q in enumerate(['Q1','Q2']):
  for b,B in enumerate(['B0','B3']):
   for p,P in enumerate(['P0','P1']):
    conf=f'{Q}_{B}_{P}'
    for c in list(range(12))+[-1]:
     case='POOLED' if c==-1 else f'D{c+1:02d}'
     for side in ['BOTH','-1','1']:
      reps=range(500) if side=='BOTH' else range(250) if side=='-1' else range(250,500)
      ids=np.array([index(q,b,cc,r,p) for cc in (range(12) if c==-1 else [c]) for r in reps]);planned=len(ids);done=ids[raw['executed'][ids]==1];valid=done[raw['accepted'][done]==1];complete=len(done)==planned;groups[conf,case,side]=done
      d=dict(config=conf,case_id=case,side=side,planned=planned,executed=len(done),accepted=len(valid),failed=len(done)-len(valid),not_evaluated=planned-len(done),complete=complete,output_rate=len(valid)/planned,failed_rate=(len(done)-len(valid))/planned)
      for name,arr in [('state_at_bound',np.any(raw['near'][done,:4],axis=1)),('bias_at_bound',np.any(raw['near'][done,4:],axis=1)),('common_at_bound',raw['near'][done,4]),('diff_at_bound',raw['near'][done,5]),('outside_original_bias',np.any(raw['outside_original'][done],axis=1)),('initial_fallback',raw['fallback'][done]),('iteration_limit',raw['solver_status'][done]==0),('weak_fit_rank',raw['rank'][done]<6)]:d[name+'_count']=int(np.sum(arr))
      diag.append(d)
      for mi,m in enumerate(METRICS):
       for mode,sel in [('unconditional',done),('accepted_only',valid)]:
        row=dict(config=conf,case_id=case,side=side,metric=m,mode=mode,planned=planned,n=len(sel),complete_panel=complete,scope='FULL_REGISTERED_GROUP' if complete else 'EXECUTED_SUBSET_ONLY')
        for label,f in [('P50',.5),('P90',.9),('P95',.95),('P99',.99)]:row[label]=quantile(err[sel,mi],f)
        row['maximum']=float(np.max(err[sel,mi])) if len(sel) else 'NOT_EVALUATED';stats.append(row);lookup[conf,case,side,m,mode]=row
      for mi,limit in [(1,.10),(2,.5),(4,.20),(5,5.)]:d[METRICS[mi]+'_large_success_count']=int(np.sum(err[valid,mi]>limit));d[METRICS[mi]+'_large_or_failed_fraction']=(d[METRICS[mi]+'_large_success_count']+len(done)-len(valid))/planned
 csvwrite(O/'STATISTICS.csv',stats);csvwrite(O/'OUTPUT_AND_BOUNDARY_DIAGNOSTICS.csv',diag)
 agg=[]
 for conf in sorted({x['config'] for x in stats}):
  for m in METRICS:
   rr=[lookup[conf,f'D{i+1:02d}','BOTH',m,'unconditional'] for i in range(12)];complete=all(x['complete_panel'] for x in rr);r=dict(config=conf,metric=m,complete_panel=complete)
   for label in ['P50','P90','P95','P99']:r['worst_scene_'+label]=max(x[label] for x in rr) if complete else 'NOT_EVALUATED'
   med=[x['P50'] for x in rr];r['mean_of_12_scene_P50']=float(np.mean(med)) if complete else 'NOT_EVALUATED';r['nearest_rank_median_of_12_scene_P50']=quantile(med,.5) if complete else 'NOT_EVALUATED';r['pooled_P50']=lookup[conf,'POOLED','BOTH',m,'unconditional']['P50'];r['pooled_P95']=lookup[conf,'POOLED','BOTH',m,'unconditional']['P95'];agg.append(r)
 csvwrite(O/'AGGREGATES.csv',agg)
 pairs=[];comparisons=[]
 for Q in ['Q1','Q2']:
  for P in ['P0','P1']:comparisons.append(('BASELINE',f'{Q}_B0_{P}',f'{Q}_B3_{P}'))
  for B in ['B0','B3']:comparisons.append(('SUPPORT',f'{Q}_{B}_P0',f'{Q}_{B}_P1'))
 for typ,a,b in comparisons:
  for c in [f'D{i+1:02d}' for i in range(12)]+['POOLED']:
   for side in ['BOTH','-1','1']:
    aa=groups[a,c,side];bb=groups[b,c,side];ia={((i//1000)%12,(i//2)%500):i for i in aa};ib={((i//1000)%12,(i//2)%500):i for i in bb};keys=sorted(ia.keys()&ib.keys());plan=6000 if c=='POOLED' else 500;plan=plan if side=='BOTH' else plan//2
    for mi,m in enumerate(METRICS):
     x=np.array([err[ia[k],mi] for k in keys]);y=np.array([err[ib[k],mi] for k in keys]);ff=np.isfinite(x)&np.isfinite(y);r=dict(comparison=typ,baseline=a,candidate=b,case_id=c,side=side,metric=m,planned=plan,paired_executed=len(keys),unpaired_not_evaluated=plan-len(keys),better=int(np.sum(y<x)),worse=int(np.sum(y>x)),equal=int(np.sum(y==x)),both_failed=int(np.sum(~np.isfinite(x)&~np.isfinite(y))),baseline_failed_candidate_valid=int(np.sum(~np.isfinite(x)&np.isfinite(y))),candidate_failed_baseline_valid=int(np.sum(np.isfinite(x)&~np.isfinite(y))),both_finite=int(np.sum(ff)),finite_difference_P50=quantile(y[ff]-x[ff],.5),finite_difference_P95=quantile(y[ff]-x[ff],.95));ratio=ff&(x>0);r['finite_ratio_P50']=quantile(y[ratio]/x[ratio],.5);r['ratio_n']=int(np.sum(ratio));pairs.append(r)
 csvwrite(O/'PAIRED_COMPARISONS.csv',pairs)
 # Complete 48000 status inventory, including nonexecuted cases, without fabricated errors.
 csvwrite(O/'CONFIGURATION_EXECUTION_STATUS.csv',[dict(configuration_id=index(q,b,c,r,p),quality=Q,geometry=B,case_id=f'D{c+1:02d}',replicate=r,side=-1 if r<250 else 1,support=P,status='NOT_EVALUATED' if raw[index(q,b,c,r,p)]['executed']==0 else 'ACCEPTED' if raw[index(q,b,c,r,p)]['accepted'] else 'SOLVER_FAILED') for q,Q in enumerate(['Q1','Q2']) for b,B in enumerate(['B0','B3']) for c in range(12) for r in range(500) for p,P in enumerate(['P0','P1'])])
 dump(O/'MAIN_COMPLETION.json',dict(planned=48000,executed=int(raw['executed'].sum()),accepted=int(raw['accepted'].sum()),not_evaluated=48000-int(raw['executed'].sum()),statistics_rows=len(stats),paired_rows=len(pairs),scope='DUAL_NODE_FOUR_PARAMETER_CONDITIONAL_DEMO',depth='NOT_EVALUATED'),True)
