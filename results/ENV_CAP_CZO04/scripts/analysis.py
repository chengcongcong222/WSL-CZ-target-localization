"""One formal deterministic information/candidate study; no MC or receiver extraction."""
from model import *
import time

def main():
 start=time.perf_counter();cfg=load(O/'DESIGN_FREEZE.json');g=grid();ss=states();r,high=curves(2);rl,low=curves(1);r0,coarse=curves(0)
 assert np.array_equal(r,rl) and np.array_equal(r,r0)
 rows=[]
 for p in range(9):
  for k in range(len(r)):
   rows.append(dict(profile=p,env=['JW','JA','ANNUAL'][p//3],source_z_m=[150,200,250][p%3],receiver_z_m=200,frequency_hz=200,range_m=r[k],I0_level_db=coarse[p,k] if np.isfinite(coarse[p,k]) else '',I1_level_db=low[p,k] if np.isfinite(low[p,k]) else '',I2_level_db=high[p,k] if np.isfinite(high[p,k]) else '',I2_missing=not np.isfinite(high[p,k])))
 table(O/'SINGLE_DEPTH_INPUT_LEVELS.csv',rows)
 mesh=[]
 for p in range(9):
  for lo,hi in [(10000,120000),(40000,66000),(45000,60000)]:
   keep=(r>=lo)&(r<=hi);good=keep & np.isfinite(high[p]) & np.isfinite(low[p]);d=high[p,good]-low[p,good]
   mesh.append(dict(profile=p,range_low_m=lo,range_high_m=hi,count=int(keep.sum()),finite_pair=int(good.sum()),unresolved=int((keep & ~good).sum()),centered_RMS_db=float(np.std(d)),max_absolute_db=float(abs(d).max()),median_absolute_db=float(np.median(abs(d)))))
 table(O/'FULL_RANGE_POINT_RECEIVER_NUMERICS.csv',mesh)
 ws={rho:weight(rho) for rho in [0,.9]};local=[];supports=[];resource=[];bsummary=[];control=[];sbar=np.deg2rad(.1);bcut=float(chi2.ppf(.975,121))
 for si,scene in enumerate(ss):
  h=np.array(scene['h']);tb,tr=geom(h);tb=tb[0];tr=tr[0];bg,rg=geom(g);bc=np.sum(((bg-tb)/sbar)**2,axis=1);ids=np.flatnonzero(bc<=bcut);gb=g[ids];ranges=rg[ids]
  assert np.all((rg>=40000)&(rg<=66000))
  bsummary.append(dict(scene=scene['id'],grid_candidates=len(g),bearing_retained=len(ids),range_low_km=float(gb[:,0].min()) if len(ids) else '',range_high_km=float(gb[:,0].max()) if len(ids) else '',unsigned_bearing=True,side_branches=int(np.any(gb[:,1]<0))+int(np.any(gb[:,1]>=0)),truth_on_grid=bool(scene['historical_grid']),truth_injected=False))
  pack={'bearing_ids':ids,'bearing_cost':bc[ids],'truth_range_m':tr,'truth_unsigned_bearing_rad':tb,'state':h}
  truth={2:np.array([np.interp(tr,r,y) for y in high]),1:np.array([np.interp(tr,r,y) for y in low])}
  pred={2:np.stack([np.interp(ranges,r,y) for y in high],axis=1),1:np.stack([np.interp(ranges,r,y) for y in low],axis=1)} if len(ids) else {i:np.empty((0,9,121)) for i in [1,2]}
  costs={};bestlabels={}
  for rho,wh in ws.items():
   for c in ['C0','C1','C2']:
    for lev in [2,1]:
     tx=project(truth[lev],wh,c);px=project(pred[lev],wh,c).reshape(-1,121)
     dist=np.maximum(np.sum(px*px,axis=1)[:,None]+np.sum(tx*tx,axis=1)[None,:]-2*px@tx.T,0).reshape(len(ids),9,9)
     cost=dist.min(axis=1).T;best=dist.argmin(axis=1).T
     key=f'I{lev}_rho{rho}_{c}';pack[key]=cost;pack[key+'_best_profile']=best;costs[(lev,rho,c)]=cost;bestlabels[(lev,rho,c)]=best
    for tp in range(9):
     for sigma in [.1,.25,.5]:
      cut=0 if c=='C2' else float(chi2.ppf(.975,120 if c=='C0' else 119))*sigma*sigma
      ah=costs[(2,rho,c)][tp]<=cut+1e-10;al=costs[(1,rho,c)][tp]<=cut+1e-10
      accepted=gb[ah];n=len(accepted);basewidth=float(np.ptp(gb[:,0])) if len(gb) else None
      width=float(np.ptp(accepted[:,0])) if n else None
      # Both representations are diagnostics; unresolved model disagreement is retained explicitly.
      union=gb[ah|al];inter=gb[ah&al]
      supports.append(dict(scene=scene['id'],truth_profile=tp,env=['JW','JA','ANNUAL'][tp//3],source_depth_m=[150,200,250][tp%3],condition=c,rho=rho,sigma_db=sigma,bearing_candidates=len(ids),I2_retained=n,I1_retained=int(al.sum()),disagreement=int((ah!=al).sum()),union_retained=len(union),intersection_retained=len(inter),range_low_km=float(accepted[:,0].min()) if n else '',range_high_km=float(accepted[:,0].max()) if n else '',range_width_km=width if n else '',bearing_range_width_km=basewidth if basewidth is not None else '',range_components=json.dumps(np.unique(accepted[:,0]).tolist()) if n else '[]',negative_theta=int((accepted[:,1]<0).sum()),nonnegative_theta=int((accepted[:,1]>=0).sum()),empty_semantics='FINITE_GRID_EMPTY_NOT_PHYSICAL_EXCLUSION' if not n else 'FINITE_POINTS_ONLY',continuous_domain='UNRESOLVED',depth_environment='MINIMUM_OVER_ALL_NINE_LABELS',C2_identity=bool(np.array_equal(ah,np.ones(len(ids),bool))) if c=='C2' else 'NA'))
     wrong=abs(gb[:,0]-h[0])>=1 if len(ids) else np.zeros(0,bool)
     if wrong.any() and c!='C2':
      ds=costs[(2,rho,c)][tp];kk=np.flatnonzero(wrong)[np.argmin(ds[wrong])];pr=int(bestlabels[(2,rho,c)][tp,kk]);gap=float(np.sqrt(ds[kk]/121));diff=truth[2][tp]-pred[2][kk,pr];coeff=np.linalg.lstsq(wh@V,wh@diff,rcond=None)[0]
      resource.append(dict(scene=scene['id'],truth_profile=tp,rho=rho,condition=c,closest_wrong_range_km=float(gb[kk,0]),profile_of_wrong=pr,whitened_residual_RMS_db=gap,deterministic_half_gap_RMS_db=gap/2,feature_sigma_screening_limit_db=float(np.sqrt(ds[kk]/chi2.ppf(.975,120 if c=='C0' else 119))),fitted_drift_over_1200s_db=float(2*coeff[1]),observed_radial_excursion_m=float(np.ptp(tr)),observation_time_s=1200,engineering_requirement='NOT_ESTABLISHED; model-conditional separation scale only'))
  for tp in range(9):
   for rho,wh in ws.items():
    jbs=[];jls=[];jls2=[]
    for lev,cur in [(2,high),(1,low)]:
     jb,jl=jac(h,r,cur[tp]);_,jl2=jac(h,r,cur[tp],2);jbs.append(jb);jls.append(jl);jls2.append(jl2)
    for c in ['C0','C1','C2']:
     ac=[project(j.T,wh,c).T/.25 for j in jls];aa=[project(j.T,wh,c).T/.25 for j in jls2];bj=jbs[0]/sbar
     ii=[info(np.vstack([bj,j])) for j in ac];base=info(bj);den=max(np.linalg.norm(ac[0]),1e-30);jd=float(np.linalg.norm(ac[0]-ac[1])/den);st=float(np.linalg.norm(ac[0]-aa[0])/den)
     er=abs(ii[0]['range_effective_per_km2']-ii[1]['range_effective_per_km2'])/max(ii[0]['range_effective_per_km2'],1e-30)
     delta=truth[2][tp]-truth[1][tp];rawp=project(delta,wh,c);nr=float(np.linalg.norm(rawp)/np.sqrt(121))
     regular=bool(np.min(abs(tb))>1e-8)
     usable=regular and jd<=.1 and er<=.2 and st<=.02 and nr<=.25/4
     local.append(dict(scene=scene['id'],truth_profile=tp,rho=rho,condition=c,sigma_db=.25,bearing_range_info=base['range_effective_per_km2'],I2_joint_range_info=ii[0]['range_effective_per_km2'],I1_joint_range_info=ii[1]['range_effective_per_km2'],relative_J_grid_difference=jd,relative_effective_info_grid_difference=er,relative_step_difference=st,whitened_template_grid_RMS_db=nr,numerically_resolved=usable,bearing_regular=regular,I2_rank=ii[0]['rank'],I2_condition=ii[0]['condition'] if ii[0]['condition'] else '',weak_direction=json.dumps(ii[0]['weak_direction']),Fisher=json.dumps((np.vstack([bj,ac[0]]).T@np.vstack([bj,ac[0]])).tolist()),local_scope='FIXED_DISCRETE_PROFILE_CONDITIONAL_UPPER_BOUND',continuous_depth_information='NOT_EVALUATED'))
     raw=jls[0].T@wh.T@wh@jls[0];pp=ac[0]*.25;dpi=float(np.linalg.eigvalsh(raw-pp.T@pp).min())
     control.append(dict(scene=scene['id'],truth_profile=tp,rho=rho,condition=c,DPI_min_eigenvalue=dpi,DPI_PASS=dpi>=-1e-7,C2_exact_zero=bool(np.array_equal(ac[0],np.zeros_like(ac[0]))) if c=='C2' else True))
  np.savez_compressed(W/(scene['id']+'_SCORES.npz'),**pack)
  print(scene['id'],'bearing',len(ids),flush=True)
 table(O/'BEARING_FINITE_GRID_BASELINE.csv',bsummary);table(O/'LOCAL_EFFECTIVE_INFORMATION.csv',local);table(O/'FINITE_CANDIDATE_SUPPORT.csv',supports);table(O/'RESOURCE_SEPARATION_SCALES.csv',resource);table(O/'ACTUAL_DATA_INFORMATION_CONTROLS.csv',control)
 assert all(x['DPI_PASS'] and x['C2_exact_zero'] for x in control)
 primary=[s for s in supports if s['rho']==0 and s['sigma_db']==.25 and s['condition']=='C1'];shr=[s for s in primary if s['I2_retained']>0 and s['bearing_range_width_km']!='' and s['range_width_km']<s['bearing_range_width_km']]
 resolved=[a for a in local if a['condition']=='C1' and a['rho']==0 and a['numerically_resolved']]
 ii0={(a['scene'],a['truth_profile'],a['rho']):a for a in local if a['condition']=='C0'}
 erosion=[]
 for a in local:
  if a['condition']=='C1':
   b=ii0[(a['scene'],a['truth_profile'],a['rho'])];x0=max(b['I2_joint_range_info']-b['bearing_range_info'],0);x1=max(a['I2_joint_range_info']-a['bearing_range_info'],0)
   erosion.append(dict(scene=a['scene'],truth_profile=a['truth_profile'],rho=a['rho'],C0_increment=x0,C1_increment=x1,retained_fraction=x1/x0 if x0>1e-15 else '',both_numerically_resolved=a['numerically_resolved'] and b['numerically_resolved']))
 table(O/'SOURCE_DRIFT_INFORMATION_LOSS.csv',erosion)
 result=dict(stage='ENV_CAP_CZO04',analysis_completed=True,candidate_scope='114576 historical finite points; 36 fixed horizontal scenarios x 9 generating labels; no continuous exclusion certificate',formal_analyses=1,new_MC=0,new_solvers=0,new_receive_records=0,primary_C1_sigma025_iid_cases=len(primary),primary_nonempty_range_shrink_cases=len(shr),primary_local_numerically_resolved=len(resolved),primary_local_cases=324,source_drift_loss_cases=sum(x['C1_increment']<.9*x['C0_increment'] for x in erosion),C2_zero_exact=True,OBSERVABLE_EXTRACTION_NOT_ESTABLISHED=True,continuous_source_depth='NOT_EVALUATED',environment_scope='ONE_WOA_FAMILY_THREE_DISCRETE_CLIMATOLOGIES; ideal boundary no absorption',R4='0%',elapsed_s=time.perf_counter()-start)
 dump(O/'ANALYSIS_EXIT.json',dict(completed=True,**result));dump(O/'CZO04_ANALYSIS_SUMMARY.json',result)
if __name__=='__main__':main()
