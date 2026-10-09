"""One frozen HLA-H1 execution, no propagation executables."""
import argparse, json, time, os, traceback
import numpy as np, pandas as pd, psutil
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits
import hla_h1_core as c
START=time.monotonic()
LIMIT=4*3600
PEAK=0
def guard():
    global PEAK
    rss=psutil.Process().memory_info().rss;PEAK=max(PEAK,rss)
    if time.monotonic()-START>LIMIT:raise RuntimeError('HARD_WALL_BUDGET')
    if rss>8*1024**3:raise RuntimeError('HARD_MEMORY_BUDGET')
def csv(name,rows):pd.DataFrame(rows).to_csv(c.OUT/name,index=False,float_format='%.17g')
def features_cache(m,states):
    path=c.LOCAL/'FULL_GRID_M0.npy'
    if path.exists():raise RuntimeError('One execution only: cache already exists')
    a=np.lib.format.open_memmap(path,mode='w+',dtype='float64',shape=(len(states),21,357))
    for start in range(0,len(states),256):
        guard();a[start:start+256]=c.project(m.levels(states[start:start+256]),'M0')
        if start%8192==0:print('cache',start,len(states),flush=True)
    a.flush()
    c.write_json(c.OUT/'LOCAL_CACHE_MANIFEST.json',dict(path=str(path),sha256=c.sha(path),shape=list(a.shape),bytes=path.stat().st_size,input_model_only=True))
    return a
def m1_from_m0(p):
    p=np.asarray(p)
    a=p[...,:180].reshape(*p.shape[:-1],3,60)
    b=p[...,180:].reshape(*p.shape[:-1],3,59)
    return np.concatenate([np.einsum('af,...ft->...at',c.UF.T,a).reshape(*p.shape[:-1],120),np.einsum('af,...ft->...at',c.UF.T,b).reshape(*p.shape[:-1],118)],axis=-1)
def distances(p,y):
    return np.maximum(np.square(p).sum(axis=-1)[...,None]+np.square(y).sum(axis=-1)-2*p@y.T,0)
def seed_noise(seed,g,rep):
    b=np.random.default_rng(np.random.SeedSequence([seed,g,rep,0])).standard_normal(121)
    l=np.random.default_rng(np.random.SeedSequence([seed,g,rep,1])).standard_normal((3,121))
    return b,l
def metrics(meta,method,states,scores,depths,true,truth_id,truth_score,unresolved=0):
    ix=np.flatnonzero(scores<=c.CUT[method]);accepted=states[ix];n=len(ix)
    best=int(np.argmin(scores));point=states[best];err=c.errors(point,true)
    r=dict(meta,method=method,count=n,empty=n==0,truth_retained=bool(truth_id is not None and truth_id in ix),truth_joint_retained=bool(truth_score<=c.CUT[method]),best_score=float(scores[best]),best_depth=float(depths[best]),unresolved=unresolved,best_scope='GLOBAL_IF_ACCEPTED_ELSE_SAFE_POOL_ONLY')
    for j,ax in enumerate(c.AXES):
        r[f'point_{ax}']=float(point[j]);r[f'point_error_{ax}']=float(err[j])
        r[f'min_{ax}']=float(accepted[:,j].min()) if n else np.nan
        r[f'max_{ax}']=float(accepted[:,j].max()) if n else np.nan
        r[f'span_{ax}']=float(np.ptp(accepted[:,j])) if n else np.nan
        r[f'tags_{ax}']=len(np.unique(accepted[:,j])) if n else 0
        r[f'worst_error_{ax}']=float(c.errors(accepted,true)[:,j].max()) if n else np.nan
    r['point_rel_r']=float(err[0]/true[0]);r['point_rel_v']=float(err[2]/true[2])
    ranges=np.unique(accepted[:,0]) if n else np.array([])
    groups=np.split(ranges,np.flatnonzero(np.diff(ranges)>1.000001)+1) if n else []
    r['range_branches']=len(groups)
    r['range_intervals']=json.dumps([[float(x[0]),float(x[-1])] for x in groups])
    r['range_interval_total_km']=sum(float(x[-1]-x[0]) for x in groups)
    r['range_branch_gap_km']=max(np.diff(ranges),default=0.) if n else np.nan
    r['worst_rel_r']=r['worst_error_r_km']/true[0] if n else np.nan
    r['worst_rel_v']=r['worst_error_v_mps']/true[2] if n else np.nan
    return r,ix
def main_matrix(m,states,pred,cache):
    rows=[];source=c.sources();obsdir=c.OUT/'observations';canddir=c.OUT/'accepted'
    obsdir.mkdir();canddir.mkdir()
    for g,true in enumerate(c.GEOMS):
        guard()
        innovations=[seed_noise(2026100901,g,i) for i in range(32)]
        bearings=c.geometry(true)[0][0]+np.radians(.1)*np.array([i[0] for i in innovations])
        truelevel=m.levels(true,True)[0,int((c.ZS[g]-150)/5)]
        levels=[];configs=[]
        for sigma in (.1,.25,.5):
            for si in range(4):
                for rep in range(32):
                    levels.append(truelevel+source[si]+sigma*innovations[rep][1])
                    configs.append((sigma,si,rep))
        levels=np.array(levels);np.savez_compressed(obsdir/f'G{g+1:02}.npz',bearings=bearings,levels=levels,configs=np.array(configs),bearing_innovations=np.array([i[0] for i in innovations]),level_innovations=np.array([i[1] for i in innovations]))
        bc=np.column_stack([c.bearing_cost(pred,o) for o in bearings])
        ids=np.flatnonzero((bc<=max(c.CUT.values())).any(axis=1))
        p0=np.array(cache[ids]);p1=m1_from_m0(p0)
        truepos=np.flatnonzero(np.all(np.isclose(states[ids],true,atol=1e-9),axis=1))
        obs0=c.project(levels,'M0');obs1=c.project(levels,'M1')
        oracle_levels=np.array([y-source[si] for y,(_,si,_) in zip(levels,configs)])
        obso=c.project(oracle_levels,'M0')
        print('main geometry',g+1,'safe pool',len(ids),flush=True)
        # Chunk observations to limit allocation and permit budget checkpoints.
        archive={}
        for method,pp,yy in [('M0',p0,obs0),('M1',p1,obs1),('ORACLE',p0,obso)]:
            for begin in range(0,len(configs),32):
                guard();end=min(begin+32,len(configs))
                ac=distances(pp.reshape(-1,pp.shape[-1]),yy[begin:end]).reshape(len(ids),21,end-begin)
                for local,ci in enumerate(range(begin,end)):
                    sigma,si,rep=configs[ci];scores_z=bc[ids,rep,None]+ac[:,:,local]/sigma**2
                    # 0.001 dB validated-scale buffer, NOT a global certificate.
                    radius=np.sqrt(363)*.001/sigma
                    lower=np.sqrt(np.maximum(ac[:,:,local],0))/sigma
                    # Any depth whose cache score can straddle the threshold is directly rebuilt.
                    near=((np.maximum(lower-radius,0)**2+bc[ids,rep,None])<=c.CUT[method])&((lower+radius)**2+bc[ids,rep,None]>=c.CUT[method])
                    nearids=np.flatnonzero(near.any(axis=1))
                    target=oracle_levels[ci] if method=='ORACLE' else levels[ci]
                    for k in range(0,len(nearids),16):
                        ii=nearids[k:k+16];exact=c.project(m.levels(states[ids[ii]],True)-target,'M1' if method=='M1' else 'M0')
                        scores_z[ii]=bc[ids[ii],rep,None]+np.square(exact).sum(axis=-1)/sigma**2
                    zz=scores_z.argmin(axis=1);sc=scores_z[np.arange(len(ids)),zz]
                    # Verify selected best with all shared depths (no direct polishing).
                    bi=int(sc.argmin());exact=c.project(m.levels(states[ids[bi]],True)-target,'M1' if method=='M1' else 'M0')
                    exactz=bc[ids[bi],rep]+np.square(exact).sum(axis=-1)/sigma**2
                    sc[bi]=exactz.min();zz[bi]=exactz.argmin()
                    tproj=c.project(truelevel-target,'M1' if method=='M1' else 'M0')
                    ts=float(c.bearing_cost(c.geometry(true)[0][0],bearings[rep])+np.square(tproj).sum()/sigma**2)
                    meta=dict(geometry=f'G{g+1:02}',sigma=sigma,source=f'S{si}',replicate=rep,config=ci,n_pool=len(ids),direct_boundary_states=len(nearids))
                    rec,keep=metrics(meta,method,states[ids],sc,c.PROFILE[zz],true,int(truepos[0]) if len(truepos) else None,ts)
                    rec['truth_joint_score']=ts;rows.append(rec)
                    archive[f'{method}_{ci}_ids']=ids[keep].astype(np.int32)
                    archive[f'{method}_{ci}_scores']=sc[keep]
                    archive[f'{method}_{ci}_depths']=c.PROFILE[zz[keep]]
        for ci,(sigma,si,rep) in enumerate(configs):
            meta=dict(geometry=f'G{g+1:02}',sigma=sigma,source=f'S{si}',replicate=rep,config=ci,n_pool=len(states),direct_boundary_states=0)
            tid=int(np.flatnonzero(np.all(np.isclose(states,true,atol=1e-9),axis=1))[0])
            rec,keep=metrics(meta,'BEARING',states,bc[:,rep],np.full(len(states),np.nan),true,tid,float(bc[tid,rep]))
            rows.append(rec);archive[f'BEARING_{ci}_ids']=keep.astype(np.int32);archive[f'BEARING_{ci}_scores']=bc[keep,rep]
        np.savez_compressed(canddir/f'G{g+1:02}.npz',**archive)
        csv('CONFIGURATION_RESULTS.csv',rows)
        c.write_json(c.OUT/'EXECUTION_CHECKPOINT.json',dict(main_completed=len(rows),offgrid_completed=0,elapsed=time.monotonic()-START))
    return rows
def generate_offgrid(m):
    rng=np.random.default_rng(2026100902);truths=[]
    while len(truths)<8:
        s=rng.uniform([47,-1.5,1.4,-10],[59,1.5,2.8,10])
        if np.min(np.abs((s-c.ORIGIN)/c.STEP-np.round((s-c.ORIGIN)/c.STEP))*c.STEP)<1e-6:continue
        truths.append(s)
    truths=np.array(truths)
    csv('OFFGRID_TRUTH_EVALUATION_ONLY.csv',[dict(panel=i,**dict(zip(c.AXES,s)),z_label=int([180,200,220][i%3])) for i,s in enumerate(truths)])
    c.write_json(c.OUT/'OFFGRID_TRUTH_FREEZE.json',dict(sha256=c.sha(c.OUT/'OFFGRID_TRUTH_EVALUATION_ONLY.csv'),before_observation_generation=True))
    bearings=[];levels=[]
    for g,s in enumerate(truths):
        z=[180,200,220][g%3];b,l=c.geometry(s)[0][0],m.levels(s,True)[0,int((z-150)/5)]
        for rep in range(4):
            eb,el=seed_noise(2026100903,g,rep)
            bearings.append(b+np.radians(.1)*eb);levels.append(l+c.sources()[2]+.25*el)
    np.savez_compressed(c.OUT/'observations/OFFGRID.npz',bearings=np.array(bearings),levels=np.array(levels))
    return truths
def offgrid_estimator(m,states,pred,cache):
    # No truth files or oracle input. Only observations and physical candidate model.
    data=np.load(c.OUT/'observations/OFFGRID.npz');bb=data['bearings'];ll=data['levels'];n=len(bb)
    obs={'M0':c.project(ll,'M0'),'M1':c.project(ll,'M1')}
    picks={method:[[] for _ in range(n)] for method in ('BEARING','M0','M1')}
    # Exhaustive registered finite coarse seed pool, in bounded blocks.
    for begin in range(0,len(states),256):
        guard();end=min(begin+256,len(states));p0=np.array(cache[begin:end]);p1=m1_from_m0(p0)
        bc=np.column_stack([c.bearing_cost(pred[begin:end],o) for o in bb])
        for method,pp in [('BEARING',None),('M0',p0),('M1',p1)]:
            if method=='BEARING':sc=bc[:,None,:];zs=np.zeros((end-begin,n),int)
            else:
                ac=distances(pp.reshape(-1,pp.shape[-1]),obs[method]).reshape(end-begin,21,n)
                sc=bc[:,None,:]+ac/.25**2
            # Keep two h,z candidates per range, with deterministic normalized dedup.
            stratum=(np.arange(begin,end)//(21*11*31)).astype(int)
            for r in np.unique(stratum):
                ii=np.flatnonzero(stratum==r)
                for ci in range(n):
                    costs=sc[ii,:,ci].ravel();order=np.argsort(costs,kind='stable')[:4]
                    items=[(float(costs[o]),int(begin+ii[o//sc.shape[1]]),int(o%sc.shape[1]),int(r)) for o in order]
                    old=picks[method][ci]
                    keep=[v for v in old if v[3]!=r]
                    selected=[]
                    for v in sorted([v for v in old if v[3]==r]+items):
                        if all(v[1:3]!=q[1:3] for q in selected):selected.append(v)
                        if len(selected)==2:break
                    picks[method][ci]=keep+selected
        if begin%16384==0:print('offgrid seed score',begin,len(states),flush=True)
    logs=[];exports={};solutions=[]
    for method in ('BEARING','M0','M1'):
        for ci in range(n):
            guard();candidates=[]
            seeds=sorted(picks[method][ci],key=lambda v:(v[3],v[0],v[1],v[2]))
            for run,(_,node,zidx,stratum) in enumerate(seeds):
                guard();s0=states[node]
                def residual(u):
                    state=c.LOW+(c.HIGH-c.LOW)*u
                    br=c.wrap(c.geometry(state)[0][0]-bb[ci])/np.radians(.1)
                    if method=='BEARING':return br
                    level=m.levels(state)[0,zidx]
                    return np.r_[br,(c.project(level,method)-obs[method][ci])/.25]
                fit=least_squares(residual,np.clip((s0-c.LOW)/(c.HIGH-c.LOW),1e-12,1-1e-12),bounds=(np.zeros(4),np.ones(4)),max_nfev=150,ftol=1e-9,xtol=1e-9,gtol=1e-9)
                state=c.LOW+(c.HIGH-c.LOW)*fit.x
                # All-depth direct verification for every exported terminal, not only aliases.
                bcost=float(c.bearing_cost(c.geometry(state)[0][0],bb[ci]))
                if method=='BEARING':score=bcost;zlabel=np.nan;original=score;spline=score
                else:
                    exact=c.project(m.levels(state,True)[0]-ll[ci],method)
                    cost=bcost+np.square(exact).sum(axis=-1)/.25**2
                    score=float(cost.min());zlabel=float(c.PROFILE[cost.argmin()]);original=float(cost[zidx])
                    spline=float(np.square(fit.fun).sum())
                candidates.append(np.r_[state,score,zlabel,original])
                if method!='BEARING' and zlabel!=c.PROFILE[zidx]:
                    candidates.append(np.r_[state,original,c.PROFILE[zidx],original])
                logs.append(dict(method=method,config=ci,panel=ci//4,replicate=ci%4,start=run,stratum=stratum,seed_id=node,seed_depth=float(c.PROFILE[zidx]) if method!='BEARING' else np.nan,nfev=fit.nfev,success=fit.success,status=fit.status,bound_touch=bool(np.any(np.abs(state-c.LOW)<1e-7)|np.any(np.abs(state-c.HIGH)<1e-7)),spline_score=spline,original_depth_direct_score=original,reprofiled_direct_score=score))
            a=np.array(candidates);exports[f'{method}_{ci}']=a
            best=a[a[:,4].argmin()] if len(a) else np.full(7,np.inf)
            solutions.append(dict(method=method,config=ci,panel=ci//4,replicate=ci%4,**dict(zip(c.AXES,best[:4])),score=float(best[4]),depth=best[5],accepted_count=int((a[:,4]<=c.CUT[method]).sum()),failed=not np.isfinite(best[4]),starts=len(seeds)))
            csv('OFFGRID_OPTIMIZER_LOG.csv',logs);csv('OFFGRID_ESTIMATES.csv',solutions)
            np.savez_compressed(c.OUT/'OFFGRID_EXPORTED_CANDIDATES.npz',**exports)
            print('offgrid',method,ci+1,n,'best',best[4],flush=True)
    return solutions
def execute():
    c.OUT.mkdir(exist_ok=True);c.LOCAL.mkdir(exist_ok=True,parents=True)
    freeze=json.loads((c.OUT/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))
    for name,expected in freeze['sha256'].items():
        if c.sha(c.ROOT/name)!=expected:raise RuntimeError('Freeze hash mismatch '+name)
    if (c.OUT/'EXECUTION_MANIFEST.json').exists():raise RuntimeError('Execution already recorded')
    # Remote verification done by entry shell before this process; local A required.
    import subprocess
    head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0]
    if head!=remote:raise RuntimeError('A local/remote mismatch')
    c.write_json(c.OUT/'DESIGN_COMMIT.json',dict(design_sha=head,remote_main_at_start=remote))
    m=c.Model();states=c.grid();pred=c.geometry(states)[0]
    csv('MAIN_TRUTH_EVALUATION_ONLY.csv',[dict(geometry=f'G{i+1:02}',**dict(zip(c.AXES,s)),z_label=int(c.ZS[i])) for i,s in enumerate(c.GEOMS)])
    np.save(c.OUT/'GRID_STATES.npy',states)
    cache=features_cache(m,states)
    main_matrix(m,states,pred,cache)
    generate_offgrid(m)
    offgrid_estimator(m,states,pred,cache)
    c.write_json(c.OUT/'EXECUTION_MANIFEST.json',dict(status='COMPLETE',design_sha=head,main_method_configurations=9216,offgrid_method_configurations=96,elapsed_seconds=time.monotonic()-START,peak_rss_bytes=PEAK,solver_calls=0,workers=1,BLAS_threads=4))
if __name__=='__main__':
    try:
        with threadpool_limits(limits=4):execute()
    except Exception as e:
        c.write_json(c.OUT/'EXECUTION_MANIFEST.json',dict(status='PARTIAL_EXECUTION',reason=str(e),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-START,peak_rss_bytes=PEAK,solver_calls=0))
        raise
