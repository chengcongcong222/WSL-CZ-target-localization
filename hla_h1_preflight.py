"""Fixed preflight only: no scientific pilot observations."""
import time, json
import numpy as np, pandas as pd
from threadpoolctl import threadpool_limits
import hla_h1_core as c
from scipy.linalg import null_space
def run():
    c.OUT.mkdir(exist_ok=True,parents=True);c.LOCAL.mkdir(exist_ok=True,parents=True)
    t=time.perf_counter();m=c.Model();checks=[]
    def check(name,error,tol):
        checks.append(dict(check=name,error=float(error),tolerance=float(tol),PASS=bool(error<=tol)))
    rng=np.random.default_rng(2026100911)
    s=rng.uniform(c.LOW,c.HIGH,(128,4));d=rng.integers(0,21,128);ti=rng.integers(0,121,128)
    a=m.levels(s,True);b=m.levels(s)
    delta=np.abs(a-b)
    check('128_full_trajectory_cache_level_max',delta.max(),.001)
    numeric=dict(max_db=float(delta.max()),P50_db=float(np.median(delta)),P95_db=float(np.quantile(delta,.95)),minimum_level_db=float(a.min()),count=int(delta.size),zero_or_clipped=bool((a<=-600).any()))
    y=rng.normal(size=(3,121))
    changed=y.copy()
    for sl in c.SECTIONS:
        n=changed[:,sl].shape[-1];changed[:,sl]+=rng.normal(size=(3,1))+rng.normal(size=(1,n))
    check('arbitrary_source_invariance',np.max(np.abs(c.project(y,'M1')-c.project(changed,'M1'))),1e-11)
    check('M0_dimension',abs(c.project(y,'M0').size-357),0)
    check('M1_dimension',abs(c.project(y,'M1').size-238),0)
    check('nested_projection',max(0,np.linalg.norm(c.project(y,'M1'))-np.linalg.norm(c.project(y,'M0'))),1e-11)
    for wi,sl in enumerate(c.SECTIONS):
        e=y[:,sl];n=e.shape[1]
        # Explicit nuisance least squares, gauge remove last temporal column.
        X=np.concatenate([np.kron(np.eye(3),np.ones((n,1))),np.tile(np.eye(n)[:,:-1],(3,1))],axis=1)
        r=e.ravel()-X@np.linalg.lstsq(X,e.ravel(),rcond=None)[0]
        double=e-e.mean(axis=0,keepdims=True)-e.mean(axis=1,keepdims=True)+e.mean()
        check(f'explicit_LS_{wi}',np.max(np.abs(r-double.ravel())),1e-11)
        Q=np.linalg.qr(X,mode='complete')[0][:,X.shape[1]:]
        check(f'independent_QR_norm_{wi}',abs(np.square(Q.T@e.ravel()).sum()-np.square(double).sum()),1e-10)
        check(f'single_frequency_zero_{wi}',np.max(np.abs(e[:1]-e[:1].mean(axis=0,keepdims=True))),0)
    # Independent Cartesian/direct modal/projected scores at all 6+32 fixed candidates.
    ss=np.r_[c.GEOMS,rng.uniform(c.LOW,c.HIGH,(32,4))]
    independent_max=0.;score_max=0.
    levels=m.levels(ss,True)
    for i,state in enumerate(ss):
        r,th,v,p=state;tt=c.TIMES;theta=np.radians(th);psi=np.radians(p)
        px=np.where(tt<=600,2*tt,1200+2*(tt-600)*np.cos(np.radians(15)))
        py=np.where(tt<=600,0,2*(tt-600)*np.sin(np.radians(15)))
        xy=np.array([r*1000*np.cos(theta),r*1000*np.sin(theta)])[:,None]+np.array([v*np.cos(psi),v*np.sin(psi)])[:,None]*tt-np.array([px,py])
        rr=np.sqrt(np.einsum('ij,ij->j',xy,xy))
        check(f'Cartesian_{i}',np.max(np.abs(rr-c.geometry(state)[1][0])),1e-9)
        rebuilt=np.empty((21,3,121))
        for fi,f in enumerate(c.FREQS):
            mod=m.mods[f];k=mod['k'].real.astype(float);alpha=-mod['k'].imag.astype(float)
            for zi,z in enumerate(c.PROFILE):
                iz=np.argmin(np.abs(mod['depths']-z));ir=np.argmin(np.abs(mod['depths']-200))
                factors=np.sqrt(2*np.pi/(k[:,None]*rr))*np.exp(-1j*k[:,None]*rr-alpha[:,None]*rr-1j*np.pi/4)
                pfield=np.sum(mod['phi'][iz,:,None]*mod['phi'][ir,:,None]*factors,axis=0)
                rebuilt[zi,fi]=20*np.log10(np.maximum(abs(pfield),1e-30))
        independent_max=max(independent_max,float(np.max(np.abs(rebuilt-levels[i]))))
        for method in ('M0','M1'):
            e=rebuilt-y
            norm=0.
            for sl in c.SECTIONS:
                v=e[...,sl];v=v-v.mean(axis=-1,keepdims=True)
                if method=='M1':v=v-v.mean(axis=-2,keepdims=True)
                norm+=np.square(v).sum(axis=(-2,-1))
            score_max=max(score_max,float(np.max(np.abs(norm-np.square(c.project(levels[i]-y,method)).sum(axis=-1)))))
    check('independent_modal_level',independent_max,1e-8)
    check('independent_projected_SSE',score_max,1e-6)
    # Full saved S7 candidate catalog: source depths and survivor sets.
    old=pd.read_csv(c.ROOT/'results/R3_RC23_CLOSEDLOOP/R3_RC23_TURN_FREQ_SUBSET_FIX/SUBSET_CANDIDATE_SCORES_FIXED.csv')
    old=old[old['subset'].astype(str)=='201+235+283']
    grid=c.grid();hist=[]
    for z,rows in old.groupby('z_true_m'):
        ids=rows.node_id.to_numpy(int);truth=m.levels([50,0,2,5],True)[0,int((z-150)/5)]
        js=[]
        for start in range(0,len(ids),32):
            es=m.levels(grid[ids[start:start+32]],True)-truth
            v=c.project(es,'M0')
            js.extend(np.sqrt(np.square(v).sum(axis=-1).min(axis=-1)/363))
        js=np.array(js);check(f'S7_all_scores_z{z}',np.max(np.abs(js-rows.J.to_numpy())),1e-8)
        kept=ids[js<=js.min()+.5]
        check(f'S7_survivors_z{z}',0 if set(kept)=={39389,39390} else 1,0)
        hist.append(dict(z=z,count=len(ids),survivors=kept.tolist(),max_error=float(np.max(abs(js-rows.J.to_numpy())))))
    if not hist:raise RuntimeError('Missing S7 rows')
    # Cost batch random independent candidates: no performance selection.
    batch=rng.uniform(c.LOW,c.HIGH,(128,4))
    bt=time.perf_counter();features=m.levels(batch);p0=c.project(features,'M0');p1=c.project(features,'M1')
    cost=time.perf_counter()-bt
    pd.DataFrame(checks).to_csv(c.OUT/'PREFLIGHT_CHECKS.csv',index=False)
    result=dict(PASS=all(x['PASS'] for x in checks),scientific_observations_generated=False,checks=len(checks),numeric=numeric,historical=hist,depth_mapping=m.mapping,cost_batch=dict(states=128,all_depths=21,seconds=cost,seconds_per_state=cost/128),elapsed_seconds=time.perf_counter()-t)
    c.write_json(c.OUT/'PREFLIGHT.json',result)
    print(json.dumps(result,indent=2))
    if not result['PASS']:raise RuntimeError('Preflight failed')
if __name__=='__main__':
    with threadpool_limits(limits=4):run()
