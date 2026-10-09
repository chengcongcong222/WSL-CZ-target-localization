"""Complete finite-grid cold acceptance audit, without new observations or optimizers."""
import time,json
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
import hla_h1_core as c
def run():
    started=time.monotonic()
    manifest=json.loads((c.OUT/'EXECUTION_MANIFEST.json').read_text())
    if manifest['status']!='COMPLETE':raise RuntimeError('Full audit requires completed matrix')
    m=c.Model();states=np.load(c.OUT/'GRID_STATES.npy');pred=c.geometry(states)[0]
    rows=[];maximum_level=0.;tested=0
    for g in range(6):
        obs=np.load(c.OUT/f'observations/G{g+1:02}.npz');saved=np.load(c.OUT/f'accepted/G{g+1:02}.npz')
        bc=np.column_stack([np.square(c.wrap(pred-b)/np.radians(.1)).sum(axis=1) for b in obs['bearings']])
        ids=np.flatnonzero((bc<=max(c.CUT.values())).any(axis=1))
        direct=np.empty((len(ids),21,3,121))
        for begin in range(0,len(ids),32):
            if manifest['elapsed_seconds']+time.monotonic()-started>14400:raise RuntimeError('TOTAL_BUDGET_REACHED')
            direct[begin:begin+32]=m.levels(states[ids[begin:begin+32]],True)
        tested+=len(ids)*21
        # Independent redundant-coordinate centering. Frobenius norm equals
        # orthogonal-coordinate SSE, without assuming redundant entries independent.
        def centered(a,method):
            out=np.array(a,copy=True)
            for sl in (slice(0,61),slice(61,121)):
                q=out[...,sl];q-=q.mean(axis=-1,keepdims=True)
                if method=='M1':q-=q.mean(axis=-2,keepdims=True)
            return out.reshape(*out.shape[:-2],-1)
        for method in ('M0','M1','ORACLE'):
            pp=centered(direct,'M1' if method=='M1' else 'M0').reshape(-1,363)
            for begin in range(0,384,32):
                end=begin+32
                yy=obs['levels'][begin:end].copy()
                if method=='ORACLE':
                    yy-=np.array([c.sources()[int(src)] for _,src,_ in obs['configs'][begin:end]])
                y=centered(yy,'M1' if method=='M1' else 'M0')
                sse=np.maximum(np.square(pp).sum(axis=1)[:,None]+np.square(y).sum(axis=1)[None,:]-2*pp@y.T,0).reshape(len(ids),21,32)
                for j,ci in enumerate(range(begin,end)):
                    sigma,src,rep=obs['configs'][ci]
                    ss=sse[:,:,j]/sigma**2+bc[ids,int(rep),None]
                    z=ss.argmin(axis=1);score=ss[np.arange(len(ids)),z]
                    accepted=ids[score<=c.CUT[method]]
                    previous=saved[f'{method}_{ci}_ids']
                    false_reject=np.setdiff1d(accepted,previous);false_accept=np.setdiff1d(previous,accepted)
                    if len(previous):
                        indices=np.searchsorted(ids,previous);delta=float(np.max(abs(score[indices]-saved[f'{method}_{ci}_scores'])))
                    else:delta=0.
                    rows.append(dict(geometry=f'G{g+1:02}',config=ci,method=method,sigma=sigma,source=f'S{int(src)}',replicate=int(rep),safe_horizontal_points=len(ids),direct_depth_states=len(ids)*21,direct_accepted=len(accepted),saved_accepted=len(previous),false_rejected=len(false_reject),false_accepted=len(false_accept),max_saved_accepted_score_difference=delta,PASS=bool(not len(false_reject) and not len(false_accept))))
        pd.DataFrame(rows).to_csv(c.OUT/'FULL_FINITE_GRID_COLD_AUDIT.csv',index=False)
        print('full audit',g+1,'states',len(ids),'PASS',all(r['PASS'] for r in rows),flush=True)
    c.write_json(c.OUT/'FULL_FINITE_GRID_COLD_AUDIT.json',dict(PASS=all(r['PASS'] for r in rows),checks=len(rows),direct_h_depth_states=tested,FALSE_REJECTED=sum(r['false_rejected'] for r in rows),FALSE_ACCEPTED=sum(r['false_accepted'] for r in rows),max_saved_accepted_score_difference=max(r['max_saved_accepted_score_difference'] for r in rows),elapsed_seconds=time.monotonic()-started,scope='Every potentially accepted original-grid point and all 21 depths, for all 6912 acoustic method configurations; independent explicit centering and direct modal fields. Outside pool bearing lower bound exceeds maximum total threshold. Not continuous certification.'))
if __name__=='__main__':
    with threadpool_limits(limits=4):run()
