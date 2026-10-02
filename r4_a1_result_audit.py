"""Reconstruct R4 pilot outputs from stored observations and catalogs.

This local consistency audit is not the research lead's independent Gate audit.
"""
import json
import numpy as np
import pandas as pd
import r4_a1_offgrid_bearing as r4


def audit():
    out = r4.OUT
    states = pd.read_csv(out/'GRID_STATE_TABLE.csv.gz')[list(r4.AXES)].to_numpy()
    panel = pd.read_csv(out/'OFFGRID_TRUTH_PANEL.csv').set_index('panel_id')
    cases = pd.read_csv(out/'CASE_LEVEL_RESULTS.csv').set_index('case_id')
    catalog = pd.read_csv(out/'CANDIDATE_SCORE_CATALOG.csv.gz')
    survivors = pd.read_csv(out/'SURVIVOR_NODES.csv.gz')
    clouds = np.load(out/'RC2_ACCEPTED_CLOUDS.npz')
    obs = np.load(out/'PILOT_OBSERVATIONS.npz')
    records = []

    def check(name, condition, detail=''):
        records.append({'check':name,'pass':bool(condition),'detail':detail})
        if not condition:
            raise AssertionError(name+': '+detail)

    # Independent geometry, expanded-SSE-free circular residuals for all nodes.
    t = obs['times_s']
    dt = np.clip(t-600,0,None)
    px = 2*np.minimum(t,600)+2*dt*np.cos(np.pi/12)
    py = 2*dt*np.sin(np.pi/12)
    def bearings(s):
        r,theta,v,psi = s.T
        x = 1000*r[:,None]*np.cos(np.radians(theta))[:,None]+v[:,None]*t*np.cos(np.radians(psi))[:,None]-px
        y = 1000*r[:,None]*np.sin(np.radians(theta))[:,None]+v[:,None]*t*np.sin(np.radians(psi))[:,None]-py
        return np.arctan2(y,x)

    check('observation_case_order',np.array_equal(obs['case_ids'],clouds['case_ids']))
    for ci,cid in enumerate(obs['case_ids']):
        row = cases.loc[cid]
        truth = panel.loc[row.panel_id]
        noisy = bearings(truth[list(r4.AXES)].to_numpy(float)[None,:])[0]
        pi = list(obs['panel_ids']).index(row.panel_id)
        eps = np.random.default_rng(np.random.SeedSequence([int(row.seed),pi])).standard_normal(len(t))
        noisy += np.radians(row.sigma_deg)*eps
        check(cid+'_continuous_observation',np.allclose(noisy,obs['bearing_rad'][ci],rtol=0,atol=1e-14))
        costs = np.empty(len(states))
        for start in range(0,len(states),2048):
            residual = bearings(states[start:start+2048])-obs['bearing_rad'][ci]
            costs[start:start+2048] = np.square(np.arctan2(np.sin(residual),np.cos(residual))).sum(axis=1)
        cutoff = costs.min()+13.3*np.radians(row.sigma_deg)**2
        ids = np.flatnonzero(costs<=cutoff)
        lo,hi = clouds['offsets'][ci:ci+2]
        saved_ids = clouds['node_ids'][lo:hi]
        check(cid+'_RC2_exact_ids',np.array_equal(ids,saved_ids),str(len(ids)))
        check(cid+'_RC2_costs',np.allclose(costs[ids],clouds['rc2_costs'][lo:hi],rtol=0,atol=1e-12))
        js = catalog[catalog.panel_id==row.panel_id].set_index('node_id').loc[ids,'J'].to_numpy()
        minimum = js.min()
        top = int(ids[np.argmin(js)])
        kept = ids[js<=minimum+.5]
        stored_kept = survivors[survivors.case_id==cid].node_id.to_numpy()
        check(cid+'_RC3_exact_outputs',top==row.top1_node_id and np.array_equal(kept,stored_kept) and abs(minimum-row.Jmin)<1e-12)
        tv = truth[list(r4.AXES)].to_numpy(float)
        def err(s):
            e = np.abs(s-tv)
            for a in (1,3):
                delta = np.radians(s[:,a]-tv[a])
                e[:,a] = np.degrees(np.abs(np.arctan2(np.sin(delta),np.cos(delta))))
            e[:,0] /= tv[0]
            e[:,2] /= tv[2]
            return e
        floor = err(states).min(axis=0)
        top_error = err(states[[top]])[0]
        worst = err(states[kept]).max(axis=0)
        names = ('rel_r','abs_theta_deg','rel_v','abs_psi_deg')
        vals = []
        for a,name in enumerate(names):
            vals.extend([abs(row['top1_'+name]-top_error[a]),abs(row['floor_'+name]-floor[a]),
                         abs(row['excess_'+name]-(top_error[a]-floor[a])),abs(row['survivor_worst_'+name]-worst[a])])
        check(cid+'_errors_floor_excess_envelope',max(vals)<1e-10)

    frozen = json.loads((out/'DESIGN_FREEZE_MANIFEST.json').read_text(encoding='utf-8'))
    for name,digest in frozen['sha256'].items():
        check('frozen_'+name,r4.sha(out/name)==digest)
    for row in pd.read_csv(out/'INPUT_SOURCE_MANIFEST.csv').itertuples():
        path = r4.ROOT/row.path
        normalized = r4.hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest() if path.suffix=='.csv' else r4.sha(path)
        check('input_'+row.path,normalized==row.sha256_lf,'LF-normalized CSV hash; exact bytes for modal files')
    axis = pd.read_csv(out/'RC2_AXIS_CASE_LEVEL.csv')
    summary = pd.read_csv(out/'RC2_BOUNDARY_SUMMARY.csv')
    for row in summary.itertuples():
        sub = axis[axis.sigma_deg==row.sigma_deg]
        check('RC2_summary_'+str(row.sigma_deg),
              abs(sub.truth_cell_vertex_in_RC2.mean()-row.truth_cell_vertex_retention_fraction)<1e-12
              and abs(sub.nearest_grid_oracle_in_RC2.mean()-row.nearest_grid_oracle_retention_fraction)<1e-12)
    # Recompute a small forward sample, then profile RMS with independent loops.
    models = {f:r4.parse_mod(r4.MODES/f'zgrid_f{f}.mod') for f in r4.FREQS}
    for pi,pid in enumerate(obs['panel_ids']):
        truth = panel.loc[pid]
        tl = r4.direct_features(models,r4.geometry(truth[list(r4.AXES)].to_numpy(float))[1],[truth.z_true_m])[0,0]
        check(pid+'_exact_truth_TL',np.allclose(tl,obs['relative_tl'][pi],rtol=0,atol=1e-12))
    sample = catalog[catalog.panel_id==obs['panel_ids'][0]].iloc[[0,len(catalog)//18,-1]]
    features = r4.direct_features(models,r4.geometry(states[sample.node_id.to_numpy()])[1])
    for i,row in enumerate(sample.itertuples()):
        residual = features[i]-obs['relative_tl'][0]
        rms = [float(np.linalg.norm(residual[z].ravel())/np.sqrt(residual[z].size)) for z in range(21)]
        check('sample_modal_score_'+str(row.node_id),abs(min(rms)-row.J)<1e-10 and r4.PROFILE[np.argmin(rms)]==row.z_star_label_m)
    return records


if __name__=='__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=2):
        results = audit()
    pd.DataFrame(results).to_csv(r4.OUT/'RESULT_RECONSTRUCTION_AUDIT.csv',index=False)
    print(f'{len(results)} reconstruction checks passed; scientific Gate audit pending.')
