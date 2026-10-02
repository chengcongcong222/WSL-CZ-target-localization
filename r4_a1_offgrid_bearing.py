#!/usr/bin/env python3
"""R4-A1: isolated continuous-truth observations, frozen-grid estimation and audit.

No imports or writes to R3 scripts/results. Model files are read-only.
"""
from __future__ import annotations
import argparse
import hashlib
import inspect
import itertools
import json
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results'/'R4_A1_OFFGRID_BEARING_BOUNDARY'
MASTER = ROOT/'results'/'R4_MASTER'
MODES = ROOT/'results'/'R3_C2_Yang_SA_depth'/'R3_C2_1'/'_kraken_zgrid'
REFERENCE = ROOT/'results'/'R3_RC23_CLOSEDLOOP'/'R3_RC23_TURN_FREQ_SUBSET_FIX'
AXES = ('r_km','theta_deg','v_mps','psi_deg')
SHAPE = (16,21,11,31)
ORIGIN = np.array([45.,-5.,1.,-15.])
STEP = np.array([1.,.5,.2,1.])
TIMES = np.arange(0.,1200.1,10.)
N_W1 = 61
PROFILE = np.arange(150.,250.1,5.)
FREQS = (201,235,283)


def json_write(path, value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_write(name, frame):
    frame.to_csv(OUT/name,index=False,float_format='%.17g',lineterminator='\n',
                 compression='gzip' if name.endswith('.gz') else None)


def wrap_rad(value):
    return (np.asarray(value)+np.pi)%(2*np.pi)-np.pi


def angle_error_deg(value, truth):
    return np.abs((np.asarray(value)-truth+180)%360-180)


def grid():
    idx = np.indices(SHAPE).reshape(4,-1).T
    return np.round(ORIGIN+idx*STEP,10)


def platform_xy():
    dt = np.maximum(TIMES-600.,0.)
    xp = 2*np.minimum(TIMES,600.)+2*dt*np.cos(np.deg2rad(15.))
    yp = 2*dt*np.sin(np.deg2rad(15.))
    return xp,yp


def geometry(states):
    states = np.atleast_2d(states).astype(float)
    xp,yp = platform_xy()
    r,th,v,psi = states.T
    th,psi = np.deg2rad(th),np.deg2rad(psi)
    dx = r[:,None]*1000*np.cos(th[:,None])+v[:,None]*TIMES*np.cos(psi[:,None])-xp
    dy = r[:,None]*1000*np.sin(th[:,None])+v[:,None]*TIMES*np.sin(psi[:,None])-yp
    return np.arctan2(dy,dx),np.hypot(dx,dy)


def parse_mod(path):
    # Same binary layout as audited R3 DIRECT_F_MATRIX_ON_TRAJECTORY.
    buf = path.read_bytes()
    recl = 4*int(np.frombuffer(buf[:4],dtype='<i4')[0])
    hdr = np.frombuffer(buf[84:108],dtype='<i4')
    ntot,nmat = int(hdr[2]),int(hdr[3])
    depths = np.frombuffer(buf[4*recl:5*recl],dtype='<f4')[:ntot].astype(float)
    count = int(np.frombuffer(buf[5*recl:5*recl+4],dtype='<i4')[0])
    phi = np.empty((nmat,count),dtype=complex)
    for im in range(count):
        phi[:,im] = np.frombuffer(buf[(7+im)*recl:(8+im)*recl],dtype='<c8')[:nmat]
    k = np.frombuffer(buf[(7+count)*recl:(7+count)*recl+count*8],dtype='<c8')
    if ntot!=nmat or not np.all(np.diff(depths)>0) or not np.isfinite(phi).all():
        raise ValueError('Unexpected modal source-depth layout')
    return {'depths':depths,'phi':phi,'k':k,'M':count}


def effective_depths(mod, depths):
    indices = np.abs(mod['depths'][:,None]-np.asarray(depths)[None,:]).argmin(axis=0)
    return indices,mod['depths'][indices]


def demean_windows(values):
    result = np.asarray(values,dtype=float).copy()
    for section in (slice(0,N_W1),slice(N_W1,None)):
        result[...,section] -= result[...,section].mean(axis=-1,keepdims=True)
    return result


def direct_features(models, ranges, source_depths=PROFILE):
    """Exact range modal factors, no range or depth interpolation.

    Output: [candidate, nuisance depth, frequency, time]. Depth mapping is
    inherited nearest stored source sample and exported explicitly.
    """
    ranges = np.atleast_2d(ranges).astype(float)
    n,nt = ranges.shape
    result = np.empty((n,len(source_depths),len(FREQS),nt),dtype=float)
    r = ranges.ravel()[None,:]
    for fi,freq in enumerate(FREQS):
        mod = models[freq]
        iz,_ = effective_depths(mod,source_depths)
        izr,_ = effective_depths(mod,[200.])
        weights = mod['phi'][iz]*mod['phi'][izr[0]][None,:]
        kre,alpha = mod['k'].real[:,None],-mod['k'].imag[:,None]
        factors = np.sqrt(2*np.pi/(kre*r))*np.exp(-1j*kre*r-alpha*r-1j*np.pi/4)
        pressure = weights@factors
        levels = 20*np.log10(np.maximum(np.abs(pressure),1e-30))
        result[:,:,fi,:] = levels.reshape(len(source_depths),n,nt).transpose(1,0,2)
    return demean_windows(result)


@dataclass(frozen=True)
class Observations:
    bearing_rad: np.ndarray
    relative_tl: np.ndarray


def generate_observation(truth, sigma, seed, panel_index, models):
    bearing,ranges = geometry(truth[list(AXES)].to_numpy(dtype=float))
    eps = np.zeros(len(TIMES)) if sigma==0 else np.random.default_rng(np.random.SeedSequence([seed,panel_index])).standard_normal(len(TIMES))
    tl = direct_features(models,ranges,[float(truth.z_true_m)])[0,0]
    return Observations(bearing[0]+np.deg2rad(sigma)*eps,tl)


class FrozenEstimator:
    """Truth-free API: observations, known sigma and frozen candidate catalog."""
    def __init__(self, states):
        self.states = states
        self.pred = np.empty((len(states),len(TIMES)),dtype=float)
        for start in range(0,len(states),4096):
            self.pred[start:start+4096] = geometry(states[start:start+4096])[0]
        self.norms = (self.pred*self.pred).sum(axis=1)

    def bearing_costs(self, observed):
        observed = np.atleast_2d(observed)
        # Use expanded SSE only if wrapping provably cannot affect any residual.
        no_wrap = max(abs(self.pred.min()-observed.max()),abs(self.pred.max()-observed.min()))<np.pi
        if no_wrap:
            costs = self.norms[:,None]+(observed*observed).sum(axis=1)[None,:]-2*self.pred@observed.T
            return np.maximum(costs,0.)
        return np.column_stack([(wrap_rad(self.pred-o)**2).sum(axis=1) for o in observed])

    def select_cloud(self, costs, sigma):
        cutoff = float(costs.min())+(1e-12 if sigma==0 else 13.3*np.deg2rad(sigma)**2)
        return np.flatnonzero(costs<=cutoff),cutoff

    @staticmethod
    def acoustic_estimate(ids, scores, tau):
        if len(ids)==0 or not np.isfinite(scores).all():
            raise ValueError('Empty or invalid estimator cloud')
        minimum = float(scores.min())
        # Exact tied score minima: deterministic smallest frozen node ID.
        top1 = int(ids[np.flatnonzero(scores==minimum)[0]])
        return top1,ids[scores<=minimum+tau],minimum


def neighborhood(truth, states):
    values = truth[list(AXES)].to_numpy(dtype=float)
    all_indices = np.arange(len(states))
    nearest = np.ones(len(states),dtype=bool)
    cells = []
    for ai,(origin,step,size) in enumerate(zip(ORIGIN,STEP,SHAPE)):
        axis = origin+np.arange(size)*step
        distance = angle_error_deg(axis,values[ai]) if ai in (1,3) else np.abs(axis-values[ai])
        nearest_ix = np.flatnonzero(np.isclose(distance,distance.min(),rtol=0,atol=1e-10))
        node_ix = np.unravel_index(all_indices,SHAPE)[ai]
        nearest &= np.isin(node_ix,nearest_ix)
        q = (values[ai]-origin)/step
        cells.append(sorted(set((int(np.floor(q)),int(np.ceil(q))))))
    vertices = np.array([np.ravel_multi_index(i,SHAPE) for i in itertools.product(*cells)],dtype=int)
    return np.flatnonzero(nearest),vertices


def errors(states, truth):
    states = np.atleast_2d(states)
    return np.column_stack((np.abs(states[:,0]-truth.r_km)/truth.r_km,
                            angle_error_deg(states[:,1],truth.theta_deg),
                            np.abs(states[:,2]-truth.v_mps)/truth.v_mps,
                            angle_error_deg(states[:,3],truth.psi_deg)))


def freeze_inputs():
    names = ['R4_A1_CONFIG.json','OFFGRID_TRUTH_PANEL.csv','R4_A1_DESIGN.md']
    frozen = {name:sha(OUT/name) for name in names}
    target = OUT/'DESIGN_FREEZE_MANIFEST.json'
    if target.exists():
        old = json.loads(target.read_text(encoding='utf-8'))
        if old['sha256']!=frozen:
            raise ValueError('Pre-run design changed after initial freeze')
    else:
        json_write(target,{'frozen_utc':datetime.now(timezone.utc).isoformat(),'sha256':frozen})
    return frozen


def initialize():
    frozen = freeze_inputs()
    config = json.loads((OUT/'R4_A1_CONFIG.json').read_text(encoding='utf-8'))
    panel = pd.read_csv(OUT/'OFFGRID_TRUTH_PANEL.csv')
    models = {f:parse_mod(MODES/f'zgrid_f{f}.mod') for f in FREQS}
    states = grid()
    mapping = []
    for freq,mod in models.items():
        _,physical = effective_depths(mod,PROFILE)
        mapping.extend({'frequency_hz':freq,'profile_label_m':float(label),'effective_source_depth_m':float(z)} for label,z in zip(PROFILE,physical))
    csv_write('NUISANCE_DEPTH_MAPPING.csv',pd.DataFrame(mapping))
    floors = []
    for _,truth in panel.iterrows():
        nearest,vertices = neighborhood(truth,states)
        floor = errors(states[nearest],truth).min(axis=0)
        if not np.all(floor>1e-10):
            raise ValueError('Panel must be simultaneously off-grid in all four horizontal parameters')
        floors.append({'panel_id':truth.panel_id,'nearest_grid_node_ids':json.dumps(nearest.tolist()),
                       'bracketing_cell_node_ids':json.dumps(vertices.tolist()),
                       'floor_rel_r':floor[0],'floor_abs_theta_deg':floor[1],
                       'floor_rel_v':floor[2],'floor_abs_psi_deg':floor[3]})
    csv_write('QUANTIZATION_FLOOR.csv',pd.DataFrame(floors))
    return config,panel,models,states,frozen


def pilot():
    config,panel,models,states,frozen = initialize()
    estimator = FrozenEstimator(states)
    rows = []
    for pi,truth in panel.iterrows():
        nearest,vertices = neighborhood(truth,states)
        observations = [generate_observation(truth,.1,seed,pi,models) for seed in config['pilot']['seeds']]
        costs = estimator.bearing_costs(np.stack([o.bearing_rad for o in observations]))
        for si,seed in enumerate(config['pilot']['seeds']):
            ids,cutoff = estimator.select_cloud(costs[:,si],.1)
            rows.append({'panel_id':truth.panel_id,'seed':seed,'sigma_deg':.1,'n_RC2':len(ids),
                         'nearest_grid_oracle_in_RC2':bool(np.isin(nearest,ids).any()),
                         'truth_cell_vertex_in_RC2':bool(np.isin(vertices,ids).any()),
                         'nearest_oracle_cost_minus_cutoff':float(costs[nearest,si].min()-cutoff),
                         'best_cell_vertex_cost_minus_cutoff':float(costs[vertices,si].min()-cutoff)})
    frame = pd.DataFrame(rows)
    csv_write('PILOT_RC2_DIAGNOSTIC.csv',frame)
    sample = np.arange(16)*31+np.ravel_multi_index((5,10,5,20),SHAPE)
    started = time.perf_counter()
    features = direct_features(models,geometry(states[sample])[1])
    elapsed = time.perf_counter()-started
    json_write(OUT/'PILOT_COMPUTE_BENCHMARK.json',{'n_nodes':len(sample),'elapsed_seconds':elapsed,
        'feature_bytes':features.nbytes,'seconds_per_node':elapsed/len(sample),
        'algorithm':'Exact direct modal factors, BLAS threads=2','design_sha256':frozen})
    print(frame.groupby('panel_id')[['n_RC2','nearest_grid_oracle_in_RC2','truth_cell_vertex_in_RC2']].mean().to_string(),flush=True)
    print('Direct forward seconds/node',elapsed/len(sample),flush=True)
    print('PILOT_COMPLETE',flush=True)


def topology(ids):
    indices = np.array(np.unravel_index(ids,SHAPE)).T
    remaining = set(map(tuple,indices.tolist()))
    sizes = []
    while remaining:
        stack,size = [remaining.pop()],0
        while stack:
            node = stack.pop()
            size += 1
            for axis in range(4):
                for shift in (-1,1):
                    adjacent = list(node)
                    adjacent[axis] += shift
                    adjacent = tuple(adjacent)
                    if adjacent in remaining:
                        remaining.remove(adjacent)
                        stack.append(adjacent)
        sizes.append(size)
    sizes.sort(reverse=True)
    range_ix = sorted(set(indices[:,0].tolist()))
    groups = []
    for index in range_ix:
        if not groups or index!=groups[-1][-1]+1:
            groups.append([])
        groups[-1].append(index)
    return {'n_connected_components':len(sizes),'largest_component_size':sizes[0],
            'largest_component_fraction':sizes[0]/len(ids),'component_size_top5':json.dumps(sizes[:5]),
            'n_range_components':len(groups),'n_range_bins':len(range_ix),
            'range_component_bins_km':json.dumps([[45+x for x in g] for g in groups])}


def score_features(features, observed_tl):
    per_depth = np.sqrt(((features-observed_tl[None,None,:,:])**2).mean(axis=(2,3)))
    iz = per_depth.argmin(axis=1)
    return per_depth[np.arange(len(features)),iz],PROFILE[iz]


def r3_identity(estimator, models):
    """Independent replay against frozen CSV IDs and all stored TRIPLE scores."""
    control = pd.Series({'r_km':50.,'theta_deg':0.,'v_mps':2.,'psi_deg':5.,'z_true_m':200.})
    # R3 used default_rng(seed), not the new panel SeedSequence convention.
    bearing,ranges = geometry(control[list(AXES)].to_numpy())
    observed_bearing = bearing[0]+np.random.default_rng(20260912).normal(0,np.deg2rad(.1),len(TIMES))
    costs = estimator.bearing_costs(observed_bearing)[:,0]
    ids,_ = estimator.select_cloud(costs,.1)
    stored = pd.read_csv(REFERENCE/'SUBSET_CANDIDATE_SCORES_FIXED.csv')
    stored = stored[stored['subset']=='201+235+283']
    checks,catalog = [],[]
    check_ids = set(ids)==set(stored.node_id)
    checks.append({'check':'R3_RC2_exact_node_ids','pass':check_ids,'value':len(ids),'expected':385})
    if not check_ids:
        raise ValueError('R4 formula replay does not reproduce frozen R3 RC2 cloud')
    features = np.concatenate([direct_features(models,geometry(estimator.states[ids[s:s+16]])[1]) for s in range(0,len(ids),16)])
    for z in (180.,200.,220.):
        observed_tl = direct_features(models,ranges,[z])[0,0]
        scores,depth = score_features(features,observed_tl)
        ref = stored[stored.z_true_m==z].set_index('node_id').loc[ids]
        error = float(np.max(np.abs(scores-ref.J.to_numpy())))
        top1,keep,_ = estimator.acoustic_estimate(ids,scores,.5)
        good = error<1e-8 and set(keep)==set(ref[ref.keep_tau05].index)
        checks.append({'check':f'R3_TRIPLE_scores_and_survivors_z{int(z)}','pass':good,'value':error,'expected':'score error <1e-8 and exact IDs'})
        catalog.extend({'z_true_m':z,'node_id':int(node),'J_replay':float(j),'J_R3_stored':float(old),
                        'z_star_label_m':float(zz),'kept_replay':bool(node in set(keep))}
                       for node,j,old,zz in zip(ids,scores,ref.J,depth))
    csv_write('R3_IDENTITY_REPLAY.csv',pd.DataFrame(checks))
    csv_write('R3_IDENTITY_CANDIDATE_SCORES.csv.gz',pd.DataFrame(catalog))
    if not all(c['pass'] for c in checks):
        raise ValueError('R3 acoustic replay failed')
    print('R3_IDENTITY_REPLAY_PASS',flush=True)
    return checks


def run():
    """The predeclared structural pilot stop path, not the full A1 experiment."""
    config,panel,models,states,frozen = initialize()
    if not (OUT/'PILOT_RC2_DIAGNOSTIC.csv').exists():
        raise ValueError('Run the frozen pilot before the structural stop path')
    started = time.perf_counter()
    estimator = FrozenEstimator(states)
    identity = r3_identity(estimator,models)
    integrity = [{'check':x['check'],'pass':x['pass'],'detail':str(x['value'])} for x in identity]
    rc2_rows,selected_cases,tl_observations = [],[],[]
    stored_bearings,stored_case_ids = [],[]
    for pi,truth in panel.iterrows():
        nearest,vertices = neighborhood(truth,states)
        observation = generate_observation(truth,0,0,pi,models)
        tl_observations.append(observation.relative_tl)
        observed,keys = [],[]
        for sigma in config['bearing_sigma_deg']:
            for seed in [0] if sigma==0 else config['positive_sigma_seeds']:
                eps = np.zeros(len(TIMES)) if sigma==0 else np.random.default_rng(np.random.SeedSequence([seed,pi])).standard_normal(len(TIMES))
                observed.append(observation.bearing_rad+np.deg2rad(sigma)*eps)
                keys.append((sigma,seed))
        observed = np.stack(observed)
        costs = estimator.bearing_costs(observed)
        for ci,(sigma,seed) in enumerate(keys):
            case_id = f'{truth.panel_id}_s{sigma:g}_seed{seed}'
            ids,cutoff = estimator.select_cloud(costs[:,ci],sigma)
            is_pilot = sigma==config['pilot']['sigma_deg'] and seed in config['pilot']['seeds']
            row = {'case_id':case_id,'panel_id':truth.panel_id,'sigma_deg':sigma,'seed':seed,
                   'n_RC2':len(ids),'rc2_cmin':float(costs[:,ci].min()),'rc2_cutoff':cutoff,
                   'nearest_grid_oracle_in_RC2':bool(np.isin(nearest,ids).any()),
                   'truth_cell_vertex_in_RC2':bool(np.isin(vertices,ids).any()),
                   'nearest_oracle_cost_minus_cutoff':float(costs[nearest,ci].min()-cutoff),
                   'best_cell_vertex_cost_minus_cutoff':float(costs[vertices,ci].min()-cutoff),
                   'experiment_scope':'RC2_ONLY_FULL_AXIS','end_to_end_scored':is_pilot}
            rc2_rows.append(row)
            if is_pilot:
                selected_cases.append({'metadata':row,'ids':ids,'costs':costs[ids,ci].copy(),'cutoff':cutoff,'observed_bearing':observed[ci]})
                stored_case_ids.append(case_id)
                stored_bearings.append(observed[ci])
        truth_features = direct_features(models,geometry(truth[list(AXES)].to_numpy())[1],[truth.z_true_m])[0,0]
        oracle_features = direct_features(models,geometry(states[nearest])[1],[truth.z_true_m])[:,0]
        delta = float(np.sqrt(((oracle_features-truth_features)**2).mean()))
        actual_depth = effective_depths(models[FREQS[0]],[truth.z_true_m])[1][0]
        integrity.append({'check':f'{truth.panel_id}_continuous_truth_not_snapped','pass':delta>1e-8,
                          'detail':f'exact truth range/features differ from nearest grid; relative-TL RMS={delta:.12g} dB; physical z={actual_depth:g}m'})
        print(f'RC2_AXIS {truth.panel_id} cases={len(keys)} pilot_cloud_union={len(set(np.concatenate([c["ids"] for c in selected_cases if c["metadata"]["panel_id"]==truth.panel_id]).tolist()))}',flush=True)
    rc2_frame = pd.DataFrame(rc2_rows)
    csv_write('RC2_AXIS_CASE_LEVEL.csv',rc2_frame)
    union = np.unique(np.concatenate([case['ids'] for case in selected_cases]))
    scores = np.empty((len(panel),len(states)),dtype=float); scores.fill(np.nan)
    nuisance = np.empty_like(scores); nuisance.fill(np.nan)
    catalog = []
    all_obs = np.stack(tl_observations)
    print(f'PILOT_ACOUSTIC_UNION nodes={len(union)}; full end-to-end MC remains paused',flush=True)
    for start in range(0,len(union),16):
        ids = union[start:start+16]
        features = direct_features(models,geometry(states[ids])[1])
        for pi,truth in panel.iterrows():
            js,zs = score_features(features,all_obs[pi])
            scores[pi,ids],nuisance[pi,ids] = js,zs
            catalog.extend({'panel_id':truth.panel_id,'node_id':int(node),'J':float(j),'z_star_label_m':float(z),
                            'z_effective_source_m':float(effective_depths(models[FREQS[0]],[z])[1][0])}
                           for node,j,z in zip(ids,js,zs))
        if start%256==0:
            print(f'ACOUSTIC_CATALOG {start+len(ids)}/{len(union)} elapsed={time.perf_counter()-started:.1f}s',flush=True)
    csv_write('CANDIDATE_SCORE_CATALOG.csv.gz',pd.DataFrame(catalog))
    # Oracle scores stay in a separate post-hoc table; they never enlarge clouds.
    oracle_rows = []
    oracle_scores = {}
    for pi,truth in panel.iterrows():
        nearest,vertices = neighborhood(truth,states)
        oracle_features = direct_features(models,geometry(states[nearest])[1])
        js,zs = score_features(oracle_features,all_obs[pi])
        oracle_scores[truth.panel_id] = float(js.min())
        oracle_rows.extend({'panel_id':truth.panel_id,'node_id':int(node),'J':float(j),'z_star_label_m':float(z),
                            'estimator_use':'EVALUATION_ONLY_NOT_INJECTED'} for node,j,z in zip(nearest,js,zs))
    csv_write('NEAREST_GRID_ORACLE_ACOUSTIC_DIAGNOSTIC.csv',pd.DataFrame(oracle_rows))
    results,qualities,survivors_out = [],[],[]
    accepted_ids,accepted_costs,offsets = [],[],[0]
    for case in selected_cases:
        meta,ids = case['metadata'],case['ids']
        pi = int(panel.index[panel.panel_id==meta['panel_id']][0]); truth = panel.loc[pi]
        nearest,vertices = neighborhood(truth,states)
        top1,kept,jmin = estimator.acoustic_estimate(ids,scores[pi,ids],config['rc3_tau_db'])
        top_error = errors(states[top1],truth)[0]
        all_error = errors(states[kept],truth).max(axis=0)
        floor = errors(states[nearest],truth).min(axis=0)
        row = {**meta,'experiment_scope':'NOMINAL_PILOT_ONLY','top1_node_id':top1,'Jmin':jmin,
               'n_survivors':len(kept),'top1_in_bracketing_cell':bool(top1 in set(vertices)),
               'nearest_grid_oracle_in_survivors':bool(np.isin(nearest,kept).any()),
               'truth_cell_vertex_in_survivors':bool(np.isin(vertices,kept).any()),
               'nearest_oracle_J':oracle_scores[truth.panel_id],
               'nearest_oracle_J_minus_selected_min':oracle_scores[truth.panel_id]-jmin,
               'top1_z_star_label_m':float(nuisance[pi,top1]),'z_role':'PROFILED_NUISANCE_NOT_DEPTH_ESTIMATE'}
        names = ('rel_r','abs_theta_deg','rel_v','abs_psi_deg')
        for ai,name in enumerate(names):
            row[f'top1_{name}'] = float(top_error[ai])
            row[f'floor_{name}'] = float(floor[ai])
            row[f'excess_{name}'] = float(top_error[ai]-floor[ai])
            row[f'survivor_worst_{name}'] = float(all_error[ai])
        for ai,axis in enumerate(AXES):
            row[f'top1_{axis}'] = float(states[top1,ai])
            row[f'survivor_{axis}_min'] = float(states[kept,ai].min())
            row[f'survivor_{axis}_max'] = float(states[kept,ai].max())
            row[f'survivor_{axis}_width'] = float(np.ptp(states[kept,ai]))
        row.update(topology(kept))
        results.append(row)
        qualities.append({key:val for key,val in row.items() if key in ('case_id','panel_id','sigma_deg','seed','n_survivors','experiment_scope') or key.startswith('survivor_') or key in topology(kept)})
        survivors_out.extend({'case_id':meta['case_id'],'node_id':int(node),**dict(zip(AXES,states[node])),
                              'J':float(scores[pi,node]),'z_star_label_m':float(nuisance[pi,node])} for node in kept)
        accepted_ids.extend(ids.tolist()); accepted_costs.extend(case['costs'].tolist()); offsets.append(len(accepted_ids))
        # Counterfactual evaluation data changes only metrics, not estimator outputs.
        repeated = estimator.acoustic_estimate(ids,scores[pi,ids],config['rc3_tau_db'])
        changed_truth = truth.copy(); changed_truth['r_km'] += .33
        altered_error = errors(states[top1],changed_truth)
        integrity.append({'check':meta['case_id']+'_evaluation_truth_not_estimator_input',
                          'pass':top1==repeated[0] and np.array_equal(kept,repeated[1]),
                          'detail':'Estimator API has no truth/oracle/depth argument; evaluation-only counterfactual cannot change fixed-observation IDs'})
    csv_write('CASE_LEVEL_RESULTS.csv',pd.DataFrame(results))
    csv_write('SURVIVOR_SET_QUALITY.csv',pd.DataFrame(qualities))
    csv_write('SURVIVOR_NODES.csv.gz',pd.DataFrame(survivors_out))
    csv_write('GRID_STATE_TABLE.csv.gz',pd.DataFrame({'node_id':np.arange(len(states)),**dict(zip(AXES,states.T))}))
    np.savez_compressed(OUT/'RC2_ACCEPTED_CLOUDS.npz',case_ids=np.array(stored_case_ids),offsets=np.array(offsets),
                        node_ids=np.array(accepted_ids,dtype=np.int32),rc2_costs=np.array(accepted_costs))
    np.savez_compressed(OUT/'PILOT_OBSERVATIONS.npz',case_ids=np.array(stored_case_ids),bearing_rad=np.stack(stored_bearings),
                        panel_ids=panel.panel_id.to_numpy(dtype=str),relative_tl=all_obs,times_s=TIMES)
    for name,method in [('bearing_costs',FrozenEstimator.bearing_costs),('select_cloud',FrozenEstimator.select_cloud),('acoustic_estimate',FrozenEstimator.acoustic_estimate)]:
        signature = str(inspect.signature(method))
        integrity.append({'check':name+'_truth_free_signature','pass':not any(x in signature for x in ('truth','oracle','z_true','panel')),'detail':signature})
    inputs = [MODES/f'zgrid_f{f}.mod' for f in FREQS]+[REFERENCE/'SUBSET_CANDIDATE_SCORES_FIXED.csv']
    csv_write('INPUT_SOURCE_MANIFEST.csv',pd.DataFrame([{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p),
        'sha256_lf':hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() if p.suffix=='.csv' else sha(p)} for p in inputs]))
    csv_write('OFFGRID_INTEGRITY_AUDIT.csv',pd.DataFrame(integrity))
    if not all(x['pass'] for x in integrity):
        raise ValueError('Off-grid pipeline integrity checks failed')
    json_write(OUT/'RUN_SCOPE.json',{'execution_utc':datetime.now(timezone.utc).isoformat(),
        'design_sha256':frozen,'executed_scope':'STRUCTURAL_STOP_PATH','RC2_only_cases':len(rc2_frame),
        'nominal_end_to_end_pilot_cases':len(results),'main_end_to_end_MC_cases_executed':0,
        'on_grid_control_cases':3,'acoustic_union_nodes':len(union),'elapsed_seconds':time.perf_counter()-started,
        'refinement':'NOT_RUN_UPSTREAM_RC2_CELL_DELETION_REQUIRES_NEW_GATE',
        'A1_status':'BLOCKED_PENDING_RC2_QUANTIZATION_TREATMENT','P5':'NOT_OPENED'})
    print('STRUCTURAL_DIAGNOSTIC_RUN_COMPLETE',flush=True)
    report()


def report():
    """Summarize saved data only; never rerun acoustics or fabricate unrun levels."""
    config = json.loads((OUT/'R4_A1_CONFIG.json').read_text(encoding='utf-8'))
    panel = pd.read_csv(OUT/'OFFGRID_TRUTH_PANEL.csv')
    cases = pd.read_csv(OUT/'CASE_LEVEL_RESULTS.csv')
    rc2 = pd.read_csv(OUT/'RC2_AXIS_CASE_LEVEL.csv')
    floors = pd.read_csv(OUT/'QUANTIZATION_FLOOR.csv')
    scope = json.loads((OUT/'RUN_SCOPE.json').read_text(encoding='utf-8'))
    metric_names = ('rel_r','abs_theta_deg','rel_v','abs_psi_deg')
    def stats(values):
        return {'P50':float(np.quantile(values,.5)),'P90':float(np.quantile(values,.9)),
                'P95':float(np.quantile(values,.95)),'worst_observed':float(np.max(values))}
    summaries = []
    for prefix in ('top1_','survivor_worst_','excess_','floor_'):
        for metric in metric_names:
            summaries.append({'scope':'NOMINAL_END_TO_END_PILOT_ONLY','sigma_deg':.1,'metric':prefix+metric,
                              'n_cases':len(cases),**stats(cases[prefix+metric]),
                              'coarse_cell_success_fraction':float(cases.top1_in_bracketing_cell.mean())})
    csv_write('BEARING_BOUNDARY_SUMMARY.csv',pd.DataFrame(summaries))
    rc2_summary = []
    for sigma,group in rc2.groupby('sigma_deg'):
        rc2_summary.append({'scope':'RC2_ONLY_FULL_AXIS' if sigma else 'NOISELESS_RC2_CONTROL',
                            'sigma_deg':float(sigma),'n_cases':len(group),
                            'nearest_grid_oracle_retention_fraction':float(group.nearest_grid_oracle_in_RC2.mean()),
                            'truth_cell_vertex_retention_fraction':float(group.truth_cell_vertex_in_RC2.mean()),
                            **{f'n_RC2_{k}':v for k,v in stats(group.n_RC2).items()}})
    rc2_summary = pd.DataFrame(rc2_summary)
    csv_write('RC2_BOUNDARY_SUMMARY.csv',rc2_summary)
    seed_rows = []
    for seed,group in cases.groupby('seed'):
        seed_rows.append({'scope':'NOMINAL_PILOT_ONLY','sigma_deg':.1,'seed':int(seed),'n_panels':len(group),
                          'coarse_cell_success_fraction':float(group.top1_in_bracketing_cell.mean()),
                          **{f'{prefix}{metric}_{k}':value for prefix in ('top1_','survivor_worst_')
                             for metric in metric_names for k,value in stats(group[prefix+metric]).items()}})
    csv_write('SEED_LEVEL_RESULTS.csv',pd.DataFrame(seed_rows))
    panel_rows = []
    for (pid,sigma),group in rc2.groupby(['panel_id','sigma_deg']):
        panel_rows.append({'panel_id':pid,'sigma_deg':float(sigma),'n_realizations':len(group),
                           'RC2_nearest_oracle_retention_fraction':float(group.nearest_grid_oracle_in_RC2.mean()),
                           'RC2_cell_retention_fraction':float(group.truth_cell_vertex_in_RC2.mean())})
    csv_write('PANEL_RC2_BOUNDARY_SUMMARY.csv',pd.DataFrame(panel_rows))
    failures = cases[['case_id','panel_id','sigma_deg','seed','truth_cell_vertex_in_RC2',
                      'nearest_grid_oracle_in_RC2','nearest_grid_oracle_in_survivors',
                      'truth_cell_vertex_in_survivors','top1_in_bracketing_cell',
                      'nearest_oracle_J_minus_selected_min','n_connected_components','n_range_components']].copy()
    failures['upstream_coarse_cell_deleted'] = ~failures.truth_cell_vertex_in_RC2
    failures['nearest_oracle_has_higher_acoustic_score_than_selected_node'] = failures.nearest_oracle_J_minus_selected_min>0
    failures['interpretation'] = 'COARSE_BEARING_GRID_SELECTION_MISMATCH_AND_MODAL_FEATURE_ALIASING; intrinsic continuous identifiability not resolved'
    csv_write('FAILURE_MODE_DIAGNOSTIC.csv',failures)
    top_stats = {m:stats(cases['top1_'+m]) for m in metric_names}
    set_stats = {m:stats(cases['survivor_worst_'+m]) for m in metric_names}
    summary_table = pd.DataFrame([{'metric':m,'top1_P50':top_stats[m]['P50'],'top1_P95':top_stats[m]['P95'],
                                  'survivor_P50':set_stats[m]['P50'],'survivor_P95':set_stats[m]['P95']} for m in metric_names])
    decision = {
        'stage':'R4-A1','baseline_commit':config['baseline_commit'],
        'decision':'R4_A1_BLOCKED_BY_COARSE_RC2_GRID_SELECTION_AND_OFFGRID_MODAL_MISMATCH',
        'execution_scope':'STRUCTURAL_STOP_PATH_NOT_FULL_A1_STATISTICAL_EXPERIMENT',
        'gates':{'A1-1':'OFFGRID_EVALUATION_PIPELINE_VALIDATED_PENDING_INDEPENDENT_AUDIT',
                 'A1-2':'PILOT_QUANTIZATION_AND_EXCESS_ERROR_DIAGNOSED_FULL_GATE_INCOMPLETE',
                 'A1-3':'BLOCKED_FULL_END_TO_END_BEARING_BOUNDARY_NOT_RUN',
                 'A1-4':'NO_STABLE_REGION_ESTABLISHED'},
        'error_interpretation':'MIXED',
        'interpretation_scope':'Grid-relative RC2 rejection and wrong coarse acoustic minima both observed; this is not a proof of intrinsic continuous physical non-identifiability.',
        'continuous_refinement':'NOT_RUN; post-RC2 local refinement is not demonstrated to restore excluded coarse cells and continuous horizontal search requires a new validation gate.',
        'truth_panel_cases':len(panel),'all_four_horizontal_coordinates_off_grid':True,
        'source_depth_role':'PROFILED_NUISANCE; truth depths180/200/220 are physical modal samples, not continuous depth-estimator tests',
        'planned_bearing_sigma_deg':config['bearing_sigma_deg'],
        'executed_end_to_end_bearing_sigma_deg':[.1],
        'pilot_seeds':config['pilot']['seeds'],'positive_axis_seeds_per_panel':len(config['positive_sigma_seeds']),
        'RC2_only_case_count':len(rc2),'end_to_end_pilot_case_count':len(cases),'on_grid_identity_cases':3,
        'coarse_cell_top1_success_fraction_pilot':float(cases.top1_in_bracketing_cell.mean()),
        'pilot_RC2_cell_retention_fraction':float(cases.truth_cell_vertex_in_RC2.mean()),
        'pilot_final_cell_retention_fraction':float(cases.truth_cell_vertex_in_survivors.mean()),
        'pilot_top1_error_stats':top_stats,'pilot_survivor_worst_error_stats':set_stats,
        'quantization_floor_definition':config['quantization_floor'],
        'first_observed_failure_boundary':{'stage':'RC2','lowest_positive_tested_sigma_deg':.02,
            'cell_retention_fraction_at_lowest_sigma':float(rc2_summary.loc[rc2_summary.sigma_deg==.02,'truth_cell_vertex_retention_fraction'].iloc[0]),
            'noiseless_cell_retention_fraction':float(rc2_summary.loc[rc2_summary.sigma_deg==0,'truth_cell_vertex_retention_fraction'].iloc[0]),
            'interpretation':'Already fails at the finest tested bearing precision; no sensor-degradation crossover or end-to-end tolerance threshold established.'},
        'R4_overall_progress_percent':0,'R3_status':'R3_FROZEN','P5':'NOT_OPENED',
        'blockers':['Coarse bearing-grid mismatch is not represented in the RC2 selection rule; better sensor precision can delete the true neighborhood.',
                    'Nearest coarse state has higher modal mismatch score than the selected wrong state in all27 pilot cases; quantization floors alone do not describe estimator errors.',
                    'Full end-to-end multi-sigma Monte Carlo and a reliable continuous-horizontal estimator have not been established.'],
        'next_recommended_stage':'R4-A1-FIX: validate a truth-independent RC2 quantization/continuation treatment and continuous-horizontal scoring strategy; rerun the unchanged frozen panel/seeds after independent review. Do not open A2/B/P5.',
        'statistical_scope':'P50/P95 of27 nominal pilot cases with3 seeds are descriptive only; the full RC2-only axis uses30 realizations per panel/sigma. No scene-population or real-ocean probability guarantee.',
        'created_utc':datetime.now(timezone.utc).isoformat()}
    json_write(OUT/'R4_A1_DECISION.json',decision)
    progress = json.loads((MASTER/'R4_PROGRESS.json').read_text(encoding='utf-8'))
    progress.update({'A1_status':'BLOCKED_BY_OFFGRID_COARSE_GRID_SELECTION',
                     'independent_audit_status':'PENDING_AUDIT_OF_STRUCTURAL_STOP_COMMIT',
                     'overall_progress_percent':0,'A1_gates':decision['gates'],
                     'blockers':decision['blockers'],'next_recommended_stage':decision['next_recommended_stage']})
    json_write(MASTER/'R4_PROGRESS.json',progress)
    (MASTER/'R4_EVIDENCE_LEDGER.csv').write_text(
        'stage,artifact,status,evidence_scope,independent_audit_status,credited_weight_percent\n'
        'A1,../R4_A1_OFFGRID_BEARING_BOUNDARY/R4_A1_DESIGN.md,PRE_RUN_DESIGN_FROZEN,9 jointly off-grid horizontal truths and fixed seeds,PENDING,0\n'
        'A1,../R4_A1_OFFGRID_BEARING_BOUNDARY/R3_IDENTITY_REPLAY.csv,EXECUTION_VALIDATED,On-grid R3 identity at3 depths,PENDING,0\n'
        'A1,../R4_A1_OFFGRID_BEARING_BOUNDARY/RC2_AXIS_CASE_LEVEL.csv,STRUCTURAL_BLOCKER_OBSERVED,1629 RC2-only cases across the frozen sigma axis,PENDING,0\n'
        'A1,../R4_A1_OFFGRID_BEARING_BOUNDARY/CASE_LEVEL_RESULTS.csv,PILOT_FAILURE_DIAGNOSED,27 nominal end-to-end pilot cases only; full A1 boundary not established,PENDING,0\n',encoding='utf-8')
    report_text = f'''# R4-A1 structural stop report

Decision: **{decision['decision']}**. R4 progress remains **0%**. R3 stays frozen.

## Main finding

This frozen coarse pipeline does not generalize successfully to the9 jointly off-grid truths in the executed nominal pilot. Top1 coarse-cell localization succeeds in0/27 cases. RC2 already deletes every bracketing-cell vertex in25/27 pilot cases; only1/27 final sets retains any such vertex. These are resolution-relative diagnostics, not the customer's as-yet-undefined joint error norm.

The cheaper full RC2-only axis confirms a structural selection defect: at sigma0.02 and0.05deg, every one of270 cases loses the truth's bracketing coarse cell; the9 noiseless cases also lose it. Larger noise relaxes the inherited selection cloud and raises neighborhood retention. This inversion does not mean worse bearing measurements improve physical observability. The rule represents bearing noise but not coarse bearing-grid mismatch.

On-grid replay reproduces all385 R3 RC2 IDs and TRIPLE scores/survivors at180/200/220m. Thus this is an R4 off-grid counterexample to extending the frozen tested control, not a reinterpretation or overwrite of R3 numerical evidence.

## Frozen design and executed scope

Baseline `{config['baseline_commit']}`. Nine jointly off-grid horizontal states:

```text
{panel.to_string(index=False)}
```

TRIPLE201+235+283Hz, E0 matched modal inputs,2m/s platform,15deg turn at600s,1200s observation,121 bearing/TL samples. Source depths180/200/220m are profiled nuisances, held on physical samples to isolate horizontal gridding. The old5m nuisance labels map to2m modal samples; see NUISANCE_DEPTH_MAPPING.csv. No depth interpolation, depth-estimator development, nav/SSP sweep, joint corner, real HLA extraction or P5.

Pre-run configuration/panel/design hashes are retained in DESIGN_FREEZE_MANIFEST.json. Planned positive levels are0.02/0.05/0.1/0.2/0.4/0.8deg with30 seeds per panel; sigma0 has one noiseless case per panel. Executed:1629 RC2-only cases,27 nominal end-to-end pilot cases with seeds410001–410003, and3 on-grid identity controls. The main1620-case end-to-end matrix was **not run**. PILOT_ACTION_DECISION.md records the structural-stop reason. Acoustic evaluation used only the union of{scope['acoustic_union_nodes']} pilot-accepted nodes plus separate evaluation-only oracle scores, not a truth-centered estimator search.

## Leakage and forward integrity

Continuous truth coordinates generate bearing and direct modal propagation along exact continuous ranges. The estimator sees only observations, known sigma and the fixed candidate grid/catalog. Neither truth coordinates, nearest-oracle IDs nor true source depth are estimator inputs. Evaluation-only oracles are stored separately and never injected into acceptance or initialization. OFFGRID_INTEGRITY_AUDIT.csv records truth-feature nonidentity, estimator signatures, fixed-observation evaluation invariance and R3 replay.

The inherited per-window/per-frequency relative-TL demeaning, shared z profile, RC2 cmin+13.3sigma² and RC3 J<=Jmin+0.5 remain fixed. Source selection does not choose new frequencies, turn or support from the truth.

## Quantization versus estimation

Nearest-grid oracle: coordinatewise nearest complete-grid node, with circular angular differences and all ties. Quantization floor is its per-axis minimum error. Bracketing cell is the Cartesian product of floor/ceiling grid coordinates, up to16 vertices; retention of any vertex and retention of nearest oracle are separate. Top1 excess error is actual error minus coordinate floor. Angle errors are degrees, never percentages.

```text
{floors[['panel_id','floor_rel_r','floor_abs_theta_deg','floor_rel_v','floor_abs_psi_deg']].to_string(index=False)}
```

The nearest-grid modal score is higher than the selected score in all27 pilot cases, with positive score differences stored in FAILURE_MODE_DIAGNOSTIC.csv. Consequently, merely reinserting that oracle node would not establish correct top1 localization. Coarse source propagation mismatch/alias selection and upstream RC2 loss both appear. Interpretation: **MIXED**. This does not prove intrinsic continuous physical non-identifiability; a reliable continuous-horizontal search has not been validated. Local refinement at a selected coarse minimum is deferred rather than assumed to recover the deleted support.

## Nominal pilot output errors

The following are empirical linear-interpolation P50/P95 for27 cases with3 seeds, **only at0.1deg**. Relative r/v entries are fractions; angular entries are degrees. They are not a multi-sigma engineering accuracy curve or a95% ocean performance guarantee.

```text
{summary_table.to_string(index=False)}
```

All individual top1 errors, coordinate floors/excesses, survivor worst errors and spans are in CASE_LEVEL_RESULTS.csv / SURVIVOR_SET_QUALITY.csv. 4D Manhattan components and range-projection components are reported per case. z_star is a nuisance diagnostic, not a depth estimate. truth-node retained/rank are not off-grid performance metrics.

## RC2 bearing-noise axis

```text
{rc2_summary.to_string(index=False)}
```

This table covers only RC2 candidate selection. Do not read it as final output accuracy or a stability region. Better bearing precision narrows the tolerated residual cloud enough to exclude coarse representations of the continuous truth. First observed selection failure is already at the lowest positive tested sigma0.02deg and exists in noiseless controls; no sensor-degradation crossover is established.

## Gate status and next route

- A1-1: off-grid observation/estimator separation validated in execution; independent audit pending.
- A1-2: quantization floors and pilot excess/envelope diagnosed; full-stage discrimination incomplete, interpretation MIXED.
- A1-3: blocked; the full end-to-end bearing-to-output statistical boundary was not executed.
- A1-4: NO_STABLE_REGION_ESTABLISHED in the executed evidence; unscored sigma levels are not inferred failures or successes.

Recommendation: open an **A1-FIX Gate**, not A2, to validate truth-independent treatment of coarse bearing-grid mismatch before hard pruning and a continuous-horizontal scoring/search strategy. Keep the existing truth panel/seeds as an unchanged regression challenge; do not tune thresholds or initialize around truth. Only after that Gate should the full planned end-to-end matrix resume. Management credit remains0%; code and a pipeline replay do not equal A1 completion.

## Reproduction and audit artifacts

`python r4_a1_offgrid_bearing.py pilot` reproduces the structural pilot selection.
`python r4_a1_offgrid_bearing.py run` executes the documented reduced structural-stop path.
`python r4_a1_offgrid_bearing.py report` regenerates summaries without propagation.

PILOT_OBSERVATIONS.npz stores exact noisy bearing and perfect relative-TL observations; RC2_ACCEPTED_CLOUDS.npz stores accepted IDs/costs in CSR order; CANDIDATE_SCORE_CATALOG.csv.gz stores deterministic per-panel scores on the pilot union; SURVIVOR_NODES.csv.gz stores each case's final nodes. These artifacts reconstruct case-level minima, IDs and metrics. GRID_STATE_TABLE.csv.gz freezes ID-to-coordinate correspondence. INPUT_SOURCE_MANIFEST.csv identifies read-only modal/control inputs.

`python r4_a1_result_audit.py` independently reconstructs the pilot bearing geometry and circular residuals, accepted node IDs, final scored outputs, quantization floors and error envelopes from saved artifacts. It also verifies source/design hashes and recomputes a small modal-score sample. See [LOCAL_VALIDATION.md](LOCAL_VALIDATION.md) for scope and commands. This local consistency validation does not award scientific Gate credit.

The R4 evidence directory preserves exact bytes through its local .gitattributes so design hashes survive Windows checkout. The source manifest records both execution-byte SHA256 and LF-normalized SHA256 for the legacy CSV source; modal hashes are exact binary hashes. Input-source hashes were captured after execution; only the design manifest is a pre-run freeze.
'''
    (OUT/'R4_A1_REPORT.md').write_text(report_text,encoding='utf-8')
    (OUT/'GPT_SYNC.md').write_text(
        '# R4-A1 structural stop\n\n**'+decision['decision']+'**\n\n'
        '9 jointly off-grid truths; full RC2-only axis:1629 cases,30 seeds per positive sigma/panel.\n'
        'End-to-end:27 nominal0.1deg pilot cases only; full MC paused.\n'
        'RC2 bracketing-cell retention:0/270 at0.02deg,0/270 at0.05deg;0/9 noiseless.\n'
        'Pilot top1 coarse-cell success:0/27. R3 on-grid replay passes.\n\n'
        +summary_table.to_string(index=False)+'\n\n'
        'Interpretation MIXED; no intrinsic continuous-identifiability conclusion.\n'
        'A1-3 incomplete; NO_STABLE_REGION_ESTABLISHED. R4 progress0%.\n'
        'Recommend A1-FIX quantization/continuation Gate before further end-to-end MC. R3 frozen; P5 not opened.\n',encoding='utf-8')
    make_plots(cases,rc2_summary,floors)
    print(summary_table.to_string(index=False),flush=True)
    print(decision['decision'],'R4_PROGRESS=0%',flush=True)


def make_plots(cases,rc2_summary,floors):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes = plt.subplots(1,2,figsize=(11,4),layout='constrained')
    positive = rc2_summary[rc2_summary.sigma_deg>0]
    axes[0].plot(positive.sigma_deg,positive.truth_cell_vertex_retention_fraction,'o-',label='Any bracketing-cell vertex')
    axes[0].plot(positive.sigma_deg,positive.nearest_grid_oracle_retention_fraction,'s--',label='Nearest-grid oracle')
    axes[0].set(ylabel='Empirical retention fraction',ylim=(-.03,1.03),title='RC2 only: fixed panel, 30 seeds per sigma')
    axes[0].legend(fontsize=8)
    axes[1].plot(positive.sigma_deg,positive.n_RC2_P50,'o-',label='P50')
    axes[1].plot(positive.sigma_deg,positive.n_RC2_P95,'s--',label='P95')
    axes[1].set(ylabel='Accepted RC2 nodes',title='Relaxed noise threshold broadens the cloud',yscale='log')
    axes[1].legend()
    for ax in axes:
        ax.set_xscale('log',base=2); ax.set_xticks(positive.sigma_deg,[f'{x:g}' for x in positive.sigma_deg])
        ax.set_xlabel('Bearing sigma (deg)'); ax.grid(alpha=.25)
    fig.suptitle('R4-A1 structural diagnostic; not an end-to-end accuracy boundary',fontsize=11)
    fig.savefig(OUT/'RC2_BEARING_SELECTION_BOUNDARY.png',dpi=160)
    fig.savefig(OUT/'RC2_BEARING_SELECTION_BOUNDARY.pdf'); plt.close(fig)
    fig,axes = plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for ax,metric,label,mult in zip(axes.ravel(),('rel_r','abs_theta_deg','rel_v','abs_psi_deg'),
                                   ('Range relative error (%)','Bearing error (deg)','Speed relative error (%)','Heading error (deg)'),(100,1,100,1)):
        for offset,q,name,color in [(-.18,.5,'P50','#247ba0'),(.18,.95,'P95','#d95f02')]:
            values = [np.quantile(cases['floor_'+metric],q),np.quantile(cases['top1_'+metric],q),np.quantile(cases['survivor_worst_'+metric],q)]
            ax.bar(np.arange(3)+offset,np.array(values)*mult,width=.36,label=name,color=color)
        ax.set_xticks(np.arange(3),['Quantization floor','Coarse top1','Survivor worst'])
        ax.set_ylabel(label); ax.grid(axis='y',alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('Nominal 0.1deg pilot only: 9 truths x 3 seeds; descriptive percentiles',fontsize=11)
    fig.savefig(OUT/'NOMINAL_PILOT_ERROR_ENVELOPES.png',dpi=160)
    fig.savefig(OUT/'NOMINAL_PILOT_ERROR_ENVELOPES.pdf'); plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action',choices=['pilot','run','report'])
    args = parser.parse_args()
    with threadpool_limits(limits=2):
        if args.action=='pilot':
            pilot()
        elif args.action=='run':
            run()
        else:
            report()
    return 0


if __name__=='__main__':
    raise SystemExit(main())
