#!/usr/bin/env python3
"""Audit frozen candidate CSVs. No propagation, simulation or score generation."""
from __future__ import annotations
import hashlib
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results' / 'R3_FINAL_PARAMETER_CONTRACTION'
BASE = ROOT / 'results' / 'R3_RC23_CLOSEDLOOP'
BASELINE = 'e2debd48bf11dfbb3a16c1461b808affe56b869d'
N0 = 114576
AXES = ('r', 'theta', 'v', 'psi')
ORIGIN = np.array([45., -5., 1., -15.])
STEPS = np.array([1., .5, .2, 1.])
SHAPE = (16, 21, 11, 31)
TRUTH = np.array([50., 0., 2., 5.])
DEPTHS = (180., 200., 220.)
TAU = .5  # Existing frozen threshold; never fitted here.
SUBSETS = {'S6': '235', 'S7': '201+235+283', 'S8': '201+235+283+338'}
LABELS = {'S0': 'INITIAL_GRID', 'S1': 'RC2_W1_600S', 'S2': 'RC2_RC3_W1_MAIN235',
          'S3': 'RC2_W1W2_STRAIGHT_1200S', 'S4': 'RC2_RC3_W1W2_MAIN235',
          'S5': 'RC2_TURN15', 'S6': 'TURN15_PLUS_235',
          'S7': 'TURN15_PLUS_TRIPLE', 'S8': 'TURN15_PLUS_FOUR'}
PARENTS = {'S1': 'S0', 'S2': 'S1', 'S3': 'S0', 'S4': 'S3',
           'S5': 'S0', 'S6': 'S5', 'S7': 'S5', 'S8': 'S5'}
PENDING = 'PARAMETER_CONTRACTION_CORE_RESULT_CONFIRMED_PENDING_STRUCTURAL_AUDIT_FIX'

def js(value):
    return json.dumps(value, ensure_ascii=False, separators=(',', ':'))

def boolean(series):
    if series.dtype == bool:
        return series
    result = series.astype(str).str.lower().map({'true': True, 'false': False})
    if result.isna().any():
        raise ValueError('Invalid stored boolean')
    return result.astype(bool)

def grid_indices(frame):
    raw = (frame[list(AXES)].to_numpy(dtype=float) - ORIGIN) / STEPS
    idx = np.rint(raw).astype(int)
    if not np.allclose(raw, idx, rtol=0, atol=1e-8):
        raise ValueError('Candidate is off the frozen grid')
    if not ((idx >= 0) & (idx < np.array(SHAPE))).all():
        raise ValueError('Candidate is outside the frozen grid')
    return idx

def normalize(frame, score=None, reference=None):
    frame = frame.copy()
    aliases = {'r': ('r', 'r0_km', 'r0'), 'theta': ('theta', 'theta0_deg', 'theta0'),
               'v': ('v', 'v_mps'), 'psi': ('psi', 'psi_deg'),
               'z_true': ('z_true', 'z_true_m', 'ztrue'), 'z_star': ('z_star', 'z_star_m')}
    for dest, options in aliases.items():
        found = next((c for c in options if c in frame), None)
        if found:
            frame[dest] = frame[found]
    if 'theta' not in frame:
        if reference is None:
            raise ValueError('Missing theta and node reference')
        frame['theta'] = frame.node_id.map(reference.set_index('node_id').theta)
        frame['theta_source'] = 'NODE_ID_JOIN_TO_VERIFIED_S5_CSV'
    else:
        frame['theta_source'] = 'SOURCE_CSV'
    if frame[list(AXES)].isna().any().any():
        raise ValueError('Incomplete candidate coordinates')
    idx = grid_indices(frame)
    ids = np.ravel_multi_index(idx.T, SHAPE)
    if not np.array_equal(frame.node_id.to_numpy(dtype=float), ids):
        raise ValueError('node_id disagrees with the stored 4D state')
    frame[list(AXES)] = np.round(ORIGIN + idx * STEPS, 10)
    frame['node_id'] = ids
    frame['is_truth'] = (frame[list(AXES)].to_numpy() == TRUTH).all(axis=1)
    for flag in ('is_true', 'is_truth_state'):
        if flag in frame and not np.array_equal(boolean(frame[flag]), frame.is_truth):
            raise ValueError(f'Stored {flag} disagrees with full truth state')
    if score:
        frame['J'] = frame[score]
        if not np.isfinite(frame.J).all():
            raise ValueError('Nonfinite stored scores')
    return frame

def survivors(frame):
    if frame.node_id.duplicated().any() or int(frame.is_truth.sum()) != 1:
        raise ValueError('Scored cloud must have unique nodes and one full truth state')
    keep = frame.J <= frame.J.min() + TAU
    for flag in ('keep', 'keep_tau05'):
        if flag in frame and not np.array_equal(boolean(frame[flag]), keep):
            raise ValueError(f'Stored {flag} disagrees with frozen tau=0.5 rule')
    return frame.loc[keep].copy()

def truth_rank(scored, kept):
    """Competition rank over the entire parent scored cloud; exact serialized ties.

    Better: J < J_truth. Tied: J == J_truth, including truth. No new tolerance.
    """
    truth = scored.loc[scored.is_truth]
    if len(truth) != 1:
        raise ValueError('Missing or duplicate full truth state')
    jt = float(truth.J.iloc[0])
    better = int((scored.J < jt).sum())
    tied = int((scored.J == jt).sum())
    best = scored.loc[scored.J == scored.J.min()]
    return {'truth_retained': bool(kept.is_truth.any()), 'truth_J': jt,
            'n_scored_candidates': len(scored), 'n_strictly_better': better,
            'n_tied_with_truth': tied, 'truth_rank_min': better + 1,
            'truth_rank_max': better + tied, 'truth_rank': better + 1,
            'truth_rank_scope': 'FULL_PARENT_SCORED_CLOUD',
            'best_range': js(sorted(best.r.unique().tolist())),
            'top1_truth_exact_in_matched_synthetic_control': bool(best.is_truth.all())}

def metrics(frame):
    if frame.empty or frame.node_id.duplicated().any():
        raise ValueError('Empty or duplicate survivor cloud')
    idx, n = grid_indices(frame), len(frame)
    result = {'n_candidates': n, 'fraction_of_initial': n/N0, 'joint_grid_fraction': n/N0,
              'contraction_from_initial': 1-n/N0, 'truth_retained': bool(frame.is_truth.any())}
    bbox_product = occupied_product = 1
    for i, axis in enumerate(AXES):
        vals = frame[axis].to_numpy()
        bins = sorted(np.unique(vals).tolist())
        bbox = int(idx[:, i].max()-idx[:, i].min()+1)
        bbox_product *= bbox
        occupied_product *= len(bins)
        error = float(np.abs(vals-TRUTH[i]).max())
        result.update({f'{axis}_min': float(vals.min()), f'{axis}_max': float(vals.max()),
                       f'{axis}_width': float(np.round(vals.max()-vals.min(), 10)),
                       f'n_occupied_{axis}_bins': len(bins), f'n_bbox_{axis}_bins': bbox,
                       f'occupied_{axis}_bins': js(bins), f'max_abs_err_{axis}': error})
        if TRUTH[i] != 0:
            result[f'max_rel_err_{axis}'] = error/abs(TRUTH[i])
    result.update({'N_bbox': bbox_product, 'occupied_product': occupied_product,
                   'bounding_box_grid_fraction': bbox_product/N0,
                   'occupied_product_grid_fraction': occupied_product/N0,
                   'correlation_sparsity': n/bbox_product,
                   'occupied_product_sparsity': n/occupied_product})
    return result

def component_metrics(frame):
    remaining = set(map(tuple, grid_indices(frame).tolist()))
    sizes = []
    while remaining:
        queue, size = [remaining.pop()], 0
        while queue:
            node = queue.pop()
            size += 1
            for axis in range(4):
                for shift in (-1, 1):
                    neighbor = list(node)
                    neighbor[axis] += shift
                    neighbor = tuple(neighbor)
                    if neighbor in remaining:
                        remaining.remove(neighbor)
                        queue.append(neighbor)
        sizes.append(size)
    sizes.sort(reverse=True)
    ranges = []
    for ix in sorted(set(grid_indices(frame)[:, 0].tolist())):
        value = float(ORIGIN[0]+ix*STEPS[0])
        if not ranges or value != ranges[-1][-1]+STEPS[0]:
            ranges.append([])
        ranges[-1].append(value)
    if sum(sizes) != len(frame):
        raise ValueError('Components do not partition candidates')
    return {'n_connected_components': len(sizes), 'largest_component_size': sizes[0],
            'largest_component_fraction': sizes[0]/len(frame), 'component_size_top5': js(sizes[:5]),
            'n_range_components': len(ranges), 'surviving_range_bins': js([v for g in ranges for v in g]),
            'range_component_bins': js(ranges)}

def final_bound(rows):
    """Only S7/S8 at all three frozen depths; strict <10%, with theta separate."""
    final = rows[rows.stage.isin(['S7', 'S8'])]
    expected = set(itertools.product(('S7', 'S8'), DEPTHS))
    if set(zip(final.stage, final.z_true)) != expected or len(final) != 6:
        raise ValueError('Final diagnostic requires S7/S8 at all three depths')
    passed = bool((final[['max_rel_err_r', 'max_rel_err_v', 'max_rel_err_psi']] < .10).all().all())
    # Theta relative error is undefined because truth=0; report the absolute error.
    return 'ESTABLISHED' if passed else 'NOT_ESTABLISHED'

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat()
    decision_path = OUT/'R3_FINAL_PARAMETER_CONTRACTION_DECISION.json'
    decision_path.write_text(json.dumps({'decision': PENDING, 'baseline_commit': BASELINE, 'created_utc': now}, indent=2)+'\n', encoding='utf-8')
    manifest, checks = [], []
    def check(name, condition, detail):
        checks.append({'check': name, 'pass': bool(condition), 'detail': detail})
        if not condition:
            raise ValueError(f'Integrity check failed: {name}: {detail}')
    def read(directory, name):
        path = BASE/directory/name
        data, frame = path.read_bytes(), pd.read_csv(path)
        manifest.append({'path': path.relative_to(ROOT).as_posix(), 'rows': len(frame),
                         'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()})
        return frame
    a, b, c = ('R3_RC23_SINGLE_WINDOW_SET_CONTRACTION', 'R3_RC23_TWO_WINDOW_TEMPORAL_CONTRACTION', 'R3_RC23_TURN_FREQ_SUBSET_FIX')
    s1 = normalize(read(a, 'RC2_ACCEPTED_CANDIDATES.csv'))
    s3 = normalize(read(b, 'RC2_CUMULATIVE_1200_ACCEPTED.csv'))
    score2 = normalize(read(a, 'RC3_PROFILED_SCORE_ALL_ACCEPTED.csv'), 'J_RC3')
    score4 = normalize(read(b, 'RC3_TWO_WINDOW_PROFILED_SCORES.csv'), 'J12')
    turn = normalize(read(c, 'SUBSET_CANDIDATE_SCORES_FIXED.csv'), 'J')
    stored = normalize(read(c, 'SUBSET_SURVIVORS_TAU05.csv'), 'J')
    s5 = turn[(turn['subset']=='235') & (turn.z_true==200)].copy()
    for (subset, z), group in turn.groupby(['subset', 'z_true']):
        check(f'shared_S5_cloud_{subset}_{z}', set(group.node_id)==set(s5.node_id), 'Exact node IDs')
    idx = np.indices(SHAPE).reshape(4, -1).T
    s0 = pd.DataFrame(np.round(ORIGIN+idx*STEPS, 10), columns=AXES)
    s0['node_id'] = np.arange(len(s0))
    s0 = normalize(s0)
    check('initial_grid_size', len(s0)==N0, 'Frozen 16 x 21 x 11 x 31 grid')
    clouds = {('S0', None): s0, ('S1', None): s1, ('S3', None): s3, ('S5', None): s5}
    scored_clouds = {}
    for stage, source in [('S2', score2), ('S4', score4)]:
        for z in DEPTHS:
            scored = source[(source.config=='MAIN') & (source.z_true==z)].copy()
            check(f'{stage}_{z}_parent_cloud', set(scored.node_id)==set(clouds[(PARENTS[stage], None)].node_id), 'All parent nodes have scores')
            scored_clouds[(stage, z)], clouds[(stage, z)] = scored, survivors(scored)
    for stage, subset in SUBSETS.items():
        for z in DEPTHS:
            scored = turn[(turn['subset']==subset) & (turn.z_true==z)].copy()
            kept = survivors(scored)
            old = stored[(stored['subset']==subset) & (stored.z_true==z)]
            check(f'{stage}_{z}_survivor_ids', set(kept.node_id)==set(old.node_id) and len(kept)==len(old), 'Frozen tau=0.5 IDs agree with survivor CSV')
            cols = ['J', 'z_star']
            check(f'{stage}_{z}_survivor_values', kept.set_index('node_id').sort_index()[cols].equals(old.set_index('node_id').sort_index()[cols]), 'Scores and profiled depths agree exactly')
            scored_clouds[(stage, z)], clouds[(stage, z)] = scored, kept
    rows = []
    for (stage, z), frame in sorted(clouds.items()):
        row = {'stage': stage, 'stage_label': LABELS[stage], 'z_true': z, 'parent_stage': PARENTS.get(stage, ''), **metrics(frame)}
        if (stage, z) in scored_clouds:
            row.update(truth_rank(scored_clouds[(stage, z)], frame))
        row['contraction_from_previous'] = None
        if stage in ('S1', 'S2', 'S4', 'S6', 'S7', 'S8'):
            row['contraction_from_previous'] = 1-len(frame)/len(clouds[(PARENTS[stage], None)])
        row.update(component_metrics(frame))
        rows.append(row)
    all_rows = pd.DataFrame(rows)
    depth = all_rows[all_rows.z_true.notna()].copy()
    check('three_depth_tables_complete', len(depth)==15, 'S2/S4/S6/S7/S8 x 180/200/220m')
    check('connected_components_complete', all_rows.n_connected_components.notna().all(), 'All 19 stage/depth clouds partitioned using 4D Manhattan adjacency')
    check('bbox_corrected', bool((all_rows.N_bbox>=all_rows.occupied_product).all()), 'BBox uses integer min/max spans')
    s6 = all_rows[(all_rows.stage=='S6') & (all_rows.z_true==200)].iloc[0]
    check('S6_range_bbox_gap', s6.n_occupied_r_bins==15 and s6.n_bbox_r_bins==16, '15 occupied bins within 16-bin bounding span')
    check('truth_rank_audited', depth.truth_rank_min.notna().all(), 'All parent scores, exact serialized ties')
    edges = [('A','S0','S1','SEQUENTIAL_UPDATE'), ('A','S1','S2','SEQUENTIAL_UPDATE'),
             ('B','S0','S3','RECOMPUTED_BRANCH'), ('B','S1','S3','RECOMPUTED_BRANCH'),
             ('B','S3','S4','SEQUENTIAL_UPDATE'), ('C','S0','S5','RECOMPUTED_BRANCH'),
             ('C','S1','S5','RECOMPUTED_BRANCH')]
    edges += [('C','S5',s,'PARALLEL_CONFIGURATION') for s in SUBSETS]
    graph = []
    for branch, parent, child, relation in edges:
        for z in DEPTHS if child in ('S2','S4','S6','S7','S8') else (None,):
            p, q = set(clouds[(parent,None)].node_id), set(clouds[(child,z)].node_id)
            graph.append({'branch': branch, 'parent_stage': parent, 'child_stage': child, 'z_true': z,
                          'parent_n': len(p), 'child_n': len(q), 'child_id_subset_of_parent': q<=p,
                          'n_child_not_in_parent': len(q-p), 'count_reduction_relative_to_parent': 1-len(q)/len(p), 'relation_type': relation})
    check('branch_semantics_fixed', all(g['child_id_subset_of_parent'] for g in graph), 'Every declared edge compared using actual node IDs')
    overlaps = []
    for z in DEPTHS:
        for left, right in itertools.combinations(SUBSETS,2):
            p, q = set(clouds[(left,z)].node_id), set(clouds[(right,z)].node_id)
            overlaps.append({'z_true': z, 'left_stage': left, 'right_stage': right, 'n_intersection': len(p&q),
                             'n_union': len(p|q), 'Jaccard': len(p&q)/len(p|q), 'left_subset_of_right': p<=q,
                             'right_subset_of_left': q<=p, 'relation_type': 'PARALLEL_CONFIGURATION'})
    final_states, final_cases = [], []
    for stage in ('S7','S8'):
        for z in DEPTHS:
            scored, kept = scored_clouds[(stage,z)], clouds[(stage,z)]
            rank = truth_rank(scored, kept)
            for _, r in kept.sort_values('node_id').iterrows():
                final_states.append({'stage': stage, 'subset': SUBSETS[stage], 'z_true': z, 'node_id': int(r.node_id),
                                     **{a: float(r[a]) for a in AXES}, 'J': float(r.J), 'z_star': float(r.z_star),
                                     'is_truth': bool(r.is_truth), 'truth_rank': rank['truth_rank'],
                                     'candidate_rank_min': int((scored.J<r.J).sum())+1})
            final_cases.append({'stage': stage, 'subset': SUBSETS[stage], 'z_true': z,
                                'states': kept[list(AXES)].sort_values(list(AXES)).to_dict('records')})
    stress_rows = []
    specs = [('amplitude','R3_RC23_TURN_AMPLITUDE_ROBUSTNESS','AMPLITUDE_CANDIDATE_SCORES.csv',['config','shape','mode','A','z_true'],'A'),
             ('tracked_freq_drift','R3_RC23_TURN_TRACKED_FREQ_ROBUSTNESS_FIX2','TRACKED_FREQ_CANDIDATE_SCORES_FIXED.csv',['config','D','z_true'],'D'),
             ('ssp_mismatch','R3_RC23_TURN_SSP_MISMATCH','SSP_MISMATCH_CANDIDATE_SCORES.csv',['config','env_truth','z_true'],'env_truth')]
    for kind, directory, name, keys, levelcol in specs:
        source = normalize(read(directory,name), 'J', s5)
        for key, scored in source.groupby(keys,sort=True):
            case = dict(zip(keys,key))
            check(f'stress_cloud_{kind}_{key}', set(scored.node_id)==set(s5.node_id), 'Full S5 cloud verified; missing theta joined by node ID')
            kept = survivors(scored)
            rank = truth_rank(scored,kept)
            row = {'stress_type': kind, 'level': case[levelcol],
                   'level_unit': {'amplitude':'dB','tracked_freq_drift':'fraction','ssp_mismatch':'environment'}[kind],
                   'config': case['config'], 'shape': case.get('shape',''), 'mode': case.get('mode',''),
                   'z_true': case['z_true'], 'n_survivors': len(kept), **metrics(kept), **rank,
                   'theta_source': scored.theta_source.iloc[0]}
            row['strict_single_bin'] = rank['truth_retained'] and rank['truth_rank']==1 and sorted(kept.r.unique())==[50.]
            stress_rows.append(row)
    stress = pd.DataFrame(stress_rows)
    check('stress_four_dimensional_quality_complete', len(stress)==144 and stress[[f'{a}_width' for a in AXES]].notna().all().all(), '96 amplitude + 30 tracked-frequency FIX2 + 18 SSP cases including controls')
    references = [
        ('amplitude', 'R3_RC23_TURN_AMPLITUDE_CLOSEOUT_FIX', 'AMPLITUDE_CASE_METRICS_FINAL.csv',
         ['config','shape','mode','level','z_true'], {'A_rms_db':'level','z_true_m':'z_true','true_rank':'truth_rank','surviving_r_width_km':'r_width'}),
        ('tracked_freq_drift', 'R3_RC23_TURN_TRACKED_FREQ_ROBUSTNESS_FIX2', 'TRACKED_FREQ_ANCHOR_CASES_FIXED.csv',
         ['config','level','z_true'], {'D':'level','z_true_m':'z_true','true_rank':'truth_rank','r_width_km':'r_width'}),
        ('ssp_mismatch', 'R3_RC23_TURN_SSP_MISMATCH', 'SSP_MISMATCH_ANCHOR_CASES.csv',
         ['config','level','z_true'], {'env_truth':'level','z_true_m':'z_true','true_rank':'truth_rank','r_width_km':'r_width'}),
    ]
    regressions = []
    for kind, directory, name, keys, rename in references:
        old = read(directory, name).rename(columns=rename)
        new = stress[stress.stress_type==kind].copy()
        old['level'], new['level'] = old.level.astype(str), new.level.astype(str)
        compared = new.merge(old, on=keys, suffixes=('_new','_old'), validate='one_to_one')
        check(f'{kind}_regression_coverage', len(compared)==len(new)==len(old), 'Every previous case matched by full case key')
        for field in ('truth_retained','truth_rank','strict_single_bin','r_width'):
            match = compared[field+'_new']==compared[field+'_old']
            check(f'{kind}_regression_{field}', bool(match.all()), f'{int(match.sum())}/{len(match)} cases agree')
            regressions.append({'stress_type':kind,'metric':field,'n_cases':len(match),'n_matches':int(match.sum()),'pass':bool(match.all())})
    summaries = []
    errorcols = ('max_rel_err_r','max_abs_err_theta','max_rel_err_v','max_rel_err_psi')
    for (kind,level,config), group in stress.groupby(['stress_type','level','config'],sort=False):
        summaries.append({'stress_type':kind,'level':level,'config':config,'n_cases':len(group),
                          'strict_count':int(group.strict_single_bin.sum()),'rank1_count':int((group.truth_rank==1).sum()),
                          'truth_retained_count':int(group.truth_retained.sum()),
                          **{f'max_{a}_width':float(group[f'{a}_width'].max()) for a in AXES},
                          **{f:float(group[f].max()) for f in errorcols}})
    summary = pd.DataFrame(summaries)
    for source in manifest:
        check(f"source_unchanged_{source['path']}", hashlib.sha256((ROOT/source['path']).read_bytes()).hexdigest()==source['sha256'], 'Read-only source SHA256 unchanged')
    bound = final_bound(depth)
    final = depth[depth.stage.isin(['S7','S8'])]
    remaining = [a for a in AXES if any(case['states'][0][a]!=s[a] for case in final_cases for s in case['states'])]
    top1 = bool(final.top1_truth_exact_in_matched_synthetic_control.all())
    worst = float(final[['max_rel_err_r','max_rel_err_v','max_rel_err_psi']].max().max())
    decision = {'stage':'R3_FINAL_PARAMETER_CONTRACTION_AUDIT_FIX','baseline_commit':BASELINE,
                'decision':'PARAMETER_CONTRACTION_AUDIT_COMPLETE','created_utc':now,
                'overall_label':'R3_CLOSEOUT_PROVISIONAL_PENDING_GPT_FINAL_AUDIT',
                'closure_checks':{name:True for name in ('branch_semantics_fixed','three_depth_tables_complete','connected_components_complete','bbox_corrected','truth_rank_audited','stress_four_dimensional_quality_complete')},
                'final_nominal_candidate_set':final_cases,'final_remaining_dimensions':remaining,
                'final_nominal_universal_lt10pct_set_bound':bound,'max_final_survivor_relative_diagnostic':worst,
                'top1_truth_exact_in_matched_synthetic_control':top1,
                'top1_truth_status':'UNIQUE_FULL_TRUTH_STATE_IS_BEST_IN_ALL_SIX_MATCHED_CONTROL_CASES' if top1 else 'NOT_UNIVERSALLY_EXACT',
                'diagnostic_label':'DIAGNOSTIC_ONLY_NOT_REQUIREMENT_VERIFICATION',
                'diagnostic_definition':'r/v/psi relative errors strictly <0.10; theta absolute error reported separately because truth=0; no new theta acceptance threshold; diagnostic only, not requirement verification',
                'truth_rank_definition':'1 + count(J < J_truth); ties use exact serialized J == J_truth; full scored parent cloud',
                'RC1_CONTRACTION_RATIO':'NOT_IDENTIFIABLE_FROM_CURRENT_EXPERIMENT',
                'exclusion_scope':'within the frozen 114576-node RC1-conditioned discrete search grid',
                'stress_candidate_case_count':len(stress),'physics_recomputed':False,
                'stress_scope':'All existing amplitude A=0/0.25/0.5/1dB, tracked-frequency FIX2 D=0/0.0025/0.005/0.01/0.02, SSP E0/E1/E2 candidate-score cases; amplitude A=2 has no candidate-level scores in this source and is outside this four-dimensional audit'}
    def csv(name, frame):
        frame.to_csv(OUT/name,index=False,lineterminator='\n')
    for names, frame in [(['PARAMETER_CONTRACTION_GRAPH.csv'],pd.DataFrame(graph)),
                          (['PARALLEL_CONFIGURATION_SET_OVERLAP.csv'],pd.DataFrame(overlaps)),
                          (['PARAMETER_CONTRACTION_BY_DEPTH_FIXED.csv','PARAMETER_CONTRACTION_BY_DEPTH.csv'],depth),
                          (['PARAMETER_CONTRACTION_STAGE_TABLE.csv'],all_rows[all_rows.z_true.isna() | (all_rows.z_true==200)]),
                          (['CANDIDATE_MULTIMODALITY_AUDIT.csv','JOINT_OCCUPANCY_AND_COMPONENTS.csv'],all_rows),
                          (['FINAL_SURVIVOR_STATE_TABLE_FIXED.csv','FINAL_SURVIVOR_STATE_TABLE.csv'],pd.DataFrame(final_states)),
                          (['STRESS_FINAL_SET_QUALITY_CASES.csv'],stress),
                          (['STRESS_FINAL_SET_QUALITY_SUMMARY.csv','STRESS_FINAL_SET_QUALITY.csv'],summary),
                          (['AUDIT_SOURCE_MANIFEST.csv'],pd.DataFrame(manifest)),
                          (['STRESS_METRIC_REGRESSION.csv'],pd.DataFrame(regressions)),
                          (['AUDIT_INTEGRITY_CHECKS.csv'],pd.DataFrame(checks))]:
        for name in names:
            csv(name,frame)
    diagnostic = depth.copy()
    diagnostic['diagnostic_role'] = np.where(diagnostic.stage.isin(['S7','S8']),'FINAL_NOMINAL','INTERMEDIATE')
    diagnostic['diagnostic_label'] = 'DIAGNOSTIC_ONLY_NOT_REQUIREMENT_VERIFICATION'
    csv('TRUTH_RELATIVE_SET_ERROR_DIAGNOSTIC.csv',diagnostic)
    report = f'''# R3 PARAMETER-CONTRACTION AUDIT-FIX

Baseline: `{BASELINE}`. Generated UTC: {now}.

`PARAMETER_CONTRACTION_AUDIT_COMPLETE` applies to this data audit. Overall status remains
`R3_CLOSEOUT_PROVISIONAL_PENDING_GPT_FINAL_AUDIT`; P5 is not opened.

## Branch semantics

```mermaid
flowchart LR
  S0["S0: 114576 RC1-conditioned grid"] --> S1["S1: 1326 W1 RC2"]
  S1 --> S2["S2: W1 + 235"]
  S0 --> S3["S3: 495 straight two-window RC2"]
  S1 -. recomputed branch .-> S3
  S3 --> S4["S4: straight two-window + 235"]
  S0 --> S5["S5: 385 turn15 RC2"]
  S1 -. recomputed branch .-> S5
  S5 --> S6["S6: 235"]
  S5 --> S7["S7: TRIPLE"]
  S5 --> S8["S8: FOUR"]
```

S6/S7/S8 are parallel configurations on the same S5 cloud. Set overlaps do not
establish serial filtering. Every graph edge uses actual node IDs. S7/S8 both
use S5 for contraction_from_previous. The nine-row overview uses z=200m for
display only; the authoritative depth table has all 15 scored stage/depth cases.

```text
{depth[['stage','z_true','n_candidates','r_width','theta_width','v_width','psi_width','truth_rank_min','n_tied_with_truth']].to_string(index=False)}
```

## Multimodality and occupancy

4D adjacency is Manhattan distance 1 on integer grid coordinates, with frozen
steps (1km, 0.5deg, 0.2m/s, 1deg), without boundary wrapping. Range components
are separate contiguous runs of the range projection. A connected range
projection alone does not prove 4D connectivity. Component sizes partition every cloud.

```text
{all_rows[['stage','z_true','n_candidates','n_connected_components','largest_component_size','n_range_components']].to_string(index=False)}
```

Bounding spans include missing internal bins. S6/z=200m has 15 occupied range bins,
16 bbox range bins, N_bbox={int(s6.N_bbox)}, bbox fraction={s6.bounding_box_grid_fraction:.12g},
correlation_sparsity={s6.correlation_sparsity:.12g}. Occupied-product uses separate fields.

## Final candidate quality

All six S7/S8 cases are generated from raw scores and cross-checked against the
existing survivor CSV. The final per-case states are in the decision JSON; the
candidate-level table includes J, profiled z, full-state truth flags and ranks.
First case (all others are separately audited):

```json
{json.dumps(final_cases[0],ensure_ascii=False,indent=2)}
```

Remaining dimensions: {js(remaining)}. Excluded fraction={float(final.contraction_from_initial.min()):.12g}
(about {float(final.contraction_from_initial.min())*100:.6f}%) **within the frozen 114576-node RC1-conditioned discrete search grid**.
`RC1_CONTRACTION_RATIO = NOT_IDENTIFIABLE_FROM_CURRENT_EXPERIMENT`.

`final_nominal_universal_lt10pct_set_bound={bound}` uses only S7/S8 x all three depths.
Worst full-survivor relative diagnostic={worst:.1%}; theta uses absolute error (truth=0).
The relative r/v/psi diagnostic uses strict <10%; theta absolute error is reported
separately, with no new acceptance threshold. A relative theta bound is undefined. This is
`DIAGNOSTIC_ONLY_NOT_REQUIREMENT_VERIFICATION`.
`top1_truth_exact_in_matched_synthetic_control={top1}` is a separate score result.
Truth rank counts all parent-cloud scores; better means J<J_truth; ties include
truth and use exact serialized equality, without a new numerical tolerance.

## Stress quality

144 cases: 96 amplitude, 30 tracked-frequency FIX2, 18 SSP, including controls.
Each case reports four parameter spans, occupied bins, worst truth errors,
full-state truth retention/rank and global-best range bins. Missing theta columns
are joined by node_id to verified S5 coordinates; r/v/psi and the frozen flattened
grid ID are checked. Range-only truth surrogates are not used.
Amplitude A=2dB has only a prior summary, no candidate-level scores in this input,
and is outside this four-dimensional audit. No absent cases are simulated.

```text
{summary.to_string(index=False)}
```

## Reproducibility

Run `python r3_final_parameter_contraction_fix.py`; the old entry point forwards here.
STRESS_METRIC_REGRESSION.csv verifies previous range widths, truth ranks, retention
and strict range-anchor metrics for every stress case (144/144 each).
AUDIT_SOURCE_MANIFEST.csv records hashes and row counts. AUDIT_INTEGRITY_CHECKS.csv
contains {len(checks)} passing checks. Source SHA256 values were verified unchanged.
Legacy filenames contain corrected tables to supersede the linear-chain and bbox errors.
No propagation, physical scenario, threshold, truth, seed or turn was changed.
'''
    (OUT/'R3_FINAL_PARAMETER_CONTRACTION_REPORT.md').write_text(report,encoding='utf-8')
    (OUT/'RC1_SUPPORT_ROLE.md').write_text('# RC1 support role\n\nRC1_CONTRACTION_RATIO = NOT_IDENTIFIABLE_FROM_CURRENT_EXPERIMENT\n\nAll exclusion fractions apply within the frozen 114576-node RC1-conditioned discrete search grid.\nRange support is 45:1:60km; no wider pre-RC1 prior was tested.\n',encoding='utf-8')
    (OUT/'GPT_SYNC.md').write_text(f'# R3 parameter-contraction audit-fix\n\nBaseline: `{BASELINE}`.\n\n**PARAMETER_CONTRACTION_AUDIT_COMPLETE**; {len(checks)} checks pass.\n\nFinal states: {js(final_cases[0]["states"])}. All six matched S7/S8 cases audited.\nRemaining dimensions: {js(remaining)}; full-set relative maximum: {worst:.1%}.\nfinal_nominal_universal_lt10pct_set_bound={bound}; top1 truth exact={top1}.\n\nDIAGNOSTIC_ONLY_NOT_REQUIREMENT_VERIFICATION. RC1 contraction is not identifiable.\nExclusion applies within the frozen 114576-node RC1-conditioned discrete search grid.\nOverall: R3_CLOSEOUT_PROVISIONAL_PENDING_GPT_FINAL_AUDIT. No new physics; P5 not opened.\n',encoding='utf-8')
    decision_path.write_text(json.dumps(decision,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    closeout_path = ROOT/'results'/'R3_FINAL_CLOSEOUT'/'R3_FINAL_DECISION.json'
    closeout = json.loads(closeout_path.read_text(encoding='utf-8'))
    closeout.update({'overall_label':decision['overall_label'],'parameter_contraction_audit':decision['decision'],
                     'note':'Structural data-only parameter-contraction audit complete; awaiting GPT final audit. Top1 truth and full-survivor diagnostics are separate.',
                     'parameter_contraction_decision':'../R3_FINAL_PARAMETER_CONTRACTION/R3_FINAL_PARAMETER_CONTRACTION_DECISION.json'})
    closeout_path.write_text(json.dumps(closeout,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(decision['decision'],f'checks={len(checks)}',f'stress_cases={len(stress)}')
    print('final_remaining_dimensions',remaining,'full_set_bound',bound,'top1_truth_exact',top1)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
