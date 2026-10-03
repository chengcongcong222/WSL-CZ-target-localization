"""Offline non-oracle projection and independent unreleased-harness checks."""
import hashlib
import inspect
import json
import subprocess
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import r4_a1_search_tractability as t
import r4_a1_search_tractability_development as d
import r4_a1_search_tractability_gates as g

SOURCE=Path('results/R4_A1_FIX2_ACOUSTIC_COVERAGE/SEARCH_BUDGET_CONVERGENCE.csv')
ALLOWED_SOURCE_COLUMNS=['case_id','panel_id','origin_group','sigma_deg','seed','cutoff']


def prepare_manifest():
    d.verify_infrastructure()
    source=pd.read_csv(t.ROOT/SOURCE,usecols=ALLOWED_SOURCE_COLUMNS,float_precision='round_trip')
    selected=source[source.sigma_deg==.1].drop_duplicates().sort_values('case_id',kind='stable').reset_index(drop=True)
    with np.load(t.ROOT/d.OBSERVATION_PATH,allow_pickle=False) as archive:
        ids=archive['case_ids'].astype(str).tolist()
        if len(ids)!=len(set(ids)):raise ValueError('Duplicate archived observation identifier')
        positions={cid:i for i,cid in enumerate(ids)}
        selected['observation_index']=[positions[cid] for cid in selected.case_id]
        selected=selected.rename(columns={'cutoff':'bearing_cutoff'})[d.MANIFEST_COLUMNS]
        d.validate_manifest(selected);d.map_observations(selected,archive)
    t.save('DEVELOPMENT_CASE_MANIFEST.csv',selected)
    t.json_write('DEVELOPMENT_MANIFEST_PROVENANCE.json',dict(parent_pre_run_commit=d.PARENT,
        source_path=SOURCE.as_posix(),source_allowed_columns=ALLOWED_SOURCE_COLUMNS,
        selection_rule='exactly sigma_deg==0.1; identical non-oracle duplicates removed; no recovery/error/score/state column read',
        observation_path=d.OBSERVATION_PATH.as_posix(),mapping_rule='archive case_ids identifier -> explicit observation_index; row-order assumptions rejected',
        source_sha256=hashlib.sha256((t.ROOT/SOURCE).read_bytes()).hexdigest(),
        observation_sha256=hashlib.sha256((t.ROOT/d.OBSERVATION_PATH).read_bytes()).hexdigest(),
        n_unique_case_ids=21,noisy_development_runs=0,development_released=False))
    return selected


def audit_harness(check):
    freeze=json.loads((t.OUT/'DEVELOPMENT_EXECUTION_FREEZE.json').read_text())
    check('Harness parent and release boundary',freeze['parent_pre_run_commit']==d.PARENT and freeze['development_released']==False and freeze['noisy_development_runs']==0 and not d.DEVELOPMENT_RELEASED,'no scientific release')
    check('Accepted infrastructure SHA constants',freeze['budget_sha256']==d.BUDGET_SHA and freeze['core_objective_sha256']==d.CORE_SHA,'accepted immutable infrastructure')
    for path,digest in freeze['sha256'].items():
        check('Harness frozen bytes '+path,hashlib.sha256((t.ROOT/path).read_bytes()).hexdigest()==digest,'SHA256')
    d.verify_infrastructure()
    frame=pd.read_csv(t.OUT/'DEVELOPMENT_CASE_MANIFEST.csv',float_precision='round_trip');d.validate_manifest(frame)
    check('Exactly 21 nominal manifest cases',len(frame)==frame.case_id.nunique()==21 and (frame.sigma_deg==.1).all(),'no noiseless or duplicates')
    check('Manifest contains only permitted columns',set(frame.columns)==set(d.MANIFEST_COLUMNS),'no privileged metrics or state columns')
    # Independent reconstruction reads only the same permitted source fields.
    projection=pd.read_csv(t.ROOT/SOURCE,usecols=ALLOWED_SOURCE_COLUMNS,float_precision='round_trip')
    projection=projection[projection.sigma_deg==.1].drop_duplicates().set_index('case_id')
    check('Frozen nominal identifiers preserved',set(frame.case_id)==set(projection.index) and len(projection)==21,'not selected by solver recovery')
    for row in frame.itertuples():
        p=projection.loc[row.case_id]
        check('Non-oracle projection '+row.case_id,row.seed==p.seed and row.sigma_deg==p.sigma_deg and row.bearing_cutoff==p.cutoff and row.panel_id==p.panel_id and row.origin_group==p.origin_group,'source identifiers and bearing cutoff only')
    with np.load(t.ROOT/d.OBSERVATION_PATH,allow_pickle=False) as archive:
        ids=archive['case_ids'].astype(str).tolist()
        for row in frame.itertuples():
            check('Explicit observation mapping '+row.case_id,ids.count(row.case_id)==1 and ids[row.observation_index]==row.case_id,'case_ids mapping')
        observations=d.map_observations(frame,archive)
        check('Selected observations finite and readonly',len(observations)==21 and all(not x.bearing_rad.flags.writeable and not x.relative_tl.flags.writeable for x in observations.values()),'mapping only; no objective evaluated')
    check('21 inherited nuisance branches',len(t.L.PROFILE)==21 and len(set(t.L.PROFILE))==21 and freeze['depth_labels_m']==t.L.PROFILE.tolist(),'no depth estimator')
    check('Witness semantics',freeze['candidate_kind']==g.KIND and freeze['scientific_semantics']=='threshold-hit / discrete witness-cluster convergence; no stationarity certificate','no minima relabeling')
    forbidden=['truth','errors','neighborhood','oracle','basin_width','basin_direction']
    for fn in [d.run_branch,d.run_case,d.run_family,d.extract_witnesses,d.rebuild_witnesses,g.raw_case_convergence,g.dual_case_agreement]:
        source=inspect.getsource(fn).lower()
        check('Search/Gate static exclusion '+fn.__name__,not any(word in source or word in str(inspect.signature(fn)).lower() for word in forbidden),'signature/source')
    try:d.execute_development()
    except SystemExit as error:stopped=str(error)=='DEVELOPMENT_NOT_RELEASED: research-lead audit required'
    else:stopped=False
    check('Unreleased real entry hard stop',stopped,'before loading runtime inputs/models')
    check('Scientific outputs absent',not any((t.OUT/name).exists() for name in ['DEVELOPMENT_SHGO_RESULTS.csv','DEVELOPMENT_DIRECT_RESULTS.csv','SHGO_WITNESS_CATALOG.csv','DIRECT_WITNESS_CATALOG.csv','SHGO_MINIMA_CATALOG.csv','DIRECT_MINIMA_CATALOG.csv','FRESH_HOLDOUT_PANEL.csv']),'no real run and no manufactured scientific catalog')
    fixtures=json.loads((t.OUT/'WITNESS_GATE_REFERENCE_FIXTURES.json').read_text())
    for fixture in fixtures['cases']:
        left=pd.DataFrame(fixture['T2_catalog']);right=pd.DataFrame(fixture['T3_catalog'])
        raw=g.raw_case_convergence(fixture['T2_result'],fixture['T3_result'],left,right)
        dual=g.dual_case_agreement(fixture['T2_result'],fixture['T3_result'],left,right)
        check('Independent synthetic semantic reference '+fixture['id'],raw['raw_T2_T3_converged']==fixture['expected_raw'] and dual['dual_solver_subthreshold_witness_cluster_agreement']==fixture['expected_dual'],'fixed expected booleans; no acoustic/noisy case run')


def freeze_execution():
    d.verify_infrastructure()
    if d.DEVELOPMENT_RELEASED:raise RuntimeError('This checkpoint must remain unreleased')
    parent_method=subprocess.check_output(['git','show',d.PARENT+':results/R4_A1_SEARCH_TRACTABILITY_AUDIT/METHOD_FREEZE.json'])
    archived=t.OUT/'ACCEPTED_PRE_RUN_METHOD_FREEZE.json'
    if archived.exists():assert archived.read_bytes()==parent_method
    else:archived.write_bytes(parent_method)
    paths=['r4_a1_search_tractability.py','r4_a1_search_tractability_development.py','r4_a1_search_tractability_gates.py',
        'r4_a1_search_tractability_harness_audit.py','r4_a1_search_tractability_audit.py','tests/test_r4_a1_search_tractability.py','.gitattributes']
    paths += [(t.OUT/name).relative_to(t.ROOT).as_posix() for name in ['DEVELOPMENT_CASE_MANIFEST.csv','DEVELOPMENT_MANIFEST_PROVENANCE.json','WITNESS_GATE_REFERENCE_FIXTURES.json','WITNESS_SEMANTICS.md','TRACTABILITY_BUDGET_FREEZE.json','TRACTABILITY_METHOD.md','A1_SEARCH_TRACTABILITY_AUDIT_DESIGN.md','DEVELOPMENT_EXECUTION_HARNESS_DESIGN.md','ACCEPTED_PRE_RUN_METHOD_FREEZE.json']]
    paths += [d.OBSERVATION_PATH.as_posix()]
    hashes={p:hashlib.sha256((t.ROOT/p).read_bytes()).hexdigest() for p in paths}
    t.json_write('DEVELOPMENT_EXECUTION_FREEZE.json',dict(parent_pre_run_commit=d.PARENT,frozen_utc=datetime.now(timezone.utc).isoformat(),
        checkpoint='DEVELOPMENT_EXECUTION_HARNESS_FREEZE',budget_sha256=d.BUDGET_SHA,core_objective_sha256=d.CORE_SHA,
        development_driver_sha256=hashes['r4_a1_search_tractability_development.py'],gate_logic_sha256=hashes['r4_a1_search_tractability_gates.py'],
        noisy_development_runs=0,development_released=False,case_manifest_count=21,depth_labels_m=t.L.PROFILE.tolist(),candidate_kind=g.KIND,
        scientific_semantics='threshold-hit / discrete witness-cluster convergence; no stationarity certificate',sha256=hashes))
    method=json.loads(parent_method)
    method['sha256'].update(hashes)
    method['sha256']['results/R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_EXECUTION_FREEZE.json']=hashlib.sha256((t.OUT/'DEVELOPMENT_EXECUTION_FREEZE.json').read_bytes()).hexdigest()
    method.update(checkpoint='DEVELOPMENT_EXECUTION_HARNESS_FREEZE_REQUIRES_SECOND_INDEPENDENT_AUDIT',parent_pre_run_commit=d.PARENT,
        frozen_utc=datetime.now(timezone.utc).isoformat(),development_enabled_in_code=False,noisy_development_runs=0)
    t.json_write('METHOD_FREEZE.json',method)
