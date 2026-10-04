"""Read-only reconstruction of saved B1 metrics/Gates; append-only closeout artifacts."""
from pathlib import Path
from collections import defaultdict
import csv
import hashlib
import json
import math
import subprocess
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY'

def read(name):
    with (OUT/name).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def js(name): return json.loads((OUT/name).read_text(encoding='utf-8-sig'))
def write(name,value): (OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def verify_saved():
    d=js('B1_DESIGN_FREEZE.json'); decision=js('B1_DECISION.json'); validation=js('B1_VALIDATION.json')
    groups=defaultdict(list)
    for r in read('B1_DEPTH_PROFILES.csv'): groups[r['case_id'],r['conditioning_id']].append(r)
    metrics={(r['case_id'],r['conditioning_id']):r for r in read('B1_PROFILE_METRICS.csv')}
    cases={c['case_id']:c for c in d['cases']}; checks=[]; tolerance=validation['numerical_reconstruction_tolerance_db']
    def check(name,ok):
        checks.append({'check':name,'pass':bool(ok)})
        if not ok: raise RuntimeError('Saved evidence reconstruction failed: '+name)
    expected={(cid,c['conditioning_id']) for cid,conditions in d['horizontal_conditioning'].items() for c in conditions}
    check('all_profile_keys',set(groups)==set(metrics)==expected)
    rebuilt=[]
    for key,rows in groups.items():
        rows.sort(key=lambda r:float(r['z_label_m'])); z=[float(r['z_label_m']) for r in rows]; j=[float(r['J_z_db']) for r in rows]
        check(str(key)+'labels',z==list(range(150,251,5)))
        c=cases[key[0]]; m=metrics[key]; i=z.index(c['z_true_m']); best=min(range(21),key=lambda k:j[k]); low=j[best]
        truth_rank=1+sum(x<j[i] for x in j); tie=sum(x<=low+tolerance for x in j)
        second=min(j[k] for k in range(21) if k!=i)-j[i]
        neighbor=min(j[k] for k in (i-1,i+1) if 0<=k<21)-j[i]
        curvature=j[i-1]-2*j[i]+j[i+1] if 0<i<20 else None
        local=sum(j[k]<j[k-1] and j[k]<j[k+1] for k in range(1,20))
        boundary=int(j[0]<j[1])+int(j[-1]<j[-2]); error=abs(z[best]-c['z_true_m'])
        width=None; component=[]
        if max(j)-low>tolerance:
            allowed=[(x-low)/(max(j)-low)<=.5 for x in j]; a=b=best
            while a>0 and allowed[a-1]: a-=1
            while b<20 and allowed[b+1]: b+=1
            width=z[b]-z[a]; component=z[a:b+1]
        numbers={'z_hat_grid_m':z[best],'true_depth_rank':truth_rank,'min_tie_count':tie,'J_true_db':j[i],'J_min_db':low,
                 'DeltaJ_second_db':second,'DeltaJ_neighbor_db':neighbor,'K_z_db':curvature,
                 'strict_interior_local_minima':local,'boundary_local_minima':boundary,
                 'dynamic_range_db':max(j)-low,'W50_z_m':width,'absolute_depth_error_m':error}
        for field,value in numbers.items():
            check(str(key)+field,m[field]=='' if value is None else math.isclose(float(m[field]),value,rel_tol=1e-12,abs_tol=1e-12))
        check(str(key)+'W50_component',json.loads(m['W50_component_labels'])==component)
        condition=next(item for item in d['horizontal_conditioning'][key[0]] if item['conditioning_id']==key[1])
        for r in rows:
            check(str(key)+'state_'+r['z_label_m'],json.loads(r['horizontal_state'])==condition['horizontal_state'])
            mapping=json.loads(r['effective_source_depth_m'])
            for f in c['frequencies_hz']:
                depth=d['modal_depth_mapping'][str(f)]['source_label_mapping'][str(int(float(r['z_label_m'])))]
                check(str(key)+'mapping_'+str(f)+'_'+r['z_label_m'],mapping[str(f)]==depth)
        rebuilt.append({'case_id':key[0],'type':condition['conditioning_type'],'unique':z[best]==c['z_true_m'] and tie==1,
                        'positive_second':second>tolerance,'positive_neighbor':neighbor>tolerance,'curvature':curvature>tolerance,
                        'error':error,'rank':truth_rank})
    h0=[r for r in rebuilt if r['type']=='H0']; h1=[r for r in rebuilt if r['type']=='H1']
    h0pass=all(r['unique'] and r['positive_second'] and r['positive_neighbor'] and r['curvature'] for r in h0)
    h1pass=all(r['error']<=10 and r['rank']<=3 for r in h1)
    label='B1_CONDITIONAL_DEPTH_IDENTIFIABILITY_ESTABLISHED' if h0pass and h1pass else 'B1_DEPTH_IDENTIFIABLE_ONLY_UNDER_TIGHT_HORIZONTAL_CONDITIONING' if h0pass else 'B1_CONDITIONAL_DEPTH_IDENTIFIABILITY_NOT_ESTABLISHED'
    check('independent_scientific_decision',label==decision['scientific_decision'])
    totals={'H0_unique_true_minima':sum(r['unique'] for r in h0),'H0_positive_second_best_margin':sum(r['positive_second'] for r in h0),
            'H0_positive_neighbor_margin':sum(r['positive_neighbor'] for r in h0),'H0_positive_curvature':sum(r['curvature'] for r in h0),
            'H1_total':len(h1),'H1_exact_depth':sum(r['error']==0 for r in h1),'H1_within_5m':sum(r['error']<=5 for r in h1),
            'H1_within_10m':sum(r['error']<=10 for r in h1),'H1_true_rank_le3':sum(r['rank']<=3 for r in h1),'worst_H1_depth_error_m':max(r['error'] for r in h1)}
    for field,value in totals.items(): check(field,decision[field]==value)
    replay=read('B1_NUMERICAL_RECONSTRUCTION.csv'); check('replay_count',len(replay)==1134)
    for row in replay:
        check('replay:'+row['case_id']+'/'+row['conditioning_id']+'/'+row['z_label_m'],
              max(abs(float(row['J_primary_db'])-float(row['J_repeat_db'])),abs(float(row['J_primary_db'])-float(row['J_independent_db'])))<=d['reconstruction']['maximum_allowed_absolute_difference_db'])
    for path,digest in d['bindings'].items(): check('frozen:'+path,sha(ROOT/path)==digest)
    for path,digest in d['protected_historical_files'].items(): check('protected:'+path,sha(ROOT/path)==digest)
    for name in ('R4_PLAN.md','R4_EVIDENCE_LEDGER.csv'):
        path='results/R4_MASTER/'+name
        old=subprocess.check_output(['git','show',d['parent_sha']+':'+path],cwd=ROOT).replace(b'\r\n',b'\n')
        check('append_only:'+name,(ROOT/path).read_bytes().replace(b'\r\n',b'\n').startswith(old))
    check('no_execution_invalid_marker',not (OUT/'B1_EXECUTION_INVALID.json').exists())
    return checks

def finalize():
    if (OUT/'B1_RESULTS_FREEZE.json').exists(): raise RuntimeError('Results already sealed')
    checks=verify_saved(); d=js('B1_DESIGN_FREEZE.json'); result=js('B1_DECISION.json'); v=js('B1_VALIDATION.json')
    testpath=Path('D:/ProjectStorage/WSL-CZ/b1-20261005/tests.xml'); root=ET.parse(testpath).getroot(); suites=list(root.iter('testsuite'))
    tests=sum(int(s.attrib['tests']) for s in suites); failed=sum(int(s.attrib['failures'])+int(s.attrib['errors']) for s in suites)
    assert tests>0 and failed==0
    (OUT/'B1_TESTS.xml').write_bytes(testpath.read_bytes()); write('B1_TESTS.json',{'passed':tests,'failed':failed})
    # Standalone scientific figure: all 54 saved profiles, no extra physical evaluations.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data=read('B1_DEPTH_PROFILES.csv'); fig,axes=plt.subplots(3,2,figsize=(12,12),sharex=True)
    for ax,case in zip(axes.T.ravel(),d['cases']):
        cid=case['case_id']; own=[r for r in data if r['case_id']==cid]
        for condition in d['horizontal_conditioning'][cid]:
            group=sorted([r for r in own if r['conditioning_id']==condition['conditioning_id']],key=lambda r:float(r['z_label_m']))
            ax.plot([float(r['z_label_m']) for r in group],[float(r['J_z_db']) for r in group],
                    linewidth=2.5 if condition['conditioning_type']=='H0' else .9,alpha=1 if condition['conditioning_type']=='H0' else .65,
                    label=condition['conditioning_id'])
        ax.axvline(case['z_true_m'],color='black',linestyle=':',linewidth=1)
        ax.set_title(cid+' ('+case['historical_stage']+')'); ax.set_ylabel('RMS relative-level error J (dB)'); ax.grid(alpha=.2); ax.legend(fontsize=7,ncol=3)
    for ax in axes[-1]: ax.set_xlabel('Source-depth label (m)')
    fig.suptitle('B1: matched conditional depth profiles; H1 = inherited one-axis grid step')
    fig.tight_layout(); fig.savefig(OUT/'B1_DEPTH_PROFILE_CURVES.svg'); fig.savefig(OUT/'B1_DEPTH_PROFILE_CURVES.png',dpi=150); plt.close(fig)
    h0=read('B1_H0_SUMMARY.csv'); h1=read('B1_H1_SUMMARY.csv'); metrics=read('B1_PROFILE_METRICS.csv')
    table='| Case | H0 z-hat | Second margin (dB) | Neighbor margin (dB) | Curvature (dB) | Interior minima | W50 (m) | H1 <=10 m | H1 rank<=3 | Worst H1 (m) |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n'
    for a,b in zip(h0,h1):
        table+='| '+a['case_id']+' | '+a['z_hat_grid_m']+' | '+format(float(a['DeltaJ_second_db']),'.6g')+' | '+format(float(a['DeltaJ_neighbor_db']),'.6g')+' | '+format(float(a['K_z_db']),'.6g')+' | '+a['strict_interior_local_minima']+' | '+a['W50_z_m']+' | '+b['within_10m']+'/8 | '+b['true_rank_le3']+'/8 | '+b['worst_error_m']+' |\n'
    axis_table='| H1 axis | Profiles | Exact | <=5 m | <=10 m | Rank<=3 | Worst error (m) |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for axis in ('r_km','theta_deg','v_mps','psi_deg'):
        rows=[r for r in metrics if r['conditioning_type']=='H1' and r['perturbed_axis']==axis]
        errors=[float(r['absolute_depth_error_m']) for r in rows]
        axis_table+='| '+axis+' | '+str(len(rows))+' | '+str(sum(e==0 for e in errors))+' | '+str(sum(e<=5 for e in errors))+' | '+str(sum(e<=10 for e in errors))+' | '+str(sum(int(r['true_depth_rank'])<=3 for r in rows))+' | '+str(max(errors))+' |\n'
    report=f"""# R4 B1 conditional depth profile identifiability

Decision: **{result['scientific_decision']}**. Execution is valid; research-lead audit is pending. Confirmed R4 progress remains 0%; proposed B1/overall progress is {result['proposed_R4_overall_progress_percent']}%.

The accepted historical six-case definition is S7 (201/235/283 Hz) x 180/200/220 m plus S8 (201/235/283/338 Hz) x the same depths. They share one horizontal scenario and noiseless matched acoustic observation rule; they are not six independent scenes or Monte Carlo replicates. The three-frequency primary route is S7; S8 retains the historical four-line companion without converting it to a new three-line score. All six are included in every all-case Gate.

Design was committed and pushed first at {result['design_freeze_sha']}, parent {d['parent_sha']}. Cases were selected by the complete accepted final_nominal_candidate_set, never by depth success. All twelve final survivor rows (psi=4 and psi=5) are separately cold reconstructed in B1_HISTORICAL_SCORE_REPLAY.csv. The old artifact stores profiled minima, not complete J(z) or a standalone acoustic observation array: B1_OBSERVATIONS.csv explicitly reconstructs the accepted generation rule after design push. No new signal/noise assumption is introduced.

H0 uses exact generation/reference horizontal state as an authorized controlled input. H1 uses all eight legal +/- inherited grid steps (1 km,0.5 degree,0.2 m/s,1 degree), one axis at a time. This describes the historically inherited local grid panel, not a promise that these errors are realistic post-acquisition errors or a demonstrated continuous tolerance radius. Score is the accepted per-frequency/per-window demeaned RMS over 121 times, without bearing or Q2 representation. Receiver label is 200 m; nearest stored source/receiver modal depths and all 21 depth labels are saved, with no interpolation.

{table}

{axis_table}

H0 unique true-depth minima: {result['H0_unique_true_minima']}/6. Positive second/neighbor margins: {result['H0_positive_second_best_margin']}/6 and {result['H0_positive_neighbor_margin']}/6. Positive curvature: {result['H0_positive_curvature']}/6. H1 exact/within5/within10/rank<=3: {result['H1_exact_depth']}/{result['H1_within_5m']}/{result['H1_within_10m']}/{result['H1_true_rank_le3']} of 48. Worst H1 error: {result['worst_H1_depth_error_m']} m. Other strict local minima and boundary minima are preserved; a unique global truth minimum does not imply one local valley. Full profiles/metrics are exported without dropping any row.

W50 is the span between terminal labels of the contiguous component J_norm<=0.5 containing the deterministic global minimum. It is not a confidence interval; a width of 0 m means one retained depth label, not zero uncertainty. No spline, parabola, continuous depth or off-grid z is used. The 5 m nuisance labels may map to nearby stored physical samples; numerical mapping is retained explicitly per frequency.

Numerical tolerance was computed by the preregistered rule max(128*eps*max(1,max_abs_J),10*max_repeat_or_independent_residual), not from a depth margin. Result: {v['numerical_reconstruction_tolerance_db']:.12g} dB. Cold repeat maximum: {v['maximum_repeat_difference_db']:.12g} dB; independent score/observation maximum: {v['maximum_independent_difference_including_observations_db']:.12g} dB. Independent parser uses struct, geometry uses scalar math, propagation uses ordered compensated mode sums, and score uses math.fsum. This supplies a numerical implementation check of the same matched physical model, not a second physical model validation. {tests} unit tests pass; runtime validation has {v['checks_passed']} passing checks, including {v['historical_files_preserved']} historical byte hashes. Saved metric/Gate reconstruction is separately recorded.

The CZ certified global acquisition branch is CLOSED: WINDOW_SHAPE_Q2 + N0 + whole-cell conservative partition over 45-60 km cold-start 4D support did not establish closure under registered budgets. Its mathematical lower-bound mechanism remains historical evidence. No C4/F4, grid/tau/enclosure fix, budget increase or 21-case CZ development is run. A1 remains BLOCKED / NOT_COMPLETED, 0/15%.

Maximum claim is ORACLE/EXACT-HORIZONTAL-CONDITIONED DEPTH MECHANISM in synthetic matched E0, on-grid source labels, with H1 outcomes explicitly bounded by the fixed panel. A CONDITIONAL result supports the exact-conditioned mechanism and identifies horizontal/depth coupling; it does not establish an untested tighter continuous neighborhood. It gives no cold-start joint inversion, end-to-end depth performance, environmental robustness, real ocean/UUV observation or final accuracy guarantee. B2 is NOT_OPENED; A2/A3/A4/B3/B4/C and P5 are NOT_OPENED. Stop after execution commit/push and wait independent audit.

See [all saved profile curves](B1_DEPTH_PROFILE_CURVES.svg), [profiles](B1_DEPTH_PROFILES.csv), [metrics](B1_PROFILE_METRICS.csv), [decision](B1_DECISION.json), [validation](B1_VALIDATION.json).
"""
    (OUT/'B1_REPORT.md').write_text(report,encoding='utf-8')
    (OUT/'GPT_SYNC.md').write_text(report,encoding='utf-8')
    master=ROOT/'results/R4_MASTER'
    plan=master/'R4_PLAN.md'; ledger=master/'R4_EVIDENCE_LEDGER.csv'
    with plan.open('a',encoding='utf-8',newline='') as file:
        file.write('\n\n## R4 B1 conditional-depth execution checkpoint\n\n'+
                   f"{result['scientific_decision']}; PENDING_RESEARCH_LEAD_AUDIT. Design push {result['design_freeze_sha']} preceded all B1 depth evaluations. Complete six historical configurations (S7 three-frequency and S8 four-frequency, three depths) are preserved. H0 unique truth minima {result['H0_unique_true_minima']}/6; H1 within10 {result['H1_within_10m']}/48, rank<=3 {result['H1_true_rank_le3']}/48, worst error {result['worst_H1_depth_error_m']} m. Proposed progress {result['proposed_R4_overall_progress_percent']}%; confirmed R4 remains 0%. No full 10% credit is awarded before lead audit.\n\n"+
                   'CZ certified global-search route CLOSED; A1 BLOCKED / NOT_COMPLETED (0/15%), acquisition implementation not established. Conditional depth mechanism is separate from cold start and end-to-end performance. [B1 report](../R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY/B1_REPORT.md) and [decision](../R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY/B1_DECISION.json). B2 and A2/A3/A4/B3/B4/C/P5 NOT_OPENED. Stop after commit/push.\n')
    with ledger.open('a',encoding='utf-8',newline='') as file:
        writer=csv.writer(file,lineterminator='\n'); writer.writerow(['B1','../R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY/B1_DECISION.json',result['scientific_decision'],'54 conditional on-grid profiles; S7 triple and historical S8 four-line companions','PENDING_RESEARCH_LEAD_AUDIT',0])
    checks.extend(verify_saved())
    with (OUT/'B1_SAVED_EVIDENCE_RECONSTRUCTION.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=['check','pass'],lineterminator='\n'); writer.writeheader(); writer.writerows(checks)
    write('B1_SAVED_EVIDENCE_RECONSTRUCTION.json',{'checks_passed':len(checks),'failed':0,'tests_passed':tests,'numeric_profiles_rerun':0,'independent_saved_metrics_gate_reconstruction':True})
    paths=sorted(p for p in OUT.iterdir() if p.is_file() and p.name!='B1_RESULTS_FREEZE.json')
    bindings={p.relative_to(ROOT).as_posix():sha(p) for p in paths}
    for name in ('R4_PLAN.md','R4_EVIDENCE_LEDGER.csv'): bindings[(master/name).relative_to(ROOT).as_posix()]=sha(master/name)
    write('B1_RESULTS_FREEZE.json',{'design_sha':result['design_freeze_sha'],'design_file_sha256':sha(OUT/'B1_DESIGN_FREEZE.json'),'decision':result['scientific_decision'],'bindings':bindings,'immutable_after_commit':True})
    print(json.dumps({'saved_reconstruction_checks':len(checks),'tests':tests,'sealed_artifacts':len(bindings),'scientific_decision':result['scientific_decision']}))

if __name__=='__main__':
    if '--verify' in __import__('sys').argv:
        seal=js('B1_RESULTS_FREEZE.json')
        for path,digest in seal['bindings'].items(): assert sha(ROOT/path)==digest,path
        checks=verify_saved(); print(json.dumps({'sealed':'PASS','reconstruction_checks':len(checks),'failed':0}))
    else: finalize()
