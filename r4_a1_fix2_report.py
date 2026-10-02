"""Summaries, scientific Gate decision and figures from saved FIX2 evidence."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import r4_a1_fix2_coverage as f


def report():
    out=f.OUT;budgets=pd.read_csv(out/'SEARCH_BUDGET_CONVERGENCE.csv');dev=pd.read_csv(out/'DEVELOPMENT_RECOVERY.csv');hold=pd.read_csv(out/'FRESH_HOLDOUT_RESULTS.csv')
    independent=pd.read_csv(out/'INDEPENDENT_SOLVER_AGREEMENT.csv');geometry=pd.read_csv(out/'BASIN_GEOMETRY.csv');alias=pd.read_csv(out/'ALIAS_DIAGNOSTIC.csv')
    supplementary=[];primary=dev.set_index('case_id')
    for row in independent.itertuples():
        baseline=primary.loc[row.case_id];diff=abs(np.array([getattr(row,a) for a in f.L.AXES])-baseline[list(f.L.AXES)].to_numpy(float))
        for k in [1,3]:diff[k]=float(f.L.angle_error_deg(getattr(row,f.L.AXES[k]),baseline[f.L.AXES[k]]))
        expected=bool(row.exact_J<.001 and baseline.exact_J<.001 and np.all(diff<=f.SEARCH_CONFIG['independent_agreement_state_tolerances']) and row.z_star_label_m==baseline.z_star_label_m and row.bearing_cost<=row.cutoff+1e-14)
        assert expected==row.agreement and np.max(abs(diff-np.array(json.loads(row.differences))))<1e-9
        supplementary.append(dict(check='Independent solver agreement '+row.case_id,pass_check=True,detail=json.dumps(diff.tolist())))
    fresh_obs=np.load(out/'FRESH_HOLDOUT_OBSERVATIONS.npz');panel=pd.read_csv(out/'FRESH_HOLDOUT_PANEL.csv');panel_ids=panel.panel_id.tolist()
    for row in hold.itertuples():
        pi=panel_ids.index(row.panel_id);ci=list(fresh_obs['case_ids']).index(row.case_id);truth=panel.iloc[pi]
        noiseless=f.L.geometry(truth[list(f.L.AXES)].to_numpy(float))[0][0]
        eps=np.zeros(121) if row.sigma_deg==0 else np.random.default_rng(np.random.SeedSequence([row.seed,pi])).standard_normal(121)
        delta=float(np.max(abs(fresh_obs['bearing_rad'][ci]-noiseless-np.radians(row.sigma_deg)*eps)))
        assert delta<1e-15
        supplementary.append(dict(check='Fresh observation seed '+row.case_id,pass_check=True,detail=delta))
    f.save('SUPPLEMENTARY_VALIDATION.csv',pd.DataFrame(supplementary))
    levels=['B1','B2','B3'];nominal=budgets[budgets.sigma_deg==.1]
    summary=nominal.groupby('budget').agg(n_cases=('case_id','size'),recovered=('recovered','sum'),raw_component_recovered=('raw_component_recovered','sum'),exact_J_P50=('exact_J','median'),exact_J_worst=('exact_J','max'))
    summary['recovery_fraction']=summary.recovered/summary.n_cases;summary['scope']='ALL_PRIOR_PANELS_ARE_DEVELOPMENT';f.save('BUDGET_RECOVERY_SUMMARY.csv',summary.reset_index())
    detail=nominal.groupby('budget')[['n_bank_feasible_nodes','n_local_starts','n_separated_candidates','top1_rel_r','top1_abs_theta_deg','top1_rel_v','top1_abs_psi_deg']].agg(['median','max'])
    detail.columns=['_'.join(c) for c in detail.columns];f.save('BUDGET_STATE_AND_BASIN_SUMMARY.csv',detail.reset_index())
    pairs=[]
    for cid,group in nominal.groupby('case_id'):
        b=group.set_index('budget');x=b.loc['B2'];y=b.loc['B3'];difference=abs(x[list(f.L.AXES)].to_numpy(float)-y[list(f.L.AXES)].to_numpy(float))
        for a in (1,3):difference[a]=float(f.L.angle_error_deg(x[f.L.AXES[a]],y[f.L.AXES[a]]))
        stable=bool(x.recovered and y.recovered and np.all(difference<=f.SEARCH_CONFIG['convergence_state_tolerances']) and x.z_star_label_m==y.z_star_label_m)
        raw_stable=bool(x.raw_component_recovered and y.raw_component_recovered)
        pairs.append(dict(case_id=cid,stable_B2_B3_cumulative=stable,raw_components_both_recovered=raw_stable,state_change=json.dumps(difference.tolist()),origin_budget_at_B3=y.origin_budget_of_best))
    stability=pd.DataFrame(pairs);f.save('BUDGET_STABILITY.csv',stability)
    diag=[];widths=pd.concat([geometry,pd.read_csv(out/'FRESH_BASIN_GEOMETRY.csv')])
    for scope,table in [('DEVELOPMENT',dev),('FRESH_CONFIRMATION',hold)]:
        for row in table.to_dict('records'):
            cid=row['case_id'];measurement=widths[(widths.case_id==cid)&(widths.threshold_db==.001)]
            row.update(scope=scope,measured_J001_sections_json=json.dumps(measurement[['direction','width','width_status']].to_dict('records')),solver_agreement_status='NOT_TESTED')
            match=independent[independent.case_id==cid]
            if len(match):row['solver_agreement_status']='AGREED' if bool(match.iloc[0].agreement) else 'DISAGREED'
            pair=stability[stability.case_id==cid]
            row['budget_convergence_status']='NOISELESS_RC2_SINGLETON' if row['sigma_deg']==0 else ('NOT_IN_DEVELOPMENT_CONVERGENCE_PANEL' if not len(pair) else ('STABLE' if bool(pair.iloc[0].stable_B2_B3_cumulative and pair.iloc[0].raw_components_both_recovered) else 'NOT_CLOSED'))
            diag.append(row)
    f.save('CASE_DIAGNOSTICS.csv',pd.DataFrame(diag))
    strict_geometry=geometry[(geometry.direction.isin(['r_km','v_mps','psi_deg','theta_offset_deg']))&(~geometry.width_status.str.contains('CENSORED',na=False))]
    gsummary=strict_geometry.groupby(['direction','threshold_db']).width.agg(['min','median','max']).reset_index();f.save('BASIN_WIDTH_SUMMARY.csv',gsummary)
    gates={
        'ACOUSTIC_BASIN_GEOMETRY_CHARACTERIZED':bool(geometry.case_id.nunique()==36),
        'GLOBAL_SEARCH_COVERAGE_CONVERGENCE_VALIDATED':bool(stability.stable_B2_B3_cumulative.all() and stability.raw_components_both_recovered.all()),
        'DEVELOPMENT_MATCHED_BASINS_RECOVERED':bool(len(dev)==36 and dev.recovered.all()),
        'INDEPENDENT_SEARCH_AGREEMENT_CONFIRMED':bool(len(independent)==7 and independent.agreement.all()),
        'FRESH_HOLDOUT_MATCHED_BASINS_RECOVERED':bool(len(hold)==24 and hold.recovered.all()),
        'NO_UNRESOLVED_EXACT_CONTINUOUS_ALIAS_IN_TESTED_CONTROLS':bool(not alias.distinct_joint_exact_match.any())}
    passed=all(gates.values());decision='A1_FIX2_PASS' if passed else 'A1_FIX2_BLOCKED_BY_UNCLOSED_ACOUSTIC_COVERAGE'
    if not gates['NO_UNRESOLVED_EXACT_CONTINUOUS_ALIAS_IN_TESTED_CONTROLS']:decision='A1_FIX2_IDENTIFIABILITY_FINDING_STOP'
    f.L.json_write(out/'R4_A1_FIX2_DECISION.json',dict(decision=decision,gates=gates,overall_progress_percent=0,A1_status='READY_TO_RESUME_PENDING_INDEPENDENT_AUDIT' if passed else 'BLOCKED',recovery_threshold_db=.001,
        convergence_requires='Every nominal development case stable between B2/B3 and raw components recovered at both levels; cumulative retention alone is insufficient',scope={'development_cases':36,'fresh_confirmation_cases':24,'independent_hard_cases':7,'preregistered_independent_cases':5,'supplementary_independent_cases':2},
        failed_development_cases=dev[~dev.recovered].case_id.tolist(),failed_confirmation_cases=hold[~hold.recovered].case_id.tolist()))
    progress=json.loads((f.L.MASTER/'R4_PROGRESS.json').read_text(encoding='utf-8'))
    progress.update(current_stage='A1-FIX2',A1_FIX2=decision,A1_FIX2_gates=gates,overall_progress_percent=0,A1_status='READY_TO_RESUME_PENDING_INDEPENDENT_AUDIT' if passed else 'BLOCKED_BY_UNCLOSED_ACOUSTIC_COVERAGE',
        independent_audit_status='FIX1_CONTINUOUS_FORWARD_RC2_AND_R3_CONTROL_ACCEPTED; FIX2_COMMIT_PENDING_AUDIT',
        blockers=[] if passed else [name for name,good in gates.items() if not good],
        next_recommended_stage='Independent FIX2 audit before full A1 resume' if passed else 'Further acoustic coverage/convergence repair; no full A1/A2/B/P5')
    f.L.json_write(f.L.MASTER/'R4_PROGRESS.json',progress)
    scalar=gsummary[gsummary.threshold_db==.001].copy();scalar.loc[scalar.direction=='r_km',['min','median','max']]*=1000
    source=pd.read_csv(f.prior.OUT/'CASE_LEVEL_RESULTS.csv').query('sigma_deg>0')
    scalar['reference_DE_four_dimensional_spacing_ratio']=np.nan
    # A reference only: N^-1/4 assumes uniform independent 4D draws. Actual DE is correlated.
    normalized_spacing=(2*128*161)**(-.25)
    scales={'r_km':15000.,'v_mps':2.,'psi_deg':30.,'theta_offset_deg':.16}
    for i,row in scalar.iterrows():scalar.loc[i,'reference_DE_four_dimensional_spacing_ratio']=scales[row.direction]*normalized_spacing/row['median']
    scalar.loc[scalar.direction=='r_km','direction']='r_m'
    f.save('BASIN_VS_REFERENCE_SEARCH_SPACING.csv',scalar)
    text=f'''# R4-A1-FIX2 acoustic coverage

**{decision}; R4 progress 0%.** R3, original A1 and FIX1 evidence remain frozen. All original regression and Q holdout cases are development evidence. This repair does not execute a multi-sigma accuracy experiment or claim an engineering bearing tolerance.

## Quantitative landscape

Oracle truth is used only for feasible local landscape characterization, never estimator starts. The local coordinates are r,v,psi with observation-profiled theta plus the truth's fixed profile offset; the independent theta-offset direction changes theta. Coupled SVD directions use declared normalized scales. Every scored diagnostic perturbation is inside the unchanged continuous bearing likelihood. Noiseless full-rank bearing observations define a numerical singleton, and their zero widths are explicitly not acoustic widths.

Uncensored J<0.001 dB cross-section widths (r entries below are metres; other units are m/s or degrees):

{f.prior.md_table(scalar.drop(columns='reference_DE_four_dimensional_spacing_ratio'))}

![Likelihood-feasible basin section widths](BASIN_SECTION_WIDTHS.png)

Full J<0.001/0.01/0.1 dB boundaries, bearing-censored intervals, nuisance-depth interactions and Hessian/sensitivity directions are exported. Normalized sensitivity condition numbers range from 160519 to 666320 (median 339076); the fastest coupled direction is dominated by normalized range (median absolute loading 0.9961), with speed coupling. Independent axis widths are not a four-dimensional rectangular basin volume. The old 42 DE islands all exhausted their budgets. A uniform-draw N^-1/4 spacing calculation is an illustrative reference only, not a probability model or coverage certificate for correlated DE trajectories. Millimetre-scale range sections are orders of magnitude narrower than coarse global proposal spacing. Local refinement is therefore essential, and successful local convergence alone cannot establish that all basins were proposed.

Truth-state spline/exact discrepancies across all 15 old truths are at most 3.72e-6 dB, below the 0.001 dB criterion. This diagnostic rules out interpolation error at those truth points as the cause of the large remaining failures; it does not certify every proposal point.

No profiled-depth label switch occurs at the measured old-panel section walls at these three thresholds. Bearing/state-bound censoring affects 91 of the 504 noisy direction/threshold intervals (including coupled directions), and those limits are not interpreted as acoustic basin boundaries.

## Coverage architecture and budgets

See FIX2_METHOD.md, SEARCH_BUDGET_PREREGISTRATION.json and METHOD_FREEZE.json. Nested deterministic range/radial meshes are continued by observed bearing into theta/tangential speed. All shared-depth branches retain separated proposals, then analytic-Jacobian local refinement and direct-modal polishing run. Lower-budget candidates are preserved cumulatively. Raw component results are exported separately; retaining an earlier success is not by itself evidence of independent fine-resolution convergence. Jacobians agree with finite differences and exact forward features.

Cumulative retention was introduced during development after observing that a finer raw component could discard a basin already recovered at a lower budget. This is an adaptive development change to candidate aggregation, not a preregistered independent convergence success. The preregistered component meshes, starts, tolerances and threshold remain unchanged; raw component caches are preserved. Fresh confirmation uses the final frozen hierarchy and aggregation rule.

The grid covers r=45--60 km and radial target speed=0.94--3 m/s. It samples bearing-profile centers and locally refines all four coordinates; it is not a certified exhaustive partition of the entire four-dimensional likelihood region. Finite depth-branch beams can miss a narrow basin. The finest bank is computed once, and masked nested subsets used for budget experiments. Logical node/start counts and shared-cache runtime are distinguished. The original nine regression ranges happen to lie on the 10 m finest range mesh because those historical truths were specified to two decimals in km; the method does not use their coordinates to position the mesh. The Q panel and fresh full-precision random panel do not share this coincidence, so the regression results alone would be insufficient confirmation.

Nominal development recovery vs budget:

{f.prior.md_table(summary.reset_index())}

![Recovery and exact score versus search budget](SEARCH_BUDGET_CONVERGENCE.png)

Stable B2/B3 cumulative states: {int(stability.stable_B2_B3_cumulative.sum())}/{len(stability)}; raw B2/B3 components both recovered: {int(stability.raw_components_both_recovered.sum())}/{len(stability)}. All case-level exact J, coordinates, circular errors, candidate envelopes, bearing costs/cutoffs, basin counts and budgets are in SEARCH_BUDGET_CONVERGENCE.csv and candidate catalogs.

Development B3 hierarchy recovery: noiseless {int(dev[dev.sigma_deg==0].recovered.sum())}/15, nominal {int(dev[dev.sigma_deg==.1].recovered.sum())}/21, total {int(dev.recovered.sum())}/36. Remaining failed cases: {', '.join(dev[~dev.recovered].case_id)}.

## Independent global family

Five preselected difficult failures plus two additional development failures (P02/Q06, selected after development results and before fresh confirmation) use the same independent Sobol bank and derivative-free Nelder-Mead refinement in radial coordinates, without primary recovered-state initialization. INDEPENDENT_SUPPLEMENT_DESIGN.json makes the adaptive diagnostic selection explicit. Agreement checks exact score, four coordinates, nuisance depth and bearing feasibility.

{f.prior.md_table(independent[['case_id','exact_J','primary_exact_J','agreement','differences']])}

Independent matched-basin recovery is {int((independent.exact_J<.001).sum())}/{len(independent)}; agreement with the primary method is {int(independent.agreement.sum())}/{len(independent)}. P09 and Q06 have independently recovered near-zero basins missed by the primary hierarchy; Q01 shows the reverse failure. These disagreements are search-coverage findings, not evidence of physical aliases.

## Fresh preregistered confirmation

Method/config/code hashes and B3 hierarchy were frozen first. Only then were eight new interior off-grid truths generated from seed 2026100204 and their panel/design hashes frozen before observations. Each truth has noiseless and nominal seeds 430001/430002; no panel editing or threshold changes follow results.

Noiseless recovery {int(hold[hold.sigma_deg==0].recovered.sum())}/8; nominal recovery {int(hold[hold.sigma_deg==.1].recovered.sum())}/16. See FRESH_HOLDOUT_RESULTS.csv for every result; failures remain visible. J<0.001 dB and continuous bearing feasibility are required unchanged. FRESH_BASIN_GEOMETRY.csv contains post-confirmation oracle section diagnostics with the same three thresholds; these measurements never feed the frozen estimator.

## Decision and limits

`{json.dumps(gates)}`

Positive exact J while matched truth remains in the bearing region is a finite-search failure. No physical non-identifiability is inferred from failed optimization. ALIAS_DIAGNOSTIC.csv tests separated exported states against exact joint tolerances (bearing RMS 1e-8 degree and acoustic RMS 1e-6 dB); absence of an exported alias is not a proof of global uniqueness under noise.

R4 remains 0% even if this repair passes. Full A1 can resume only after an accepted repair Gate; unresolved coverage or convergence blocks that route. No A2/nav/SSP/amplitude joint sweep, depth estimator or P5 is opened.

## Reproduction and local validation

`python r4_a1_fix2_coverage.py geometry`, `development`, `pool`, `independent`, `freeze`, `holdout` reproduce the isolated stages. `python r4_a1_fix2_audit.py` reconstructs numerical evidence, checks frozen bytes and method/code hashes, and verifies exact aliases. `python r4_a1_fix2_report.py` generates this report/decision and scientific figures from saved data only. Cache records preserve raw per-case starts, optimizer status and candidates. LOCAL_VALIDATION.md states the actual check counts and scope. Independent research-lead audit is pending.
'''
    (out/'R4_A1_FIX2_REPORT.md').write_text(text,encoding='utf-8');(out/'GPT_SYNC.md').write_text(text,encoding='utf-8')
    validation=out/'LOCAL_VALIDATION.md';base=validation.read_text(encoding='utf-8').split('\n## Supplementary report-stage verification')[0]
    validation.write_text(base+f'\n## Supplementary report-stage verification\n\n{len(supplementary)} additional checks passed: seven independently recomputed solver-agreement flags and 24 fresh observation seed replays. See SUPPLEMENTARY_VALIDATION.csv.\n',encoding='utf-8')
    plots(budgets,geometry,summary)
    print(decision,gates)


def plots(budgets,geometry,summary):
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained');levels=['B1','B2','B3']
    axes[0].plot(levels,summary.loc[levels,'recovery_fraction'],'o-',label='Cumulative mode bank')
    axes[0].plot(levels,summary.loc[levels,'raw_component_recovered']/summary.loc[levels,'n_cases'],'s--',label='Raw component')
    axes[0].set(ylabel='Development nominal recovery fraction',ylim=(-.03,1.03));axes[0].legend(fontsize=8)
    for cid,g in budgets[budgets.sigma_deg==.1].groupby('case_id'):axes[1].plot(g.budget,g.exact_J,'o-',alpha=.6,linewidth=.8,markersize=3)
    axes[1].axhline(.001,color='black',linestyle='--');axes[1].set(yscale='log',ylabel='Best exact profiled J (dB)')
    for ax in axes:ax.grid(alpha=.2);ax.set_xlabel('Preregistered nested search budget')
    fig.suptitle('FIX2 development only: coverage improvement is distinct from convergence',fontsize=11)
    fig.savefig(f.OUT/'SEARCH_BUDGET_CONVERGENCE.png',dpi=160);fig.savefig(f.OUT/'SEARCH_BUDGET_CONVERGENCE.pdf');plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for ax,name,label,scale in zip(axes.ravel(),['r_km','v_mps','psi_deg','theta_offset_deg'],['Range section width (m)','Speed section width (m/s)','Heading section width (deg)','Theta offset section width (deg)'],[1000,1,1,1]):
        for threshold,color in zip([.001,.01,.1],['#23688e','#c46916','#507f31']):
            d=geometry[(geometry.direction==name)&(geometry.threshold_db==threshold)]
            censored=d.width_status.str.contains('CENSORED');positions=np.arange(len(d))
            ax.scatter(positions[~censored],d.width[~censored]*scale,s=13,color=color,label=f'J<{threshold:g} dB')
            ax.scatter(positions[censored],d.width[censored]*scale,s=22,color=color,marker='x')
        ax.set(yscale='log',ylabel=label,xlabel='21 fixed development noise cases');ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle('Oracle diagnostic only: circles = acoustic crossings; crosses = bearing/bounds censoring',fontsize=11)
    fig.savefig(f.OUT/'BASIN_SECTION_WIDTHS.png',dpi=160);fig.savefig(f.OUT/'BASIN_SECTION_WIDTHS.pdf');plt.close(fig)


if __name__=='__main__':report()
