"""Close out RC2 pilot with all input gaps and unresolved supports retained."""
import json
import r4_rc2_support as P
E=P.E;G=E.G;O=P.O
def main():
    P.verify();d=E.read(O/"RC2_SUPPORT_DECISION.json");v=E.read(O/"VALIDATION.json");f=E.read(O/"DESIGN_FREEZE.json")
    diag=G.rows(O/"BRANCH_AND_BUDGET_DIAGNOSTICS.csv");truth=G.rows(O/"TRUTH_COVERAGE_EVALUATION.csv")
    summary={}
    for name in ["range_km","speed_mps","theta_deg","heading_deg"]:
        values=[float(r[name+"_width"]) for r in diag if r[name+"_width"]!="EMPTY"]
        summary[name+"_width_min"]=min(values) if values else "EMPTY"
        summary[name+"_width_max"]=max(values) if values else "EMPTY"
    summary["budget_exhausted_cases"]=sum(r["budget_exhausted"]=="True" for r in diag)
    summary["retained_boxes"]=sum(int(r["retained_boxes"]) for r in diag)
    summary["truth_compatible"]=sum(r["truth_compatible_event"]=="True" for r in truth)
    summary["truth_retained"]=sum(r["algorithm_kept_truth"]=="True" for r in truth)
    d.update(summary=summary,design_SHA=E.read(O/"EXECUTION_STARTED.json")["design_SHA"],execution_SHA="COMMIT_CONTAINING_THIS_REPORT")
    E.dump(O/"RC2_SUPPORT_DECISION.json",d)
    relation="""# H3 finite25-node diagnostic relationship\n\nH3-G1's finite horizontal nodes are centered on H01/H06/H12 generating truths with fixed angular/heading coordinates. They are neither observations nor admissible RC2 seeds. They were not injected into this constructor.\n\nThis pilot's available U0 observations are instead P/Q historical cases. Preferred H01/H06/H12 raw bearing sequences and reported dual-node tracks were never stored. Hence SAME_OBSERVATION_H3_GRID_NESTING=NOT_EVALUATED; no cross-panel membership counts are presented as H3 coverage. Full H3-C1e depth survival and continuous depth coverage remain NOT_EVALUATED.\n\nThe registered physical prior spans15km in range and2m/s in speed. Pilot output widths are reported percase in BRANCH_AND_BUDGET_DIAGNOSTICS.csv. A broad envelope means future acoustic/depth-set construction must handle wide multi-state supports, rather than the finite truth-centered grid. Width comparison is only a scale diagnostic; it does not certify any depth interval or discard any H3 node.\n"""
    (O/"H3_FINITE_GRID_RELATION.md").write_text(relation,encoding="utf-8")
    text=f"""# RC2 observation-defined multibranch outer support pilot

Parent: {P.PARENT}
Design SHA: {d['design_SHA']}
Execution SHA: commit containing this report;final remote record written after push.

## Architecture-specific decisions
{json.dumps(d['architecture_classifications'],indent=2)}

U0 uses all21 pre-existing noisy P/Q case arrays; not H01/H06/H12. U dual-node is input-incomplete:draws,truths and point estimates are present but raw bearings/reported node positions are absent. No old scene(truth,...) calls or new observation draws filled the gap. Same-scene single/auxiliary comparison is unavailable. Hardware resources remain separate.

## Conservative support and statistical scope
Physical prior is r45..60km,theta-5..5deg,v1..3m/s,heading-15..15deg;exact prescribed1200s platform track,zero navigation/system bias and0.1deg121epoch Gaussian noise are inherited only for U0. Real navigation/systematic/temporal correlation and HLA extractability are not certified.

The compatibility set uses residual squared sum<=195sigma^2 and every residual<=5sigma plus registered model rounding slack. A reflection union retains possible left/right branches as a conservative expansion of the stored signed-angle model. The conditional simultaneous coverage lower bound over all21 fixed cases is {f['confidence_rule']['simultaneous_conditional_coverage_lower']:.10f}, from a self-contained Gaussian-integral Chernoff bound and marginal union bound, not these21 empirical successes. Cross-case independence is unnecessary;within-case iid Gaussian covariance is required by the archived contract.

The interval algorithm receives observations,time,platform only. Directed rounded rational Taylor/Machin enclosures bound Cartesian trajectories and bearings. Only certified incompatibility is discarded. Independent Decimal rejection replay and closed dyadic-tree reconstruction certify exported union inclusion. Finite corners/control tests are supplementary to the geometric proof.

## Widths and budget
{json.dumps(summary,indent=2)}

All budget-unresolved and minimum-resolution leaves are included in SUPPORT_BOXES.csv. These are outer boxes, not proven feasible states. Disconnected branches cannot be lost by the union inclusion argument, but exact connected-component counts remain unresolved. Sign-sector labels are not a global modal count. Exhausting the subdivision budget does not turn a broad outer union into localization failure or physical non-observability.

Finite evaluation truth was read only after all support construction. Truth compatibility,actual retention and any compatible-but-deleted event are separate columns. Small-panel retention is not a sea-trial coverage certificate.

## Independent validation and H3 boundary
{v['checks']} checks,{v['PASS']} PASS,{v['FAIL']} FAIL;Decimal rejection certificates {v['Decimal_rejection_certificates']}. Review scope:full partition coverage and every discarded-box inclusion certificate,analytic co-linear/turn/auxiliary geometry controls and supplemental corner tests.

H3 finite25-node comparison cannot be certified across different panels. Those truth-centered nodes were not candidate inputs. H3-G0 conditional local information remains accepted;G1 unresolved,G2E compressed extraction unreliable,P0/P1 unchanged,C1e NOT_EVALUATED. No acoustic/depth or joint5D result was generated.

## Stop
New KRAKEN/FIELD/MC/noise draws/recordings=0;R4=0%. No automatic SSP,matching,time-domain,depth-MC or joint metric stage. U needs separately authorized archiving of actual bearings and reported tracks before a matching support trial. U0 broad/unfinished supports require a separately reviewed tightening or support/depth-set design, never truth supplementation or deletion of unsettled branches.
"""
    (O/"GPT_SYNC.md").write_text(text,encoding="utf-8")
    master=E.ROOT/"results/R4_MASTER"
    with (master/"R4_PLAN.md").open("a",encoding="utf-8") as fh:fh.write("\n\n## RC2 observation-defined outer support\n\n"+str(d["architecture_classifications"])+". HistoricalP/Q fallback only;U raw observation gap;R4=0%;STOP.\n")
    import csv
    with (master/"R4_EVIDENCE_LEDGER.csv").open("a",encoding="utf-8",newline="") as fh:
        csv.writer(fh).writerow(["R4_RC2_OBSERVATION_DEFINED_SUPPORT_PILOT","../R4_RC2_OBSERVATION_SUPPORT/RC2_SUPPORT_DECISION.json",";".join(d["architecture_classifications"].values()),"ARCHIVED_U0_PQ_CONDITIONAL_OUTER_UNION;DUAL_OBSERVATIONS_MISSING;NO_DEPTH_CREDIT","INDEPENDENT_REVIEW_"+d["independent_review"]+";RESEARCH_LEAD_AUDIT_PENDING;STOP",0])
    E.dump(master/"R4_RC2_SUPPORT_PROGRESS.json",{"classifications":d["architecture_classifications"],"R4_percent":0,"STOP":True})
    E.dump(O/"OUTPUT_MANIFEST.json",{"files":[{"path":str(p.relative_to(O)).replace("\\","/"),"sha256":E.sha(p)} for p in sorted(O.rglob("*")) if p.is_file() and p.name!="OUTPUT_MANIFEST.json"]})
    print(json.dumps(d,indent=2))
if __name__=="__main__":main()
