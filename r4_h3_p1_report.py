"""Finalize frozen stage record; no further science."""
from collections import Counter
import json
import r4_h3_p1 as P
E=P.E;G=P.G;O=P.O
def rows(n):return G.rows(O/n) if (O/n).exists() else []
def main():
    P.verify();d=E.read(O/"H3_P1_DECISION.json");v=E.read(O/"VALIDATION.json");start=E.read(O/"EXECUTION_STARTED.json")
    checks=rows("PREREGISTERED_CHECKS.csv");loc=rows("SSP_SECANT_LOCALITY.csv");stable=rows("RELATIVE_RESPONSE_NUMERICAL_STABILITY.csv")
    modes=rows("MODE_BOUNDARY_AND_CONTINUITY.csv")
    failed=[r for r in checks if r["PASS"]!="True"]
    summary={"failed_checks_by_type":dict(Counter(r["check"].split(":")[0] for r in failed)),
        "two_step_C1b_relative_max":max((float(r["two_step_C1b_direction_relative"]) for r in loc),default="NOT_EVALUATED"),
        "forward_backward_C1b_relative_max":max((float(r["forward_backward_C1b_relative"]) for r in loc),default="NOT_EVALUATED"),
        "source_removed_phase_max_rad":max((float(r["max_source_removed_principal_phase_rad"]) for r in loc),default="NOT_EVALUATED"),
        "CLOW_margin_min":min((float(r["CLOW_margin_min"]) for r in modes),default="NOT_EVALUATED"),
        "CHIGH_margin_min":min((float(r["CHIGH_margin_min"]) for r in modes),default="NOT_EVALUATED")}
    for kind in ["GRID","WINDOW"]:
        subset=[r for r in stable if r["kind"]==kind]
        summary[kind+"_raw_relative_max"]=max((float(r["raw_field_relative"]) for r in subset),default="NOT_EVALUATED")
        summary[kind+"_source_response_relative_max"]=max((float(r["source_removed_field_relative"]) for r in subset),default="NOT_EVALUATED")
        summary[kind+"_response_noise1_ratio_max"]=max((float(r["response_noise1_ratio_"+n]) for r in subset for n in P.NOISE),default="NOT_EVALUATED")
    d.update(summary=summary,design_SHA=start["design_SHA"],execution_SHA="COMMIT_CONTAINING_THIS_REPORT")
    E.dump(O/"H3_P1_DECISION.json",d)
    text=f"""# H3-P1 unified modal window and local SSP response

Parent: {P.PARENT}
Design SHA: {start['design_SHA']}
Execution SHA: commit containing this report; final verified local/remote SHA stored after push.

## Result and scope
{chr(10).join('- '+s for s in d['classifications'])}

This is a frozen three-frequency numerical screen, not a source-depth/SSP identifiability or depth-accuracy claim. OriginalP0 remains H3_P0_NUMERICAL_OR_FORMULA_INCOMPLETE. True source-depth derivative accuracy remains ACCURACY_UNCERTIFIED. C1e-P0 depth survival NOT_EVALUATED; complete RC2 support NOT_ESTABLISHED; R4=0%.

## Provider
New KRAKEN calls {d['new_KRAKEN']}/60; generated files {d['generated_modal_files']}/60. Every window includes its own new nominal baseline. No old1500m/s field is used as a new-window nominal. Frequencies150/200/250Hz; grids80001/160001; CLOW1490/1480m/s; CHIGH1800m/s; water-speed shifts0,±0.01,±0.02m/s. All other frozen environment fields and13 exact source/receive sample depths are unchanged.

CHIGH has not received an independent truncation certification. Agreement between the two lower windows is evidence only within this tested contract, not a universal propagation completeness theorem. The SSP step sizes are locality probes, not a measured ocean error range.

## Frozen diagnostics
{json.dumps(summary,indent=2)}

All three H01/H06/H12 scenarios and B0/B1/B2 resources, two meshes, two windows and both inherited noise assumptions are retained. Full raw-field, source-removed spatial-projector and fixed-gain-profiled response comparisons appear in RELATIVE_RESPONSE_NUMERICAL_STABILITY.csv. Secants use only the new nominal baseline and are compared at both steps; no source-depth D1 enters them.

C1b removes unknown complex source per frequency/time and complex element/frequency gain fixed across all snapshots. Principal phase is measured on element/reference source-invariant ratios; it is not a certified unwrap or infinitesimal derivative. Low-pressure locations and absolute-noise-floor sensitivity remain in the per-case tables.

## Complete modal sets
MODE_BOUNDARY_AND_CONTINUITY.csv records counts, phase speeds and boundary margins. MODE_MATCHING.csv records sampled-shape Hungarian assignment, correlations and apparent unmatched indices. Matching based on13 depth samples is not an exact identity theorem. Count differences alone do not fail a Gate. No modes are discarded.

WINDOW_FIELD_CONTRIBUTION.csv rebuilds unmatched and low-correlation subsets while retaining the complete full field for all scientific comparisons. Assignment is only diagnostic; SSP secants are differences of full fields. Full-window comparisons, both grids, both secant steps and nuisance-profiled reproducibility jointly determine admission.

## Independent review
{v['checks']} checks; {v['PASS']} PASS; {v['FAIL']} FAIL. Cold pressure relative maximum {v['cold_pressure_relative_max']}. Independent binary parser, Cartesian geometry and termwise modal summation reconstruct every saved field, with a separate nuisance chart and pivoted QR for response and secant checks. Scientific Gate failures are distinct from implementation reconstruction failures.

## Stop
No FIELD calls, MC, receive recordings or continuous depth localization. No automatic lower window, smaller SSP step, new frequency/grid or original±1m/s stress trial. If numerical support is accepted, only a separately designed full-band research-lead review is suggested. Otherwise stop the current propagation derivative route and prioritize a separate RC2-horizontal-support/depth-set design decision. No next stage is opened here.
"""
    (O/"GPT_SYNC.md").write_text(text,encoding="utf-8")
    master=E.ROOT/"results/R4_MASTER"
    with (master/"R4_PLAN.md").open("a",encoding="utf-8") as f:
        f.write("\n\n## H3-P1 unified window local SSP response\n\n"+"; ".join(d["classifications"])+". OriginalP0 unchanged;C1e-P0 NOT_EVALUATED;R4=0%;STOP;no next stage opened.\n")
    import csv
    with (master/"R4_EVIDENCE_LEDGER.csv").open("a",encoding="utf-8",newline="") as f:
        csv.writer(f).writerow(["R4_H3_P1_UNIFIED_MODAL_WINDOW_AND_LOCAL_SSP_RESPONSE","../R4_H3_P1_UNIFIED_WINDOW_SSP/H3_P1_DECISION.json",";".join(d["classifications"]),"THREE_FREQUENCY_UNIFIED_WINDOW_LOCAL_SSP_NUMERICAL_SCREEN_ONLY;NO_SOURCE_DEPTH_SURVIVAL","INDEPENDENT_REVIEW_"+d["independent_review"]+";RESEARCH_LEAD_AUDIT_PENDING;STOP",0])
    E.dump(master/"R4_H3_P1_PROGRESS.json",{"classifications":d["classifications"],"R4_percent":0,"STOP":True})
    E.dump(O/"OUTPUT_MANIFEST.json",{"files":[{"path":str(p.relative_to(O)).replace("\\","/"),"sha256":E.sha(p)} for p in sorted(O.rglob("*")) if p.is_file() and p.name!="OUTPUT_MANIFEST.json"]})
    print(json.dumps(d,indent=2))
if __name__=="__main__":main()
