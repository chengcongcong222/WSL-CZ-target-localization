"""Write review findings after independent reconstruction; no new experiment."""
import json
import numpy as np
import r4_h3_p0_derivative_review as R
E=R.E;G=R.G;O=R.O
def main():
    R.verify();v=E.read(O/"VALIDATION.json");assert v["FAIL"]==0
    d=E.read(O/"REVIEW_DECISION.json");a=G.rows(O/"TRUNCATION_LAW_DIAGNOSTIC.csv")
    ratios=np.array([float(r["fitted_h2_ratio"]) for r in a])
    residuals=np.array([float(r["relative_h2_vector_residual"]) for r in a])
    info=G.rows(O/"EFFECTIVE_INFORMATION_SENSITIVITY.csv")
    grid=G.rows(O/"GRID_NUMERICAL_DIAGNOSTIC.csv")
    canc=G.rows(O/"CANCELLATION_AND_LOW_AMPLITUDE.csv")
    env=G.rows(O/"SSP_ENVIRONMENT_WINDOW_AUDIT.csv")
    marg=G.rows(O/"NOMINAL_MODAL_WINDOW_MARGINS.csv")
    start=E.read(O/"EXECUTION_STARTED.json")
    changes=[float(r["relative_to_D1"]) for r in info if r["method"]=="R1"]
    negative=[r for r in env if r["offset"]=="-1"]
    d.update(independent_validation="PASS",old_2pct_gate_failures_retained=20,
        fitted_h2_ratio_min=float(ratios.min()),fitted_h2_ratio_median=float(np.median(ratios)),fitted_h2_ratio_max=float(ratios.max()),
        h2_vector_residual_max=float(residuals.max()),R1_over_D1_information_min=min(changes),R1_over_D1_information_max=max(changes),
        grid_all_methods_information_max=max(float(r["effective_information_grid_relative"]) for r in grid),
        nominal_mode_min_CLOW_margin=min(float(r["lower_margin"]) for r in marg),
        nominal_modes_within_1mps_CLOW=sum(r["within_1mps_lower_boundary"]=="True" for r in marg),
        nominal_mode_entries=len(marg),SSP_env_files_reviewed=len(negative)*2,design_SHA=start["design_SHA"],
        execution_SHA="COMMIT_CONTAINING_THIS_REPORT")
    E.dump(O/"REVIEW_DECISION.json",d)
    text=f"""# H3-P0 saved depth derivative and modal-window review

Parent: {R.PARENT}
Design SHA: {start['design_SHA']}
Execution SHA: commit containing this report; final local/remote SHA recorded separately after push.

## Frozen outcome
{chr(10).join('- '+s for s in d['classifications'])}

Original P0 remains NUMERICAL_OR_FORMULA_INCOMPLETE. Its source-column 2% Gate still has20/36 failures;126/126 original G0 replay is not revoked. SSP depth information NOT_EVALUATED; R4=0%.

## Three separate diagnostics
All five derivatives and all ten pairwise comparisons were retained across H01/H06/H12, B0/B1/B2, both meshes and both noise models. Effective information includes both registered sigma levels.

Stencil: projected h2 fitted ratio range {ratios.min():.8g}–{ratios.max():.8g}, median {np.median(ratios):.8g}; expected4 under a smooth asymptotic central-difference expansion. Maximum h2 vector residual {residuals.max():.8g}. This is a truncation-pattern diagnostic, not a true-error bound.

R1/R2 maximum projected-direction difference {d['max_R1_R2_direction_relative']:.8g}; maximum effective-information difference {d['max_R1_R2_effective_information_relative']:.8g}. R1 effective information / oldD1 ranges {min(changes):.8g}–{max(changes):.8g}.

Grid: high-order same-method maximum nuisance-profiled direction difference {d['max_high_order_grid_direction_relative']:.8g}; maximum effective-information difference {d['max_high_order_grid_information_relative']:.8g}. Complete same-method grid comparisons appear in GRID_NUMERICAL_DIAGNOSTIC.csv.

Effective information is a conditional local tangent diagnostic. It is not achieved depth accuracy, an SSP-robust bound, a finite-sample estimator result or a global branch certificate. Null directions receive INF scale, never pseudoinverse zero variance.

## Sampling, cancellation and weighting
All source and receiver depths are exact saved samples. Source200m lies at a C-linear SSP interpolation knot; the left/right sound-speed slopes differ. This does not by itself prove derivative inaccuracy, but the smoothness needed for an unqualified fourth-order remainder bound has not been certified. Mode shapes are serialized complex float32, and R1/R2 share samples. No independent analytic derivative or rigorous remainder bound is available.

CANCELLATION_AND_LOW_AMPLITUDE.csv retains per-case minimum pressure, modal cancellation and stencil cancellation for each method. Modal cancellation measures sum of absolute contributions / magnitude of their complex sum; stencil cancellation measures coefficient-weighted pressure magnitude / derivative magnitude. A large ratio identifies sensitivity, not independently measured roundoff error. ABSOLUTE_FLOOR weights are sqrt(2)*|pressure|/REF, so weak-pressure observations contribute less than under RELATIVE whitening; these are frozen design assumptions, not measured SNR.

## Modal-window risk
All52 prepared positive/negative SSP environments and26 nominal environments were read; no perturbed modes were generated. Negative profiles have minimum1499m/s while CLOW=1500m/s; positive minimum1501m/s. CHIGH=1800m/s exceeds the sampled water speeds, but this alone is not a mode-convergence certificate. The narrowest nominal lower-window margin is {d['nominal_mode_min_CLOW_margin']:.8g}m/s; {d['nominal_modes_within_1mps_CLOW']} nominal entries lie within1m/s across both grids/frequencies.

The [official KRAKEN manual](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html) documents C-linear interpolation and exclusion of slower modes by nonzero CLOW. Therefore the present contract lacks a certificate that the ±1m/s SSP parameter difference uses a comparable physically relevant mode family. This is a design risk, not observed loss or predicted perturbed counts.

Any future contract must preregister a common window covering every SSP offset with boundary margin, or independently bound the contribution of excluded modes. It must check boundary/count/contribution convergence and rebuild the nominal provider using that same window. New and old nominal/perturbed providers must not be mixed. No new window is executed or authorized here.

## Independent verification and stop
Independent binary parsing, termwise modal summation, a separate nuisance chart and pivoted QR rebuilt {v['checks']} checks: {v['PASS']} PASS, {v['FAIL']} FAIL. Max pressure discrepancy {v['cold_pressure_relative_max']:.8g}; derivative {v['cold_derivative_relative_max']:.8g}; effective information {v['cold_information_relative_max']:.8g}. The original20 depth-column failures remain reconstructed.

New KRAKEN/FIELD/MC/audio=0. Internal consistency does not certify true derivative accuracy. The current local Fisher environment route stops at ACCURACY_UNCERTIFIED; a redesigned propagation contract is a separate future decision. No automatic SSP solve or source-depth expansion follows. This review contributes no negative physical evidence about SSP-eroded depth information, which has not been evaluated.
"""
    (O/"GPT_SYNC.md").write_text(text,encoding="utf-8")
    (O/"MODAL_WINDOW_RISK.md").write_text(text[text.index("## Modal-window risk"):text.index("## Independent verification")],encoding="utf-8")
    master=E.ROOT/"results/R4_MASTER"
    with (master/"R4_PLAN.md").open("a",encoding="utf-8") as f:
        f.write("\n\n## H3-P0 saved depth and modal-window review\n\n"+", ".join(d["classifications"])+". OriginalP0 unchanged; SSP information NOT_EVALUATED; R4=0%; independent checks PASS. STOP; no further numeric stage opened.\n")
    with (master/"R4_EVIDENCE_LEDGER.csv").open("a",encoding="utf-8",newline="") as f:
        import csv
        csv.writer(f).writerow(["R4_H3_P0_DEPTH_DERIVATIVE_AND_MODAL_WINDOW_REVIEW","../R4_H3_P0_DERIVATIVE_REVIEW/REVIEW_DECISION.json",";".join(d["classifications"]),"SAVED_STENCIL_INTERNAL_DIAGNOSTIC;NO_TRUE_ACCURACY_CERTIFICATE;SSP_NOT_EVALUATED","INDEPENDENT_RECONSTRUCTION_PASS;RESEARCH_LEAD_AUDIT_PENDING;STOP",0])
    E.dump(master/"R4_H3_P0_DERIVATIVE_REVIEW_PROGRESS.json",{"stage":"R4_H3_P0_DEPTH_DERIVATIVE_AND_MODAL_WINDOW_REVIEW","statuses":d["classifications"],"R4_percent":0,"STOP":True})
    E.dump(O/"OUTPUT_MANIFEST.json",{"files":[{"path":str(p.relative_to(O)).replace("\\","/"),"sha256":E.sha(p)} for p in sorted(O.rglob("*")) if p.is_file() and p.name!="OUTPUT_MANIFEST.json"]})
    print(json.dumps(d,indent=2))
if __name__=="__main__":main()
