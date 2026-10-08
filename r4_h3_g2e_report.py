"""Frozen reporting/classification; no propagation or statistical inference."""
import csv,json
from collections import defaultdict
import numpy as np
import r4_h3_g2e as E
import r4_h3_g2e_interface as I

def report():
    o=E.OUT
    rows=list(csv.DictReader((o/"EXTRACTION_BY_CELL.csv").open(encoding="utf-8")))
    stability=list(csv.DictReader((o/"NUMERICAL_EXTRACTION_STABILITY.csv").open(encoding="utf-8")))
    groups=defaultdict(list)
    for r in rows:
        key=tuple(r[x] for x in ["scene","resource","pack","K","noise_model","sigma"])
        groups[key].append(r)
    summary=[]; bias=[]
    for key,rr in groups.items():
        scene,res,pack,ks,noise,sigma=key
        errors=[float(x["after_max_projector_error"]) for x in rr if x["after_max_projector_error"]]
        oracle=[float(x["oracle_max_projector_error"]) for x in rr if x["oracle_max_projector_error"]]
        summary.append(dict(zip(["scene","resource","pack","K","noise_model","sigma"],key))|{
            "registered_realizations":16,
            "accepted_quality_realizations":sum(x["accepted_quality"]=="True" for x in rr),
            "failed_blocks":sum(int(x["failed_blocks"]) for x in rr),
            "blocks":sum(int(x["blocks"]) for x in rr),
            "mean_max_error_detected_only":float(np.mean(errors)) if errors else None,
            "mean_max_oracle_error_detected_only":float(np.mean(oracle)) if oracle else None,
            "failures_not_removed":True})
        ri=list(E.CH).index(res); ki=E.KS.index(int(ks)); ch=E.CH[res]
        p=E.fields(scene,160001)
        for t in range(3):
            for f in E.PACK[pack]:
                projections=[]
                for r in rr:
                    base=f'{scene}_{noise}_s{float(sigma):.2f}_rep{int(r["realization"]):02d}'
                    with np.load(o/"features"/f"{base}_n160001.npz") as z:
                        if z["detected"][ri,ki,t,f]:
                            projections.append(I.projector(z["u"][ri,ki,t,f,:len(ch)]))
                if projections:
                    mean=np.mean(projections,axis=0)
                    b=float(np.linalg.norm(mean-I.projector(p[t,ch,f]),"fro"))
                    v=float(np.mean([np.linalg.norm(x-mean,"fro")**2 for x in projections]))
                else:
                    b=v=None
                bias.append(dict(zip(["scene","resource","pack","K","noise_model","sigma"],key))|{
                    "template_s":[0,600,1200][t],"frequency_hz":int(E.F[f]),
                    "registered_realizations":16,"detected":len(projections),
                    "failed":16-len(projections),"conditional_projector_bias":b,
                    "conditional_projector_variance":v,
                    "scope":"CONDITIONAL_ON_DETECTION;ALL_FAILURE_DENOMINATORS_RETAINED"})
    E.write_csv(o/"EXTRACTION_SUMMARY.csv",summary)
    E.write_csv(o/"NORMALIZED_RESPONSE_BIAS_VARIANCE.csv",bias)
    failure_rate=1-sum(r["accepted_quality"]=="True" for r in rows)/5184
    unstable=sum(int(r["detection_disagreements"])>0 or
                 (bool(r["max_accepted_projector_grid_difference"]) and
                  float(r["max_accepted_projector_grid_difference"])>.05) for r in stability)
    validation=json.loads((o/"VALIDATION.json").read_text(encoding="utf-8"))
    if not validation["PASS"]:
        decision="H3_G2E_STATISTICAL_IMPLEMENTATION_INVALID"
    elif failure_rate>.05 or unstable/5184>.05:
        decision="H3_G2E_COMPRESSED_RESPONSE_EXTRACTION_UNRELIABLE"
    else:
        decision="H3_G2E_CONDITIONAL_FREQUENCY_STATISTIC_ESTABLISHED"
    d={"stage":"R4_H3_G2E_CONDITIONAL_FINITE_SAMPLE_EXTRACTION",
       "classification":decision,"registered_cells":5184,
       "quality_accepted_cells":sum(r["accepted_quality"]=="True" for r in rows),
       "quality_failure_fraction":failure_rate,"grid_unstable_cells":int(unstable),
       "joint_statistical_chain":"PASS" if validation["PASS"] else "INVALID",
       "normalized_response_Fisher_efficiency":"NOT_EVALUATED",
       "C1b_interpretation":"GAIN_DISTORTED_SOURCE_SCALAR_INVARIANT_RESPONSE;NOT_CALIBRATED_FIELD",
       "unknown_source_power":"NO_TRUTH_POWER_USED_BY_ESTIMATOR;FINITE_RANDOM_POWER_AFFECTS_DETECTION",
       "H3_G0":"UNCHANGED","H3_G1_gaps":"UNCHANGED","C1g_free_modal_zero":"RETAINED",
       "C1g_constrained_propagation":"NOT_ESTABLISHED","C1e_environment_pose":"NOT_ESTABLISHED",
       "actual_UUV_source_occupancy":"NOT_ESTABLISHED","TIME_DOMAIN_RECEIVED_SIGNAL_EXTRACTION":"NOT_OPENED",
       "FULL_RC2_DEPTH_SET":"NOT_ESTABLISHED","R4_percent":0,
       "next_stage":"NOT_AUTHORIZED;STOP","scientific_scope":"CONDITIONAL_STATIC_FREQUENCY_STATISTICS_ONLY"}
    E.write_json(o/"H3_G2E_DECISION.json",d)
    meta=json.loads((o/"EXECUTION_METADATA.json").read_text(encoding="utf-8"))
    likes=list(csv.DictReader((o/"RAW_VS_CSD_LIKELIHOOD_AUDIT.csv").open(encoding="utf-8")))
    byk=[]
    for k in ["1","8","32"]:
        xx=[r for r in rows if r["K"]==k]
        byk.append(f'- K={k}: accepted {sum(r["accepted_quality"]=="True" for r in xx)}/{len(xx)}; '
                   f'failed template/frequency blocks {sum(int(r["failed_blocks"]) for r in xx)}; '
                   f'oracle failures {sum(int(r["oracle_failed_blocks"]) for r in xx)}.')
    matched=list(csv.DictReader((o/"RESOURCE_MATCHED_EXTRACTION.csv").open(encoding="utf-8")))
    b1=sum(r["B1_quality"]=="True" for r in matched); b2=sum(r["B2_quality"]=="True" for r in matched)
    diff=[float(r["B1_minus_B2_error"]) for r in matched if r["B1_minus_B2_error"]]
    note=f"""# H3-G2E conditional frequency-statistic synchronization

Parent: 3f01469c4a718413d7368db6c564a269c0ca5ce7.
Design commit: DESIGN_COMMIT_VERIFICATION.json.
Execution commit: the commit containing this report. Remote verification is outside the repository to avoid self-SHA recursion.

Decision: {decision}. R4 remains 0%; submitted then STOP.

The full Gaussian sample/joint-CSD chain passed {validation["checks"]} cold checks, {validation["FAIL"]} FAIL.
All {len(likes)} registered likelihood witnesses preserve every channel and selected cross-frequency block.
Maximum relative NLL difference: {max(float(r["relative_difference"]) for r in likes):.8g}.
Sample CSD is never inverted, ridge is never added and samples are never forced to rank one.
K=1 is a singular sample statistic with a positive-definite MODEL covariance.

{chr(10).join(byk)}

Overall diagnostic quality accepted {d["quality_accepted_cells"]}/5184; grid-unstable cells {unstable}.
Quality means every template/frequency block exceeds the energy gate and unit-projector error <=0.1.
This is a frozen extraction diagnostic, not a new project performance PASS.
The per-block conservative 1% Markov null bound can reject informative samples.
Failure does not establish physical information absence.
Every actual/oracle rejection, before/after error and background-estimation error remains in the records.

Resource-matched B1/B2 accepted {b1}/1728 vs {b2}/1728.
Median paired B1-minus-B2 maximum detected projector error {float(np.median(diff)) if diff else None}.
A response-stability comparison does not demonstrate depth localization.
EXTRACTION_SUMMARY.csv includes every scene/resource/K/noise/package.
NORMALIZED_RESPONSE_BIAS_VARIANCE.csv labels detected-only moments and retains all 16-realization denominators.

Noise is estimated from 256 independent background samples, using full empirical covariance.
Oracle noise is separate. Hann 249/250-Hz correlation is +1/6.
Likelihood witnesses use synthetic non-truth fixed-gain candidates and template-dependent unknown source powers;
they are statistical algebra checks, NOT physical model fits or recovered source/calibration values.
Full CSD is sufficient within the SAME zero-mean Gaussian covariance family.
No compressed-response Fisher efficiency is claimed without a matched finite-sample model.

Unknown source scalar is removed by normalization; unknown fixed channel/frequency gains remain in the observed direction.
Physical propagation response therefore is not identified without external calibration constraints.
Source/gain gauges, separable-component depth absorption, free-modal zero and C2 zero remain explicit.

New stochastic simulation: 48 independent innovation pools (3 geometries x16 realizations), each containing 3 static template draws.
192 scaled noise conditions, K prefixes, resource/frequency masks and two grids are paired derivatives.
5184 diagnostic cells are NOT 5184 independent experiments.
Saved source/noise/background arrays, locked fields and raw Y/background SHA records reconstruct every observation.
No new propagation, localization-error Monte Carlo, time-series extraction or 64-s stationarity evidence was generated.

Execution {meta["execution_seconds"]:.2f}s; execution plus cold audit {validation["execution_plus_audit_seconds"]:.2f}s.
Stored bytes at audit {validation["delivered_bytes_at_audit"]}; cap 1,000,000,000.
Independent algorithms are internal; external research-lead acceptance remains pending.

H3-G0/G1 conclusions unchanged. Actual UUV source occupancy, constrained propagation, environment/pose robustness,
full RC2 support and continuous depth performance NOT_ESTABLISHED. Time-domain extraction NOT_OPENED.
No next module is authorized.

Formula references: [SciPy spectral analysis](https://docs.scipy.org/doc/scipy/tutorial/signal.html#spectral-analysis)
and [Gaussian spectral-likelihood research](https://pmc.ncbi.nlm.nih.gov/articles/PMC10050575/).
"""
    (o/"GPT_SYNC.md").write_text(note,encoding="utf-8")
    manifest=[{"path":str(p.relative_to(E.ROOT)).replace("\\","/"),"sha256":E.sha(p),
               "bytes":p.stat().st_size} for p in sorted(o.rglob("*"))
              if p.is_file() and p.name!="DELIVERABLE_MANIFEST.json"]
    E.write_json(o/"DELIVERABLE_MANIFEST.json",{"artifacts":manifest,
        "total_bytes":sum(x["bytes"] for x in manifest)})
    master=E.ROOT/"results/R4_MASTER"
    with (master/"R4_PLAN.md").open("a",encoding="utf-8") as f:
        f.write(f"\n\n## H3-G2E conditional frequency extraction\n\n{decision}. "
                "Full joint CSD statistical interface audited; compressed response diagnostic only. "
                "All source/propagation/full horizontal-support boundaries retained. R4=0%; STOP.\n")
    with (master/"R4_EVIDENCE_LEDGER.csv").open("a",encoding="utf-8",newline="") as f:
        csv.writer(f).writerow(["R4_H3_G2E_CONDITIONAL_FINITE_SAMPLE_EXTRACTION",
            "../R4_H3_G2E_CONDITIONAL_EXTRACTION/H3_G2E_DECISION.json",decision,
            "CONDITIONAL_GAUSSIAN_STATIC_FREQUENCY_EXTRACTION;NO_DEPTH_ESTIMATION",
            "PENDING_RESEARCH_LEAD_AUDIT;STOP;NO_CREDIT",0])
    E.write_json(master/"R4_H3_G2E_CONDITIONAL_EXTRACTION_PROGRESS.json",d)
    print(json.dumps(d,indent=2))
if __name__=="__main__":
    report()
