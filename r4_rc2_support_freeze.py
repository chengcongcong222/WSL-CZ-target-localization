"""Freeze archived fallback inputs and conservative support contract before one run."""
import math
import numpy as np
import r4_rc2_support as P
import r4_a1_offgrid_bearing as L
E=P.E;O=P.O
def main():
    assert E.G.git("rev-parse","HEAD")==P.PARENT;assert not O.exists();O.mkdir()
    with np.load(P.SOURCE/"CASE_OBSERVATIONS.npz") as z:
        selected=[str(v) for v in z["case_ids"] if "_s0.1_" in str(v)]
        t=z["times_s"].copy()
    assert len(selected)==21 and np.array_equal(t,np.arange(121)*10)
    xp,yp=L.platform_xy();platform=np.column_stack([xp,yp])
    files=[E.ROOT/n for n in ["r4_rc2_interval.py","r4_rc2_support.py","r4_rc2_support_freeze.py","r4_rc2_support_audit.py","r4_rc2_support_report.py","r4_a1_offgrid_bearing.py","r4_a1_fix_continuous.py","r4_a1_new_dynamic.py","r4_h3_p0.py","r4_h3_g0.py"]]
    files += [P.SOURCE/n for n in ["CASE_OBSERVATIONS.npz","R4_A1_FIX_CONFIG.json","HOLDOUT_TRUTH_PANEL.csv"]]
    files += [L.OUT/n for n in ["R4_A1_CONFIG.json","OFFGRID_TRUTH_PANEL.csv"]]
    for n in ["A1_NEW_DRAWS_A.npz","A1_NEW_DESIGN_FREEZE.json","A1_NEW_ERROR_MODEL.md","A1_NEW_TRUTH_PANEL.csv"]:
        files.append(E.ROOT/"results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC"/n)
    files += [E.ROOT/"results/R4_H3_G1_EXTRACTION_SUPPORT_REVIEW"/n for n in ["SUPPORT_PROVENANCE.csv","FINITE_PROFILE_SUMMARY.csv"]]
    files += [E.ROOT/"results/R4_H3_G0_FINITE_VERTICAL_RESPONSE/FINITE_HORIZONTAL_LABEL_DIAGNOSTIC.csv"]
    binding=[{"path":str(p.relative_to(E.ROOT)).replace("\\","/"),"sha256":E.sha(p)} for p in files]
    risk=21*math.exp(-7)+2*21*121*math.exp(-12.5)
    E.dump(O/"DESIGN_FREEZE.json",{
        "stage":"R4_RC2_OBSERVATION_DEFINED_SUPPORT_PILOT","parent_SHA":P.PARENT,"case_ids":selected,
        "U0_source":"R4_A1_FIX_CONTINUOUS_SEARCH/CASE_OBSERVATIONS.npz:bearing_rad;ALL21saved noisy cases;no new draws",
        "U0_fallback":"P/Q panel;preferred H01/H06/H12 raw observations absent",
        "U_status":"RC2_OBSERVATION_OR_NOISE_CONTRACT_INCOMPLETE;raw bearings and reported nodes not archived;do not regenerate from truths/draws",
        "state_order":["r_km","theta_deg","v_mps","psi_deg"],"low":P.LOW.tolist(),"high":P.HIGH.tolist(),
        "physical_prior":"inherited RC2 bounded near-forward domain;not full ocean/global headings",
        "platform_xy":platform.tolist(),"times_s":t.tolist(),"navigation_and_bias":"EXACT_FROZEN_U0_TRACK;ZERO_NAV_ERROR_ZERO_SYSTEM_BIAS;NOT_DUAL_NODE_ERROR_CONTRACT",
        "noise":"archived0.1deg marginal Gaussian,independent121epoch vector percase;no independence across cases required for union bound",
        "confidence_rule":{"simultaneous_all_case_failure_upper":risk,"simultaneous_conditional_coverage_lower":1-risk,
            "SSE_normalized_upper":195,"per_epoch_absolute_sigma_upper":5,
            "statistical_proof":"Direct Gaussian-square MGF Chernoff at lambda1/6:exp(-195/6)*(3/2)^(121/2)<=exp(-7),using ln(3/2)<=5/12 from alternating log series. Gaussian marginal union:2*21*121exp(-25/2). Sum failure budgets;no cross-case independence",
            "correlated_actual_measurements":"NOT_CERTIFIED_BY_THIS_U0_CONTRACT",
            "model_rounding_slack_rad":1e-10,
            "left_right":"minimum residual to observed and reflected observed;conservative superset of archived signed support;arbitrary sign branch retained;not real HLA extractability certification",
            "periodicity":"declared physical +x sector keeps prediction and obs inside(-1/3,1/3)rad;no +/-pi seam inside this domain;full-domain angular prior not expanded"},
        "interval_certificate":{"arithmetic":"IEEE754binary64 each primitive outward nextafter;coefficients exact rational conversion brackets",
            "pi":"rational Machin16atan(1/5)-4atan(1/239),35term alternating bounds",
            "sin_cos":"Taylor13/12 with exact rational(1/3)^15/15!, (1/3)^14/14! remainders",
            "atan":"Taylor17 with rational(1/3)^19/19 remainder;dx strictlypositive;ratio abs<1/3",
            "inclusion_proof":"monotone endpoint sin,cos extrema;interval Cartesian trajectory;monotone atan endpoints;distance to union is lower bound;sum squares rounded DOWN",
            "rejection":"only certified lower SSE>195+1e-6 OR certified minimum epoch distance>5sigma+1e-6sigma;all other boxes kept or split",
            "independent":"Decimal50digit polynomial/rational remainder reconstruction of EVERYreject certificate and dyadic partition proof;analytic geometry controls;finite control corners are supplementary not the proof"},
        "subdivision":"FIFO breadthfirst;split max(width/registered_resolution);lowest index wins ties;exact midpoint shared by bothclosedchildren;no truth ordering",
        "resolution":P.RESOLUTION.tolist(),"budget":{"evaluations_percase":P.CAP,"wall_percase_s":30,"primary_total_s":660,
            "independent_review_s":900,"overall_s":1800,"delivery_bytes":128000000,"new_KRAKEN":0,"new_FIELD":0,"new_MC":0,"new_bearing_draws":0},
        "unresolved":"allqueue boxes and smallest-resolution leaves remain in outer union;never relabel undecided boxes feasible",
        "branch_records":"initialtheta sign sectors plus crosseszero;not exact connectedcomponent enumeration;outer union covers all possible branches",
        "classification":{"A":"certificates PASS,allcasesnonempty,range width<=1km,speedwidth<=0.2mps;diagnostic tightness only not depth precision",
            "B":"certificatesPASS,eachcase at least one dimension shrinks10% vs physical prior, but range/speed stillbroad",
            "C":"no useful10% contraction or inclusion/partition certificates fail;inclusionfailure also IMPLEMENTATION_INVALID",
            "D":"archive observation/noise contract missing;UalwaysD;no synthetic substitution",
            "budget":"finite-budget unsettled boxes alone do not invalidate a certified outerunion;reported distinctly from convergence"},
        "truth_policy":"constructor has ONLY observations/time/platform;truthfiles read after allsupports complete;finite coverage evaluation only",
        "H3_relation":"H3H01/H06/H12 25nodes aretruth-centered;cannot evaluate same-observation nesting against P/Qfallback;read-only provenance/scale comparison;no candidate injection",
        "H3_C1e":"NOT_EVALUATED","R4_percent":0,"stop":"CommitB,push,remoteverify,STOP;no automatically opened stage","bindings":binding})
    (O/"OBSERVATION_AND_ERROR_CONTRACT.md").write_text(
        "# Observation and error contract\n\n"
        "U0 uses all21 saved0.1deg P/Q bearing arrays,121epochs. Archived relative_TL is not used. Exact platform motion and zero navigation/system bias are inherited, not extended to physical hardware. Physical domain r45..60km,theta-5..5deg,v1..3m/s,heading-15..15deg is a restrictive inherited prior.\n\n"
        "Preferred H01/H06/H12 raw bearings are absent. Dual-node archives contain random draws and truths but not received bearings/reported positions. Reconstructing via scene(truth,...) is forbidden here. U is input-incomplete. No same-scene single/dual quantitative comparison or engineering coverage claim is possible.\n\n"
        "A conservative reflection union retains unresolved left/right sectors even though the stored generator produced signed bearings. It never adds information or removes signed-compatible states. Gaussian coverage is conditional on the archived121-epoch iid vector law, exact navigation, constant motion and physical domain. Cross-case dependence is allowed by a union bound. Real correlated noise/navigation or array sign resolution is not certified.\n",encoding="utf-8")
    (O/"INTERVAL_GEOMETRY_PROOF.md").write_text(
        "# Geometric outer-cover proof\n\n"
        "A35term alternating rational Machin series encloses pi. Degree-to-radian and polynomial coefficients are enclosed by exact Fraction comparison with each binary64 endpoint. Every arithmetic primitive is nextafter-expanded. Sin/cosTaylor remainders bound all registered angles(abs<=15deg<1/3rad);monotonic sine and endpoint/zero cosine extrema enclose each coordinate interval.\n\n"
        "Interval products enclose r*cos(theta)+t*v*cos(psi)-platform_x and the analogous y. The full inherited physical domain has dx>0,abs(dy/dx)<1/3;division and monotone atan17Taylor endpoint/remainder bounds therefore enclose every predicted bearing. The domain lies away from the circular seam. It includes both negative and positive initialtheta. Distance from a bearing interval to the observed/reflected angle union is a LOWER bound on every state's absolute residual. After subtracting registered rounding slack, divisions,squares and sums are rounded DOWN. Only lower bounds exceeding the frozen compatibility thresholds plus security gap are rejected.\n\n"
        "Root domain is covered by two CLOSED midpoint children at every split. Rejected terminals contain no compatible states by the enclosure theorem. Every unprocessed queue box and unresolved leaf is exported. Induction on the dyadic tree proves that the exported union covers the entire compatible set, including all disconnected branches. Finite point/corner controls supplement this proof;they cannot establish it alone. Exact connected components are not enumerated.\n\n"
        "Statistical event is separate:for a121component iid standardized Gaussian vector,the Gaussian-integral MGF gives P(sumZ^2>=195)<=exp(-195/6)*(3/2)^(121/2). The alternating log series gives ln(3/2)<=1/2-1/8+1/24=5/12,so the exponent is at most-7.2916667<-7. Each marginal exceeds5sigma with probability<=2exp(-12.5). Summing over21cases and2541epochs yields conditional coverage>=1-21exp(-7)-5082exp(-12.5). Reflected compatibility only enlarges the set. No observed21-case retention is used to prove this probability statement.\n",encoding="utf-8")
    print("FROZEN21 U0 archived cases;U incomplete",len(binding),"bindings")
if __name__=="__main__":main()
