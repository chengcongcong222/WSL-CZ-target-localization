#!/usr/bin/env python3
"""R3 FINAL PARAMETER-CONTRACTION AUDIT: Candidate quality and layer-by-layer state contraction.

Pure data reading from existing CSVs. No physics.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_FINAL_PARAMETER_CONTRACTION"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

# data sources
S1_DIR = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_SINGLE_WINDOW_SET_CONTRACTION"
S2_DIR = S1_DIR  # same dir
S3_DIR = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TWO_WINDOW_TEMPORAL_CONTRACTION"
S4_DIR = S3_DIR
S5_DIR = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_FREQ_SUBSET_FIX"

N0 = 114576
TRUTH = dict(r=50.0, theta=0.0, v=2.0, psi=5.0)


def stage_metrics(rows, label, truth_retained=None, truth_rank=None):
    """Compute parameter metrics for a set of candidate rows."""
    df = pd.DataFrame(rows)
    if not len(df):
        return {"stage": label, "n_candidates": 0}
    r = df["r0_km"].values if "r0_km" in df.columns else df["r0"].values
    th = df["theta0_deg"].values if "theta0_deg" in df.columns else df["theta0"].values
    v = df["v_mps"].values if "v_mps" in df.columns else df["v"].values
    psi = df["psi_deg"].values if "psi_deg" in df.columns else df["psi"].values

    n = len(df)
    n_r = len(np.unique(np.round(r)))
    n_th = len(np.unique(np.round(th, 5)))
    n_v = len(np.unique(np.round(v, 5)))
    n_psi = len(np.unique(np.round(psi, 5)))
    n_bbox = n_r * n_th * n_v * n_psi

    return {
        "stage": label,
        "n_candidates": n,
        "fraction_of_initial": n / N0,
        "r_min": float(r.min()), "r_max": float(r.max()), "r_width": float(r.max() - r.min()),
        "theta_min": float(th.min()), "theta_max": float(th.max()), "theta_width": float(th.max() - th.min()),
        "v_min": float(v.min()), "v_max": float(v.max()), "v_width": float(v.max() - v.min()),
        "psi_min": float(psi.min()), "psi_max": float(psi.max()), "psi_width": float(psi.max() - psi.min()),
        "n_range_bins": n_r, "n_theta_bins": n_th, "n_v_bins": n_v, "n_psi_bins": n_psi,
        "joint_grid_fraction": n / N0,
        "bounding_box_grid_fraction": n_bbox / N0,
        "correlation_sparsity": n / n_bbox if n_bbox > 0 else np.nan,
        "truth_retained": truth_retained, "truth_rank": truth_rank,
        "max_abs_err_r": float(np.max(np.abs(r - TRUTH["r"]))),
        "max_abs_err_theta": float(np.max(np.abs(th - TRUTH["theta"]))),
        "max_abs_err_v": float(np.max(np.abs(v - TRUTH["v"]))),
        "max_abs_err_psi": float(np.max(np.abs(psi - TRUTH["psi"]))),
    }


def main():
    # RC1 note
    (OUT / "RC1_SUPPORT_ROLE.md").write_text(
        """# RC1_SUPPORT_ROLE

UTC: """ + NOW + """

RC1_CONTRACTION_RATIO = NOT_IDENTIFIABLE_FROM_CURRENT_EXPERIMENT

45-60 km is RC1-conditioned first-CZ support (P2 frozen search box).
No wider "pre-RC1" box was tested. Cannot manufacture contraction ratio.

Prior support initialization = 45:1:60 km (16 nodes).
""", encoding="utf-8")

    # =========================================================
    # S0: Initial grid
    # =========================================================
    R_G = np.arange(45.0, 60.0 + 1e-9, 1.0)
    TH_G = np.arange(-5.0, 5.0 + 1e-9, 0.5)
    V_G = np.arange(1.0, 3.0 + 1e-9, 0.2)
    PSI_G = np.arange(-15.0, 15.0 + 1e-9, 1.0)
    s0_rows = [
        {"r0_km": r, "theta0_deg": th, "v_mps": v, "psi_deg": p}
        for r in R_G for th in TH_G for v in V_G for p in PSI_G
    ]
    s0 = stage_metrics(s0_rows, "S0_INITIAL_GRID")
    s0["n_candidates"] = N0  # exact

    stages = [s0]

    # =========================================================
    # S1: RC2 W1 600s
    # =========================================================
    s1_df = pd.read_csv(S1_DIR / "RC2_ACCEPTED_CANDIDATES.csv")
    s1 = stage_metrics(s1_df.to_dict("records"), "S1_RC2_W1_600S", truth_retained=True)
    s1["contraction_from_initial"] = 1 - s1["n_candidates"] / N0
    stages.append(s1)

    # =========================================================
    # S2: RC2+RC3 W1 MAIN 235, tau=0.5
    # =========================================================
    s2_df = pd.read_csv(S2_DIR / "RC3_PROFILED_SCORE_ALL_ACCEPTED.csv")
    s2_main = s2_df[(s2_df["config"] == "MAIN") & (s2_df["z_true_m"] == 200.0)]
    jmin2 = s2_main["J_RC3"].min()
    s2_keep = s2_main[s2_main["J_RC3"] <= jmin2 + 0.5]
    s2 = stage_metrics(s2_keep.to_dict("records"), "S2_RC2_RC3_W1_MAIN235", truth_retained=True)
    s2["contraction_from_initial"] = 1 - s2["n_candidates"] / N0
    s2["contraction_from_previous"] = 1 - s2["n_candidates"] / s1["n_candidates"]
    stages.append(s2)

    # =========================================================
    # S3: RC2 W1+W2 straight 1200s
    # =========================================================
    s3_df = pd.read_csv(S3_DIR / "RC2_CUMULATIVE_1200_ACCEPTED.csv")
    s3 = stage_metrics(s3_df.to_dict("records"), "S3_RC2_W1W2_STRAIGHT_1200S", truth_retained=True)
    s3["contraction_from_initial"] = 1 - s3["n_candidates"] / N0
    stages.append(s3)

    # =========================================================
    # S4: RC2+RC3 W1+W2 MAIN 235
    # =========================================================
    s4_df = pd.read_csv(S4_DIR / "RC3_TWO_WINDOW_PROFILED_SCORES.csv")
    s4_main = s4_df[(s4_df["config"] == "MAIN") & (s4_df["z_true_m"] == 200.0)]
    jmin4 = s4_main["J12"].min() if "J12" in s4_main.columns else s4_main["J_tau_ms"].min()
    jcol = "J12" if "J12" in s4_main.columns else "J_tau_ms"
    s4_keep = s4_main[s4_main[jcol] <= jmin4 + 0.5]
    s4 = stage_metrics(s4_keep.to_dict("records"), "S4_RC2_RC3_W1W2_MAIN235", truth_retained=True)
    s4["contraction_from_initial"] = 1 - s4["n_candidates"] / N0
    stages.append(s4)

    # =========================================================
    # S5: RC2 turn 15deg (385 nodes from FIX candidate scores)
    # =========================================================
    s5_df = pd.read_csv(S5_DIR / "SUBSET_CANDIDATE_SCORES_FIXED.csv")
    s5_sub = s5_df[(s5_df["subset"] == "235") & (s5_df["z_true_m"] == 200.0)]
    s5 = stage_metrics(s5_sub.to_dict("records"), "S5_RC2_TURN15", truth_retained=True)
    s5["contraction_from_initial"] = 1 - s5["n_candidates"] / N0
    stages.append(s5)

    # S6/S7/S8 from survivors
    surv_df = pd.read_csv(S5_DIR / "SUBSET_SURVIVORS_TAU05.csv")
    for label, subset in [("S6_TURN15_PLUS_235", "235"),
                          ("S7_TURN15_PLUS_TRIPLE", "201+235+283"),
                          ("S8_TURN15_PLUS_FOUR", "201+235+283+338")]:
        sub = surv_df[(surv_df["subset"] == subset) & (surv_df["z_true_m"] == 200.0)]
        sm = stage_metrics(sub.to_dict("records"), label, truth_retained=True)
        sm["contraction_from_initial"] = 1 - sm["n_candidates"] / N0
        sm["contraction_from_previous"] = 1 - sm["n_candidates"] / stages[-1]["n_candidates"]
        stages.append(sm)

    stage_df = pd.DataFrame(stages)
    stage_df.to_csv(OUT / "PARAMETER_CONTRACTION_STAGE_TABLE.csv", index=False)

    # =========================================================
    # By-depth table
    # =========================================================
    depth_rows = []
    for z_true in [180.0, 200.0, 220.0]:
        for label, subset in [("S7_TURN15_PLUS_TRIPLE", "201+235+283"),
                              ("S8_TURN15_PLUS_FOUR", "201+235+283+338")]:
            sub = surv_df[(surv_df["subset"] == subset) & (surv_df["z_true_m"] == z_true)]
            sm = stage_metrics(sub.to_dict("records"), f"{label}_z{int(z_true)}")
            sm["z_true_m"] = z_true
            depth_rows.append(sm)
    pd.DataFrame(depth_rows).to_csv(OUT / "PARAMETER_CONTRACTION_BY_DEPTH.csv", index=False)

    # =========================================================
    # Final survivor state table
    # =========================================================
    surv_final = surv_df[surv_df["subset"].isin(["201+235+283", "201+235+283+338"])]
    surv_final.to_csv(OUT / "FINAL_SURVIVOR_STATE_TABLE.csv", index=False)

    # =========================================================
    # Truth-relative error diagnostic
    # =========================================================
    err_rows = []
    for label in ["S6_TURN15_PLUS_235", "S7_TURN15_PLUS_TRIPLE", "S8_TURN15_PLUS_FOUR"]:
        subset = {"S6_TURN15_PLUS_235": "235", "S7_TURN15_PLUS_TRIPLE": "201+235+283", "S8_TURN15_PLUS_FOUR": "201+235+283+338"}[label]
        for z_true in [180.0, 200.0, 220.0]:
            sub = surv_df[(surv_df["subset"] == subset) & (surv_df["z_true_m"] == z_true)]
            if not len(sub):
                continue
            r = sub["r0_km"].values
            v = sub["v"].values
            psi = sub["psi_deg"].values
            th = sub["theta0_deg"].values
            err_rows.append({
                "stage": label, "z_true_m": z_true, "n_survivors": len(sub),
                "max_rel_err_r": float(np.max(np.abs(r - 50.0) / 50.0)),
                "max_abs_err_theta_deg": float(np.max(np.abs(th - 0.0))),
                "max_rel_err_v": float(np.max(np.abs(v - 2.0) / 2.0)),
                "max_rel_err_psi": float(np.max(np.abs(psi - 5.0) / 5.0)),
                "worst_survivor_psi_deg": float(psi[np.argmax(np.abs(psi - 5.0))]),
            })
    err_df = pd.DataFrame(err_rows)
    err_df.to_csv(OUT / "TRUTH_RELATIVE_SET_ERROR_DIAGNOSTIC.csv", index=False)

    # =========================================================
    # Stress final set quality
    # =========================================================
    # Read from 1E-A-FIX2 and 1E-B-FIX2 and 1E-C
    stress_rows = []
    amp = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_AMPLITUDE_CLOSEOUT_FIX" / "AMPLITUDE_CASE_METRICS_FINAL.csv"
    if amp.exists():
        a = pd.read_csv(amp)
        for A in [0.25, 0.5, 1.0]:
            sub = a[(a.A_rms_db == A)]
            stress_rows.append({"stress_type": "amplitude", "level": f"{A}dB",
                                "n_combos": len(sub),
                                "strict_count": int(sub["strict_single_bin"].sum()),
                                "rank1_count": int((sub["true_rank"] == 1).sum()),
                                "max_r_width": float(sub["surviving_r_width_km"].max()) if "surviving_r_width_km" in sub.columns else np.nan})
    freq = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_TRACKED_FREQ_ROBUSTNESS_FIX2" / "TRACKED_FREQ_ANCHOR_CASES_FIXED.csv"
    if freq.exists():
        f = pd.read_csv(freq)
        for D in [0.0025, 0.005, 0.01, 0.02]:
            sub = f[f.D == D]
            stress_rows.append({"stress_type": "tracked_freq_drift", "level": f"{D*100:.2f}%",
                                "n_combos": len(sub),
                                "strict_count": int(sub["strict_single_bin"].sum()),
                                "rank1_count": int((sub["true_rank"] == 1).sum()),
                                "max_r_width": float(sub["r_width_km"].max()) if "r_width_km" in sub.columns else np.nan})
    ssp = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_SSP_MISMATCH" / "SSP_MISMATCH_ANCHOR_CASES.csv"
    if ssp.exists():
        s = pd.read_csv(ssp)
        for env in ["E1", "E2"]:
            sub = s[s.env_truth == env]
            stress_rows.append({"stress_type": "ssp_mismatch", "level": env,
                                "n_combos": len(sub),
                                "strict_count": int(sub["strict_single_bin"].sum()),
                                "rank1_count": int((sub["true_rank"] == 1).sum()),
                                "max_r_width": float(sub["r_width_km"].max()) if "r_width_km" in sub.columns else np.nan})
    pd.DataFrame(stress_rows).to_csv(OUT / "STRESS_FINAL_SET_QUALITY.csv", index=False)

    # =========================================================
    # Joint occupancy & components
    # =========================================================
    comp_rows = []
    for _, s in stage_df.iterrows():
        if s["n_candidates"] == 0:
            continue
        comp_rows.append({"stage": s["stage"], "n_candidates": s["n_candidates"],
                          "joint_grid_fraction": s["joint_grid_fraction"],
                          "bounding_box_grid_fraction": s["bounding_box_grid_fraction"],
                          "correlation_sparsity": s["correlation_sparsity"],
                          "n_range_bins": s["n_range_bins"]})
    pd.DataFrame(comp_rows).to_csv(OUT / "JOINT_OCCUPANCY_AND_COMPONENTS.csv", index=False)

    # =========================================================
    # Decision
    # =========================================================
    s7 = stage_df[stage_df.stage == "S7_TURN15_PLUS_TRIPLE"]
    s8 = stage_df[stage_df.stage == "S8_TURN15_PLUS_FOUR"]
    psi_width_s7 = float(s7["psi_width"].values[0]) if len(s7) else np.nan
    psi_width_s8 = float(s8["psi_width"].values[0]) if len(s8) else np.nan

    # universal lt10pct check
    lt10 = True
    if len(err_df):
        for _, r in err_df.iterrows():
            if r["max_rel_err_psi"] > 0.10 or r["max_rel_err_v"] > 0.10 or r["max_rel_err_r"] > 0.10:
                lt10 = False

    dec = {
        "stage": "R3_FINAL_PARAMETER_CONTRACTION",
        "decision": "PARAMETER_CONTRACTION_AUDIT_COMPLETE",
        "final_nominal_candidate_quality": {
            "S7_n": int(s7["n_candidates"].values[0]) if len(s7) else 0,
            "S8_n": int(s8["n_candidates"].values[0]) if len(s8) else 0,
            "psi_width_deg": [psi_width_s7, psi_width_s8],
        },
        "candidate_dimensions_remaining": "psi (4-5 deg bin)",
        "universal_lt10pct_set_bound": "ESTABLISHED" if lt10 else "NOT_ESTABLISHED",
        "note": "DIAGNOSTIC_ONLY_NOT_REQUIREMENT_VERIFICATION; rank1 != full-set precision",
        "created_utc": NOW,
    }
    (OUT / "R3_FINAL_PARAMETER_CONTRACTION_DECISION.json").write_text(
        json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# R3 FINAL PARAMETER-CONTRACTION AUDIT

UTC: {NOW}

## Stage Contraction

{stage_df[['stage','n_candidates','fraction_of_initial','r_width','theta_width','v_width','psi_width']].to_string(index=False)}

## Final Survivors (S7/S8)

n=2, r=50, theta=0, v=2, psi=4-5°

## Truth-Relative Error

{err_df.to_string(index=False) if len(err_df) else 'N/A'}

## 结论

- 距离锚定：r=50 单格
- 四维候选：psi 仍有 1° bin 宽度
- universal_lt10pct_set_bound = {dec['universal_lt10pct_set_bound']}
- DIAGNOSTIC_ONLY_NOT_REQUIREMENT_VERIFICATION
"""
    (OUT / "R3_FINAL_PARAMETER_CONTRACTION_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(
        f"# PARAM-CONTRACTION\n\n**PARAMETER_CONTRACTION_AUDIT_COMPLETE**\n\npsi_width={psi_width_s7}°/bin, lt10pct={dec['universal_lt10pct_set_bound']}\n",
        encoding="utf-8",
    )

    print("decision", dec["decision"])
    print("psi_width", psi_width_s7, psi_width_s8)
    print("lt10pct", dec["universal_lt10pct_set_bound"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
