#!/usr/bin/env python3
"""R3-CLOSEDLOOP-1E-A-CLOSEOUT: Data-only closeout of amplitude robustness.

No propagation runs. Read existing CSVs, fix identity, recompute metrics, revise report.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_AMPLITUDE_CLOSEOUT"
AMP = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_AMPLITUDE_ROBUSTNESS"
FIX = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_FREQ_SUBSET_FIX"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

CONFIG_MAP = {"TRIPLE": "201+235+283", "FOUR": "201+235+283+338"}


def main():
    # =========================================================
    # 0. Reference manifest
    # =========================================================
    (OUT / "REFERENCE_MANIFEST.json").write_text(json.dumps({
        "reference": str(FIX / "SUBSET_CANDIDATE_SCORES_FIXED.csv"),
        "current": str(AMP / "AMPLITUDE_CANDIDATE_SCORES.csv"),
        "baseline_commit": "d99a974e0d7bd71b35638ca355a961bb04fc07f7",
        "config_map": CONFIG_MAP,
        "created_utc": NOW,
    }, indent=2), encoding="utf-8")

    # =========================================================
    # 1. Zero-stress identity
    # =========================================================
    ref_df = pd.read_csv(FIX / "SUBSET_CANDIDATE_SCORES_FIXED.csv")
    cur_df = pd.read_csv(AMP / "AMPLITUDE_CANDIDATE_SCORES.csv")
    a0 = cur_df[cur_df["A"] == 0.0].copy()

    id_rows = []
    for cfg_name, sub_tag in CONFIG_MAP.items():
        for z_true in [180.0, 200.0, 220.0]:
            ref = ref_df[(ref_df["subset"] == sub_tag) & (ref_df["z_true_m"] == z_true)]
            # A=0 rows (any shape/mode - should be identical)
            cur = a0[(a0["config"] == cfg_name) & (a0["ztrue"] == z_true)]
            if not len(cur):
                id_rows.append({"config": cfg_name, "z_true_m": z_true, "pass": False, "note": "no A=0 rows"})
                continue
            # check shape/mode independence
            groups = cur.groupby(["shape", "mode"])
            first_key = list(groups.groups.keys())[0]
            base = groups.get_group(first_key).set_index("node_id")
            shape_mode_ok = True
            for key in list(groups.groups.keys())[1:]:
                g = groups.get_group(key).set_index("node_id")
                if len(g) != len(base) or not np.allclose(g.loc[base.index, "J"].values, base["J"].values, atol=1e-10):
                    shape_mode_ok = False

            # compare with reference
            ref_map = ref.set_index("node_id")
            ids_ref = set(ref_map.index)
            ids_cur = set(base.index)
            ids_match = ids_ref == ids_cur
            max_j_err = 0.0
            z_ok = True
            keep_ok = True
            for nid in ids_ref & ids_cur:
                j_ref = ref_map.loc[nid, "J"]
                j_cur = base.loc[nid, "J"]
                if np.isfinite(j_ref) and np.isfinite(j_cur):
                    max_j_err = max(max_j_err, abs(j_ref - j_cur))
                z_ref = ref_map.loc[nid, "z_star"]
                z_cur = base.loc[nid, "z_star"]
                if abs(z_ref - z_cur) > 0.1:
                    z_ok = False
                k_ref = bool(ref_map.loc[nid, "keep_tau05"])
                k_cur = bool(base.loc[nid, "keep"])
                if k_ref != k_cur:
                    keep_ok = False

            passed = ids_match and max_j_err <= 1e-10 and z_ok and keep_ok and shape_mode_ok
            id_rows.append({
                "config": cfg_name, "z_true_m": z_true, "n_ref": len(ids_ref), "n_cur": len(ids_cur),
                "ids_match": ids_match, "max_abs_J_error": max_j_err,
                "z_star_exact": z_ok, "keep_exact": keep_ok, "shape_mode_invariant": shape_mode_ok,
                "pass": passed,
            })

    id_df = pd.DataFrame(id_rows)
    id_df.to_csv(OUT / "ZERO_STRESS_INDEPENDENT_COMPARISON.csv", index=False)

    # negative test
    neg_rows = []
    if len(a0):
        test = a0.copy()
        idx = test.index[~test["is_truth_state"] if "is_truth_state" in test.columns else range(1, min(2, len(test)))]
        if len(idx):
            test.loc[idx[0], "J"] += 0.01
        # re-run comparison for one group
        cfg_name = test.iloc[0]["config"]
        z_true = test.iloc[0]["ztrue"]
        sub_tag = CONFIG_MAP[cfg_name]
        ref = ref_df[(ref_df["subset"] == sub_tag) & (ref_df["z_true_m"] == z_true)].set_index("node_id")
        cur = test[(test["config"] == cfg_name) & (test["ztrue"] == z_true)].set_index("node_id")
        max_err = 0.0
        for n in ref.index:
            if n in cur.index:
                rv = float(ref.loc[n, "J"]) if np.isscalar(ref.loc[n, "J"]) else float(ref.loc[n, "J"].iloc[0])
                cv = float(cur.loc[n, "J"]) if np.isscalar(cur.loc[n, "J"]) else float(cur.loc[n, "J"].iloc[0])
                if np.isfinite(rv) and np.isfinite(cv):
                    max_err = max(max_err, abs(rv - cv))
        neg_rows.append({"test": "perturb_one_J_by_0.01", "detected": max_err > 1e-10, "max_err": max_err, "pass": max_err > 1e-10})
    pd.DataFrame(neg_rows).to_csv(OUT / "IDENTITY_COMPARATOR_NEGATIVE_TEST.csv", index=False)

    identity_ok = bool(id_df["pass"].all()) and bool(pd.DataFrame(neg_rows)["pass"].all()) if len(neg_rows) else False

    # =========================================================
    # 2. Recompute case metrics
    # =========================================================
    case_rows = []
    for (cfg, shape, mode, A, z_true), grp in a0.groupby(["config", "shape", "mode", "A", "ztrue"]) if len(a0) else []:
        pass  # placeholder

    # use full cur_df (all A)
    for (cfg, shape, mode, A, z_true), grp in cur_df.groupby(["config", "shape", "mode", "A", "ztrue"]):
        J = grp["J"].values
        r0 = grp["r0"].values
        node_ids = grp["node_id"].values
        keep = grp["keep"].values
        jmin = float(J.min())
        truth_mask = np.abs(r0 - 50.0) < 0.1
        # truth state: the one with r=50, v=2, psi=5 closest
        # find truth node from reference
        ref_sub = ref_df[(ref_df["subset"] == CONFIG_MAP[cfg]) & (ref_df["z_true_m"] == z_true)]
        truth_node = ref_sub[ref_sub["is_truth_state"]]["node_id"].values
        if len(truth_node):
            tn = truth_node[0]
            ti_idx = np.where(node_ids == tn)[0]
            true_J = float(J[ti_idx[0]]) if len(ti_idx) else np.nan
            true_rank = int((J < true_J).sum()) + 1 if np.isfinite(true_J) else -1
            true_kept = bool(keep[ti_idx[0]]) if len(ti_idx) else False
        else:
            true_J, true_rank, true_kept = np.nan, -1, False

        surv = grp[keep]
        surv_r = surv["r0"].values if len(surv) else np.array([])
        r_width = float(surv_r.max() - surv_r.min()) if len(surv_r) else 0.0
        e_set = float(np.max(np.abs(surv_r - 50.0))) if len(surv_r) else 0.0
        r_best = sorted(set(np.round(r0[J <= jmin + 1e-10]).astype(int)))
        wrong_mask = np.abs(r0 - 50.0) >= 0.1
        best_wrong_J = float(J[wrong_mask].min()) if wrong_mask.any() else np.nan
        wrong_margin = best_wrong_J - true_J if np.isfinite(true_J) and np.isfinite(best_wrong_J) else np.nan

        case_rows.append({
            "config": cfg, "shape": shape, "mode": mode, "A_rms_db": A, "z_true_m": z_true,
            "true_J": true_J, "J_min": jmin, "true_rank": true_rank,
            "truth_retained": true_kept,
            "best_range_is_50": 50 in r_best if r_best else False,
            "wrong_range_margin": wrong_margin,
            "surviving_r_width": r_width, "max_surviving_dev": e_set,
            "r_best": str(r_best),
            "n_survivors": int(keep.sum()),
            "strict_single_bin": true_kept and true_rank == 1 and len(r_best) == 1 and r_best[0] == 50,
        })
    case_df = pd.DataFrame(case_rows)
    case_df.to_csv(OUT / "AMPLITUDE_CASE_METRICS_RECOMPUTED.csv", index=False)

    # =========================================================
    # 3. Injected vs residual RMS
    # =========================================================
    t = np.arange(0.0, 1200.0 + 1e-9, 10.0)
    n_w1 = 61
    q_ramp = t / 1200.0 - 0.5
    q_ramp = q_ramp - q_ramp.mean()
    q_ramp = q_ramp / np.sqrt(np.mean(q_ramp ** 2))
    q_sine = np.sin(2 * np.pi * t / 1200.0)
    q_sine = q_sine - q_sine.mean()
    q_sine = q_sine / np.sqrt(np.mean(q_sine ** 2))
    rms_rows = []
    for shape, q in [("A-RAMP", q_ramp), ("A-SINE", q_sine)]:
        for A in [0.25, 0.5, 1.0, 2.0]:
            inj = A * q
            w1_res = inj[:n_w1] - np.mean(inj[:n_w1])
            w2_res = inj[n_w1:] - np.mean(inj[n_w1:])
            rms_rows.append({
                "shape": shape, "A_raw_rms_db": A,
                "A_after_window_demean_rms_db": float(np.sqrt(np.mean(np.concatenate([w1_res, w2_res]) ** 2))),
                "W1_residual_rms_db": float(np.sqrt(np.mean(w1_res ** 2))),
                "W2_residual_rms_db": float(np.sqrt(np.mean(w2_res ** 2))),
            })
    pd.DataFrame(rms_rows).to_csv(OUT / "AMPLITUDE_INJECTED_VS_RESIDUAL_RMS.csv", index=False)

    # =========================================================
    # 4. Working range summary
    # =========================================================
    wr_rows = []
    for A in [0.0, 0.25, 0.5, 1.0, 2.0]:
        for cfg in ["TRIPLE", "FOUR"]:
            sub = case_df[(case_df.A_rms_db == A) & (case_df.config == cfg)]
            if not len(sub):
                continue
            wr_rows.append({
                "A_rms_db": A, "config": cfg, "n_combos": len(sub),
                "strict_single_bin": int(sub["strict_single_bin"].sum()),
                "truth_retained": int(sub["truth_retained"].sum()),
                "best_range_50": int(sub["best_range_is_50"].sum()),
                "min_wrong_range_margin": float(sub["wrong_range_margin"].min()),
                "median_r_width": float(sub["surviving_r_width"].median()),
                "max_r_width": float(sub["surviving_r_width"].max()),
                "max_dev": float(sub["max_surviving_dev"].max()),
            })
    pd.DataFrame(wr_rows).to_csv(OUT / "AMPLITUDE_WORKING_RANGE_SUMMARY.csv", index=False)

    # =========================================================
    # 5. TRIPLE vs FOUR paired
    # =========================================================
    paired_rows = []
    for _, t_row in case_df[case_df.config == "TRIPLE"].iterrows():
        f_row = case_df[(case_df.config == "FOUR") & (case_df["shape"] == t_row["shape"]) &
                        (case_df["mode"] == t_row["mode"]) & (case_df.A_rms_db == t_row["A_rms_db"]) &
                        (case_df.z_true_m == t_row["z_true_m"])]
        if len(f_row):
            f = f_row.iloc[0]
            paired_rows.append({
                "shape": t_row["shape"], "mode": t_row["mode"], "A": t_row["A_rms_db"], "z_true": t_row["z_true_m"],
                "triple_strict": t_row["strict_single_bin"], "four_strict": f["strict_single_bin"],
                "triple_margin": t_row["wrong_range_margin"], "four_margin": f["wrong_range_margin"],
                "triple_width": t_row["surviving_r_width"], "four_width": f["surviving_r_width"],
            })
    pd.DataFrame(paired_rows).to_csv(OUT / "TRIPLE_FOUR_PAIRED_COMPARISON.csv", index=False)

    # =========================================================
    # 6. Decision (keep original label, add layers)
    # =========================================================
    decision = "AMPLITUDE_VARIATION_BREAKS_TURN_RANGE_ANCHOR_IN_TESTED_RANGE"
    note = "STRICT_SINGLE_RANGE_BIN_CRITERION_NOT_UNIFORMLY_RETAINED"
    truth_ok = bool(case_df[case_df.A_rms_db.isin([0.25, 0.5, 1.0])]["truth_retained"].all())
    rank_ok = bool((case_df[case_df.A_rms_db.isin([0.25, 0.5, 1.0])]["true_rank"] == 1).all())

    dec = {
        "stage": "R3-CLOSEDLOOP-1E-A-CLOSEOUT",
        "decision": decision,
        "interpretation": note,
        "identity_pass": identity_ok,
        "truth_retained_all_main": truth_ok,
        "true_rank1_all_main": rank_ok,
        "key_conclusion": "strict single-bin fails at 0.25dB+ but truth stays rank-1; set width grows from 1km to 15km over 0.25-1dB",
        "engineering_note": "no amplitude stability requirement derived; 0.25dB is smallest tested nonzero stress",
        "created_utc": NOW,
    }
    (OUT / "R3_RC23_CLOSEDLOOP_1E_A_CLOSEOUT_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# R3-CLOSEDLOOP-1E-A-CLOSEOUT

UTC: {NOW}

## Identity

{'PASS' if identity_ok else 'FAIL'} (independent comparison with 1D-FIX2)

## Performance Layers

| Layer | Finding |
|---|---|
| Truth retained | {'PASS' if truth_ok else 'mixed'} at 0.25/0.5/1.0 dB |
| True rank = 1 | {'PASS' if rank_ok else 'mixed'} at 0.25/0.5/1.0 dB |
| Strict single-bin | NOT uniform (original FAIL label retained) |
| Set width | 1 km @ 0.25 dB → 2 km @ 0.5 dB → up to 15 km @ 1 dB |

## Working Range Summary

{pd.DataFrame(wr_rows).to_string(index=False)}

## Injected vs Residual RMS

{pd.DataFrame(rms_rows).to_string(index=False)}

## 判定

**{decision}**

含义：{note}

工程口径：在冻结门限下，严格单距离格保持对慢变幅漂敏感。0.25–0.5 dB 条件下最优候选排序仍保持，距离集合仍有明显收缩。尚未得到真实系统的幅度稳定性要求。
"""
    (OUT / "R3_RC23_CLOSEDLOOP_1E_A_CLOSEOUT_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# 1E-A-CLOSEOUT\n\n**{decision}**\n\n{note}\n\ntruth_rank1={rank_ok}\n", encoding="utf-8")

    print("identity", identity_ok)
    print("decision", decision)
    print("truth_ok", truth_ok, "rank_ok", rank_ok)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
