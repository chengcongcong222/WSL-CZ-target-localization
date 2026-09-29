#!/usr/bin/env python3
"""RANGE-UB-2B-FIX0A: Xu2024 secondary slope-family identity audit at 50 km.

No new BELLHOP. Read existing Radon maps, build theory envelope, check family identity.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "RANGE_UB" / "E_STD_VLA_SLOPE_FAMILY_ID"
FIX0 = ROOT / "results" / "RANGE_UB" / "E_STD_VLA_DELAY_DEPTH_RIDGE"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

# E-STD frozen parameters
H = 5000.0  # m
R = 50000.0  # m
ZS_LIST = [180.0, 200.0, 220.0]
C_MIN, C_MAX = 1500.0, 1635.2  # m/s from E-STD SSP
NBEAMS = [1601, 3201]
DK = 0.002  # ms/m local-max check step


def theory_envelope():
    """Enumerate k32/k42/k43 from Xu2024 Eq.(10),(11),(13) with n_i=+-1."""
    rows = []
    for zs in ZS_LIST:
        for c in np.linspace(C_MIN, C_MAX, 20):
            for n2 in [-1, 1]:
                for n3 in [-1, 1]:
                    # k32 = [-6H - (n3+n2)*zs] / (cR)  (Eq.10)
                    k32 = (-6 * H - (n3 + n2) * zs) / (c * R) * 1000  # ms/m
                    rows.append({"family": "k32", "zs": zs, "c": c, "n2": n2, "n3": n3, "k_ms_per_m": k32})
                for n4 in [-1, 1]:
                    # k42 = [2H + (n4-n2)*zs] / (cR)  (Eq.11)
                    k42 = (2 * H + (n4 - n2) * zs) / (c * R) * 1000
                    rows.append({"family": "k42", "zs": zs, "c": c, "n2": n2, "n4": n4, "k_ms_per_m": k42})
                    # k43 = [8H + (n4+n3)*zs] / (cR)  (from tau43 = tau41 - tau31)
                    k43 = (8 * H + (n4 + n3) * zs) / (c * R) * 1000
                    rows.append({"family": "k43", "zs": zs, "c": c, "n3": n3, "n4": n4, "k_ms_per_m": k43})
    return pd.DataFrame(rows)


def load_radon_map(nb, zs):
    fp = FIX0 / f"RADON_SLOPE_MAP_nb{nb}_zs{int(zs)}.csv"
    if not fp.exists():
        return None
    df = pd.read_csv(fp)
    return df


def slope_profile(df):
    """P(k) = max_b S(k,b)."""
    prof = df.groupby("k_ms_per_m")["score"].max().reset_index()
    prof = prof.sort_values("k_ms_per_m").reset_index(drop=True)
    return prof


def find_family_peak(prof, k_lo, k_hi, full_df=None):
    """Find local max of P(k) within [k_lo, k_hi]."""
    mask = (prof["k_ms_per_m"] >= k_lo) & (prof["k_ms_per_m"] <= k_hi)
    sub = prof[mask]
    if len(sub) == 0:
        return None
    ks = sub["k_ms_per_m"].values
    ps = sub["score"].values
    best_i = int(np.argmax(ps))
    k_star = ks[best_i]
    p_star = ps[best_i]
    is_local = True
    if best_i > 0 and ps[best_i] <= ps[best_i - 1]:
        is_local = False
    if best_i < len(ps) - 1 and ps[best_i] < ps[best_i + 1]:
        is_local = False
    if not is_local:
        return None
    b_star = np.nan
    if full_df is not None:
        rows_at_k = full_df[full_df["k_ms_per_m"] == k_star]
        if len(rows_at_k) > 0:
            b_star = float(rows_at_k.loc[rows_at_k["score"].idxmax(), "b_ms"])
    return {"k_star": float(k_star), "b_star": float(b_star), "score": float(p_star), "is_local_max": True}


def main():
    config = {
        "stage": "RANGE-UB-2B-FIX0A", "created_utc": NOW,
        "baseline_commit": "b15bc03aae38b4aa5603d8a58711de0aae20b4be",
        "goal": "Xu2024 secondary slope-family identity audit at 50km",
        "no_new_bellhop": True,
    }
    (OUT / "RANGE_UB_2B_FIX0A_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    # scope correction
    (OUT / "FIX0_SCOPE_CORRECTION.md").write_text(
        """# FIX0_SCOPE_CORRECTION

UTC: """ + NOW + """

FIX0 原判 E_STD_VLA_SECONDARY_DELAY_RIDGES_CONFIRMED 收窄为:
GENERIC_VLA_DELAY_DEPTH_RIDGES_STABLE_AT_50KM

原因: convergence gate triplet != g-vector triplet。
g-vector 标记 SUPERSEDED_BY_ARBITRARY_RIDGE_ORDERING。

本轮做 theory-family identity audit。
""", encoding="utf-8")

    # =========================================================
    # 1. Theory envelope
    # =========================================================
    env_df = theory_envelope()
    env_summary = env_df.groupby("family")["k_ms_per_m"].agg(["min", "max"]).reset_index()
    env_summary.columns = ["family", "k_min_ms_per_m", "k_max_ms_per_m"]
    env_summary.to_csv(OUT / "E_STD_XU2024_SLOPE_ENVELOPE.csv", index=False)
    env_dict = {r["family"]: (r["k_min_ms_per_m"], r["k_max_ms_per_m"]) for _, r in env_summary.iterrows()}

    # =========================================================
    # 2. Family peak detection in Radon maps
    # =========================================================
    peak_rows = []
    conv_rows = []
    g_rows = []
    top_audit_rows = []

    for zs in ZS_LIST:
        maps = {}
        full_maps = {}
        for nb in NBEAMS:
            df = load_radon_map(nb, zs)
            if df is not None:
                maps[nb] = slope_profile(df)
                full_maps[nb] = df

        for fam in ["k32", "k42", "k43"]:
            k_lo, k_hi = env_dict[fam]
            results = {}
            for nb in NBEAMS:
                prof = maps.get(nb)
                if prof is None:
                    results[nb] = None
                    continue
                peak = find_family_peak(prof, k_lo, k_hi, full_maps.get(nb))
                results[nb] = peak
                # global score rank
                if peak:
                    all_scores = prof["score"].values
                    rank = int(np.sum(all_scores > peak["score"])) + 1
                    gmax = float(all_scores.max())
                    peak_rows.append({
                        "zs": zs, "family": fam, "nbeams": nb,
                        "k_lo": k_lo, "k_hi": k_hi,
                        "k_star": peak["k_star"], "b_star": peak["b_star"],
                        "score": peak["score"], "global_rank": rank,
                        "score_frac": peak["score"] / gmax if gmax > 0 else 0,
                        "detected": True,
                    })
                else:
                    peak_rows.append({
                        "zs": zs, "family": fam, "nbeams": nb,
                        "k_lo": k_lo, "k_hi": k_hi,
                        "detected": False,
                    })

            # convergence
            r1, r2 = results.get(1601), results.get(3201)
            if r1 and r2:
                sign_ok = np.sign(r1["k_star"]) == np.sign(r2["k_star"])
                rel_diff = abs(r1["k_star"] - r2["k_star"]) / max(abs(r2["k_star"]), 0.02)
                b_diff = abs(r1["b_star"] - r2["b_star"])
                ok = sign_ok and rel_diff <= 0.05 and b_diff <= 50.0
                conv_rows.append({
                    "zs": zs, "family": fam,
                    "k_1601": r1["k_star"], "k_3201": r2["k_star"],
                    "rel_diff": rel_diff, "b_diff_ms": b_diff,
                    "sign_ok": sign_ok, "pass": ok,
                })
            else:
                conv_rows.append({"zs": zs, "family": fam, "pass": False, "note": "FAMILY_NOT_DETECTED"})

            # g-vector (use 1601)
            if r1:
                if fam == "k32":
                    g_rows.append({"zs": zs, "family": fam, "k": r1["k_star"], "g": -r1["k_star"] / 6})
                elif fam == "k42":
                    g_rows.append({"zs": zs, "family": fam, "k": r1["k_star"], "g": r1["k_star"] / 2})
                else:
                    g_rows.append({"zs": zs, "family": fam, "k": r1["k_star"], "g": r1["k_star"] / 8})

        # top-score vs theory audit
        prof = maps.get(1601)
        if prof is not None:
            top3 = prof.nlargest(3, "score")
            for i, (_, r) in enumerate(top3.iterrows()):
                k = r["k_ms_per_m"]
                in_k32 = env_dict["k32"][0] <= k <= env_dict["k32"][1]
                in_k42 = env_dict["k42"][0] <= k <= env_dict["k42"][1]
                in_k43 = env_dict["k43"][0] <= k <= env_dict["k43"][1]
                top_audit_rows.append({
                    "zs": zs, "rank": i + 1, "k": k, "score": r["score"],
                    "in_k32_band": in_k32, "in_k42_band": in_k42, "in_k43_band": in_k43,
                    "matches_theory": in_k32 or in_k42 or in_k43,
                })

    peak_df = pd.DataFrame(peak_rows)
    peak_df.to_csv(OUT / "SLOPE_PROFILE_FAMILY_PEAKS.csv", index=False)
    conv_df = pd.DataFrame(conv_rows)
    conv_df.to_csv(OUT / "XU2024_FAMILY_RIDGE_CONVERGENCE.csv", index=False)
    g_df = pd.DataFrame(g_rows)
    g_df.to_csv(OUT / "THEORY_FAMILY_G_VECTOR.csv", index=False)
    top_df = pd.DataFrame(top_audit_rows)
    top_df.to_csv(OUT / "TOP_SCORE_VS_THEORY_FAMILY_AUDIT.csv", index=False)

    # =========================================================
    # 3. Decision
    # =========================================================
    all_9_pass = bool(conv_df["pass"].all()) and len(conv_df) == 9
    any_pass = bool(conv_df["pass"].any())
    n_detected = int(peak_df[peak_df["detected"] == True].shape[0])

    if all_9_pass:
        decision = "XU2024_SECONDARY_SLOPE_FAMILIES_IDENTIFIED_AT_50KM"
    elif any_pass and n_detected >= 6:
        decision = "XU2024_SECONDARY_SLOPE_FAMILIES_PARTIAL_AT_50KM"
    else:
        decision = "GENERIC_RIDGES_EXIST_XU2024_FAMILY_IDENTITY_NOT_ESTABLISHED"

    why = f"all_9_pass={all_9_pass}, any_pass={any_pass}, n_detected={n_detected}/18"
    dec = {
        "stage": "RANGE-UB-2B-FIX0A", "decision": decision, "why": why,
        "note": "50km family identity only; no range discrimination yet",
        "created_utc": NOW,
    }
    (OUT / "RANGE_UB_2B_FIX0A_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# RANGE-UB-2B-FIX0A — Slope Family Identity Audit

UTC: {NOW}

基线：`b15bc03aae38b4aa5603d8a58711de0aae20b4be`

## Theory Envelope (E-STD)

{env_summary.to_string(index=False)}

## Family Peak Detection

{peak_df.to_string(index=False)}

## Family Convergence (1601 vs 3201)

{conv_df.to_string(index=False)}

## g-vector (theory family)

{g_df.to_string(index=False) if len(g_df) else 'N/A'}

## Top-score vs Theory Family

{top_df.to_string(index=False) if len(top_df) else 'N/A'}

## 判定

### `{decision}`

{why}

## 未做

45–60 km、651 候选、新 BELLHOP、噪声、HLA 融合、P5。
"""
    (OUT / "RANGE_UB_2B_FIX0A_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# RANGE-UB-2B-FIX0A\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
