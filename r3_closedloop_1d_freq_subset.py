#!/usr/bin/env python3
"""R3-CLOSEDLOOP-1D: Frequency-subset sufficiency of turn-assisted range anchor.

delta=15deg, 15 frequency subsets, range anchor definition frozen.
"""
from __future__ import annotations

import itertools
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_FREQ_SUBSET"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
A1C = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_SINGLE_TURN_RANGE_ANCHOR"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_kraken"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0
T_TURN = 600.0
T_END = 1200.0
DT_OBS = 10.0
SIGMA_DEG = 0.1
RNG_SEED = 20260912
DELTA = 15.0
TRUTH = dict(r0_m=50e3, theta0_deg=0.0, v=2.0, psi_deg=5.0)
R_MIN_M, R_MAX_M = 45e3, 60e3
R_GRID_KM = np.arange(45.0, 60.0 + 1e-9, 1.0)
TH_GRID_DEG = np.arange(-5.0, 5.0 + 1e-9, 0.5)
V_GRID = np.arange(1.0, 3.0 + 1e-9, 0.2)
PSI_GRID_DEG = np.arange(-15.0, 15.0 + 1e-9, 1.0)
N_NODES = int(len(R_GRID_KM) * len(TH_GRID_DEG) * len(V_GRID) * len(PSI_GRID_DEG))
Z_TRUE_LIST = [180.0, 200.0, 220.0]
Z_PROFILE = np.arange(150.0, 250.0 + 1e-9, 5.0)
FREQS_ALL = [201.0, 235.0, 283.0, 338.0]
ZR = 200.0
TAU = 0.5


# ---- geometry ----
def platform_states_turn(t, delta_deg, t_turn=T_TURN, u=U_PLAT):
    t = np.asarray(t, float)
    xp = np.empty_like(t)
    yp = np.empty_like(t)
    d = math.radians(delta_deg)
    c, s = math.cos(d), math.sin(d)
    if abs(delta_deg) < 1e-12:
        xp[:] = u * t
        yp[:] = 0.0
        return xp, yp
    mask = t <= t_turn
    xp[mask] = u * t[mask]
    yp[mask] = 0.0
    dt = t[~mask] - t_turn
    xp[~mask] = u * t_turn + u * dt * c
    yp[~mask] = u * dt * s
    return xp, yp


def target_states(t, r0_m, theta0_rad, v, psi_rad):
    xt = r0_m * np.cos(theta0_rad) + v * t * np.cos(psi_rad)
    yt = r0_m * np.sin(theta0_rad) + v * t * np.sin(psi_rad)
    return xt, yt


def predict_bearings_batch(t, r0, th0, v, psi, delta_deg):
    xp, yp = platform_states_turn(t, delta_deg)
    xt = r0[:, None] * np.cos(th0[:, None]) + v[:, None] * t[None, :] * np.cos(psi[:, None])
    yt = r0[:, None] * np.sin(th0[:, None]) + v[:, None] * t[None, :] * np.sin(psi[:, None])
    return np.arctan2(yt - yp[None, :], xt - xp[None, :])


def wrap_angle(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def range_traj(t, r0_m, theta0_deg, v, psi_deg, delta_deg):
    xp, yp = platform_states_turn(t, delta_deg)
    xt, yt = target_states(t, r0_m, np.deg2rad(theta0_deg), v, np.deg2rad(psi_deg))
    return np.hypot(xt - xp, yt - yp)


def chi2_thr(sigma_rad):
    return 1e-12 if sigma_rad <= 0 else 13.3 * (sigma_rad ** 2)


# ---- acoustics ----
def parse_mod(path):
    buf = path.read_bytes()
    recl = 4 * int(np.frombuffer(buf[:4], dtype="<i4")[0])
    hdr = np.frombuffer(buf[84:108], dtype="<i4")
    ntot, nmat = int(hdr[2]), int(hdr[3])
    depths = np.frombuffer(buf[4 * recl: 5 * recl], dtype="<f4")[:ntot].astype(float)
    M = int(np.frombuffer(buf[5 * recl: 5 * recl + 4], dtype="<i4")[0])
    phi = np.zeros((nmat, M), complex)
    for im in range(M):
        off = (7 + im) * recl
        chunk = np.frombuffer(buf[off: off + recl], dtype="<c8")
        take = min(nmat, chunk.size)
        phi[:take, im] = chunk[:take]
    k = np.frombuffer(buf[(7 + M) * recl: (7 + M) * recl + M * 8], dtype="<c8")
    return {"depths": depths, "M": M, "phi": phi, "k": k}


def F_matrix(mod, r):
    kre = mod["k"].real
    alpha = -mod["k"].imag
    r = np.asarray(r, float)
    amp = np.sqrt(2 * np.pi / (kre[:, None] * r[None, :]))
    phase = np.exp(-1j * kre[:, None] * r[None, :] - alpha[:, None] * r[None, :] - 1j * np.pi / 4)
    return amp * phase


def L_profile(mod, F, z_list, zr=ZR):
    izr = int(np.argmin(np.abs(mod["depths"] - zr)))
    pr = mod["phi"][izr]
    W = mod["phi"] * pr[None, :]
    izs = [int(np.argmin(np.abs(mod["depths"] - z))) for z in z_list]
    P = W[izs, :] @ F
    return 20.0 * np.log10(np.maximum(np.abs(P), 1e-30))


def demean(L):
    return np.asarray(L, float) - float(np.mean(L))


def rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def main():
    config = {
        "stage": "R3-CLOSEDLOOP-1D", "created_utc": NOW,
        "baseline_commit": "e3030060e3bc36fc35e3232e90a2dbb9ead92da9",
        "delta_deg": DELTA, "tau_db": TAU,
        "freqs": FREQS_ALL,
        "n_subsets": 15,
    }
    (OUT / "CLOSEDLOOP_1D_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    # =========================================================
    # 1. RC2 at delta=15
    # =========================================================
    t_full = np.arange(0.0, T_END + 1e-9, DT_OBS)
    n_w1 = int(np.sum(t_full <= T_TURN + 1e-9))
    t_w1, t_w2 = t_full[:n_w1], t_full[n_w1:]
    n_w2 = len(t_w2)
    sigma_rad = np.deg2rad(SIGMA_DEG)
    thr = chi2_thr(sigma_rad)

    rr, tt, vv, pp = np.meshgrid(R_GRID_KM * 1e3, np.deg2rad(TH_GRID_DEG), V_GRID, np.deg2rad(PSI_GRID_DEG), indexing="ij")
    r0_all = rr.ravel()
    th0_all = tt.ravel()
    v_all = vv.ravel()
    psi_all = pp.ravel()

    rng = np.random.default_rng(RNG_SEED)
    noise = rng.normal(0.0, sigma_rad, size=t_full.shape)
    th_true = np.arctan2(
        target_states(t_full, TRUTH["r0_m"], np.deg2rad(TRUTH["theta0_deg"]), TRUTH["v"], np.deg2rad(TRUTH["psi_deg"]))[1]
        - platform_states_turn(t_full, DELTA)[1],
        target_states(t_full, TRUTH["r0_m"], np.deg2rad(TRUTH["theta0_deg"]), TRUTH["v"], np.deg2rad(TRUTH["psi_deg"]))[0]
        - platform_states_turn(t_full, DELTA)[0],
    )
    obs = th_true + noise

    pred = predict_bearings_batch(t_full, r0_all, th0_all, v_all, psi_all, DELTA)
    cost = np.sum(wrap_angle(pred - obs[None, :]) ** 2, axis=1)
    cmin = float(np.min(cost))
    acc = cost <= (cmin + thr)
    n_rc2 = int(acc.sum())
    idx_acc = np.where(acc)[0]
    truth_idx = int(np.argmin(
        (r0_all - TRUTH["r0_m"]) ** 2 + (th0_all - np.deg2rad(TRUTH["theta0_deg"])) ** 2
        + (v_all - TRUTH["v"]) ** 2 + (psi_all - np.deg2rad(TRUTH["psi_deg"])) ** 2))

    # RC2 identity
    a1c_rc2 = pd.read_csv(A1C / "RC2_DELTA15_IDENTITY.csv") if (A1C / "RC2_DELTA15_IDENTITY.csv").exists() else pd.DataFrame()
    if len(a1c_rc2):
        n_1c = len(a1c_rc2)
    else:
        n_1c = 385  # from 1C
    rc2_ok = n_rc2 == 385 and bool(acc[truth_idx])
    pd.DataFrame([{"n_RC2": n_rc2, "expected": 385, "truth_in": bool(acc[truth_idx]), "pass": rc2_ok}]).to_csv(
        OUT / "RC2_DELTA15_IDENTITY.csv", index=False)
    if not rc2_ok:
        (OUT / "R3_RC23_CLOSEDLOOP_1D_DECISION.json").write_text(
            json.dumps({"stage": "R3-CLOSEDLOOP-1D", "decision": "CLOSEDLOOP_1D_BLOCKED_BY_RC2_IDENTITY", "created_utc": NOW}, indent=2),
            encoding="utf-8")
        print("BLOCKED_BY_RC2_IDENTITY")
        return 1

    # =========================================================
    # 2. Precompute acoustic L for all freqs, both windows
    # =========================================================
    mods = {f: parse_mod(ZGRID / f"zgrid_f{int(f)}.mod") for f in FREQS_ALL}
    n_z = Z_PROFILE.size

    # observation L per freq per z_true per window
    r_true = range_traj(t_full, TRUTH["r0_m"], TRUTH["theta0_deg"], TRUTH["v"], TRUTH["psi_deg"], DELTA)
    L_obs = {}  # (z_true, f) -> (w1_demeaned, w2_demeaned)
    for f in FREQS_ALL:
        F_t = F_matrix(mods[f], r_true)
        L_t = L_profile(mods[f], F_t, Z_TRUE_LIST)
        for i, zt in enumerate(Z_TRUE_LIST):
            L_obs[(zt, f)] = (demean(L_t[i, :n_w1]), demean(L_t[i, n_w1:]))

    # candidate L per (node, f) -> (w1, w2) arrays per z
    # precompute for accepted nodes
    L_cand = {}  # (node_idx, f) -> (L1_all_z, L2_all_z)
    for j in idx_acc:
        r_c = range_traj(t_full, r0_all[j], np.rad2deg(th0_all[j]), v_all[j], np.rad2deg(psi_all[j]), DELTA)
        r1, r2 = r_c[:n_w1], r_c[n_w1:]
        for f in FREQS_ALL:
            F1 = F_matrix(mods[f], r1)
            F2 = F_matrix(mods[f], r2)
            L1 = L_profile(mods[f], F1, Z_PROFILE)
            L2 = L_profile(mods[f], F2, Z_PROFILE)
            L_cand[(j, f)] = (L1, L2)

    # =========================================================
    # 3. Score all 15 subsets
    # =========================================================
    subsets = []
    for r in range(1, 5):
        for combo in itertools.combinations(FREQS_ALL, r):
            subsets.append(list(combo))

    score_rows = []
    surv_rows = []
    anchor_rows = []
    alias_by_freq = []

    for subset in subsets:
        n_f = len(subset)
        sub_tag = "+".join(f"{int(f)}" for f in subset)
        for z_true in Z_TRUE_LIST:
            # compute J for each accepted candidate
            J_all = np.empty(len(idx_acc))
            for i, j in enumerate(idx_acc):
                J_z = np.zeros(n_z)
                for f in subset:
                    o1, o2 = L_obs[(z_true, f)]
                    L1, L2 = L_cand[(j, f)]
                    for iz in range(n_z):
                        e1 = demean(L1[iz]) - o1
                        e2 = demean(L2[iz]) - o2
                        J_z[iz] += np.sum(e1 ** 2) + np.sum(e2 ** 2)
                J_z = np.sqrt(J_z / ((n_w1 + n_w2) * n_f))
                J_all[i] = J_z.min()
            jmin = float(J_all.min())
            order = np.argsort(J_all)
            pos_true = int(np.where(idx_acc == truth_idx)[0][0]) if truth_idx in idx_acc else -1
            true_rank = int(np.where(order == pos_true)[0][0]) + 1 if pos_true >= 0 else -1
            true_J = float(J_all[pos_true]) if pos_true >= 0 else np.nan
            best_false_i = int(order[1]) if len(order) > 1 else -1
            best_false_J = float(J_all[best_false_i]) if best_false_i >= 0 else np.nan
            best_false_r = float(r0_all[idx_acc[best_false_i]] / 1e3) if best_false_i >= 0 else np.nan

            # survivors at tau=0.5
            keep = J_all <= (jmin + TAU)
            n_surv = int(keep.sum())
            surv_r = r0_all[idx_acc][keep] / 1e3
            r_min = float(surv_r.min()) if n_surv else np.nan
            r_max = float(surv_r.max()) if n_surv else np.nan
            r_bins = len(set(np.round(surv_r).astype(int))) if n_surv else 0
            true_kept = bool(keep[pos_true]) if pos_true >= 0 else False

            # range anchor
            anchored = true_kept and true_rank == 1 and r_bins == 1 and abs(r_min - 50.0) < 0.1

            score_rows.append({
                "subset": sub_tag, "n_lines": n_f, "z_true_m": z_true,
                "n_RC2": len(idx_acc), "n_survivors_tau05": n_surv,
                "r_min_km": r_min, "r_max_km": r_max,
                "r_width_km": r_max - r_min if n_surv else np.nan,
                "r_bins_occupied": r_bins,
                "true_retained": true_kept, "true_rank": true_rank, "true_J": true_J,
                "best_false_J": best_false_J, "best_false_r_km": best_false_r,
                "margin_to_best_false": true_J - best_false_J if np.isfinite(true_J) and np.isfinite(best_false_J) else np.nan,
                "range_anchored": anchored,
            })

            # alias overlap
            for r_km in [45.0, 50.0, 56.0, 58.0, 60.0]:
                surv_at_r = bool((np.abs(surv_r - r_km) < 0.1).any()) if n_surv else False
                alias_by_freq.append({"subset": sub_tag, "z_true_m": z_true, "r_km": r_km, "survives": surv_at_r})

    score_df = pd.DataFrame(score_rows)
    score_df.to_csv(OUT / "SUBSET_PROFILED_SCORES.csv", index=False)

    # subset anchor summary
    anchor_summary = []
    for subset in subsets:
        sub_tag = "+".join(f"{int(f)}" for f in subset)
        n_f = len(subset)
        sub = score_df[score_df["subset"] == sub_tag]
        all_anchored = bool(sub["range_anchored"].all()) if len(sub) else False
        anchor_summary.append({"subset": sub_tag, "n_lines": n_f, "three_depth_range_anchor": all_anchored,
                               "n_anchored": int(sub["range_anchored"].sum())})
    anchor_df = pd.DataFrame(anchor_summary)
    anchor_df.to_csv(OUT / "SUBSET_ANCHOR_SUMMARY.csv", index=False)

    # =========================================================
    # 4. Identity gates vs 1C
    # =========================================================
    # 235Hz identity: compare with 1C MAIN
    one_c_main = pd.DataFrame([{"subset": "235", "note": "compare with 1C MAIN delta=15 tau=0.5"}])
    one_c_main.to_csv(OUT / "ONE_C_MAIN_IDENTITY.csv", index=False)
    # FOUR identity
    four_sub = score_df[score_df["subset"] == "201+235+283+338"]
    four_anchored = bool(four_sub["range_anchored"].all()) if len(four_sub) else False
    pd.DataFrame([{"subset": "201+235+283+338", "four_identity": four_anchored,
                   "expected": "all three depths r=50 only"}]).to_csv(OUT / "ONE_C_FOUR_IDENTITY.csv", index=False)

    # =========================================================
    # 5. Frequency complementarity
    # =========================================================
    comp_rows = []
    for f in FREQS_ALL:
        f_str = f"{int(f)}"
        subsets_with_f = [s for s in subsets if f in s]
        anchored_with_f = sum(1 for s in subsets_with_f if anchor_df[anchor_df["subset"] == "+".join(f"{int(x)}" for x in s)]["three_depth_range_anchor"].any())
        comp_rows.append({"freq_hz": f, "n_subsets_containing": len(subsets_with_f),
                          "n_anchored_subsets": anchored_with_f})
    pd.DataFrame(comp_rows).to_csv(OUT / "FREQUENCY_COMPLEMENTARITY.csv", index=False)
    pd.DataFrame(alias_by_freq).to_csv(OUT / "RANGE_ALIAS_OVERLAP_BY_FREQUENCY.csv", index=False)

    # =========================================================
    # 6. Decision
    # =========================================================
    singles = anchor_df[anchor_df["n_lines"] == 1]
    doubles = anchor_df[anchor_df["n_lines"] == 2]
    triples = anchor_df[anchor_df["n_lines"] == 3]
    four = anchor_df[anchor_df["n_lines"] == 4]

    any_single = bool(singles["three_depth_range_anchor"].any())
    all_single = bool(singles["three_depth_range_anchor"].all())
    any_double = bool(doubles["three_depth_range_anchor"].any())
    all_double = bool(doubles["three_depth_range_anchor"].all())
    any_triple = bool(triples["three_depth_range_anchor"].any())
    four_ok = bool(four["three_depth_range_anchor"].any()) if len(four) else False

    if not four_ok:
        decision = "FOUR_LINE_ANCHOR_NOT_REPRODUCED"
    elif all_single:
        decision = "ANY_SINGLE_LINE_TURN_RANGE_ANCHOR"
    elif any_single:
        decision = "SOME_SINGLE_LINE_TURN_RANGE_ANCHOR"
    elif all_double:
        decision = "ANY_TWO_OF_FOUR_TURN_RANGE_ANCHOR"
    elif any_double:
        decision = "SOME_TWO_LINE_TURN_RANGE_ANCHOR"
    elif any_triple:
        decision = "SOME_THREE_LINE_TURN_RANGE_ANCHOR"
    else:
        decision = "FOUR_LINE_REQUIRED_IN_TESTED_CASE"

    why = f"any_single={any_single}, all_single={all_single}, any_double={any_double}, all_double={all_double}, any_triple={any_triple}, four={four_ok}"
    dec = {"stage": "R3-CLOSEDLOOP-1D", "decision": decision, "why": why, "created_utc": NOW}
    (OUT / "R3_RC23_CLOSEDLOOP_1D_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# R3-CLOSEDLOOP-1D — Frequency-Subset Sufficiency of Turn-Assisted Range Anchor

UTC: {NOW}

基线：`e3030060e3bc36fc35e3232e90a2dbb9ead92da9`

## RC2 Identity

n_RC2 = {n_rc2} (expected 385), truth_in = {bool(acc[truth_idx])}

## Subset Anchor Summary

{anchor_df.to_string(index=False)}

## 判定

### `{decision}`

{why}

## Frequency Complementarity

{pd.DataFrame(comp_rows).to_string(index=False)}

## 未做

SSP、幅漂、频漂、噪声、P5。
"""
    (OUT / "R3_RC23_CLOSEDLOOP_1D_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# R3-CLOSEDLOOP-1D\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
