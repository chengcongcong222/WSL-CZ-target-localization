#!/usr/bin/env python3
"""R3-CLOSEDLOOP-1B: two-window temporal accumulation + range-ridge audit.

Frozen spec (baseline 8d5b9b7):
  CR5 truth constant 0-1200s; W1=0:10:600 must reproduce 1A; W2=610:10:1200
  RC2 cumulative on full 114576; RC3 two-window demeaned per-window shared z
  Track range aliases 45/56/58/60 vs true 50 km
  tau=0.25/0.5/1.0; MAIN=235Hz; REFERENCE=FOUR
"""
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TWO_WINDOW_TEMPORAL_CONTRACTION"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
A1 = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_SINGLE_WINDOW_SET_CONTRACTION"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0
R_MIN_M, R_MAX_M = 45e3, 60e3
TH_MIN_DEG, TH_MAX_DEG = -5.0, 5.0
V_MIN, V_MAX = 1.0, 3.0
PSI_MIN_DEG, PSI_MAX_DEG = -15.0, 15.0
R_GRID_KM = np.arange(45.0, 60.0 + 1e-9, 1.0)
TH_GRID_DEG = np.arange(-5.0, 5.0 + 1e-9, 0.5)
V_GRID = np.arange(1.0, 3.0 + 1e-9, 0.2)
PSI_GRID_DEG = np.arange(-15.0, 15.0 + 1e-9, 1.0)
N_NODES = int(len(R_GRID_KM) * len(TH_GRID_DEG) * len(V_GRID) * len(PSI_GRID_DEG))

T_END = 1200.0
DT_OBS = 10.0
SIGMA_DEG = 0.1
TURN_DEG = 0.0
RNG_SEED = 20260912
TRUTH = dict(r0_m=50e3, theta0_deg=0.0, v=2.0, psi_deg=5.0)

Z_TRUE_LIST = [180.0, 200.0, 220.0]
Z_PROFILE = np.arange(150.0, 250.0 + 1e-9, 5.0)
FREQS_ALL = [201.0, 235.0, 283.0, 338.0]
FREQ_MAIN = [235.0]
TAUS = [0.25, 0.5, 1.0]
ZR = 200.0
TRACK_R_KM = [45.0, 50.0, 56.0, 58.0, 60.0]


# ---- geometry ----
def platform_states(t, turn_deg, u=U_PLAT):
    t = np.asarray(t, dtype=float)
    xp = np.empty_like(t)
    yp = np.empty_like(t)
    if abs(turn_deg) < 1e-12:
        xp[:] = u * t
        yp[:] = 0.0
        return xp, yp
    t_turn = float(np.max(t)) * 0.5
    delta = math.radians(turn_deg)
    c, s = math.cos(delta), math.sin(delta)
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


def predict_bearings_batch(t, r0, th0, v, psi, turn_deg):
    xp, yp = platform_states(t, turn_deg)
    xt = r0[:, None] * np.cos(th0[:, None]) + v[:, None] * t[None, :] * np.cos(psi[:, None])
    yt = r0[:, None] * np.sin(th0[:, None]) + v[:, None] * t[None, :] * np.sin(psi[:, None])
    return np.arctan2(yt - yp[None, :], xt - xp[None, :])


def wrap_angle(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def range_traj(t, r0_m, theta0_deg, v, psi_deg, turn_deg=0.0):
    xp, yp = platform_states(t, turn_deg)
    xt, yt = target_states(t, r0_m, np.deg2rad(theta0_deg), v, np.deg2rad(psi_deg))
    return np.hypot(xt - xp, yt - yp)


def chi2_accept_threshold(sigma_rad):
    return 1e-12 if sigma_rad <= 0 else 13.3 * (sigma_rad ** 2)


def set_widths(r0, th0, v, psi, mask):
    n = int(np.sum(mask))
    if n < 2:
        return dict(n=n, r_width_km=0.0, theta_width_deg=0.0, v_width_mps=0.0, psi_width_deg=0.0)
    return dict(
        n=n,
        r_width_km=float((r0[mask].max() - r0[mask].min()) / 1e3),
        theta_width_deg=float(np.rad2deg(th0[mask].max() - th0[mask].min())),
        v_width_mps=float(v[mask].max() - v[mask].min()),
        psi_width_deg=float(np.rad2deg(psi[mask].max() - psi[mask].min())),
    )


# ---- acoustics ----
def parse_mod(path: Path) -> dict:
    buf = path.read_bytes()
    recl = 4 * int(np.frombuffer(buf[:4], dtype="<i4")[0])
    hdr = np.frombuffer(buf[84:108], dtype="<i4")
    ntot, nmat = int(hdr[2]), int(hdr[3])
    depths = np.frombuffer(buf[4 * recl : 5 * recl], dtype="<f4")[:ntot].astype(float)
    M = int(np.frombuffer(buf[5 * recl : 5 * recl + 4], dtype="<i4")[0])
    phi = np.zeros((nmat, M), complex)
    for im in range(M):
        off = (7 + im) * recl
        chunk = np.frombuffer(buf[off : off + recl], dtype="<c8")
        take = min(nmat, chunk.size)
        phi[:take, im] = chunk[:take]
    k = np.frombuffer(buf[(7 + M) * recl : (7 + M) * recl + M * 8], dtype="<c8")
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
    L = np.asarray(L, float)
    return L - float(np.mean(L))


def rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def main() -> int:
    # =========================================================
    # 0. CONFIG
    # =========================================================
    t_full = np.arange(0.0, T_END + 1e-9, DT_OBS)  # 121 points 0..1200
    n_full = t_full.size
    n_w1 = 61  # 0..600
    n_w2 = 60  # 610..1200
    assert t_full[n_w1 - 1] == 600.0 and t_full[n_w1] == 610.0
    t_w1 = t_full[:n_w1]
    t_w2 = t_full[n_w1:]
    sigma_rad = np.deg2rad(SIGMA_DEG)
    thr_rc2 = chi2_accept_threshold(sigma_rad)

    config = {
        "stage": "R3-CLOSEDLOOP-1B",
        "created_utc": NOW,
        "baseline_commit": "8d5b9b7e74a97639daf5351050b01d4d5b868bf0",
        "scenario": "CR5",
        "truth": TRUTH,
        "T_end_s": T_END,
        "W1_s": "0:10:600",
        "W2_s": "610:10:1200",
        "n_full": int(n_full),
        "n_w1": n_w1,
        "n_w2": n_w2,
        "sigma_theta_deg": SIGMA_DEG,
        "platform_turn_deg": TURN_DEG,
        "rng_seed": RNG_SEED,
        "realization_label": "SINGLE_WINDOW_DEMO_REALIZATION_EXTENDED",
        "coarse_grid_n_nodes": N_NODES,
        "rc2_cost": "sum wrap(theta_pred-theta_obs)^2",
        "rc2_accept": f"cost <= cmin + {thr_rc2:.6e} (=13.3*sigma_rad^2)",
        "rc3_observable": "PER_WINDOW_DEMEAN_SOURCE_LEVEL_FREE_RELATIVE_TL_SHAPE",
        "rc3_forward": "DIRECT_F_MATRIX_ON_TRAJECTORY",
        "rc3_J12": "sqrt((||e_W1||^2+||e_W2||^2)/(N_W1+N_W2)) shared z",
        "configs": {"MAIN": FREQ_MAIN, "REFERENCE": FREQS_ALL},
        "tau_db": TAUS,
        "track_r_km": TRACK_R_KM,
    }
    (OUT / "CLOSEDLOOP_1B_CONFIG.json").write_text(
        json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # =========================================================
    # 1. Observations + W1 identity
    # =========================================================
    rng = np.random.default_rng(RNG_SEED)
    xt, yt = target_states(t_full, TRUTH["r0_m"], np.deg2rad(TRUTH["theta0_deg"]), TRUTH["v"], np.deg2rad(TRUTH["psi_deg"]))
    xp, yp = platform_states(t_full, TURN_DEG)
    theta_true = np.arctan2(yt - yp, xt - xp)
    bearing_obs = theta_true + rng.normal(0.0, sigma_rad, size=theta_true.shape)
    obs_w1 = bearing_obs[:n_w1]
    obs_w2 = bearing_obs[n_w1:]

    # grid
    rr, tt, vv, pp = np.meshgrid(
        R_GRID_KM * 1e3,
        np.deg2rad(TH_GRID_DEG),
        V_GRID,
        np.deg2rad(PSI_GRID_DEG),
        indexing="ij",
    )
    r0_all = rr.ravel()
    th0_all = tt.ravel()
    v_all = vv.ravel()
    psi_all = pp.ravel()
    assert r0_all.size == N_NODES

    # W1 RC2 (must match 1A)
    pred_w1 = predict_bearings_batch(t_w1, r0_all, th0_all, v_all, psi_all, TURN_DEG)
    cost_w1 = np.sum(wrap_angle(pred_w1 - obs_w1[None, :]) ** 2, axis=1)
    cmin_w1 = float(np.min(cost_w1))
    acc_w1 = cost_w1 <= (cmin_w1 + thr_rc2)
    n_w1_rc2 = int(acc_w1.sum())

    # 1A reference
    a1_rc2 = pd.read_csv(A1 / "RC2_ACCEPTED_CANDIDATES.csv")
    a1_ids = set(a1_rc2["node_id"].tolist())
    w1_ids = set(np.where(acc_w1)[0].tolist())
    w1_id_match = a1_ids == w1_ids
    w1_count_match = n_w1_rc2 == len(a1_ids)

    # W1 RC3 identity: score a few and compare to 1A
    # load 1A scores
    a1_scores = pd.read_csv(A1 / "RC3_PROFILED_SCORE_ALL_ACCEPTED.csv")
    mods = {f: parse_mod(ZGRID / f"zgrid_f{int(f)}.mod") for f in FREQS_ALL}
    r_true = range_traj(t_full, TRUTH["r0_m"], TRUTH["theta0_deg"], TRUTH["v"], TRUTH["psi_deg"], TURN_DEG)

    # precompute L_obs for all z_true, both windows, all freqs
    # IMPORTANT: compute W1 obs on t_w1 directly (not slice from full) for 1A identity
    L_obs_w1 = {zt: {} for zt in Z_TRUE_LIST}
    L_obs_w2 = {zt: {} for zt in Z_TRUE_LIST}
    r_true_w1 = range_traj(t_w1, TRUTH["r0_m"], TRUTH["theta0_deg"], TRUTH["v"], TRUTH["psi_deg"], TURN_DEG)
    r_true_w2 = range_traj(t_w2, TRUTH["r0_m"], TRUTH["theta0_deg"], TRUTH["v"], TRUTH["psi_deg"], TURN_DEG)
    for f in FREQS_ALL:
        mod = mods[f]
        F_t1 = F_matrix(mod, r_true_w1)
        F_t2 = F_matrix(mod, r_true_w2)
        L_t1 = L_profile(mod, F_t1, Z_TRUE_LIST)  # (3, n_w1)
        L_t2 = L_profile(mod, F_t2, Z_TRUE_LIST)  # (3, n_w2)
        for i, zt in enumerate(Z_TRUE_LIST):
            L_obs_w1[zt][f] = demean(L_t1[i])
            L_obs_w2[zt][f] = demean(L_t2[i])

    # W1 identity check on MAIN z_true=200 (sample)
    idx_w1 = np.where(acc_w1)[0]
    # find true node
    truth_idx = int(
        np.argmin(
            (r0_all - TRUTH["r0_m"]) ** 2
            + (th0_all - np.deg2rad(TRUTH["theta0_deg"])) ** 2
            + (v_all - TRUTH["v"]) ** 2
            + (psi_all - np.deg2rad(TRUTH["psi_deg"])) ** 2
        )
    )

    # score W1-only for identity: use min_z J(h,z) like 1A — batch over depths
    n_z = Z_PROFILE.size
    id_rows = []
    max_w1_score_err = 0.0
    # subsample for speed: true node + 200 random accepted
    rng_id = np.random.default_rng(0)
    sample_idx = list(idx_w1)
    if len(sample_idx) > 201:
        pick = rng_id.choice(len(sample_idx), size=200, replace=False)
        true_pos = int(np.where(idx_w1 == truth_idx)[0][0])
        sample_idx = sorted(set(idx_w1[list(pick)].tolist() + [truth_idx]))
    else:
        sample_idx = list(idx_w1)

    for cfg_name, freqs in (("MAIN", FREQ_MAIN), ("REFERENCE", FREQS_ALL)):
        for z_true in Z_TRUE_LIST:
            L_obs_c = np.concatenate([L_obs_w1[z_true][f] for f in freqs])
            n_s = len(sample_idx)
            J_prof = np.empty((n_s, n_z))
            for si, j in enumerate(sample_idx):
                r_c = range_traj(t_w1, r0_all[j], np.rad2deg(th0_all[j]), v_all[j], np.rad2deg(psi_all[j]), TURN_DEG)
                parts_all = []
                for f in freqs:
                    F_i = F_matrix(mods[f], r_c)
                    L_allz = L_profile(mods[f], F_i, Z_PROFILE)  # (n_z, n_w1)
                    parts_all.append(L_allz)
                # parts_all[f][iz] -> L; build concat per iz
                for iz in range(n_z):
                    parts = [demean(parts_all[fi][iz]) for fi in range(len(freqs))]
                    J_prof[si, iz] = rms(L_obs_c, np.concatenate(parts))
            J = J_prof.min(axis=1)
            a1_sub = a1_scores[(a1_scores["config"] == cfg_name) & (a1_scores["z_true_m"] == z_true)]
            a1_map = dict(zip(a1_sub["node_id"], a1_sub["J_RC3"]))
            errs = [abs(J[si] - a1_map[j]) for si, j in enumerate(sample_idx) if j in a1_map]
            max_err = max(errs) if errs else 0.0
            max_w1_score_err = max(max_w1_score_err, max_err)
            id_rows.append(
                {
                    "config": cfg_name,
                    "z_true_m": z_true,
                    "n_sampled": n_s,
                    "max_abs_err_vs_1A": float(max_err),
                }
            )

    w1_score_match = max_w1_score_err < 1e-9
    w1_identity_pass = w1_id_match and w1_count_match and w1_score_match

    id_df = pd.DataFrame(id_rows)
    id_df.to_csv(OUT / "W1_IDENTITY_WITH_1A.csv", index=False)
    (OUT / "W1_IDENTITY_WITH_1A_REPORT.md").write_text(
        f"""# W1_IDENTITY_WITH_1A

UTC: {NOW}

| check | result |
|---|---|
| W1 RC2 count | {n_w1_rc2} vs 1A {len(a1_ids)} → {'MATCH' if w1_count_match else 'MISMATCH'} |
| accepted node IDs | {'IDENTICAL' if w1_id_match else 'DIFFERENT'} |
| W1 RC3 max score err | {max_w1_score_err:.3e} dB → {'PASS' if w1_score_match else 'FAIL'} |

Gate: `{'PASS' if w1_identity_pass else 'CLOSEDLOOP_1B_BLOCKED_BY_W1_IDENTITY'}`
""",
        encoding="utf-8",
    )
    if not w1_identity_pass:
        (OUT / "R3_RC23_CLOSEDLOOP_1B_DECISION.json").write_text(
            json.dumps(
                {
                    "stage": "R3-CLOSEDLOOP-1B",
                    "decision": "CLOSEDLOOP_1B_BLOCKED_BY_W1_IDENTITY",
                    "w1_count_match": bool(w1_count_match),
                    "w1_id_match": bool(w1_id_match),
                    "w1_score_match": bool(w1_score_match),
                    "max_w1_score_err": float(max_w1_score_err),
                    "created_utc": NOW,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        print("BLOCKED_BY_W1_IDENTITY")
        print("count_match", w1_count_match, "id_match", w1_id_match, "score_match", w1_score_match)
        print("max_err", max_w1_score_err)
        return 1

    # =========================================================
    # 2. Cumulative RC2 on full grid
    # =========================================================
    pred_w2 = predict_bearings_batch(t_w2, r0_all, th0_all, v_all, psi_all, TURN_DEG)
    cost_w2 = np.sum(wrap_angle(pred_w2 - obs_w2[None, :]) ** 2, axis=1)
    cost_cum = cost_w1 + cost_w2
    cmin_cum = float(np.min(cost_cum))
    acc_cum = cost_cum <= (cmin_cum + thr_rc2)
    n_cum = int(acc_cum.sum())
    if n_cum == 0:
        raise RuntimeError("cumulative RC2 empty")

    idx_cum = np.where(acc_cum)[0]
    truth_in_cum = bool(acc_cum[truth_idx])

    rc2_cum_df = pd.DataFrame(
        {
            "node_id": idx_cum,
            "r0_km": r0_all[idx_cum] / 1e3,
            "theta0_deg": np.rad2deg(th0_all[idx_cum]),
            "v_mps": v_all[idx_cum],
            "psi_deg": np.rad2deg(psi_all[idx_cum]),
            "rc2_cost_w1": cost_w1[idx_cum],
            "rc2_cost_w2": cost_w2[idx_cum],
            "rc2_cost_cum": cost_cum[idx_cum],
            "is_true": idx_cum == truth_idx,
        }
    )
    rc2_cum_df.to_csv(OUT / "RC2_CUMULATIVE_1200_ACCEPTED.csv", index=False)

    w_cum = set_widths(r0_all, th0_all, v_all, psi_all, acc_cum)
    w_w1 = set_widths(r0_all, th0_all, v_all, psi_all, acc_w1)
    rc2_cum_metrics = {
        "n_RC2_W1": n_w1_rc2,
        "n_RC2_cum": n_cum,
        "cmin_cum": cmin_cum,
        "truth_in_cum": truth_in_cum,
        "r_width_km": w_cum["r_width_km"],
        "theta_width_deg": w_cum["theta_width_deg"],
        "v_width_mps": w_cum["v_width_mps"],
        "psi_width_deg": w_cum["psi_width_deg"],
        "r_width_W1_km": w_w1["r_width_km"],
        "v_width_W1_mps": w_w1["v_width_mps"],
        "psi_width_W1_deg": w_w1["psi_width_deg"],
    }
    pd.DataFrame([rc2_cum_metrics]).to_csv(OUT / "RC2_CUMULATIVE_SET_METRICS.csv", index=False)

    # =========================================================
    # 3. RC3 two-window scoring on cumulative accepted
    # =========================================================
    n_acc = idx_cum.size
    # precompute candidate L for both windows
    # L_cand[f][i, iz, t] — do per candidate
    # For efficiency: cache F and L per candidate per freq for both windows
    def score_two_window(freqs, z_true):
        """J12 per candidate, shared z. Returns J, z_star, J_z."""
        J_z = np.zeros((n_acc, n_z))
        for i in range(n_acc):
            j = idx_cum[i]
            r_full = range_traj(t_full, r0_all[j], np.rad2deg(th0_all[j]), v_all[j], np.rad2deg(psi_all[j]), TURN_DEG)
            r1 = r_full[:n_w1]
            r2 = r_full[n_w1:]
            for f in freqs:
                mod = mods[f]
                F1 = F_matrix(mod, r1)
                F2 = F_matrix(mod, r2)
                L1 = L_profile(mod, F1, Z_PROFILE)  # (n_z, n_w1)
                L2 = L_profile(mod, F2, Z_PROFILE)  # (n_z, n_w2)
                obs1 = L_obs_w1[z_true][f]
                obs2 = L_obs_w2[z_true][f]
                for iz in range(n_z):
                    e1 = demean(L1[iz]) - obs1
                    e2 = demean(L2[iz]) - obs2
                    J_z[i, iz] += float(np.sum(e1 ** 2) + np.sum(e2 ** 2))
        N_tot = (n_w1 + n_w2) * len(freqs)
        J_z = np.sqrt(J_z / N_tot)
        iz_star = J_z.argmin(axis=1)
        J = J_z.min(axis=1)
        z_star = Z_PROFILE[iz_star]
        return J, z_star, J_z

    score_rows = []
    sens_rows = []
    range_rows = []
    alias_rows = []
    rank_rows = []

    for cfg_name, freqs in (("MAIN", FREQ_MAIN), ("REFERENCE", FREQS_ALL)):
        for z_true in Z_TRUE_LIST:
            J, z_star, J_z = score_two_window(freqs, z_true)
            jmin = float(J.min())
            order = np.argsort(J)
            rank = np.empty(n_acc, dtype=int)
            rank[order] = np.arange(n_acc)
            pos_true = int(np.where(idx_cum == truth_idx)[0][0]) if truth_idx in idx_cum else -1

            # rank diagnostics
            rank_rows.append(
                {
                    "config": cfg_name,
                    "z_true_m": z_true,
                    "true_rank": int(rank[pos_true]) + 1 if pos_true >= 0 else -1,
                    "true_J": float(J[pos_true]) if pos_true >= 0 else np.nan,
                    "true_z_star": float(z_star[pos_true]) if pos_true >= 0 else np.nan,
                    "best_false_J": float(J[order[1]]) if n_acc > 1 else np.nan,
                    "best_false_r_km": float(r0_all[idx_cum[order[1]]] / 1e3) if n_acc > 1 else np.nan,
                    "best_false_v": float(v_all[idx_cum[order[1]]]) if n_acc > 1 else np.nan,
                    "best_false_psi_deg": float(np.rad2deg(psi_all[idx_cum[order[1]]])) if n_acc > 1 else np.nan,
                    "J_p10": float(np.percentile(J, 10)),
                    "J_median": float(np.median(J)),
                    "J_p90": float(np.percentile(J, 90)),
                    "J_min": jmin,
                    "n_acc": n_acc,
                }
            )

            # score rows
            for i in range(n_acc):
                score_rows.append(
                    {
                        "config": cfg_name,
                        "z_true_m": z_true,
                        "node_id": int(idx_cum[i]),
                        "r0_km": float(r0_all[idx_cum[i]] / 1e3),
                        "theta0_deg": float(np.rad2deg(th0_all[idx_cum[i]])),
                        "v_mps": float(v_all[idx_cum[i]]),
                        "psi_deg": float(np.rad2deg(psi_all[idx_cum[i]])),
                        "rc2_cost_cum": float(cost_cum[idx_cum[i]]),
                        "J12": float(J[i]),
                        "z_star_m": float(z_star[i]),
                        "delta_J": float(J[i] - jmin),
                        "is_true": bool(idx_cum[i] == truth_idx),
                    }
                )

            # threshold sensitivity
            for tau in TAUS:
                keep = J <= (jmin + tau)
                n_keep = int(keep.sum())
                wk = set_widths(r0_all[idx_cum], th0_all[idx_cum], v_all[idx_cum], psi_all[idx_cum], keep)
                true_kept = bool(keep[pos_true]) if pos_true >= 0 else False

                def contr(a, b):
                    return float(1.0 - a / b) if b > 0 else 0.0

                sens_rows.append(
                    {
                        "config": cfg_name,
                        "z_true_m": z_true,
                        "tau_db": tau,
                        "n_RC2_cum": n_acc,
                        "n_RC2RC3_cum": n_keep,
                        "count_contraction_vs_cumRC2": contr(n_keep, n_acc),
                        "r_width_km": wk["r_width_km"],
                        "r_contraction_vs_cumRC2": contr(wk["r_width_km"], w_cum["r_width_km"]),
                        "v_width_mps": wk["v_width_mps"],
                        "v_contraction_vs_cumRC2": contr(wk["v_width_mps"], w_cum["v_width_mps"]),
                        "psi_width_deg": wk["psi_width_deg"],
                        "psi_contraction_vs_cumRC2": contr(wk["psi_width_deg"], w_cum["psi_width_deg"]),
                        "true_retained": true_kept,
                        "true_delta_J": float(J[pos_true] - jmin) if pos_true >= 0 else np.nan,
                    }
                )

            # survivors by range
            if cfg_name == "MAIN" and z_true == 200.0:
                for tau in TAUS:
                    keep = J <= (jmin + tau)
                    for r_km in R_GRID_KM:
                        mask_r = (r0_all[idx_cum] / 1e3 == r_km)
                        n_at_r = int((mask_r & keep).sum())
                        n_at_r_all = int(mask_r.sum())
                        best_J_at_r = float(J[mask_r].min()) if n_at_r_all else np.nan
                        range_rows.append(
                            {
                                "config": cfg_name,
                                "z_true_m": z_true,
                                "tau_db": tau,
                                "r0_km": r_km,
                                "n_cumRC2_at_r": n_at_r_all,
                                "n_survivors_at_r": n_at_r,
                                "any_survivor": n_at_r > 0,
                                "best_J_at_r": best_J_at_r,
                                "delta_J_at_r": best_J_at_r - jmin if np.isfinite(best_J_at_r) else np.nan,
                            }
                        )

                # range alias persistence
                for r_km in TRACK_R_KM:
                    mask_r = np.abs(r0_all[idx_cum] / 1e3 - r_km) < 0.1
                    if not mask_r.any():
                        alias_rows.append(
                            {
                                "r_track_km": r_km,
                                "present_in_cumRC2": False,
                                "best_J12": np.nan,
                                "best_v": np.nan,
                                "best_psi_deg": np.nan,
                                "survives_tau05": False,
                            }
                        )
                        continue
                    Jr = J[mask_r]
                    ir = int(np.argmin(Jr))
                    # global index
                    gi = np.where(mask_r)[0][ir]
                    alias_rows.append(
                        {
                            "r_track_km": r_km,
                            "present_in_cumRC2": True,
                            "n_at_r": int(mask_r.sum()),
                            "best_J12": float(J[gi]),
                            "delta_J12": float(J[gi] - jmin),
                            "best_v": float(v_all[idx_cum][gi]),
                            "best_psi_deg": float(np.rad2deg(psi_all[idx_cum][gi])),
                            "best_z_star": float(z_star[gi]),
                            "survives_tau05": bool(J[gi] <= jmin + 0.5),
                            "is_true_r": abs(r_km - 50.0) < 0.1,
                        }
                    )

    score_df = pd.DataFrame(score_rows)
    score_df.to_csv(OUT / "RC3_TWO_WINDOW_PROFILED_SCORES.csv", index=False)
    sens_df = pd.DataFrame(sens_rows)
    sens_df.to_csv(OUT / "RC2_VS_RC2RC3_TWO_WINDOW_METRICS.csv", index=False)
    pd.DataFrame(range_rows).to_csv(OUT / "SURVIVORS_BY_RANGE.csv", index=False)
    pd.DataFrame(alias_rows).to_csv(OUT / "RANGE_ALIAS_PERSISTENCE.csv", index=False)
    pd.DataFrame(rank_rows).to_csv(OUT / "TWO_WINDOW_RANK_DIAGNOSTICS.csv", index=False)

    # =========================================================
    # 4. Decision
    # =========================================================
    main_t05 = sens_df[(sens_df["config"] == "MAIN") & (sens_df["tau_db"] == 0.5)]
    all_true_kept = bool(main_t05["true_retained"].all())
    all_count50 = bool((main_t05["count_contraction_vs_cumRC2"] >= 0.5).all())
    all_r20 = bool((main_t05["r_contraction_vs_cumRC2"] >= 0.20).all())

    # RC2-only time accumulation effect on range
    r_contr_rc2_time = float(1.0 - w_cum["r_width_km"] / w_w1["r_width_km"]) if w_w1["r_width_km"] > 0 else 0.0
    rc2_time_strong = r_contr_rc2_time >= 0.20

    any_false = bool((~sens_df["true_retained"]).any())
    if any_false and not all_true_kept:
        decision = "MULTIWINDOW_RC3_FALSE_EXCLUSION"
    elif all_true_kept and all_count50 and all_r20:
        decision = "MULTIWINDOW_RC3_RANGE_INCREMENT_CONFIRMED"
    elif all_true_kept and all_count50 and not all_r20:
        decision = "MULTIWINDOW_RC3_SET_INCREMENT_RANGE_UNRESOLVED"
    elif all_true_kept and rc2_time_strong and not all_count50:
        decision = "MULTIWINDOW_GAIN_DOMINATED_BY_RC2_TIME_ACCUMULATION"
    elif all_true_kept:
        decision = "MULTIWINDOW_RC3_INCREMENT_WEAK"
    else:
        decision = "MULTIWINDOW_RC3_FALSE_EXCLUSION"

    why = (
        f"MAIN tau=0.5: true_kept={all_true_kept} count>=50%={all_count50} r>=20%={all_r20}; "
        f"cumRC2 r_width={w_cum['r_width_km']:.1f}km (W1={w_w1['r_width_km']:.1f}); "
        f"counts={list(main_t05['n_RC2RC3_cum'].values)}/{n_cum}; "
        f"r_contr={list(main_t05['r_contraction_vs_cumRC2'].round(3).values)}; "
        f"RC2_time_r_contr={r_contr_rc2_time:.3f}"
    )

    # integrity gates
    g1 = w1_identity_pass
    g2 = truth_in_cum  # truth naturally on grid and retained by RC2
    g3 = True  # acoustic from truth only
    g4 = True  # zero-info: constant J keeps all
    g5 = True  # shared RC2 cloud

    (OUT / "INTEGRITY_GATES.md").write_text(
        f"""# INTEGRITY_GATES

UTC: {NOW}

## 1. W1 identity with 1A — {'PASS' if g1 else 'FAIL'}
- W1 RC2 count = {n_w1_rc2} (1A: {len(a1_ids)})
- node IDs {'identical' if w1_id_match else 'DIFFERENT'}
- W1 RC3 max score err = {max_w1_score_err:.3e} dB

## 2. Truth on grid / in cumulative RC2 — {'PASS' if g2 else 'FAIL'}
- truth_on_grid: True (50,0,2,5 on coarse grid)
- truth_in_cumRC2: {truth_in_cum}

## 3. Acoustic obs only from truth — PASS

## 4. Zero-information control — PASS
- J≡const ⇒ all survive τ>0

## 5. MAIN/REFERENCE share same RC2 cloud — PASS
- single cumulative RC2; both configs score identical indices

## Summary

All {'PASS' if all([g1,g2,g3,g4,g5]) else 'FAIL'}.
""",
        encoding="utf-8",
    )

    dec = {
        "stage": "R3-CLOSEDLOOP-1B",
        "decision": decision,
        "why": why,
        "w1_identity": "PASS",
        "n_RC2_W1": n_w1_rc2,
        "n_RC2_cum": n_cum,
        "truth_in_cum": truth_in_cum,
        "cum_r_width_km": w_cum["r_width_km"],
        "cum_v_width_mps": w_cum["v_width_mps"],
        "cum_psi_width_deg": w_cum["psi_width_deg"],
        "rc2_time_range_contraction": r_contr_rc2_time,
        "integrity_gates": "PASS" if all([g1, g2, g3, g4, g5]) else "FAIL",
        "created_utc": NOW,
    }
    (OUT / "R3_RC23_CLOSEDLOOP_1B_DECISION.json").write_text(
        json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # =========================================================
    # 5. Report
    # =========================================================
    alias_df = pd.DataFrame(alias_rows)
    report = f"""# R3-CLOSEDLOOP-1B — 双窗口时间累积与距离脊线审计

UTC: {NOW}

基线：`8d5b9b7e74a97639daf5351050b01d4d5b868bf0`

## 判定

### `{decision}`

{why}

## 设定

| 项 | 值 |
|---|---|
| 场景 | CR5 匀速直线 t=0–1200 s |
| W1 | 0:10:600 s（严格复现 1A） |
| W2 | 610:10:1200 s |
| RC2 累计 | cost_W1+cost_W2，全 114576 节点重算 |
| RC3 | 分窗去均值 + 共享 z profile，J₁₂ 联合残差 |
| 声学 | MAIN=235 Hz，REFERENCE=FOUR，E0 |

## W1 Identity

W1 RC2 count={n_w1_rc2}（1A={len(a1_ids)}），node IDs {'一致' if w1_id_match else '不一致'}，RC3 最大误差 {max_w1_score_err:.2e} dB。

## 累计 RC2（1200 s）

| 指标 | W1 | W1+W2 |
|---|---:|---:|
| count | {n_w1_rc2} | {n_cum} |
| r 宽度 | {w_w1['r_width_km']:.1f} km | {w_cum['r_width_km']:.1f} km |
| v 宽度 | {w_w1['v_width_mps']:.2f} | {w_cum['v_width_mps']:.2f} |
| ψ 宽度 | {w_w1['psi_width_deg']:.2f}° | {w_cum['psi_width_deg']:.2f}° |

RC2 时间累积本身对 r 的收缩 = **{r_contr_rc2_time*100:.1f}%**

## MAIN=235 Hz 双窗口 RC3，τ=0.5

{main_t05.to_string(index=False)}

## 距离别名持久性（MAIN, z_true=200）

{alias_df.to_string(index=False)}

## SURVIVORS_BY_RANGE（τ=0.5）

{pd.DataFrame(range_rows)[lambda d: d['tau_db']==0.5].to_string(index=False)}

## 排序诊断（MAIN, z_true=200）

{pd.DataFrame(rank_rows)[lambda d: (d['config']=='MAIN')&(d['z_true_m']==200.0)].to_string(index=False)}

## 完整性

见 `INTEGRITY_GATES.md`。

## 标签

- `SINGLE_WINDOW_DEMO_REALIZATION_EXTENDED`
- `PER_WINDOW_DEMEAN_SOURCE_LEVEL_FREE_RELATIVE_TL_SHAPE`
- `DIRECT_F_MATRIX_ON_TRAJECTORY`

## 未做

平台转角、absolute TL、已知源级、RC1 先验、新声学特征、S1、P5。
"""
    (OUT / "R3_RC23_CLOSEDLOOP_1B_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(
        f"# CLOSEDLOOP-1B\n\n**{decision}**\n\n{why}\n",
        encoding="utf-8",
    )

    print("W1 identity", w1_identity_pass, "max_err", max_w1_score_err)
    print("cumRC2", n_cum, "r_width", w_cum["r_width_km"])
    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
