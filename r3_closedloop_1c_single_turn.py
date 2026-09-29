#!/usr/bin/env python3
"""R3-CLOSEDLOOP-1C: single small platform turn virtual baseline vs range ridge.

Baseline b8fb1ca. Turn at t=600s absolute. delta in {0,2,5,10,15} deg.
delta=0 must reproduce 1B. Paired noise (same seed for all branches).
Each branch re-runs full 114576 RC2, then RC3 on its own accepted cloud.
"""
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_SINGLE_TURN_RANGE_ANCHOR"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
A1B = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TWO_WINDOW_TEMPORAL_CONTRACTION"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0
R_MIN_M, R_MAX_M = 45e3, 60e3
R_GRID_KM = np.arange(45.0, 60.0 + 1e-9, 1.0)
TH_GRID_DEG = np.arange(-5.0, 5.0 + 1e-9, 0.5)
V_GRID = np.arange(1.0, 3.0 + 1e-9, 0.2)
PSI_GRID_DEG = np.arange(-15.0, 15.0 + 1e-9, 1.0)
N_NODES = int(len(R_GRID_KM) * len(TH_GRID_DEG) * len(V_GRID) * len(PSI_GRID_DEG))

T_TURN = 600.0
T_END = 1200.0
DT_OBS = 10.0
SIGMA_DEG = 0.1
RNG_SEED = 20260912
TRUTH = dict(r0_m=50e3, theta0_deg=0.0, v=2.0, psi_deg=5.0)

Z_TRUE_LIST = [180.0, 200.0, 220.0]
Z_PROFILE = np.arange(150.0, 250.0 + 1e-9, 5.0)
FREQS_ALL = [201.0, 235.0, 283.0, 338.0]
FREQ_MAIN = [235.0]
TAUS = [0.25, 0.5, 1.0]
ZR = 200.0
DELTAS = [0.0, 2.0, 5.0, 10.0, 15.0]
TRACK_R_KM = [45.0, 50.0, 56.0, 58.0, 60.0]


# ---- platform with absolute turn time ----
def platform_states_turn(t, delta_deg, t_turn=T_TURN, u=U_PLAT):
    t = np.asarray(t, dtype=float)
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


def predict_bearings(t, r0, th0, v, psi, delta_deg):
    xp, yp = platform_states_turn(t, delta_deg)
    xt, yt = target_states(t, r0, th0, v, psi)
    return np.arctan2(yt - yp, xt - xp)


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


def occupied_r_bins(r0, mask):
    return int(len(np.unique(np.round(r0[mask] / 1e3).astype(int))))


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
    t_full = np.arange(0.0, T_END + 1e-9, DT_OBS)
    n_w1 = int(np.sum(t_full <= T_TURN + 1e-9))  # 61
    n_w2 = t_full.size - n_w1  # 60
    t_w1 = t_full[:n_w1]
    t_w2 = t_full[n_w1:]
    sigma_rad = np.deg2rad(SIGMA_DEG)
    thr = chi2_thr(sigma_rad)

    # grid
    rr, tt, vv, pp = np.meshgrid(
        R_GRID_KM * 1e3, np.deg2rad(TH_GRID_DEG), V_GRID, np.deg2rad(PSI_GRID_DEG), indexing="ij"
    )
    r0_all = rr.ravel()
    th0_all = tt.ravel()
    v_all = vv.ravel()
    psi_all = pp.ravel()
    assert r0_all.size == N_NODES

    truth_idx = int(
        np.argmin(
            (r0_all - TRUTH["r0_m"]) ** 2
            + (th0_all - np.deg2rad(TRUTH["theta0_deg"])) ** 2
            + (v_all - TRUTH["v"]) ** 2
            + (psi_all - np.deg2rad(TRUTH["psi_deg"])) ** 2
        )
    )

    # shared noise
    rng = np.random.default_rng(RNG_SEED)
    noise = rng.normal(0.0, sigma_rad, size=t_full.shape)

    mods = {f: parse_mod(ZGRID / f"zgrid_f{int(f)}.mod") for f in FREQS_ALL}
    n_z = Z_PROFILE.size

    config = {
        "stage": "R3-CLOSEDLOOP-1C",
        "created_utc": NOW,
        "baseline_commit": "b8fb1ca5cc3b4f78b7820e7196aff15094db50d7",
        "turn_s": T_TURN,
        "T_end_s": T_END,
        "delta_deg_list": DELTAS,
        "paired_noise_seed": RNG_SEED,
        "rc2_cost": "sum wrap^2, accept cost<=cmin+13.3*sigma^2",
        "rc3": "per-window demean, shared z, DIRECT_F_MATRIX_ON_TRAJECTORY",
        "configs": {"MAIN": FREQ_MAIN, "REFERENCE": FREQS_ALL},
        "tau_db": TAUS,
    }
    (OUT / "CLOSEDLOOP_1C_CONFIG.json").write_text(
        json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # =========================================================
    # Platform turn geometry integrity
    # =========================================================
    geom_rows = []
    xp0, yp0 = platform_states_turn(t_full, 0.0)
    for delta in DELTAS:
        xp, yp = platform_states_turn(t_full, delta)
        w1_identical = bool(np.allclose(xp[:n_w1], xp0[:n_w1], atol=1e-12) and np.allclose(yp[:n_w1], yp0[:n_w1], atol=1e-12))
        # t=600 position
        i600 = n_w1 - 1
        pos600_ok = bool(abs(xp[i600] - U_PLAT * T_TURN) < 1e-9 and abs(yp[i600]) < 1e-9)
        # W2 lateral displacement
        b_perp = abs(yp[-1] - yp[i600])
        expected_b = U_PLAT * (T_END - T_TURN) * abs(math.sin(math.radians(delta)))
        b_ok = abs(b_perp - expected_b) < 1e-9
        geom_rows.append(
            {
                "delta_deg": delta,
                "W1_platform_identical_to_zero": w1_identical,
                "t600_position_ok": pos600_ok,
                "B_perp_m": b_perp,
                "B_perp_expected_m": expected_b,
                "B_perp_ok": b_ok,
                "delta_theta_geom_deg": math.degrees(math.atan2(b_perp, 50e3)),
            }
        )
    geom_df = pd.DataFrame(geom_rows)
    geom_df.to_csv(OUT / "PLATFORM_TURN_GEOMETRY_INTEGRITY.csv", index=False)

    # =========================================================
    # Zero-turn identity with 1B
    # =========================================================
    # load 1B reference
    b_scores = pd.read_csv(A1B / "RC3_TWO_WINDOW_PROFILED_SCORES.csv")
    b_cum = pd.read_csv(A1B / "RC2_CUMULATIVE_1200_ACCEPTED.csv")
    b_ids = set(b_cum["node_id"].tolist())

    # delta=0 observations
    theta_true_0 = predict_bearings(t_full, TRUTH["r0_m"], np.deg2rad(TRUTH["theta0_deg"]), TRUTH["v"], np.deg2rad(TRUTH["psi_deg"]), 0.0)
    obs_0 = theta_true_0 + noise

    # delta=0 RC2 cumulative
    pred_0 = predict_bearings_batch(t_full, r0_all, th0_all, v_all, psi_all, 0.0)
    cost_0 = np.sum(wrap_angle(pred_0 - obs_0[None, :]) ** 2, axis=1)
    cmin_0 = float(np.min(cost_0))
    acc_0 = cost_0 <= (cmin_0 + thr)
    n_cum_0 = int(acc_0.sum())
    ids_0 = set(np.where(acc_0)[0].tolist())
    zero_id_match = ids_0 == b_ids
    zero_count_match = n_cum_0 == len(b_ids)

    # W1 identity across all deltas
    w1_bear_ok = True
    for delta in DELTAS:
        th_d = predict_bearings(t_full, TRUTH["r0_m"], np.deg2rad(TRUTH["theta0_deg"]), TRUTH["v"], np.deg2rad(TRUTH["psi_deg"]), delta)
        if not np.allclose(th_d[:n_w1], theta_true_0[:n_w1], atol=1e-15):
            w1_bear_ok = False
        obs_d = th_d + noise
        if not np.allclose(obs_d[:n_w1], obs_0[:n_w1], atol=1e-15):
            w1_bear_ok = False

    # RC3 identity: sample check for delta=0
    rng_id = np.random.default_rng(0)
    idx0 = np.where(acc_0)[0]
    sample = sorted(set(idx0[rng_id.choice(len(idx0), size=min(150, len(idx0)), replace=False)].tolist() + [truth_idx]))
    max_score_err = 0.0
    for cfg_name, freqs in (("MAIN", FREQ_MAIN), ("REFERENCE", FREQS_ALL)):
        for z_true in Z_TRUE_LIST:
            # obs
            r_true_0 = range_traj(t_full, TRUTH["r0_m"], TRUTH["theta0_deg"], TRUTH["v"], TRUTH["psi_deg"], 0.0)
            obs_w1_parts, obs_w2_parts = [], []
            for f in freqs:
                F_t = F_matrix(mods[f], r_true_0)
                L_t = L_profile(mods[f], F_t, [z_true])[0]
                obs_w1_parts.append(demean(L_t[:n_w1]))
                obs_w2_parts.append(demean(L_t[n_w1:]))
            J_s = np.empty(len(sample))
            for si, j in enumerate(sample):
                r_c = range_traj(t_full, r0_all[j], np.rad2deg(th0_all[j]), v_all[j], np.rad2deg(psi_all[j]), 0.0)
                J_z = np.zeros(n_z)
                for f in freqs:
                    F_i = F_matrix(mods[f], r_c)
                    L_all = L_profile(mods[f], F_i, Z_PROFILE)
                    for iz in range(n_z):
                        e1 = demean(L_all[iz, :n_w1]) - obs_w1_parts[freqs.index(f)]
                        e2 = demean(L_all[iz, n_w1:]) - obs_w2_parts[freqs.index(f)]
                        J_z[iz] += np.sum(e1**2) + np.sum(e2**2)
                J_z = np.sqrt(J_z / ((n_w1 + n_w2) * len(freqs)))
                J_s[si] = J_z.min()
            b_sub = b_scores[(b_scores["config"] == cfg_name) & (b_scores["z_true_m"] == z_true)]
            b_map = dict(zip(b_sub["node_id"], b_sub["J12"]))
            for si, j in enumerate(sample):
                if j in b_map:
                    max_score_err = max(max_score_err, abs(J_s[si] - b_map[j]))

    zero_score_ok = max_score_err < 1e-8
    zero_pass = zero_id_match and zero_count_match and zero_score_ok and w1_bear_ok

    pd.DataFrame(
        [
            {
                "w1_bearing_identical_all_deltas": w1_bear_ok,
                "zero_count": n_cum_0,
                "zero_count_match_1B": zero_count_match,
                "zero_id_match_1B": zero_id_match,
                "zero_score_max_err": max_score_err,
                "zero_score_ok": zero_score_ok,
                "gate_pass": zero_pass,
            }
        ]
    ).to_csv(OUT / "ZERO_TURN_IDENTITY.csv", index=False)

    if not zero_pass:
        (OUT / "R3_RC23_CLOSEDLOOP_1C_DECISION.json").write_text(
            json.dumps(
                {
                    "stage": "R3-CLOSEDLOOP-1C",
                    "decision": "CLOSEDLOOP_1C_BLOCKED_BY_ZERO_TURN_IDENTITY",
                    "w1_bearing_ok": bool(w1_bear_ok),
                    "zero_count_match": bool(zero_count_match),
                    "zero_id_match": bool(zero_id_match),
                    "zero_score_ok": bool(zero_score_ok),
                    "max_score_err": float(max_score_err),
                    "created_utc": NOW,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        print("BLOCKED_BY_ZERO_TURN_IDENTITY")
        print("w1_bear", w1_bear_ok, "count", zero_count_match, "id", zero_id_match, "score", zero_score_ok, "err", max_score_err)
        return 1

    # =========================================================
    # Per-delta RC2 + RC3
    # =========================================================
    rc2_rows = []
    rc23_rows = []
    alias_rows = []
    rank_rows = []
    surv_rows = []
    base_rows = []

    for delta in DELTAS:
        th_true = predict_bearings(t_full, TRUTH["r0_m"], np.deg2rad(TRUTH["theta0_deg"]), TRUTH["v"], np.deg2rad(TRUTH["psi_deg"]), delta)
        obs = th_true + noise

        pred = predict_bearings_batch(t_full, r0_all, th0_all, v_all, psi_all, delta)
        cost = np.sum(wrap_angle(pred - obs[None, :]) ** 2, axis=1)
        cmin = float(np.min(cost))
        acc = cost <= (cmin + thr)
        n_acc = int(acc.sum())
        idx_acc = np.where(acc)[0]
        truth_in = bool(acc[truth_idx])

        # RC2-only metrics
        w = set_widths(r0_all, th0_all, v_all, psi_all, acc)
        n_r_bins = occupied_r_bins(r0_all, acc)
        rc2_rows.append(
            {
                "delta_deg": delta,
                "n_RC2": n_acc,
                "r_width_km": w["r_width_km"],
                "r_bins_occupied": n_r_bins,
                "v_width_mps": w["v_width_mps"],
                "psi_width_deg": w["psi_width_deg"],
                "theta_width_deg": w["theta_width_deg"],
                "true_retained": truth_in,
            }
        )

        # virtual baseline
        xp, yp = platform_states_turn(t_full, delta)
        b_perp = abs(yp[-1] - yp[n_w1 - 1])
        base_rows.append(
            {
                "delta_deg": delta,
                "B_perp_m": b_perp,
                "delta_theta_geom_deg": math.degrees(math.atan2(b_perp, 50e3)),
            }
        )

        # RC3 two-window
        # obs L
        r_true = range_traj(t_full, TRUTH["r0_m"], TRUTH["theta0_deg"], TRUTH["v"], TRUTH["psi_deg"], delta)
        L_obs_w1 = {}
        L_obs_w2 = {}
        for f in FREQS_ALL:
            F_t = F_matrix(mods[f], r_true)
            L_t = L_profile(mods[f], F_t, Z_TRUE_LIST)
            for i, zt in enumerate(Z_TRUE_LIST):
                L_obs_w1[(zt, f)] = demean(L_t[i, :n_w1])
                L_obs_w2[(zt, f)] = demean(L_t[i, n_w1:])

        pos_true = int(np.where(idx_acc == truth_idx)[0][0]) if truth_idx in idx_acc else -1

        for cfg_name, freqs in (("MAIN", FREQ_MAIN), ("REFERENCE", FREQS_ALL)):
            for z_true in Z_TRUE_LIST:
                n = idx_acc.size
                J_z = np.zeros((n, n_z))
                for i in range(n):
                    j = idx_acc[i]
                    r_c = range_traj(t_full, r0_all[j], np.rad2deg(th0_all[j]), v_all[j], np.rad2deg(psi_all[j]), delta)
                    r1, r2 = r_c[:n_w1], r_c[n_w1:]
                    for f in freqs:
                        F1 = F_matrix(mods[f], r1)
                        F2 = F_matrix(mods[f], r2)
                        L1 = L_profile(mods[f], F1, Z_PROFILE)
                        L2 = L_profile(mods[f], F2, Z_PROFILE)
                        o1 = L_obs_w1[(z_true, f)]
                        o2 = L_obs_w2[(z_true, f)]
                        for iz in range(n_z):
                            e1 = demean(L1[iz]) - o1
                            e2 = demean(L2[iz]) - o2
                            J_z[i, iz] += np.sum(e1**2) + np.sum(e2**2)
                J_z = np.sqrt(J_z / ((n_w1 + n_w2) * len(freqs)))
                J = J_z.min(axis=1)
                z_star = Z_PROFILE[J_z.argmin(axis=1)]
                jmin = float(J.min())
                order = np.argsort(J)
                rank = np.empty(n, dtype=int)
                rank[order] = np.arange(n)

                rank_rows.append(
                    {
                        "delta_deg": delta,
                        "config": cfg_name,
                        "z_true_m": z_true,
                        "true_rank": int(rank[pos_true]) + 1 if pos_true >= 0 else -1,
                        "true_J": float(J[pos_true]) if pos_true >= 0 else np.nan,
                        "best_false_J": float(J[order[1]]) if n > 1 else np.nan,
                        "best_false_r_km": float(r0_all[idx_acc[order[1]]] / 1e3) if n > 1 else np.nan,
                        "J_min": jmin,
                        "n_acc": n,
                    }
                )

                for tau in TAUS:
                    keep = J <= (jmin + tau)
                    wk = set_widths(r0_all[idx_acc], th0_all[idx_acc], v_all[idx_acc], psi_all[idx_acc], keep)
                    nb = occupied_r_bins(r0_all[idx_acc], keep)
                    true_kept = bool(keep[pos_true]) if pos_true >= 0 else False
                    r_min = float(r0_all[idx_acc][keep].min() / 1e3) if keep.any() else np.nan
                    r_max = float(r0_all[idx_acc][keep].max() / 1e3) if keep.any() else np.nan
                    rc23_rows.append(
                        {
                            "delta_deg": delta,
                            "config": cfg_name,
                            "z_true_m": z_true,
                            "tau_db": tau,
                            "n_RC2": n,
                            "n_RC2RC3": int(keep.sum()),
                            "r_min_km": r_min,
                            "r_max_km": r_max,
                            "r_width_km": wk["r_width_km"],
                            "r_bins_occupied": nb,
                            "v_width_mps": wk["v_width_mps"],
                            "psi_width_deg": wk["psi_width_deg"],
                            "true_retained": true_kept,
                            "true_rank": int(rank[pos_true]) + 1 if pos_true >= 0 else -1,
                        }
                    )

                # range alias + survivors for MAIN z=200
                if cfg_name == "MAIN" and z_true == 200.0:
                    for r_km in TRACK_R_KM:
                        mask_r = np.abs(r0_all[idx_acc] / 1e3 - r_km) < 0.1
                        if not mask_r.any():
                            alias_rows.append(
                                {
                                    "delta_deg": delta,
                                    "r_track_km": r_km,
                                    "present_in_RC2": False,
                                    "n_at_r": 0,
                                    "best_J12": np.nan,
                                    "best_z_star": np.nan,
                                    "best_v": np.nan,
                                    "best_psi_deg": np.nan,
                                    "survives_tau05": False,
                                }
                            )
                            continue
                        Jr = J[mask_r]
                        gi = np.where(mask_r)[0][int(np.argmin(Jr))]
                        alias_rows.append(
                            {
                                "delta_deg": delta,
                                "r_track_km": r_km,
                                "present_in_RC2": True,
                                "n_at_r": int(mask_r.sum()),
                                "best_J12": float(J[gi]),
                                "delta_J12": float(J[gi] - jmin),
                                "best_z_star": float(z_star[gi]),
                                "best_v": float(v_all[idx_acc][gi]),
                                "best_psi_deg": float(np.rad2deg(psi_all[idx_acc][gi])),
                                "survives_tau05": bool(J[gi] <= jmin + 0.5),
                                "is_true_r": abs(r_km - 50.0) < 0.1,
                            }
                        )
                    # survivors by range at tau=0.5
                    for tau in (0.5,):
                        keep = J <= (jmin + tau)
                        for r_km in R_GRID_KM:
                            mask_r = np.abs(r0_all[idx_acc] / 1e3 - r_km) < 0.1
                            surv_rows.append(
                                {
                                    "delta_deg": delta,
                                    "tau_db": tau,
                                    "r0_km": r_km,
                                    "n_survivors": int((mask_r & keep).sum()),
                                    "any_survivor": bool((mask_r & keep).any()),
                                }
                            )

    rc2_df = pd.DataFrame(rc2_rows)
    rc2_df.to_csv(OUT / "RC2_TURN_ONLY_SET_METRICS.csv", index=False)
    rc23_df = pd.DataFrame(rc23_rows)
    rc23_df.to_csv(OUT / "RC2_RC3_TURN_SET_METRICS.csv", index=False)
    pd.DataFrame(alias_rows).to_csv(OUT / "TURN_RANGE_ALIAS_PERSISTENCE.csv", index=False)
    pd.DataFrame(rank_rows).to_csv(OUT / "TURN_RANK_DIAGNOSTICS.csv", index=False)
    pd.DataFrame(surv_rows).to_csv(OUT / "SURVIVORS_BY_RANGE_AND_TURN.csv", index=False)
    pd.DataFrame(base_rows).to_csv(OUT / "VIRTUAL_BASELINE_DIAGNOSTIC.csv", index=False)

    # =========================================================
    # Decision
    # =========================================================
    main05 = rc23_df[(rc23_df["config"] == "MAIN") & (rc23_df["tau_db"] == 0.5)]
    base0 = main05[main05["delta_deg"] == 0.0].set_index("z_true_m")
    found_delta = None
    min_delta = None
    for delta in [2.0, 5.0, 10.0, 15.0]:
        sub = main05[main05["delta_deg"] == delta].set_index("z_true_m")
        if not all(zt in sub.index for zt in Z_TRUE_LIST):
            continue
        all_true = all(sub.loc[zt, "true_retained"] for zt in Z_TRUE_LIST)
        r_ok = all(
            (1.0 - sub.loc[zt, "r_width_km"] / base0.loc[zt, "r_width_km"]) >= 0.20
            if base0.loc[zt, "r_width_km"] > 0
            else True
            for zt in Z_TRUE_LIST
        )
        bins_ok = all(
            (1.0 - sub.loc[zt, "r_bins_occupied"] / base0.loc[zt, "r_bins_occupied"]) >= 0.25
            if base0.loc[zt, "r_bins_occupied"] > 0
            else True
            for zt in Z_TRUE_LIST
        )
        if all_true and r_ok and bins_ok:
            found_delta = delta
            if min_delta is None:
                min_delta = delta

    any_false = bool((~main05["true_retained"]).any())

    # auxiliary labels
    aux = []
    # RC2-only range contraction?
    rc2_base = rc2_df[rc2_df["delta_deg"] == 0.0].iloc[0]
    for delta in [2.0, 5.0, 10.0, 15.0]:
        row = rc2_df[rc2_df["delta_deg"] == delta].iloc[0]
        if rc2_base["r_width_km"] > 0 and (1.0 - row["r_width_km"] / rc2_base["r_width_km"]) >= 0.20:
            aux.append("TURN_GEOMETRY_RC2_RANGE_INCREMENT")
            break
    if found_delta is not None and "TURN_GEOMETRY_RC2_RANGE_INCREMENT" not in aux:
        aux.append("TURN_RC3_RANGE_SYNERGY")

    if any_false and not all(main05["true_retained"]):
        decision = "TURN_BRANCH_FALSE_EXCLUSION"
    elif found_delta is not None:
        decision = "SMALL_TURN_RANGE_ANCHOR_CONFIRMED"
    else:
        # check if candidates keep contracting but range doesn't
        decision = "SMALL_TURN_RANGE_RIDGE_PERSISTS"

    why = (
        f"MAIN tau=0.5: found_delta={found_delta}, min_delta={min_delta}, aux={aux}; "
        + "; ".join(
            f"δ={d}°: n={int(main05[main05['delta_deg']==d]['n_RC2RC3'].iloc[0])} "
            f"r_w={main05[main05['delta_deg']==d]['r_width_km'].iloc[0]:.1f} "
            f"bins={int(main05[main05['delta_deg']==d]['r_bins_occupied'].iloc[0])}"
            for d in DELTAS
        )
    )

    # integrity
    g1 = w1_bear_ok
    g2 = zero_pass
    g3 = bool(geom_df["W1_platform_identical_to_zero"].all())
    g4 = True

    (OUT / "INTEGRITY_GATES.md").write_text(
        f"""# INTEGRITY_GATES

UTC: {NOW}

## 1. W1 bearing identical across all deltas — {'PASS' if g1 else 'FAIL'}
## 2. delta=0 reproduces 1B — {'PASS' if g2 else 'FAIL'}
- cumRC2 count {n_cum_0} (1B: {len(b_ids)})
- node IDs {'match' if zero_id_match else 'DIFFER'}
- RC3 score max err {max_score_err:.2e} dB
## 3. W1 platform trajectory identical for all deltas — {'PASS' if g3 else 'FAIL'}
## 4. Paired noise (same seed) — PASS

All {'PASS' if all([g1,g2,g3,g4]) else 'FAIL'}.
""",
        encoding="utf-8",
    )

    dec = {
        "stage": "R3-CLOSEDLOOP-1C",
        "decision": decision,
        "why": why,
        "auxiliary_labels": aux,
        "min_tested_turn_with_range_contraction": min_delta,
        "zero_turn_identity": "PASS",
        "integrity_gates": "PASS" if all([g1, g2, g3, g4]) else "FAIL",
        "created_utc": NOW,
    }
    (OUT / "R3_RC23_CLOSEDLOOP_1C_DECISION.json").write_text(
        json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    # report
    report = f"""# R3-CLOSEDLOOP-1C — 单次小转向虚拟基线与距离脊线

UTC: {NOW}

基线：`b8fb1ca5cc3b4f78b7820e7196aff15094db50d7`

## 判定

### `{decision}`

{why}

辅助：{aux or '—'}

## 设定

| 项 | 值 |
|---|---|
| 转向时刻 | t=600 s（绝对时间） |
| 转角 δ | 0°, 2°, 5°, 10°, 15° |
| 噪声 | paired：同一 seed 20260912 |
| W1 | 0:10:600 s（所有 δ 平台轨迹严格一致） |
| W2 | 610:10:1200 s |

## δ=0 Identity

cumRC2={n_cum_0}（1B={len(b_ids)}），IDs {'一致' if zero_id_match else '不一致'}，RC3 max err {max_score_err:.2e} dB。

## RC2-only 各转角

{rc2_df.to_string(index=False)}

## RC2+RC3（MAIN, τ=0.5）

{main05.to_string(index=False)}

## 虚拟基线

{pd.DataFrame(base_rows).to_string(index=False)}

## 距离别名持久性（MAIN, z=200）

{pd.DataFrame(alias_rows).to_string(index=False)}

## 排序诊断（MAIN, z=200）

{pd.DataFrame(rank_rows)[lambda d: (d['config']=='MAIN') & (d['z_true_m']==200.0)].to_string(index=False)}

## 完整性

见 `INTEGRITY_GATES.md`。

## 未做

TDOA、RC1 先验、absolute TL、已知源级、新声学特征、S1、P5、Z型转向。
"""
    (OUT / "R3_RC23_CLOSEDLOOP_1C_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# CLOSEDLOOP-1C\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("zero_pass", zero_pass, "err", max_score_err)
    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
