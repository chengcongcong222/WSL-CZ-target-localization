#!/usr/bin/env python3
"""R3-CLOSEDLOOP-1E-B-FIX: Frequency-model integrity + real tracked-drift rerun.

5B env method: copy zgrid_f235.env, change title+freq only. Hard gates. No zero-fill.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_TRACKED_FREQ_ROBUSTNESS_FIX2"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
FIX = ROOT / "results" / "R3_RC23_CLOSEDLOOP" / "R3_RC23_TURN_FREQ_SUBSET_FIX"
K5B = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_LINE_INSTABILITY_BOUNDARY" / "_kfreq"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_kraken"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0; T_TURN = 600.0; T_END = 1200.0; DT_OBS = 10.0
SIGMA_DEG = 0.1; RNG_SEED = 20260912; DELTA = 15.0
TRUTH = dict(r0_m=50e3, theta0_deg=0.0, v=2.0, psi_deg=5.0)
R_GRID_KM = np.arange(45.0, 60.0 + 1e-9, 1.0)
TH_GRID_DEG = np.arange(-5.0, 5.0 + 1e-9, 0.5)
V_GRID = np.arange(1.0, 3.0 + 1e-9, 0.2)
PSI_GRID_DEG = np.arange(-15.0, 15.0 + 1e-9, 1.0)
N_NODES = int(len(R_GRID_KM) * len(TH_GRID_DEG) * len(V_GRID) * len(PSI_GRID_DEG))
Z_TRUE_LIST = [180.0, 200.0, 220.0]
Z_PROFILE = np.arange(150.0, 250.0 + 1e-9, 5.0)
ZR = 200.0; TAU = 0.5
CONFIGS = {"TRIPLE": [201.0, 235.0, 283.0], "FOUR": [201.0, 235.0, 283.0, 338.0]}
D_F_LIST = [0.0, 0.0025, 0.005, 0.01, 0.02]
Q_F = np.array([0.0, 0.5, 1.0, 0.5, 0.0, -0.5, -1.0, -0.5])
ALL_BASE_F = [201.0, 235.0, 283.0, 338.0]


def freq_snap(f):
    """Exact 4-decimal rounding per 5B-FIX. No grid quantization."""
    return round(float(f), 4)


def freq_ok(f0, D, q):
    """Check if exact frequency model exists."""
    fq = round(f0 * (1 + D * q), 4)
    return fq in mod_cache if 'mod_cache' in dir() else False


def find_model(fq, tol=0.001):
    """Tolerance-based model lookup to avoid float key mismatches."""
    if fq in mod_cache:
        return mod_cache[fq]
    for k, v in mod_cache.items():
        if abs(k - fq) < tol:
            return v
    return None


def write_env_freq(path: Path, freq: float, env0: Path):
    """5B method: copy zgrid env, change title + freq only."""
    lines = env0.read_text(encoding="utf-8").splitlines()
    out = []
    for i, ln in enumerate(lines):
        if i == 0:
            out.append("'PERT'")
        elif i == 1:
            out.append(f"{freq:.4f}")
        else:
            out.append(ln)
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


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


def F_mat(mod, r):
    kre, alpha = mod["k"].real, -mod["k"].imag
    r = np.asarray(r, float)
    amp = np.sqrt(2 * np.pi / (kre[:, None] * r[None, :]))
    return amp * np.exp(-1j * kre[:, None] * r[None, :] - alpha[:, None] * r[None, :] - 1j * np.pi / 4)


def L_prof(mod, F, z_list, zr=ZR):
    izr = int(np.argmin(np.abs(mod["depths"] - zr)))
    pr = mod["phi"][izr]
    W = mod["phi"] * pr[None, :]
    izs = [int(np.argmin(np.abs(mod["depths"] - z))) for z in z_list]
    return 20.0 * np.log10(np.maximum(np.abs(W[izs, :] @ F), 1e-30))


def demean(L):
    return np.asarray(L, float) - float(np.mean(L))


def platform_turn(t, delta_deg, t_turn=T_TURN, u=U_PLAT):
    t = np.asarray(t, float)
    xp, yp = np.empty_like(t), np.empty_like(t)
    d = math.radians(delta_deg)
    c, s = math.cos(d), math.sin(d)
    if abs(delta_deg) < 1e-12:
        xp[:] = u * t; yp[:] = 0.0; return xp, yp
    mask = t <= t_turn
    xp[mask] = u * t[mask]; yp[mask] = 0.0
    dt = t[~mask] - t_turn
    xp[~mask] = u * t_turn + u * dt * c; yp[~mask] = u * dt * s
    return xp, yp


def target_xy(t, r0_m, th0_rad, v, psi_rad):
    return r0_m * np.cos(th0_rad) + v * t * np.cos(psi_rad), r0_m * np.sin(th0_rad) + v * t * np.sin(psi_rad)


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def range_traj(t, r0_m, th0_deg, v, psi_deg, delta_deg):
    xp, yp = platform_turn(t, delta_deg)
    xt, yt = target_xy(t, r0_m, np.deg2rad(th0_deg), v, np.deg2rad(psi_deg))
    return np.hypot(xt - xp, yt - yp)


def main():
    config = {"stage": "R3-CLOSEDLOOP-1E-B-FIX", "created_utc": NOW,
              "baseline_commit": "8c134df8d962d7e557cdfe32ccea2402966ce181",
              "env_method": "5B write_env_freq: copy zgrid_f235.env, title+freq only"}
    (OUT / "CLOSEDLOOP_1E_B_FIX_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    (OUT / "1E_B_INVALIDATION_AUDIT.md").write_text(
        f"""# 1E_B_INVALIDATION_AUDIT

UTC: {NOW}

旧1E-B判定降级为 SUPERSEDED_INVALID_FREQUENCY_MODEL_MANIFEST。

原因:
1. 所有非标称漂移频率 model_exists=False
2. manifest 失败未停止执行
3. q!=0 时间片被静默置零
4. D>0 结果无效

保留旧输出作为失败审计。
""", encoding="utf-8")

    # =========================================================
    # 1. Build frequency models (5B method)
    # =========================================================
    env0 = ZGRID / "zgrid_f235.env"
    needed_freqs = set()
    for f0 in ALL_BASE_F:
        for D in D_F_LIST:
            for q in [-1, -0.5, 0, 0.5, 1]:
                needed_freqs.add(round(f0 * (1 + D * q), 4))
    needed_freqs = sorted(needed_freqs)

    # Exact schedule audit + exact vs snapped
    sched_rows = []
    cmp_rows = []
    t_full_sched = np.arange(0.0, T_END + 1e-9, DT_OBS)
    seg_sched = np.minimum((t_full_sched / T_END * 8).astype(int), 7)
    q_t_sched = Q_F[seg_sched]
    for f0 in ALL_BASE_F:
        for D in D_F_LIST:
            for q in [-1, -0.5, 0, 0.5, 1]:
                exact_f = round(f0 * (1 + D * q), 4)
                snapped_f = round(exact_f * 10.0) / 10.0
                sched_rows.append({"base_f": f0, "D": D, "q": q, "exact_actual_f": exact_f,
                                   "old_snapped_f": snapped_f, "abs_diff_hz": abs(exact_f - snapped_f)})
            for ti, t in enumerate(t_full_sched):
                q = q_t_sched[ti]
                exact_f = round(f0 * (1 + D * q), 4)
                sched_rows.append({"base_f": f0, "D": D, "time_s": t, "segment": seg_sched[ti],
                                   "q_t": q, "actual_f": exact_f})
    pd.DataFrame(sched_rows).to_csv(OUT / "DRIFT_SCHEDULE_5B_IDENTITY_FIXED2.csv", index=False)

    for f0 in ALL_BASE_F:
        for D in D_F_LIST:
            for q in [-1, -0.5, 0, 0.5, 1]:
                exact_f = round(f0 * (1 + D * q), 4)
                snapped_f = round(exact_f * 10.0) / 10.0
                intended_offset = f0 * D * q
                diff_frac = abs(exact_f - snapped_f) / abs(intended_offset) if abs(intended_offset) > 1e-12 else 0.0
                cmp_rows.append({"base_f": f0, "D": D, "q": q, "exact_f": exact_f,
                                 "old_f": snapped_f, "difference_hz": abs(exact_f - snapped_f),
                                 "difference_as_fraction_of_intended_offset": diff_frac})
    pd.DataFrame(cmp_rows).to_csv(OUT / "EXACT_VS_OLD_SNAPPED_FREQUENCY.csv", index=False)

    (OUT / "FREQUENCY_QUANTIZATION_AUDIT.md").write_text(
        f"""# FREQUENCY_QUANTIZATION_AUDIT

UTC: {NOW}

旧 freq_snap (0.1Hz grid) 已删除。
唯一合法: actual_f = round(f0 * (1 + D*q), 4) (5B-FIX exact).

最大量化偏差:
{pd.DataFrame(cmp_rows)["difference_hz"].max():.4f} Hz

最大频偏比例误差:
{pd.DataFrame(cmp_rows)["difference_as_fraction_of_intended_offset"].max():.3f}
""", encoding="utf-8")

    build_rows = []
    manifest_rows = []
    mod_cache = {}
    build_failed = False

    for fq in needed_freqs:
        key = f"f{fq:.4f}".replace(".", "p")
        modp = WORK / f"{key}.mod"
        is_base = fq in ALL_BASE_F

        # try 5B cache first for non-base
        reused = False
        if not is_base:
            k5b = K5B / f"{key}.mod"
            if k5b.exists():
                try:
                    mod_cache[fq] = parse_mod(k5b)
                    reused = True
                    build_rows.append({"actual_f": fq, "source": "5B_KFREQ", "reused_existing": True,
                                       "return_code": 0, "mod_exists": True, "mod_size": k5b.stat().st_size,
                                       "parse_ok": True, "n_modes": mod_cache[fq]["M"], "fallback_used": False})
                    manifest_rows.append({"base_f": fq if is_base else None, "actual_f": fq,
                                          "source_model": "5B_KFREQ", "model_exists": True,
                                          "parse_ok": True, "fallback_used": False})
                    continue
                except Exception:
                    pass

        if is_base:
            zp = ZGRID / f"zgrid_f{int(fq)}.mod"
            mod_cache[fq] = parse_mod(zp)
            build_rows.append({"actual_f": fq, "source": "ZGRID", "reused_existing": True,
                               "return_code": 0, "mod_exists": True, "mod_size": zp.stat().st_size,
                               "parse_ok": True, "n_modes": mod_cache[fq]["M"], "fallback_used": False})
            manifest_rows.append({"base_f": fq, "actual_f": fq, "source_model": "ZGRID",
                                  "model_exists": True, "parse_ok": True, "fallback_used": False})
            continue

        # generate with 5B method
        envp = WORK / f"{key}.env"
        write_env_freq(envp, fq, env0)
        try:
            proc = subprocess.run([str(AT_BIN / "kraken.exe"), key], cwd=str(WORK),
                                  capture_output=True, text=True, timeout=120)
            rc = proc.returncode
            stderr_tail = proc.stderr[-200:] if proc.stderr else ""
            stdout_tail = proc.stdout[-200:] if proc.stdout else ""
        except Exception as e:
            rc, stderr_tail, stdout_tail = -1, str(e), ""

        mod_exists = modp.exists()
        parse_ok = False
        n_modes = 0
        if mod_exists:
            try:
                m = parse_mod(modp)
                mod_cache[fq] = m
                n_modes = m["M"]
                parse_ok = n_modes > 0
            except Exception:
                parse_ok = False

        build_rows.append({"actual_f": fq, "source": "GENERATED", "reused_existing": reused,
                           "return_code": rc, "stderr_tail": stderr_tail, "stdout_tail": stdout_tail,
                           "mod_exists": mod_exists, "mod_size": modp.stat().st_size if mod_exists else 0,
                           "parse_ok": parse_ok, "n_modes": n_modes, "fallback_used": False})
        manifest_rows.append({"base_f": None, "actual_f": fq, "source_model": "GENERATED",
                              "model_exists": mod_exists, "parse_ok": parse_ok, "fallback_used": False})

        if rc != 0 or not mod_exists or not parse_ok:
            build_failed = True

    build_df = pd.DataFrame(build_rows)
    build_df.to_csv(OUT / "FREQUENCY_MODEL_BUILD_AUDIT.csv", index=False)
    man_df = pd.DataFrame(manifest_rows)
    man_df.to_csv(OUT / "FREQUENCY_MODEL_MANIFEST_FIXED.csv", index=False)

    if build_failed:
        (OUT / "R3_RC23_CLOSEDLOOP_1E_B_FIX_DECISION.json").write_text(
            json.dumps({"decision": "CLOSEDLOOP_1E_B_FIX_BLOCKED_BY_FREQUENCY_MODEL_BUILD",
                        "n_failed": int((~build_df["parse_ok"]).sum()), "created_utc": NOW}, indent=2),
            encoding="utf-8")
        print("BLOCKED_BY_FREQ_BUILD", int((~build_df["parse_ok"]).sum()), "failed")
        return 1

    # distinctness audit
    dist_rows = []
    for f0 in ALL_BASE_F:
        mod0 = mod_cache.get(f0)
        if mod0 is None:
            continue
        for fq in needed_freqs:
            if abs(fq - f0) < 1e-9:
                continue
            mod_q = mod_cache.get(fq)
            if mod_q is None:
                continue
            k0, kq = mod0["k"], mod_q["k"]
            n = min(len(k0), len(kq))
            k_diff = float(np.linalg.norm(k0[:n].real - kq[:n].real)) if n > 0 else 0.0
            dist_rows.append({"base_f": f0, "actual_f": fq, "base_M": mod0["M"], "actual_M": mod_q["M"],
                              "k_real_diff_norm": k_diff, "distinct": k_diff > 1e-12})
    pd.DataFrame(dist_rows).to_csv(OUT / "FREQUENCY_MODEL_DISTINCTNESS.csv", index=False)

    # =========================================================
    # 2. Setup RC2 + trajectories
    # =========================================================
    t_full = np.arange(0.0, T_END + 1e-9, DT_OBS)
    n_w1 = int(np.sum(t_full <= T_TURN + 1e-9))
    n_w2 = len(t_full) - n_w1
    seg = np.minimum((t_full / T_END * 8).astype(int), 7)
    q_t = Q_F[seg]
    sigma_rad = np.deg2rad(SIGMA_DEG)
    thr = 13.3 * sigma_rad ** 2

    rr, tt, vv, pp = np.meshgrid(R_GRID_KM * 1e3, np.deg2rad(TH_GRID_DEG), V_GRID, np.deg2rad(PSI_GRID_DEG), indexing="ij")
    r0_all, th0_all, v_all, psi_all = rr.ravel(), tt.ravel(), vv.ravel(), pp.ravel()

    rng = np.random.default_rng(RNG_SEED)
    noise = rng.normal(0.0, sigma_rad, size=t_full.shape)
    xp, yp = platform_turn(t_full, DELTA)
    xt, yt = target_xy(t_full, TRUTH["r0_m"], np.deg2rad(TRUTH["theta0_deg"]), TRUTH["v"], np.deg2rad(TRUTH["psi_deg"]))
    obs_bearing = np.arctan2(yt - yp, xt - xp) + noise

    pred = np.arctan2(
        r0_all[:, None] * np.sin(th0_all[:, None]) + v_all[:, None] * t_full[None, :] * np.sin(psi_all[:, None]) - yp[None, :],
        r0_all[:, None] * np.cos(th0_all[:, None]) + v_all[:, None] * t_full[None, :] * np.cos(psi_all[:, None]) - xp[None, :])
    cost = np.sum(wrap(pred - obs_bearing[None, :]) ** 2, axis=1)
    acc = cost <= (float(cost.min()) + thr)
    idx_acc = np.where(acc)[0]
    n_rc2 = len(idx_acc)
    truth_idx = int(np.argmin((r0_all - TRUTH["r0_m"]) ** 2 + (th0_all - np.deg2rad(TRUTH["theta0_deg"])) ** 2
                              + (v_all - TRUTH["v"]) ** 2 + (psi_all - np.deg2rad(TRUTH["psi_deg"])) ** 2))

    r_true = range_traj(t_full, TRUTH["r0_m"], TRUTH["theta0_deg"], TRUTH["v"], TRUTH["psi_deg"], DELTA)
    r_cand = {}
    for j in idx_acc:
        r_cand[j] = range_traj(t_full, r0_all[j], np.rad2deg(th0_all[j]), v_all[j], np.rad2deg(psi_all[j]), DELTA)

    # =========================================================
    # 3. Precompute acoustic L (all freqs, no fallback)
    # =========================================================
    n_z = Z_PROFILE.size
    L_cache = {}
    for fq in needed_freqs:
        mod = mod_cache.get(fq)
        if mod is None:
            raise RuntimeError(f"Missing model for {fq}")
        F_t = F_mat(mod, r_true)
        L_true = L_prof(mod, F_t, Z_TRUE_LIST)
        L_cache[fq] = {"true": {zt: L_true[i].copy() for i, zt in enumerate(Z_TRUE_LIST)}, "cand": {}}
        for j in idx_acc:
            F_c = F_mat(mod, r_cand[j])
            L_cache[fq]["cand"][j] = L_prof(mod, F_c, Z_PROFILE)

    # =========================================================
    # 4. D-effect gate
    # =========================================================
    d_eff_rows = []
    for f0 in ALL_BASE_F:
        L_d = {}
        for D in D_F_LIST:
            L_full = np.zeros(len(t_full))
            for qi, q in enumerate(q_t):
                fq = round(f0 * (1 + D * q), 4)
                # tolerance-based lookup
                m = None
                for k, v in L_cache.items():
                    if abs(k - fq) < 0.001:
                        m = v
                        break
                if m is None:
                    raise RuntimeError(f"Missing model {fq}")
                L_full[qi] = m["true"][200.0][qi]
            L_d[D] = L_full
        for D1 in D_F_LIST:
            for D2 in D_F_LIST:
                if D1 < D2:
                    rms_diff = float(np.sqrt(np.mean((L_d[D1] - L_d[D2]) ** 2)))
                    d_eff_rows.append({"base_f": f0, "D1": D1, "D2": D2, "rms_diff_db": rms_diff})
    d_eff_df = pd.DataFrame(d_eff_rows)
    d_eff_df.to_csv(OUT / "TRACKED_TL_D_EFFECT_AUDIT.csv", index=False)

    # D-effect gate: nonzero-D must differ from D=0
    nonzero_diffs = d_eff_df[d_eff_df.D1 == 0.0]["rms_diff_db"].values
    d_effect_ok = all(d > 1e-6 for d in nonzero_diffs) if len(nonzero_diffs) else False
    if not d_effect_ok:
        (OUT / "R3_RC23_CLOSEDLOOP_1E_B_FIX_DECISION.json").write_text(
            json.dumps({"decision": "CLOSEDLOOP_1E_B_FIX_BLOCKED_BY_D_EFFECT_GATE", "created_utc": NOW}, indent=2),
            encoding="utf-8")
        print("BLOCKED_BY_D_EFFECT_GATE")
        return 1

    # coverage gate
    cov_rows = []
    for f0 in ALL_BASE_F:
        for D in D_F_LIST:
            all_ok = True
            for q in [-1, -0.5, 0, 0.5, 1]:
                fq = round(f0 * (1 + D * q), 4)
                found = any(abs(k - fq) < 0.001 for k in L_cache.keys())
                if not found:
                    all_ok = False
            cov_rows.append({"base_f": f0, "D": D, "n_expected": len(t_full),
                             "n_filled": len(t_full) if all_ok else 0,
                             "fraction_filled": 1.0 if all_ok else 0.0})
    cov_df = pd.DataFrame(cov_rows)
    cov_df.to_csv(OUT / "TRACKED_INPUT_COVERAGE_FIXED2.csv", index=False)
    if not (cov_df["fraction_filled"] == 1.0).all():
        (OUT / "R3_RC23_CLOSEDLOOP_1E_B_FIX2_DECISION.json").write_text(
            json.dumps({"decision": "CLOSEDLOOP_1E_B_FIX2_BLOCKED_BY_INPUT_COVERAGE", "created_utc": NOW}, indent=2),
            encoding="utf-8")
        print("BLOCKED_BY_INPUT_COVERAGE")
        return 1

    # =========================================================
    # 5. Score (same as 1E-B but with hard no-fallback)
    # =========================================================
    cand_rows = []
    case_rows = []
    truth_id_rows = []

    for cfg_name in ["TRIPLE", "FOUR"]:
        freqs = CONFIGS[cfg_name]
        for D in D_F_LIST:
            for z_true in Z_TRUE_LIST:
                L_obs = {}
                for f0 in freqs:
                    L_full = np.zeros(len(t_full))
                    for qi, q in enumerate(q_t):
                        fq = round(f0 * (1 + D * q), 4)
                        m = None
                        for k, v in L_cache.items():
                            if abs(k - fq) < 0.001:
                                m = v
                                break
                        if m is None:
                            raise RuntimeError(f"Missing model {fq}")
                        L_full[qi] = m["true"][z_true][qi]
                    L_obs[f0] = (demean(L_full[:n_w1]), demean(L_full[n_w1:]))

                J_all = np.empty(n_rc2)
                z_star_all = np.empty(n_rc2)
                for i, j in enumerate(idx_acc):
                    J_z = np.zeros(n_z)
                    for f0 in freqs:
                        o1, o2 = L_obs[f0]
                        L_full_z = np.zeros((n_z, len(t_full)))
                        for qi, q in enumerate(q_t):
                            fq = round(f0 * (1 + D * q), 4)
                            m = None
                            for k, v in L_cache.items():
                                if abs(k - fq) < 0.001:
                                    m = v
                                    break
                            if m is None:
                                raise RuntimeError(f"Missing model {fq}")
                            L_full_z[:, qi] = m["cand"][j][:, qi]
                        for iz in range(n_z):
                            e1 = demean(L_full_z[iz, :n_w1]) - o1
                            e2 = demean(L_full_z[iz, n_w1:]) - o2
                            J_z[iz] += np.sum(e1 ** 2) + np.sum(e2 ** 2)
                    J_z = np.sqrt(J_z / ((n_w1 + n_w2) * len(freqs)))
                    J_all[i] = J_z.min()
                    z_star_all[i] = Z_PROFILE[J_z.argmin()]

                jmin = float(J_all.min())
                pos_true = int(np.where(idx_acc == truth_idx)[0][0])
                true_J = float(J_all[pos_true])
                true_rank = int((J_all < true_J).sum()) + 1
                true_kept = bool(true_J <= jmin + TAU)
                keep = J_all <= (jmin + TAU)

                truth_id_rows.append({"config": cfg_name, "D": D, "z_true_m": z_true,
                                      "true_J": true_J, "pass": true_J <= 1e-10})

                surv = idx_acc[keep]
                surv_r = r0_all[surv] / 1e3
                r_keep = sorted(set(int(round(x)) for x in surv_r)) if len(surv_r) else []
                r_best = sorted(set(int(round(x)) for x in r0_all[idx_acc][J_all <= jmin + 1e-10] / 1e3))
                wrong_mask = np.abs(r0_all[idx_acc] / 1e3 - 50.0) >= 0.1
                best_wrong_J = float(J_all[wrong_mask].min()) if wrong_mask.any() else np.nan
                strict = true_kept and true_rank == 1 and r_keep == [50]

                case_rows.append({
                    "config": cfg_name, "D": D, "D_pct": D * 100, "z_true_m": z_true,
                    "true_J": true_J, "J_min": jmin, "true_rank": true_rank,
                    "truth_retained": true_kept, "best_range_bins": str(r_best),
                    "best_range_unique_50": r_best == [50], "surviving_r_bins": str(r_keep),
                    "n_surviving_r_bins": len(r_keep),
                    "r_width_km": float(surv_r.max() - surv_r.min()) if len(surv_r) else 0.0,
                    "max_surviving_deviation_km": float(np.max(np.abs(surv_r - 50.0))) if len(surv_r) else 0.0,
                    "best_wrong_range_J": best_wrong_J,
                    "wrong_range_margin": best_wrong_J - true_J if np.isfinite(best_wrong_J) else np.nan,
                    "strict_single_bin": strict,
                })

                for i, j in enumerate(idx_acc):
                    cand_rows.append({"config": cfg_name, "D": D, "D_pct": D * 100, "ztrue": z_true,
                                      "node_id": int(j), "r0": float(r0_all[j] / 1e3),
                                      "v": float(v_all[j]), "psi": float(np.rad2deg(psi_all[j])),
                                      "J": float(J_all[i]), "z_star": float(z_star_all[i]),
                                      "keep": bool(keep[i]), "is_truth_state": bool(j == truth_idx)})

    cand_df = pd.DataFrame(cand_rows)
    cand_df.to_csv(OUT / "TRACKED_FREQ_CANDIDATE_SCORES_FIXED.csv", index=False)
    case_df = pd.DataFrame(case_rows)
    case_df.to_csv(OUT / "TRACKED_FREQ_ANCHOR_CASES_FIXED.csv", index=False)
    pd.DataFrame(truth_id_rows).to_csv(OUT / "TRACKED_TRUTH_MODEL_IDENTITY.csv", index=False)

    # D=0 identity
    ref = pd.read_csv(FIX / "SUBSET_CANDIDATE_SCORES_FIXED.csv")
    id_rows = []
    for cfg_name, sub_tag in [("TRIPLE", "201+235+283"), ("FOUR", "201+235+283+338")]:
        for z_true in Z_TRUE_LIST:
            cur = cand_df[(cand_df.config == cfg_name) & (cand_df.D == 0) & (cand_df.ztrue == z_true)]
            r = ref[(ref["subset"] == sub_tag) & (ref["z_true_m"] == z_true)]
            ids_c = set(cur["node_id"]); ids_r = set(r["node_id"])
            max_err = 0.0
            for nid in ids_c & ids_r:
                j_c = float(cur[cur.node_id == nid]["J"].values[0])
                j_r = float(r[r.node_id == nid]["J"].values[0])
                max_err = max(max_err, abs(j_c - j_r))
            id_rows.append({"config": cfg_name, "z_true_m": z_true, "ids_match": ids_c == ids_r,
                            "max_J_err": max_err, "pass": ids_c == ids_r and max_err <= 1e-10})
    id_df = pd.DataFrame(id_rows)
    id_df.to_csv(OUT / "ZERO_DRIFT_1D_FIX2_IDENTITY.csv", index=False)

    # working range + D dependence
    wr_rows = []
    for cfg in ["TRIPLE", "FOUR"]:
        for D in D_F_LIST:
            sub = case_df[(case_df.config == cfg) & (case_df.D == D)]
            wr_rows.append({"config": cfg, "D": D, "D_pct": D * 100,
                            "strict": int(sub["strict_single_bin"].sum()),
                            "truth_ok": int(sub["truth_retained"].sum()),
                            "rank1": int((sub["true_rank"] == 1).sum()),
                            "best_50": int(sub["best_range_unique_50"].sum()),
                            "max_width": float(sub["r_width_km"].max()),
                            "min_margin": float(sub["wrong_range_margin"].min())})
    wr_df = pd.DataFrame(wr_rows)
    wr_df.to_csv(OUT / "TRACKED_FREQ_WORKING_RANGE_SUMMARY_FIXED.csv", index=False)
    case_df.to_csv(OUT / "D_DEPENDENCE_AUDIT.csv", index=False)

    # TRIPLE vs FOUR
    cmp_rows = []
    for D in D_F_LIST:
        for z_true in Z_TRUE_LIST:
            t = case_df[(case_df.config == "TRIPLE") & (case_df.D == D) & (case_df.z_true_m == z_true)]
            f = case_df[(case_df.config == "FOUR") & (case_df.D == D) & (case_df.z_true_m == z_true)]
            if len(t) and len(f):
                cmp_rows.append({"D": D, "z_true_m": z_true,
                                 "triple_strict": bool(t.iloc[0]["strict_single_bin"]),
                                 "four_strict": bool(f.iloc[0]["strict_single_bin"])})
    pd.DataFrame(cmp_rows).to_csv(OUT / "TRIPLE_VS_FOUR_TRACKED_FREQ_ROBUSTNESS_FIXED.csv", index=False)

    # decision
    triple_all = bool(case_df[case_df.config == "TRIPLE"]["strict_single_bin"].all())
    four_all = bool(case_df[case_df.config == "FOUR"]["strict_single_bin"].all())
    all_rank = bool((case_df["true_rank"] == 1).all()) and bool(case_df["truth_retained"].all())
    all_50 = bool(case_df["best_range_unique_50"].all())

    if triple_all and four_all:
        decision = "TRIPLE_AND_FOUR_STRICT_RANGE_ANCHOR_SURVIVES_TRACKED_2PCT_FREQ_DRIFT"
    elif four_all and not triple_all:
        decision = "FOUR_REDUNDANCY_IMPROVES_TRACKED_FREQUENCY_DRIFT_ROBUSTNESS"
    elif all_rank and all_50:
        decision = "TRACKED_FREQUENCY_DRIFT_WEAKENS_STRICT_ANCHOR_BUT_PRESERVES_OPTIMAL_RANGE"
    else:
        decision = "TRACKED_FREQUENCY_DRIFT_DEGRADES_RANGE_RANKING_IN_TESTED_RANGE"

    why = f"triple_strict={triple_all}, four_strict={four_all}, rank1={all_rank}, best50={all_50}"
    dec = {"stage": "R3-CLOSEDLOOP-1E-B-FIX", "decision": decision, "why": why,
           "frequency_models": f"{len(needed_freqs)} unique freqs all built", "created_utc": NOW}
    (OUT / "R3_RC23_CLOSEDLOOP_1E_B_FIX_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# 1E-B-FIX

UTC: {NOW}

## 频率模型

{len(needed_freqs)} 个唯一频率全部构建成功（5B env 方法）。

## D-effect Gate

非零 D 与 D=0 RMS 差异 > 1e-6 dB：{'PASS' if d_effect_ok else 'FAIL'}

## 判定

### `{decision}`

{why}

## Working Range

{wr_df.to_string(index=False)}
"""
    (OUT / "R3_RC23_CLOSEDLOOP_1E_B_FIX_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# 1E-B-FIX\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
