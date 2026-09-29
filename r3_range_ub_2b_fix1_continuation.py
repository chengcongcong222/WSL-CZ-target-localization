#!/usr/bin/env python3
"""RANGE-UB-2B-FIX1: 50km anchored slope-family continuation + range info at zs=200.

DATA_DRIVEN_LOCAL_RIDGE_CONTINUATION from 50km seed. No per-range theory bands.
"""
from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "RANGE_UB" / "E_STD_VLA_RANGE_CONTINUATION_Z200"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
FIX0B = ROOT / "results" / "RANGE_UB" / "E_STD_VLA_SLOPE_FAMILY_ID_FIX"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_bellhop"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 235.0
H_DEPTH = 5000.0
ZBOX_M = 5001.0
RBOX_KM = 70.0
ZS = 200.0
ZR_VLA = np.arange(20.0, 1620.0 + 1e-9, 10.0)
R_ALL_KM = np.arange(45.0, 60.0 + 1e-9, 0.5)
NBEAMS_LIST = [1601, 3201]
TAU_BIN_MS = 5.0
TAU_MAX_MS = 3000.0
N_TAU = int(TAU_MAX_MS / TAU_BIN_MS)
K_SCAN = np.arange(-0.8, 0.8 + 1e-9, 0.002)
B_SCAN = np.arange(0.0, 3000.0 + 1e-9, 5.0)
DK = 0.002
# continuation windows (frozen)
DK_CONT = 0.03  # ms/m
DB_CONT = 100.0  # ms
# FIX0B seeds at 50km
SEED_1601 = {"k32": -0.368, "k42": 0.138, "k43": 0.486}
SEED_3201 = {"k32": -0.366, "k42": 0.136, "k43": 0.486}


def write_env(path, nbeams, tag):
    src = ZGRID / "zgrid_f235.env"
    lines = src.read_text(encoding="utf-8").splitlines()
    ssp = []
    for ln in lines[5:]:
        p = ln.split()
        if len(p) >= 2:
            try:
                ssp.append(f"  {float(p[0]):.1f} {float(p[1]):.4f} /")
            except ValueError:
                break
        else:
            break
    n_rz = len(ZR_VLA)
    n_rr = len(R_ALL_KM)
    out = [f"'{tag}'", f"{FREQ:.1f}", "1", "'CVN'", f"  51 0.0 {H_DEPTH:.1f}", *ssp,
           "'R' 0.0", "1", f"{ZS:.1f} /",
           str(n_rz), f"{ZR_VLA[0]:.1f} {ZR_VLA[-1]:.1f} /",
           str(n_rr), f"{R_ALL_KM[0]:.1f} {R_ALL_KM[-1]:.1f} /",
           "'A'", str(nbeams), "-85.0 85.0 /",
           f"0.0 {ZBOX_M:.1f} {RBOX_KM:.1f}"]
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def parse_arr(fp):
    if not fp.exists():
        return []
    lines = [l.strip() for l in fp.read_text(encoding="utf-8", errors="replace").splitlines() if l.strip()]
    if not lines:
        return []
    h = lines[0].split()
    Nsd, Nrd, Nrr = int(float(h[1])), int(float(h[2])), int(float(h[3]))
    idx = 1
    sd = [float(x) for x in lines[idx].split()[:Nsd]]; idx += 1
    rd = [float(x) for x in lines[idx].split()[:Nrd]]; idx += 1
    rr = [float(x) for x in lines[idx].split()[:Nrr]]; idx += 1
    out = []
    for isd in range(Nsd):
        if idx >= len(lines): break
        idx += 1
        for ird in range(Nrd):
            for ir in range(Nrr):
                if idx >= len(lines): break
                narr = int(float(lines[idx])); idx += 1
                for _ in range(narr):
                    if idx >= len(lines): break
                    p = lines[idx].split(); idx += 1
                    if len(p) >= 8:
                        out.append({"r_km": rr[ir]/1000 if ir < len(rr) else None,
                                    "zr": rd[ird] if ird < len(rd) else None,
                                    "amp": float(p[0]), "delay_s": float(p[2])})
    return out


def build_I1(arrs, zr_list, n_tau=N_TAU, tau_bin=TAU_BIN_MS):
    I1 = np.zeros((len(zr_list), n_tau))
    for zi, zr in enumerate(zr_list):
        sel = [a for a in arrs if abs(a["zr"] - zr) < 0.5]
        if len(sel) < 2:
            continue
        amps = np.array([a["amp"] for a in sel])
        delays = np.array([a["delay_s"] for a in sel])
        w2 = amps ** 2
        for i in range(len(delays)):
            for j in range(i + 1, len(delays)):
                dt = abs(delays[j] - delays[i]) * 1000.0
                b = int(dt / tau_bin)
                if 0 <= b < n_tau:
                    I1[zi, b] += w2[i] * w2[j]
        mx = I1[zi].max()
        if mx > 0:
            I1[zi] /= mx
    return I1


def build_I2(I1, n_tau=N_TAU):
    nz = I1.shape[0]
    I2 = np.zeros((nz, n_tau))
    for zi in range(nz):
        row = I1[zi]
        for d in range(1, n_tau):
            I2[zi, d] = np.sum(row[:n_tau - d] * row[d:])
        mx = I2[zi].max()
        if mx > 0:
            I2[zi] /= mx
    return I2


def radon_scan(I2, zr_vals, k_scan, b_scan, tau_bin=TAU_BIN_MS):
    """Vectorized Radon scan: S(k,b) = sum_zr I2(zr, k*zr+b)."""
    nz = I2.shape[0]
    n_tau = I2.shape[1]
    zr_arr = np.asarray(zr_vals, float)
    # For each k, compute tau = k*zr+b for all b and zr
    # tau_idx[zr, b] = (k*zr + b) / tau_bin
    S = np.zeros((len(k_scan), len(b_scan)))
    for ki, k in enumerate(k_scan):
        # tau for all zr and b: shape (nz, n_b)
        tau = k * zr_arr[:, None] + b_scan[None, :]
        t_idx = tau / tau_bin
        # linear interpolation
        i0 = np.floor(t_idx).astype(int)
        frac = t_idx - i0
        valid = (i0 >= 0) & (i0 < n_tau - 1)
        i0c = np.clip(i0, 0, n_tau - 2)
        vals = (1 - frac) * I2[np.arange(nz)[:, None], i0c] + frac * I2[np.arange(nz)[:, None], i0c + 1]
        vals[~valid] = 0
        S[ki] = vals.sum(axis=0)
    return S


def find_2d_local_max_in_box(S, k_scan, b_scan, k_prev, b_prev, dk=DK_CONT, db=DB_CONT):
    """Find 2D local max within continuity box around (k_prev, b_prev)."""
    k_lo, k_hi = k_prev - dk, k_prev + dk
    b_lo, b_hi = b_prev - db, b_prev + db
    ki_mask = (k_scan >= k_lo) & (k_scan <= k_hi)
    bi_mask = (b_scan >= b_lo) & (b_scan <= b_hi)
    if not ki_mask.any() or not bi_mask.any():
        return None
    ki_idx = np.where(ki_mask)[0]
    bi_idx = np.where(bi_mask)[0]
    sub = S[np.ix_(ki_idx, bi_idx)]
    # find 2D local maxima in sub
    candidates = []
    for i in range(sub.shape[0]):
        for j in range(sub.shape[1]):
            v = sub[i, j]
            is_max = True
            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    if di == 0 and dj == 0:
                        continue
                    ni, nj = i + di, j + dj
                    if 0 <= ni < sub.shape[0] and 0 <= nj < sub.shape[1]:
                        if sub[ni, nj] > v:
                            is_max = False
                            break
                if not is_max:
                    break
            if is_max and v > 0:
                ki = ki_idx[i]
                bi = bi_idx[j]
                D = ((k_scan[ki] - k_prev) / dk) ** 2 + ((b_scan[bi] - b_prev) / db) ** 2
                candidates.append({"k": float(k_scan[ki]), "b": float(b_scan[bi]),
                                   "score": float(v), "D": float(D)})
    if not candidates:
        return None
    # sort by D, then by -score
    candidates.sort(key=lambda x: (x["D"], -x["score"]))
    return candidates[0]


def main():
    config = {"stage": "RANGE-UB-2B-FIX1", "created_utc": NOW,
              "baseline_commit": "37b0f3544a86ff9d8ad2c98dd89bef70acb61c40",
              "zs_m": ZS, "ranges_km": R_ALL_KM.tolist(),
              "continuation_window": {"dk_ms_per_m": DK_CONT, "db_ms": DB_CONT},
              "seeds_1601": SEED_1601, "seeds_3201": SEED_3201,
              "observable": "PATH_LABEL_FREE_DELAY_DEPTH_OBSERVABLE"}
    (OUT / "RANGE_UB_2B_FIX1_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "FIX0B_SCOPE_LOCK.md").write_text(
        """# FIX0B_SCOPE_LOCK

UTC: """ + NOW + """

zs=200: k32/k42/k43 全部 PASS (local-max + convergence).
zs=180: k32 NOT_LOCAL_MAX, k42 LOCAL_MAX_BUT_CONV_FAIL, k43 PASS.
zs=220: k32 NOT_LOCAL_MAX, k42 NOT_LOCAL_MAX, k43 PASS.

FIX1 admission: zs=200 only.
""", encoding="utf-8")

    # =========================================================
    # 1. BELLHOP at all ranges, both beam densities
    # =========================================================
    radon_db = {}  # (nb, r_km) -> S
    for nb in NBEAMS_LIST:
        tag = f"cont_nb{nb}"
        write_env(WORK / f"{tag}.env", nb, tag)
        subprocess.run([str(AT_BIN / "bellhop.exe"), tag], cwd=str(WORK), capture_output=True, timeout=900)
        arrs = parse_arr(WORK / f"{tag}.arr")
        by_r = {}
        for a in arrs:
            rk = round((a.get("r_km") or 0) * 2) / 2
            by_r.setdefault(rk, []).append(a)
        zr_list = sorted(set(a["zr"] for a in arrs if a["zr"] is not None))
        for rk in R_ALL_KM:
            ra = by_r.get(round(rk * 2) / 2, [])
            I1 = build_I1(ra, zr_list)
            I2 = build_I2(I1)
            zr_arr = np.array(zr_list)
            S = radon_scan(I2, zr_arr, K_SCAN, B_SCAN)
            radon_db[(nb, float(rk))] = S

    # =========================================================
    # 2. Continuation tracking
    # =========================================================
    def track_family(nb, fam, seed_k):
        """Track family from 50km seed outward."""
        S50 = radon_db.get((nb, 50.0))
        if S50 is None:
            return []
        # find seed point in 50km map (b from max score at seed_k)
        ki = int(np.argmin(np.abs(K_SCAN - seed_k)))
        bi = int(np.argmax(S50[ki]))
        k_prev, b_prev = float(K_SCAN[ki]), float(B_SCAN[bi])
        results = [{"r_km": 50.0, "family": fam, "k": k_prev, "b": b_prev,
                     "score": float(S50[ki, bi]), "prev_k": seed_k, "prev_b": b_prev,
                     "delta_k": 0.0, "delta_b": 0.0, "D": 0.0, "branch_alive": True}]
        # track toward 45
        k_p, b_p = k_prev, b_prev
        for rk in sorted([r for r in R_ALL_KM if r < 50.0], reverse=True):
            S = radon_db.get((nb, float(rk)))
            if S is None:
                results.append({"r_km": rk, "family": fam, "branch_alive": False})
                continue
            peak = find_2d_local_max_in_box(S, K_SCAN, B_SCAN, k_p, b_p)
            if peak is None:
                results.append({"r_km": rk, "family": fam, "branch_alive": False, "note": "FAMILY_BRANCH_LOST"})
                break
            results.append({"r_km": rk, "family": fam, "k": peak["k"], "b": peak["b"],
                             "score": peak["score"], "prev_k": k_p, "prev_b": b_p,
                             "delta_k": peak["k"] - k_p, "delta_b": peak["b"] - b_p,
                             "D": peak["D"], "branch_alive": True})
            k_p, b_p = peak["k"], peak["b"]
        # track toward 60
        k_p, b_p = k_prev, b_prev
        for rk in sorted([r for r in R_ALL_KM if r > 50.0]):
            S = radon_db.get((nb, float(rk)))
            if S is None:
                results.append({"r_km": rk, "family": fam, "branch_alive": False})
                continue
            peak = find_2d_local_max_in_box(S, K_SCAN, B_SCAN, k_p, b_p)
            if peak is None:
                results.append({"r_km": rk, "family": fam, "branch_alive": False, "note": "FAMILY_BRANCH_LOST"})
                break
            results.append({"r_km": rk, "family": fam, "k": peak["k"], "b": peak["b"],
                             "score": peak["score"], "prev_k": k_p, "prev_b": b_p,
                             "delta_k": peak["k"] - k_p, "delta_b": peak["b"] - b_p,
                             "D": peak["D"], "branch_alive": True})
            k_p, b_p = peak["k"], peak["b"]
        return results

    cont_1601 = []
    cont_3201 = []
    for fam in ["k32", "k42", "k43"]:
        cont_1601.extend(track_family(1601, fam, SEED_1601[fam]))
        cont_3201.extend(track_family(3201, fam, SEED_3201[fam]))
    df1 = pd.DataFrame(cont_1601)
    df2 = pd.DataFrame(cont_3201)
    df1.to_csv(OUT / "FAMILY_CONTINUATION_1601.csv", index=False)
    df2.to_csv(OUT / "FAMILY_CONTINUATION_3201.csv", index=False)

    # =========================================================
    # 3. Convergence + working range
    # =========================================================
    conv_rows = []
    for rk in R_ALL_KM:
        for fam in ["k32", "k42", "k43"]:
            r1 = df1[(df1["r_km"] == rk) & (df1["family"] == fam) & (df1.get("branch_alive", False) == True)]
            r2 = df2[(df2["r_km"] == rk) & (df2["family"] == fam) & (df2.get("branch_alive", False) == True)]
            if len(r1) and len(r2):
                k1, k2 = float(r1.iloc[0]["k"]), float(r2.iloc[0]["k"])
                b1, b2 = float(r1.iloc[0]["b"]), float(r2.iloc[0]["b"])
                sign_ok = np.sign(k1) == np.sign(k2)
                rel = abs(k1 - k2) / max(abs(k2), 0.02)
                b_diff = abs(b1 - b2)
                ok = sign_ok and rel <= 0.05 and b_diff <= 50.0
                conv_rows.append({"r_km": rk, "family": fam, "k_1601": k1, "k_3201": k2,
                                  "rel_diff": rel, "b_diff_ms": b_diff, "pass": ok})
            else:
                conv_rows.append({"r_km": rk, "family": fam, "pass": False, "note": "BRANCH_LOST"})
    conv_df = pd.DataFrame(conv_rows)
    conv_df.to_csv(OUT / "FAMILY_CONTINUATION_CONVERGENCE.csv", index=False)

    # working range
    wr_rows = []
    for rk in R_ALL_KM:
        fams = {}
        for fam in ["k32", "k42", "k43"]:
            c = conv_df[(conv_df["r_km"] == rk) & (conv_df["family"] == fam)]
            fams[fam] = bool(len(c) and c.iloc[0]["pass"])
        wr_rows.append({"r_km": rk, "k32_valid": fams["k32"], "k42_valid": fams["k42"],
                        "k43_valid": fams["k43"], "all_three_valid": fams["k32"] and fams["k42"] and fams["k43"]})
    wr_df = pd.DataFrame(wr_rows)
    wr_df.to_csv(OUT / "VLA_FAMILY_WORKING_RANGE.csv", index=False)

    # =========================================================
    # 4. Post-hoc theory audit + g-vector
    # =========================================================
    H, R0 = 5000.0, 50000.0
    C_REF = 1550.0
    post_rows = []
    g_rows = []
    for rk in R_ALL_KM:
        row = wr_df[wr_df["r_km"] == rk]
        if not len(row) or not row.iloc[0]["all_three_valid"]:
            continue
        for fam in ["k32", "k42", "k43"]:
            c = conv_df[(conv_df["r_km"] == rk) & (conv_df["family"] == fam)]
            if not len(c):
                continue
            k_val = float(c.iloc[0]["k_1601"])
            # theory value at this range
            if fam == "k32":
                k_theory = -6 * H / (C_REF * rk * 1000) * 1000
                g = -k_val / 6
            elif fam == "k42":
                k_theory = 2 * H / (C_REF * rk * 1000) * 1000
                g = k_val / 2
            else:
                k_theory = 8 * H / (C_REF * rk * 1000) * 1000
                g = k_val / 8
            post_rows.append({"r_km": rk, "family": fam, "k_tracked": k_val,
                              "k_theory": k_theory, "inside_envelope": abs(k_val - k_theory) / abs(k_theory) < 0.3})
            g_rows.append({"r_km": rk, "family": fam, "k": k_val, "g": g})
    post_df = pd.DataFrame(post_rows)
    post_df.to_csv(OUT / "POSTHOC_THEORY_CONSISTENCY.csv", index=False)
    g_df = pd.DataFrame(g_rows)
    g_df.to_csv(OUT / "RANGE_G_VECTOR.csv", index=False)

    # median g per range
    g_med_rows = []
    for rk in R_ALL_KM:
        sub = g_df[g_df["r_km"] == rk]
        if len(sub) == 3:
            gs = sub["g"].values
            g_med_rows.append({"r_km": rk, "g_median": float(np.median(gs)),
                               "g_mean": float(np.mean(gs)), "g_std": float(np.std(gs)),
                               "g_cv": float(np.std(gs) / np.mean(gs)) if np.mean(gs) > 0 else np.nan,
                               "g_max_min": float(np.max(gs) / np.min(gs)) if np.min(gs) > 0 else np.nan})
    gmed_df = pd.DataFrame(g_med_rows)
    gmed_df.to_csv(OUT / "G_VS_INVERSE_RANGE.csv", index=False)

    # =========================================================
    # 5. 1/r law + range discrimination
    # =========================================================
    if len(gmed_df) >= 3:
        r_vals = gmed_df["r_km"].values
        g_vals = gmed_df["g_median"].values
        inv_r = 1.0 / r_vals
        pearson = float(np.corrcoef(inv_r, g_vals)[0, 1])
        # simple Spearman via rank correlation
        def rankdata(a):
            order = np.argsort(a)
            ranks = np.empty_like(order, dtype=float)
            ranks[order] = np.arange(len(a), dtype=float)
            return ranks
        spearman = float(np.corrcoef(rankdata(inv_r), rankdata(g_vals))[0, 1])
        # linear fit g vs 1/r
        A = np.vstack([inv_r, np.ones_like(inv_r)]).T
        coef, res, _, _ = np.linalg.lstsq(A, g_vals, rcond=None)
        g_pred = A @ coef
        ss_res = np.sum((g_vals - g_pred) ** 2)
        ss_tot = np.sum((g_vals - np.mean(g_vals)) ** 2)
        r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else 0
        Q = r_vals * g_vals
        g_stats = {"pearson": pearson, "spearman": spearman, "r2": r2,
                   "Q_mean": float(np.mean(Q)), "Q_cv": float(np.std(Q) / np.mean(Q)) if np.mean(Q) > 0 else np.nan,
                   "Q_max_min": float(np.max(Q) / np.min(Q)) if np.min(Q) > 0 else np.nan}
    else:
        g_stats = {"pearson": np.nan, "spearman": np.nan, "r2": np.nan}

    # range score at 50km
    score_rows = []
    if len(gmed_df):
        g50 = gmed_df[gmed_df["r_km"] == 50.0]
        if len(g50):
            g_true = np.array([g_df[(g_df["r_km"] == 50) & (g_df["family"] == f)]["g"].values[0]
                               for f in ["k32", "k42", "k43"]])
            for rk in R_ALL_KM:
                sub = g_df[g_df["r_km"] == rk]
                if len(sub) == 3:
                    g_r = np.array([sub[sub["family"] == f]["g"].values[0] for f in ["k32", "k42", "k43"]])
                    J = float(np.sqrt(np.mean((g_r - g_true) ** 2)))
                    score_rows.append({"r_km": rk, "J_g": J})
            score_df = pd.DataFrame(score_rows)
            score_df.to_csv(OUT / "VLA_RANGE_SCORE_50KM.csv", index=False)
            if len(score_df):
                valid = score_df[np.isfinite(score_df["J_g"])]
                if len(valid):
                    order = valid.sort_values("J_g")
                    true_rank = int(np.where(np.abs(order["r_km"].values - 50.0) < 0.1)[0][0]) + 1 if any(np.abs(order["r_km"].values - 50.0) < 0.1) else -1
                    false = valid[np.abs(valid["r_km"] - 50.0) >= 0.1]
                    nearest_false = false.loc[false["J_g"].idxmin()] if len(false) else None
                else:
                    true_rank, nearest_false = -1, None
            else:
                true_rank, nearest_false = -1, None
        else:
            true_rank, nearest_false = -1, None
            score_df = pd.DataFrame()
    else:
        true_rank, nearest_false = -1, None
        score_df = pd.DataFrame()

    # one-point range law
    law_rows = []
    if len(gmed_df):
        g50_row = gmed_df[gmed_df["r_km"] == 50.0]
        if len(g50_row):
            C50 = 50.0 * float(g50_row.iloc[0]["g_median"])
            for _, r in gmed_df.iterrows():
                r_hat = C50 / r["g_median"] if r["g_median"] > 0 else np.nan
                law_rows.append({"r_true_km": r["r_km"], "r_hat_km": r_hat,
                                 "abs_error_km": abs(r_hat - r["r_km"]) if np.isfinite(r_hat) else np.nan,
                                 "relative_error": abs(r_hat - r["r_km"]) / r["r_km"] if np.isfinite(r_hat) else np.nan})
    pd.DataFrame(law_rows).to_csv(OUT / "ONE_POINT_RANGE_LAW_DIAGNOSTIC.csv", index=False)

    # =========================================================
    # 6. Decision
    # =========================================================
    # check 49.5 and 50.5 have all three families
    near_ok = all(
        bool(wr_df[wr_df["r_km"] == r].iloc[0]["all_three_valid"]) if len(wr_df[wr_df["r_km"] == r]) else False
        for r in [49.5, 50.5]
    ) if len(wr_df) else False
    # working range span
    valid_ranges = wr_df[wr_df["all_three_valid"]]["r_km"].values if len(wr_df) else []
    if len(valid_ranges):
        # contiguous span containing 50
        in_span = sorted([r for r in valid_ranges if 45 <= r <= 60])
        # find contiguous block containing 50
        block = [r for r in in_span if 48 <= r <= 52]
        span = max(block) - min(block) if block else 0
    else:
        span = 0
    # J_g unique
    jg_unique = True
    if len(score_df):
        false = score_df[np.abs(score_df["r_km"] - 50.0) >= 0.1]
        if len(false) and (false["J_g"] < 1e-6).any():
            jg_unique = False
    spearman_ok = g_stats.get("spearman", 0) >= 0.90 if np.isfinite(g_stats.get("spearman", np.nan)) else False

    if near_ok and span >= 5.0 and jg_unique and spearman_ok:
        decision = "E_STD_VLA_SLOPE_RANGE_INFORMATION_CONFIRMED_AT_ZS200"
    elif near_ok and span < 5.0:
        decision = "E_STD_VLA_SLOPE_RANGE_INFORMATION_LOCAL_ONLY"
    elif not near_ok:
        decision = "E_STD_VLA_SLOPE_BRANCH_CONTINUATION_PARTIAL"
    else:
        decision = "E_STD_VLA_SLOPE_RANGE_INFORMATION_NOT_ESTABLISHED"

    why = f"near_ok={near_ok}, span={span:.1f}km, jg_unique={jg_unique}, spearman={g_stats.get('spearman', np.nan):.3f}"
    dec = {"stage": "RANGE-UB-2B-FIX1", "decision": decision, "why": why,
           "g_stats": g_stats, "true_rank": true_rank, "created_utc": NOW}
    (OUT / "RANGE_UB_2B_FIX1_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# RANGE-UB-2B-FIX1 — 50km Anchored Slope-Family Continuation

UTC: {NOW}

基线：`37b0f3544a86ff9d8ad2c98dd89bef70acb61c40`

## Continuation (1601)

{df1.to_string(index=False) if len(df1) else 'N/A'}

## Working Range

{wr_df.to_string(index=False) if len(wr_df) else 'N/A'}

## g-vector Statistics

{g_stats}

## Range Score at 50 km

true_rank = {true_rank}
nearest_false = {nearest_false}

## One-point Range Law

{pd.DataFrame(law_rows).to_string(index=False) if law_rows else 'N/A'}

## 判定

### `{decision}`

{why}
"""
    (OUT / "RANGE_UB_2B_FIX1_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# RANGE-UB-2B-FIX1\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
