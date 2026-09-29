#!/usr/bin/env python3
"""RANGE-UB-1C-FIX2: Early-lag normalized autocorrelation shape (no peak selection).

O2 = EARLY_LAG_NORMALIZED_AUTOCORR_SHAPE
- Complex autocorrelation with phase
- 5ms bins, 0-5s window
- Normalize: q_b = u_b / sqrt(sum u_b^2)
- Score = L2 of shape difference
- No top-K, no Hungarian, no path metadata
"""
from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "RANGE_UB" / "CZ_TDOA_EARLY_LAG_AUTOCORR"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_bellhop"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 235.0
ZR = 200.0
H_DEPTH = 5000.0
ZBOX_M = 5001.0
RBOX_KM = 70.0
Z_S_ALL = np.arange(150.0, 250.0 + 1e-9, 5.0)
R_ALL_KM = np.arange(45.0, 60.0 + 1e-9, 0.5)
Z_TRUE_LIST = [180.0, 200.0, 220.0]
TAU_BIN_MS = 5.0
TAU_MAX_S = 5.0
N_BINS = int(TAU_MAX_S * 1000 / TAU_BIN_MS)  # 1000 bins
NBEAMS_CONV = [1601, 3201]
CONV_RS = [45.0, 50.0, 56.0, 60.0]
CONV_ZS = [180.0, 200.0, 220.0]
NBEAMS_MAIN = 1601
COS_GATE = 0.99
L2_GATE = 0.15


def write_env(path, zs, nbeams, tag):
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
    out = [f"'{tag}'", f"{FREQ:.1f}", "1", "'CVN'", f"  51 0.0 {H_DEPTH:.1f}", *ssp,
           "'R' 0.0", "1", f"{zs:.1f} /", "1", f"{ZR:.1f} /",
           f"{len(R_ALL_KM)}", f"{R_ALL_KM[0]:.1f} {R_ALL_KM[-1]:.1f} /",
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
                                    "zs": sd[isd] if isd < len(sd) else None,
                                    "amp": float(p[0]), "phase_deg": float(p[1]),
                                    "delay_s": float(p[2])})
    return out


def early_lag_shape(arrs, tau_bin_ms=TAU_BIN_MS, tau_max_s=TAU_MAX_S, n_bins=N_BINS):
    """Complex autocorrelation -> |C_b| per 5ms bin in (0, tau_max] -> normalized shape q_b."""
    if len(arrs) < 2:
        return np.zeros(n_bins), 0.0
    delays = np.array([a["delay_s"] for a in arrs])
    a_cplx = np.array([a["amp"] * np.exp(1j * np.radians(a["phase_deg"])) for a in arrs])
    # bins: b = int(tau_ms / 5), only 0 < tau <= 5000ms
    C = np.zeros(n_bins, dtype=complex)
    for i in range(len(delays)):
        for j in range(len(delays)):
            tau_ms = (delays[j] - delays[i]) * 1000.0
            if 0 < tau_ms <= tau_max_s * 1000:
                b = int(tau_ms / tau_bin_ms)
                if 0 <= b < n_bins:
                    C[b] += a_cplx[j] * np.conj(a_cplx[i])
    u = np.abs(C)
    energy = float(np.sum(u ** 2))
    if energy <= 0:
        return np.zeros(n_bins), 0.0
    q = u / np.sqrt(energy)
    return q, energy


def cosine_sim(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na <= 0 or nb <= 0:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


def l2_norm_diff(a, b):
    return float(np.sqrt(np.mean((a - b) ** 2)))


def main():
    config = {
        "stage": "RANGE-UB-1C-FIX2", "created_utc": NOW,
        "baseline_commit": "312e4906407fb51c628d539ee73038cde9a0d98b",
        "observable": "EARLY_LAG_NORMALIZED_AUTOCORR_SHAPE",
        "tau_bin_ms": TAU_BIN_MS, "tau_max_s": TAU_MAX_S, "n_bins": N_BINS,
        "nbeams_conv": NBEAMS_CONV, "nbeams_main": NBEAMS_MAIN,
        "cosine_gate": COS_GATE, "l2_gate": L2_GATE,
        "window_label": "PAPER_GROUNDED_MAIN_CZ_EARLY_LAG_WINDOW",
    }
    (OUT / "RANGE_UB_1C_FIX2_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    # O1 correction note
    (OUT / "O1_REPORT_CORRECTION.md").write_text(
        """# O1_REPORT_CORRECTION

UTC: """ + NOW + """

## 纠正

旧 1C-FIX 摘要写 "O1 convergence PASS 12/12" 有误。

实际 O1 = 10/12：
- zs=180, r=56: matched_fraction=0.754 < 0.8
- zs=200, r=50: matched_fraction=0.455 < 0.8

正确标签：`O1_STRONG_PATH_MEMBERSHIP_NOT_FULLY_CONVERGED`

## 保留现象

匹配上的强路径时延本身非常稳定（RMS 0.02–0.11 ms）。
不稳定来自个别 eigenray 跨越 -20 dB 选择边界。

不改变旧 1C-FIX 主判。
""", encoding="utf-8")

    # =========================================================
    # 1. Beam density convergence (12 points, 1601 vs 3201)
    # =========================================================
    conv_q = {}  # (nb, zs, rk) -> (q, energy)
    for nb in NBEAMS_CONV:
        for zs in CONV_ZS:
            tag = f"el_nb{nb}_zs{int(zs)}"
            write_env(WORK / f"{tag}.env", zs, nb, tag)
            subprocess.run([str(AT_BIN / "bellhop.exe"), tag], cwd=str(WORK), capture_output=True, timeout=300)
            arrs = parse_arr(WORK / f"{tag}.arr")
            by_r = {}
            for a in arrs:
                rk = round((a.get("r_km") or 0) * 2) / 2
                by_r.setdefault(rk, []).append(a)
            for rk in CONV_RS:
                ra = by_r.get(round(rk * 2) / 2, [])
                q, e = early_lag_shape(ra)
                conv_q[(nb, zs, rk)] = (q, e)

    conv_rows = []
    conv_pass = True
    for zs in CONV_ZS:
        for rk in CONV_RS:
            q1, e1 = conv_q.get((1601, zs, rk), (np.zeros(N_BINS), 0))
            q2, e2 = conv_q.get((3201, zs, rk), (np.zeros(N_BINS), 0))
            cos = cosine_sim(q1, q2)
            l2 = l2_norm_diff(q1, q2)
            # max peak lag
            lag1 = float(np.argmax(q1) * TAU_BIN_MS) if q1.max() > 0 else np.nan
            lag2 = float(np.argmax(q2) * TAU_BIN_MS) if q2.max() > 0 else np.nan
            ok = cos >= COS_GATE and l2 <= L2_GATE
            if not ok:
                conv_pass = False
            conv_rows.append({"zs": zs, "r_km": rk, "cosine": cos, "l2_err": l2,
                              "peak_lag_1601_ms": lag1, "peak_lag_3201_ms": lag2,
                              "energy_1601": e1, "energy_3201": e2, "pass": ok})
    conv_df = pd.DataFrame(conv_rows)
    conv_df.to_csv(OUT / "EARLY_LAG_CONVERGENCE.csv", index=False)

    # save sample signature
    q_s, e_s = conv_q.get((1601, 200.0, 50.0), (np.zeros(N_BINS), 0))
    sig_rows = [{"bin": b, "lag_ms": b * TAU_BIN_MS, "q_b": float(q_s[b])} for b in range(N_BINS)]
    pd.DataFrame(sig_rows).to_csv(OUT / "EARLY_LAG_SIGNATURE_SAMPLE.csv", index=False)

    if not conv_pass:
        (OUT / "RANGE_UB_1C_FIX2_DECISION.json").write_text(
            json.dumps({"stage": "RANGE-UB-1C-FIX2",
                        "decision": "SINGLE_DEPTH_AUTOCORR_OBSERVABLE_NOT_NUMERICALLY_STABLE",
                        "conv_pass": False, "created_utc": NOW}, indent=2), encoding="utf-8")
        print("NOT_NUMERICALLY_STABLE")
        return 1

    # =========================================================
    # 2. Full grid at 1601
    # =========================================================
    arrival_db = {}
    shape_db = {}  # (zs, rk) -> q
    for zs in Z_S_ALL:
        tag = f"el_zs{int(zs)}"
        write_env(WORK / f"{tag}.env", zs, NBEAMS_MAIN, tag)
        subprocess.run([str(AT_BIN / "bellhop.exe"), tag], cwd=str(WORK), capture_output=True, timeout=300)
        arrs = parse_arr(WORK / f"{tag}.arr")
        by_r = {}
        for a in arrs:
            rk = round((a.get("r_km") or 0) * 2) / 2
            by_r.setdefault(rk, []).append(a)
        for rk in R_ALL_KM:
            ra = by_r.get(round(rk * 2) / 2, [])
            arrival_db[(float(zs), float(rk))] = ra
            q, _ = early_lag_shape(ra)
            shape_db[(float(zs), float(rk))] = q

    # =========================================================
    # 3. Scoring
    # =========================================================
    score_rows = []
    for z_true in Z_TRUE_LIST:
        q_true = shape_db.get((z_true, 50.0))
        if q_true is None or q_true.max() <= 0:
            continue
        for z_s in Z_S_ALL:
            for rk in R_ALL_KM:
                q_cand = shape_db.get((float(z_s), float(rk)))
                if q_cand is None or q_cand.max() <= 0:
                    score_rows.append({"z_true_m": z_true, "r_km": rk, "z_s_m": z_s,
                                       "J_shape": np.nan, "cosine": np.nan, "status": "UNAVAILABLE"})
                    continue
                J = l2_norm_diff(q_true, q_cand)
                cos = cosine_sim(q_true, q_cand)
                score_rows.append({"z_true_m": z_true, "r_km": rk, "z_s_m": z_s,
                                   "J_shape": J, "cosine": cos, "status": "OK"})

    score_df = pd.DataFrame(score_rows)
    score_df.to_csv(OUT / "EARLY_LAG_RZ_SCORE.csv", index=False)

    # =========================================================
    # 4. Profiles + local + global gap
    # =========================================================
    profile_rows = []
    local_rows = []
    gap_rows = []
    rank_rows = []

    for z_true in Z_TRUE_LIST:
        sub = score_df[(score_df["z_true_m"] == z_true) & (score_df["status"] == "OK")]
        if len(sub) < 10:
            continue
        Js = sub["J_shape"].values
        rs = sub["r_km"].values
        zs_arr = sub["z_s_m"].values
        order = np.argsort(Js)
        tm = (np.abs(rs - 50.0) < 0.1) & (np.abs(zs_arr - z_true) < 0.1)
        if tm.any():
            ti = int(np.where(tm)[0][int(np.argmin(Js[tm]))])
            true_rank = int(np.where(order == ti)[0][0]) + 1
            true_J = float(Js[ti])
        else:
            true_rank, true_J = -1, np.nan
        rank_rows.append({"z_true_m": z_true, "true_rank": true_rank, "true_J": true_J, "n_valid": len(sub)})

        for rk in R_ALL_KM:
            m = np.abs(rs - rk) < 0.1
            if m.any():
                J_star = float(Js[m].min())
                z_star = float(zs_arr[m][int(np.argmin(Js[m]))])
            else:
                J_star = z_star = np.nan
            profile_rows.append({"z_true_m": z_true, "r_km": rk, "J_shape_star": J_star, "z_star_m": z_star})
            # local 48-52
            if 48.0 <= rk <= 52.0:
                local_rows.append({"z_true_m": z_true, "r_km": rk, "J_shape_star": J_star, "z_star_m": z_star})

        # global false-range gap
        false_gaps = []
        for rk in R_ALL_KM:
            if abs(rk - 50.0) < 0.1:
                continue
            m = np.abs(rs - rk) < 0.1
            if m.any():
                J_star = float(Js[m].min())
                z_star = float(zs_arr[m][int(np.argmin(Js[m]))])
                false_gaps.append({"r_km": rk, "J_shape_star": J_star, "z_star_m": z_star})
        false_gaps.sort(key=lambda x: x["J_shape_star"])
        if false_gaps:
            g1 = false_gaps[0]
            g2 = false_gaps[1] if len(false_gaps) > 1 else {}
            gap_rows.append({
                "z_true_m": z_true,
                "nearest_false_r_km": g1["r_km"], "global_min_false_gap": g1["J_shape_star"],
                "nearest_false_z_star": g1["z_star_m"],
                "second_false_r_km": g2.get("r_km", np.nan), "second_false_J": g2.get("J_shape_star", np.nan),
            })

    prof_df = pd.DataFrame(profile_rows)
    prof_df.to_csv(OUT / "EARLY_LAG_RANGE_PROFILE.csv", index=False)
    local_df = pd.DataFrame(local_rows)
    local_df.to_csv(OUT / "LOCAL_RANGE_PROFILE.csv", index=False)
    gap_df = pd.DataFrame(gap_rows)
    gap_df.to_csv(OUT / "GLOBAL_FALSE_RANGE_GAP.csv", index=False)
    rank_df = pd.DataFrame(rank_rows)

    # =========================================================
    # 5. Decision
    # =========================================================
    def range_unique(z_true):
        prof = [p for p in profile_rows if p["z_true_m"] == z_true and np.isfinite(p.get("J_shape_star", np.nan))]
        if not prof:
            return False
        k = int(np.argmin([p["J_shape_star"] for p in prof]))
        if abs(prof[k]["r_km"] - 50.0) > 0.1:
            return False
        for p in prof:
            if abs(p["r_km"] - 50.0) > 0.1 and p["J_shape_star"] < 1e-6:
                return False
        return True

    all_unique = all(range_unique(z) for z in Z_TRUE_LIST)
    any_unique = any(range_unique(z) for z in Z_TRUE_LIST)
    all_rank1 = bool((rank_df["true_rank"] == 1).all()) if len(rank_df) else False

    if all_unique and all_rank1:
        decision = "EARLY_LAG_AUTOCORR_RANGE_INFORMATION_RETAINED"
    elif any_unique:
        decision = "EARLY_LAG_AUTOCORR_RANGE_INFORMATION_PARTIAL"
    else:
        decision = "EARLY_LAG_AUTOCORR_RANGE_INFORMATION_LOST"

    why = f"all_unique={all_unique}, any_unique={any_unique}, all_rank1={all_rank1}, conv_pass={conv_pass}"
    dec = {
        "stage": "RANGE-UB-1C-FIX2", "decision": decision, "why": why,
        "observable": "EARLY_LAG_NORMALIZED_AUTOCORR_SHAPE",
        "window": "0 < tau <= 5s, 5ms bins",
        "conv_pass": conv_pass,
        "created_utc": NOW,
    }
    (OUT / "RANGE_UB_1C_FIX2_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# RANGE-UB-1C-FIX2 — Early-lag Normalized Autocorrelation Shape

UTC: {NOW}

基线：`312e4906407fb51c628d539ee73038cde9a0d98b`

## Observable

`EARLY_LAG_NORMALIZED_AUTOCORR_SHAPE`：0 < τ ≤ 5 s，5 ms bins，复数相干累加 → |C_b| → 归一化 q_b。

不再选峰、不再用 path metadata、不再需要 Hungarian。

## Convergence Gate (1601 vs 3201)

{conv_df.to_string(index=False)}

**{'PASS' if conv_pass else 'FAIL'}**

## 判定

### `{decision}`

{why}

## 排序诊断

{rank_df.to_string(index=False)}

## 全局假距离 gap

{gap_df.to_string(index=False) if len(gap_df) else 'N/A'}

## 局部距离 profile (48–52 km)

{local_df.to_string(index=False) if len(local_df) else 'N/A'}

## 强制停止规则

本轮后不再调 top-N、不再调 threshold、不再加 beams。
若失败则关闭单深度 TDOA 支线，转 VERTICAL_APERTURE_TDOA_REFERENCE_UPPER_BOUND。
"""
    (OUT / "RANGE_UB_1C_FIX2_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# RANGE-UB-1C-FIX2\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("conv_pass", conv_pass)
    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
