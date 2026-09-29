#!/usr/bin/env python3
"""RANGE-UB-2B-FIX0: E-STD VLA delay-depth ridge admission.

No path labels. Build I1(zr,tau) delay-difference image, I2(zr,delta) secondary image,
Radon slope scan, ridge detection and convergence.
"""
from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "RANGE_UB" / "E_STD_VLA_DELAY_DEPTH_RIDGE"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_bellhop"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 235.0
H_DEPTH = 5000.0
ZBOX_M = 5001.0
RBOX_KM = 70.0
R_FIX = 50.0
Z_S_LIST = [180.0, 200.0, 220.0]
ZR_VLA = np.arange(20.0, 1620.0 + 1e-9, 10.0)  # 161
NBEAMS_LIST = [1601, 3201]
TAU_BIN_MS = 5.0
TAU_MAX_MS = 3000.0
N_TAU = int(TAU_MAX_MS / TAU_BIN_MS)  # 600
K_SCAN = np.arange(-0.8, 0.8 + 1e-9, 0.002)  # ms/m
B_SCAN = np.arange(0.0, 3000.0 + 1e-9, 5.0)  # ms


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
    n_rz = len(ZR_VLA)
    out = [f"'{tag}'", f"{FREQ:.1f}", "1", "'CVN'", f"  51 0.0 {H_DEPTH:.1f}", *ssp,
           "'R' 0.0", "1", f"{zs:.1f} /",
           str(n_rz), f"{ZR_VLA[0]:.1f} {ZR_VLA[-1]:.1f} /",
           "1", f"{R_FIX:.1f} /",
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
                        out.append({"zr": rd[ird] if ird < len(rd) else None,
                                    "amp": float(p[0]), "delay_s": float(p[2])})
    return out


def build_I1(arrs, zr_list, n_tau=N_TAU, tau_bin=TAU_BIN_MS):
    """I1(zr_idx, tau_bin) = sum |Ai|^2|Aj|^2 for |tj-ti| in bin."""
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
        # normalize per depth
        mx = I1[zi].max()
        if mx > 0:
            I1[zi] /= mx
    return I1


def build_I2(I1, n_tau=N_TAU):
    """I2(zr, delta) = autocorrelation of I1 along tau, excluding zero lag."""
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
    """S(k,b) = sum_zr I2(zr, k*zr+b) with linear interpolation."""
    nz = I2.shape[0]
    n_tau = I2.shape[1]
    S = np.zeros((len(k_scan), len(b_scan)))
    for ki, k in enumerate(k_scan):
        for bi, b in enumerate(b_scan):
            s = 0.0
            for zi in range(nz):
                tau = k * zr_vals[zi] + b
                t_idx = tau / tau_bin
                if 0 <= t_idx < n_tau - 1:
                    frac = t_idx - int(t_idx)
                    s += (1 - frac) * I2[zi, int(t_idx)] + frac * I2[zi, int(t_idx) + 1]
            S[ki, bi] = s
    return S


def find_ridges(S, k_scan, b_scan, n_ridges=3):
    """Find top ridge lines: non-maximum suppression on (k,b) map."""
    ridges = []
    S_work = S.copy()
    for _ in range(n_ridges * 3):
        ki, bi = np.unravel_index(np.argmax(S_work), S_work.shape)
        score = S_work[ki, bi]
        if score <= 0:
            break
        ridges.append({"k_ms_per_m": float(k_scan[ki]), "b_ms": float(b_scan[bi]), "score": float(score)})
        # suppress neighbourhood
        dk = max(1, int(0.05 / (k_scan[1] - k_scan[0])))
        db = max(1, int(50 / (b_scan[1] - b_scan[0])))
        k0, k1 = max(0, ki - dk), min(S_work.shape[0], ki + dk + 1)
        b0, b1 = max(0, bi - db), min(S_work.shape[1], bi + db + 1)
        S_work[k0:k1, b0:b1] = 0
    return ridges


def main():
    config = {
        "stage": "RANGE-UB-2B-FIX0", "created_utc": NOW,
        "baseline_commit": "48b3e2d77dc3b71e2498a60c351e21ae8a933e2b",
        "observable": "PATH_LABEL_FREE_DELAY_DEPTH_OBSERVABLE",
        "r_km": R_FIX, "zs_list": Z_S_LIST, "nbeams": NBEAMS_LIST,
        "tau_bin_ms": TAU_BIN_MS, "tau_max_ms": TAU_MAX_MS,
        "k_scan": "[-0.8, 0.8] step 0.002 ms/m",
        "note": "no reflection-order classification; delay-depth ridge only",
    }
    (OUT / "RANGE_UB_2B_FIX0_CONFIG.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    # 2B closure note
    (OUT / "2B_PATH_MAPPING_CLOSURE.md").write_text(
        """# 2B_PATH_MAPPING_CLOSURE

UTC: """ + NOW + """

旧 2B 的 reflection-order + angle-sign 分类无法映射 Xu2024 时间簇。
本轮改为 PATH_LABEL_FREE_DELAY_DEPTH_OBSERVABLE：
从 delay-depth 二维图直接找斜率脊线，不再手工标 C2/C3/C4。

2B 状态已修正为 RANGE_UB_2B_BLOCKED_BY_PATH_CLASS_MAPPING（可复现）。
""", encoding="utf-8")

    # =========================================================
    # 1. BELLHOP runs
    # =========================================================
    I2_db = {}  # (nb, zs) -> I2
    I1_db = {}
    for nb in NBEAMS_LIST:
        for zs in Z_S_LIST:
            tag = f"dd_nb{nb}_zs{int(zs)}"
            write_env(WORK / f"{tag}.env", zs, nb, tag)
            subprocess.run([str(AT_BIN / "bellhop.exe"), tag], cwd=str(WORK), capture_output=True, timeout=600)
            arrs = parse_arr(WORK / f"{tag}.arr")
            zr_list = sorted(set(a["zr"] for a in arrs if a["zr"] is not None))
            I1 = build_I1(arrs, zr_list)
            I2 = build_I2(I1)
            I1_db[(nb, zs)] = (I1, zr_list)
            I2_db[(nb, zs)] = (I2, zr_list)

    # summaries
    for nb in NBEAMS_LIST:
        for zs in Z_S_LIST:
            I1, zr_list = I1_db[(nb, zs)]
            I2, _ = I2_db[(nb, zs)]
            pd.DataFrame({
                "zr_m": zr_list,
                "I1_max": I1.max(axis=1),
                "I1_energy": (I1 ** 2).sum(axis=1),
            }).to_csv(OUT / f"DELAY_IMAGE_I1_SUMMARY_nb{nb}_zs{int(zs)}.csv", index=False)
            pd.DataFrame({
                "zr_m": zr_list,
                "I2_max": I2.max(axis=1),
                "I2_energy": (I2 ** 2).sum(axis=1),
            }).to_csv(OUT / f"SECONDARY_IMAGE_I2_SUMMARY_nb{nb}_zs{int(zs)}.csv", index=False)

    # =========================================================
    # 2. Radon scan at 1601
    # =========================================================
    ridge_db = {}  # (nb, zs) -> ridges
    for nb in NBEAMS_LIST:
        for zs in Z_S_LIST:
            I2, zr_list = I2_db[(nb, zs)]
            zr_arr = np.array(zr_list)
            S = radon_scan(I2, zr_arr, K_SCAN, B_SCAN)
            ridges = find_ridges(S, K_SCAN, B_SCAN, n_ridges=6)
            ridge_db[(nb, zs)] = ridges
            # save compact map summary
            pd.DataFrame({
                "k_ms_per_m": np.repeat(K_SCAN, len(B_SCAN)),
                "b_ms": np.tile(B_SCAN, len(K_SCAN)),
                "score": S.ravel(),
            }).to_csv(OUT / f"RADON_SLOPE_MAP_nb{nb}_zs{int(zs)}.csv", index=False)

    # detected ridges
    det_rows = []
    for nb in NBEAMS_LIST:
        for zs in Z_S_LIST:
            for i, r in enumerate(ridge_db.get((nb, zs), [])):
                det_rows.append({"nbeams": nb, "zs": zs, "rank": i + 1,
                                 "slope_ms_per_m": r["k_ms_per_m"], "intercept_ms": r["b_ms"],
                                 "line_score": r["score"], "sign": "+" if r["k_ms_per_m"] > 0 else "-"})
    det_df = pd.DataFrame(det_rows)
    det_df.to_csv(OUT / "DETECTED_SLOPE_RIDGES.csv", index=False)

    # =========================================================
    # 3. Convergence: 1601 vs 3201 on ridge slopes
    # =========================================================
    conv_rows = []
    conv_pass = True
    for zs in Z_S_LIST:
        r1 = [x for x in ridge_db.get((1601, zs), []) if abs(x["k_ms_per_m"]) > 0.01]
        r2 = [x for x in ridge_db.get((3201, zs), []) if abs(x["k_ms_per_m"]) > 0.01]
        neg1 = [x for x in r1 if x["k_ms_per_m"] < 0]
        pos1 = [x for x in r1 if x["k_ms_per_m"] > 0]
        neg2 = [x for x in r2 if x["k_ms_per_m"] < 0]
        pos2 = [x for x in r2 if x["k_ms_per_m"] > 0]
        has_structure = len(neg1) >= 1 and len(pos1) >= 2 and len(neg2) >= 1 and len(pos2) >= 2
        # match slopes
        def match_slopes(a_list, b_list):
            matched = []
            used = set()
            for a in a_list:
                best, best_j, best_d = None, None, 1e9
                for j, b in enumerate(b_list):
                    if j in used: continue
                    d = abs(a["k_ms_per_m"] - b["k_ms_per_m"])
                    if d < best_d:
                        best_d, best, best_j = d, b, j
                if best is not None:
                    used.add(best_j)
                    matched.append((a, best))
            return matched

        all_match = match_slopes(r1[:3], r2[:3])
        rel_diffs = []
        sign_ok = True
        for a, b in all_match:
            denom = max(abs(b["k_ms_per_m"]), 0.02)
            rel = abs(a["k_ms_per_m"] - b["k_ms_per_m"]) / denom
            rel_diffs.append(rel)
            if np.sign(a["k_ms_per_m"]) != np.sign(b["k_ms_per_m"]):
                sign_ok = False
        ok = has_structure and sign_ok and all(rd <= 0.05 for rd in rel_diffs) if rel_diffs else False
        if not ok: conv_pass = False
        conv_rows.append({
            "zs": zs, "has_1neg_2pos": has_structure,
            "n_neg_1601": len(neg1), "n_pos_1601": len(pos1),
            "n_neg_3201": len(neg2), "n_pos_3201": len(pos2),
            "max_rel_diff": max(rel_diffs) if rel_diffs else np.nan,
            "sign_ok": sign_ok, "pass": ok,
        })
    conv_df = pd.DataFrame(conv_rows)
    conv_df.to_csv(OUT / "RAY_DENSITY_RIDGE_CONVERGENCE.csv", index=False)

    if not conv_pass:
        (OUT / "RANGE_UB_2B_FIX0_DECISION.json").write_text(
            json.dumps({"stage": "RANGE-UB-2B-FIX0",
                        "decision": "E_STD_VLA_DELAY_DEPTH_RIDGE_NOT_ESTABLISHED",
                        "conv_pass": False, "created_utc": NOW}, indent=2), encoding="utf-8")
        print("RIDGE_NOT_ESTABLISHED")
        return 1

    # =========================================================
    # 4. Source-depth comparison + g-vector
    # =========================================================
    comp_rows = []
    for nb in NBEAMS_LIST:
        for zs in Z_S_LIST:
            ridges = [x for x in ridge_db.get((nb, zs), []) if abs(x["k_ms_per_m"]) > 0.01]
            neg = sorted([x for x in ridges if x["k_ms_per_m"] < 0], key=lambda x: x["k_ms_per_m"])
            pos = sorted([x for x in ridges if x["k_ms_per_m"] > 0], key=lambda x: x["k_ms_per_m"])
            if neg and len(pos) >= 2:
                k_neg = neg[0]["k_ms_per_m"]
                k_pos1 = pos[0]["k_ms_per_m"]
                k_pos2 = pos[1]["k_ms_per_m"]
                g32 = -k_neg / 6
                g42 = k_pos1 / 2
                g43 = k_pos2 / 8
                comp_rows.append({"nbeams": nb, "zs": zs,
                                  "k_neg": k_neg, "k_pos1": k_pos1, "k_pos2": k_pos2,
                                  "g32": g32, "g42": g42, "g43": g43,
                                  "g_spread": max(g32, g42, g43) / max(min(g32, g42, g43), 1e-9)})
    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(OUT / "SOURCE_DEPTH_SLOPE_COMPARISON.csv", index=False)

    # =========================================================
    # 5. Decision
    # =========================================================
    all_zs_pass = all(r["pass"] for r in conv_rows)
    any_zs_pass = any(r["pass"] for r in conv_rows)
    if all_zs_pass:
        decision = "E_STD_VLA_SECONDARY_DELAY_RIDGES_CONFIRMED"
    elif any_zs_pass:
        decision = "E_STD_VLA_DELAY_DEPTH_RIDGE_PARTIAL"
    else:
        decision = "E_STD_VLA_DELAY_DEPTH_RIDGE_NOT_ESTABLISHED"

    why = f"all_zs_pass={all_zs_pass}, any_zs_pass={any_zs_pass}"
    dec = {
        "stage": "RANGE-UB-2B-FIX0", "decision": decision, "why": why,
        "observable": "PATH_LABEL_FREE_DELAY_DEPTH_OBSERVABLE",
        "note": "50km mechanism admission only; not yet range discrimination",
        "created_utc": NOW,
    }
    (OUT / "RANGE_UB_2B_FIX0_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f"""# RANGE-UB-2B-FIX0 — Delay-Depth Ridge Admission

UTC: {NOW}

基线：`48b3e2d77dc3b71e2498a60c351e21ae8a933e2b`

## Observable

`PATH_LABEL_FREE_DELAY_DEPTH_OBSERVABLE`：I1(zr,τ) → I2(zr,δ) → Radon 斜率脊线。
无路径标签、无 C2/C3/C4 分类。

## Convergence (1601 vs 3201)

{conv_df.to_string(index=False)}

**{'PASS' if conv_pass else 'FAIL'}**

## 判定

### `{decision}`

{why}

## 检测脊线（1601）

{det_df[det_df['nbeams']==1601].to_string(index=False) if len(det_df) else 'N/A'}

## g-vector 诊断

{comp_df.to_string(index=False) if len(comp_df) else 'N/A'}

## 未做

45–60 km 扫描、651 候选、距离误差、噪声、HLA 融合、RC2/RC3、P5。
"""
    (OUT / "RANGE_UB_2B_FIX0_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "GPT_SYNC.md").write_text(f"# RANGE-UB-2B-FIX0\n\n**{decision}**\n\n{why}\n", encoding="utf-8")

    print("conv_pass", conv_pass)
    print("decision", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
