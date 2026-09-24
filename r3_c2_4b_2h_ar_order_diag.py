#!/usr/bin/env python3
"""R3-C2.4B-2H: AR order / 20dB spectrum diagnostic. No p sweep, no D(z)."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4B_2H_AR_ORDER_DIAGNOSTIC"
MOD_PATH = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2C" / "_kraken_paper" / "yang2014_f350.mod"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 350.0
R0 = 5010.0
DR = 2.5
ZR = 70.0
L_SIDE = 5
N_EL = 2 * L_SIDE + 1
P7 = 7  # PRINTED_LITERAL_CONTROL — not changed
N_OMEGA = 16384
DW = 2 * np.pi / N_OMEGA
K_MIN, K_MAX = 1.20, 1.55
TWO_PI = 2 * np.pi
PROM = 0.05
SPANS = {"sufficient_4990m": 4990.0, "insufficient_1990m": 1990.0}
ZS_LIST = [4.0, 50.0]
THETAS = [0.0, 30.0, 60.0, 90.0]


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


def ssp_c(z):
    if z <= 10:
        return 1533.0
    if z <= 40:
        return 1533.0 + (1478.0 - 1533.0) * (z - 10) / 30.0
    return 1478.0


def complex_mcov_ar(y, p):
    y = np.asarray(y, dtype=complex).ravel()
    N = y.size
    Af = np.array([[y[n - kk] for kk in range(1, p + 1)] for n in range(p, N)], dtype=complex)
    Ab = np.array([[np.conj(y[n - p + kk]) for kk in range(1, p + 1)] for n in range(p, N)], dtype=complex)
    bf = np.array([-y[n] for n in range(p, N)], dtype=complex)
    bb = np.array([-np.conj(y[n - p]) for n in range(p, N)], dtype=complex)
    A = np.vstack([Af, Ab])
    b = np.concatenate([bf, bb])
    a, *_ = np.linalg.lstsq(A, b, rcond=None)
    u = np.array([y[n] + np.dot(a, [y[n - kk] for kk in range(1, p + 1)]) for n in range(p, N)])
    return a, float(np.mean(np.abs(u) ** 2))


def ar_spectrum(a, s2, omega, power=1):
    a = np.asarray(a, dtype=complex)
    acc = np.ones_like(omega, dtype=complex)
    for kk in range(1, a.size + 1):
        acc = acc + a[kk - 1] * np.exp(-1j * omega * kk)
    return s2 / np.maximum(np.abs(acc) ** power, 1e-30)


def ar_poles(a):
    """Roots of A(z)=1+sum a_k z^{-k} → multiply z^p: z^p + a1 z^{p-1}+...+ap = 0."""
    a = np.asarray(a, dtype=complex)
    p = a.size
    coeffs = np.concatenate([[1.0 + 0j], a])  # z^p + a1 z^{p-1} + ...
    roots = np.roots(coeffs)
    return roots


def pole_to_k(root, dr):
    """angle ω = -angle(z) or +? use frozen map k=-(ω+2πm)/dr with ω=angle(z^{-1})= -angle(z).

    y~z^{-n}=e^{-j ω n} ⇒ z=e^{jω}, ω=angle(z). Peak of 1/|A| at z on unit circle.
    Frozen: ω_peak = wrap(-k dr) ⇒ k = -(ω+2πm)/dr, ω=angle(z).
    """
    omega = float(np.angle(root))
    cands = [-(omega + TWO_PI * m) / dr for m in range(-8, 9)]
    return omega, cands


def detect_local_maxima(P):
    return [i for i in range(1, P.size - 1) if P[i] >= P[i - 1] and P[i] >= P[i + 1] and P[i] > 0]


def detect_peaks_5pct(P, x):
    n = P.size
    prom = PROM * float(P.max())
    cand = detect_local_maxima(P)
    scored = []
    for i in cand:
        l = i
        while l > 0 and P[l - 1] <= P[i]:
            l -= 1
        r = i
        while r < n - 1 and P[r + 1] <= P[i]:
            r += 1
        pr = P[i] - max(P[l : i + 1].min(), P[i : r + 1].min())
        if pr >= prom:
            scored.append(i)
    scored.sort(key=lambda i: -P[i])
    keep = []
    for i in scored:
        if all(abs(i - j) > 3 for j in keep):
            keep.append(i)
    keep.sort()
    return keep, [float(x[i]) for i in keep]


def map_omega_to_k(w, dr):
    ek = DW / (2 * dr) + 1e-12
    cands = [-(w + TWO_PI * m) / dr for m in range(-10, 11)]
    return sorted({round(c, 8) for c in cands if (K_MIN - ek) <= c <= (K_MAX + ek)})


def gen_hla(mod, zs, r, d, theta_deg, cph):
    depths, phi, k = mod["depths"], mod["phi"], mod["k"]
    M = mod["M"]
    izs = int(np.argmin(np.abs(depths - zs)))
    izr = int(np.argmin(np.abs(depths - ZR)))
    th = np.deg2rad(theta_deg)
    alpha = -k.imag
    k_re = k.real
    lvals = np.arange(-L_SIDE, L_SIDE + 1)
    B = np.zeros(r.size, dtype=complex)
    for m in range(M):
        km, am = k_re[m], alpha[m]
        ps, pr = phi[izs, m], phi[izr, m]
        for i, ri in enumerate(r):
            A_m = np.sqrt(2 * np.pi) * np.exp(-1j * np.pi / 4) * ps * pr / np.sqrt(km * ri)
            s = 0j
            for l in lvals:
                s += np.exp(1j * (2 * np.pi * FREQ / cph) * l * d * np.sin(th)) * np.exp(
                    -1j * (km - 1j * am) * (ri + l * d * np.sin(th))
                )
            B[i] += A_m * (s / N_EL)
    return B


def add_noise_eq28(B, snr_db=20.0):
    """Eq.(28): SNR=10 log10(Ps/Pn)|_{r=r0}; array gain already in SNR.

    Pn = Ps / 10^(snr_db/10)
    """
    Ps = float(abs(B[0]) ** 2)
    Pn = Ps / (10 ** (snr_db / 10.0))
    rng = np.random.default_rng(0)
    n = np.sqrt(Pn / 2.0) * (rng.standard_normal(B.size) + 1j * rng.standard_normal(B.size))
    snr_meas = 10.0 * np.log10(Ps / max(float(np.mean(np.abs(n) ** 2)), 1e-300))
    return B + n, snr_meas


def main() -> int:
    mod = parse_mod(MOD_PATH)
    k_sorted = np.sort(mod["k"].real.copy())
    c_zr = ssp_c(ZR)
    d = (c_zr / FREQ) / 2
    omega = np.linspace(-np.pi, np.pi, N_OMEGA, endpoint=False)

    # ---- 1) look-angle 20dB fix audit ----
    th_rows = []
    for th in THETAS:
        n_r = int(round(4990.0 / DR)) + 1
        r = R0 + np.arange(n_r) * DR
        for zs in ZS_LIST:
            B0 = gen_hla(mod, zs, r, d, th, c_zr)
            B, snr_meas = add_noise_eq28(B0, 20.0)
            S = np.full(n_r, np.mean(np.abs(B) ** 2) ** (-0.5))
            y = B * S
            a, s2 = complex_mcov_ar(y, P7)
            P = ar_spectrum(a, s2, omega, 1)
            _, xs = detect_peaks_5pct(P, omega)
            k_mm = [map_omega_to_k(w, DR)[0] for w in xs if map_omega_to_k(w, DR)]
            th_rows.append(
                {
                    "theta_deg": th,
                    "zs": zs,
                    "SNR_AT_R0_DB": snr_meas,
                    "n_peaks_5pct": len(k_mm),
                    "ok_snr_20db": abs(snr_meas - 20.0) < 0.5,
                }
            )
    th_df = pd.DataFrame(th_rows)
    th_df.to_csv(OUT / "look_angle_20db_fixed.csv", index=False)
    spread = th_df.groupby("zs")["n_peaks_5pct"].agg(lambda s: int(s.max() - s.min()))
    material = bool((spread >= 2).any())
    look_label = (
        "PAPER_LOOK_ANGLE_CONFIG_AMBIGUITY_MATERIAL"
        if material
        else "PAPER_LOOK_ANGLE_AMBIGUITY_NONMATERIAL"
    )
    (OUT / "LOOK_ANGLE_20DB_FIX_AUDIT.md").write_text(
        f"""# LOOK_ANGLE_20DB_FIX_AUDIT

UTC: {NOW}

## Bug 修复

错误：`Pn = Ps / 10**20` ⇒ SNR=200 dB  
正确：`Pn = Ps / 10**(20/10)` = Ps/100 ⇒ SNR=20 dB

反算 `SNR_AT_R0_DB`：见 `look_angle_20db_fixed.csv`（均应 ≈20）。

## 冻结（仅在修正后成立）

**{look_label}**（θ=0/30/60/90°，4990 m，p=7，20 dB）

R0 下非材料性结论保留；**20 dB 版本以本表为准**。
""",
        encoding="utf-8",
    )

    # ---- 2-3) poles vs local maxima vs 5% peaks (p=7 fixed) ----
    pole_rows = []
    peak_rows = []
    for span_name, span in SPANS.items():
        n_r = int(round(span / DR)) + 1
        r = R0 + np.arange(n_r) * DR
        for zs in ZS_LIST:
            B0 = gen_hla(mod, zs, r, d, 0.0, c_zr)
            for snr_tag, use_noise in [("R0_NOISELESS", False), ("R1_20dB", True)]:
                if use_noise:
                    B, snr_meas = add_noise_eq28(B0, 20.0)
                else:
                    B, snr_meas = B0.copy(), np.inf
                S = np.full(n_r, np.mean(np.abs(B) ** 2) ** (-0.5))
                y = B * S
                a, s2 = complex_mcov_ar(y, P7)
                roots = ar_poles(a)
                for zroot in roots:
                    w, cands = pole_to_k(zroot, DR)
                    in_band = [c for c in cands if K_MIN <= c <= K_MAX]
                    kbest = min(cands, key=lambda c: min(abs(c - km) for km in k_sorted))
                    dmin = float(np.min(np.abs(k_sorted - kbest)))
                    pole_rows.append(
                        {
                            "span": span_name,
                            "zs": zs,
                            "snr": snr_tag,
                            "pole_re": float(zroot.real),
                            "pole_im": float(zroot.imag),
                            "radius": float(abs(zroot)),
                            "angle_omega": w,
                            "k_mapped_best": kbest,
                            "in_paper_band": bool(K_MIN <= kbest <= K_MAX),
                            "dist_to_nearest_kraken": dmin,
                            "n_k_cands_in_band": len(in_band),
                            "matches_R0_mode_within_5e-3": bool(dmin < 5e-3),
                            "SNR_AT_R0_DB": snr_meas,
                        }
                    )
                P = ar_spectrum(a, s2, omega, 1)
                n_loc = len(detect_local_maxima(P))
                idx5, xs5 = detect_peaks_5pct(P, omega)
                k5 = [map_omega_to_k(w, DR)[0] for w in xs5 if map_omega_to_k(w, DR)]
                # physical poles near unit circle in band
                n_pole_phys = sum(
                    1
                    for zroot in roots
                    if 0.7 < abs(zroot) < 1.3
                    and any(K_MIN <= c <= K_MAX for c in pole_to_k(zroot, DR)[1])
                )
                peak_rows.append(
                    {
                        "span": span_name,
                        "zs": zs,
                        "snr": snr_tag,
                        "p": P7,
                        "N_POLES": int(roots.size),
                        "N_POLES_NEAR_UNIT_IN_BAND": n_pole_phys,
                        "N_LOCAL_MAXIMA": n_loc,
                        "N_5PCT_PEAKS": len(idx5),
                        "N_5PCT_IN_BAND": len(k5),
                        "SNR_AT_R0_DB": snr_meas,
                    }
                )
    pd.DataFrame(pole_rows).to_csv(OUT / "P7_POLE_STRUCTURE_R0_VS_20DB.csv", index=False)
    pk_df = pd.DataFrame(peak_rows)
    pk_df.to_csv(OUT / "P7_LOCAL_MAXIMA_VS_PROMINENCE.csv", index=False)

    r1 = pk_df[pk_df["snr"] == "R1_20dB"]
    many_local_few_5 = bool(((r1["N_LOCAL_MAXIMA"] >= 3) & (r1["N_5PCT_PEAKS"] <= 2)).any())
    few_poles = bool((r1["N_POLES_NEAR_UNIT_IN_BAND"] <= 1).all())

    # ---- 4-5) AR order reference audit ----
    (OUT / "AR_ORDER_REFERENCE18_AUDIT.md").write_text(
        f"""# AR_ORDER_REFERENCE18_AUDIT

UTC: {NOW}

| 来源 | 含义 |
| --- | --- |
| Liang Eq.(18) 印刷 | p often (2/3)(2L+1)；L=5 ⇒ **p=7** `PRINTED_LITERAL_CONTROL` |
| Ref.[18] / modified covariance 体系 | **p ≤ (2/3)N** 为相对数据长度 N 的非奇异性/阶数**上限**，非“固定取 2N/3” |

因此 `OUR_MOVING_SAMPLE_INTERPRETATION` 改称：

**`MAX_ORDER_BOUND_SENSITIVITY`**（非“论文可能主实现”；禁止因结果好而采用）。

## 结论

Liang Fig.3/4 实际 simulation p **未由 primary + Ref.[18] 恢复**。

→ **`PAPER_AR_ORDER_NOT_RECOVERABLE`**（合法科学结论）

禁止 p sweep / 用 Fig.3 反推最佳 p。
""",
        encoding="utf-8",
    )
    (OUT / "PAPER_AR_ORDER_STATUS.md").write_text(
        f"""# PAPER_AR_ORDER_STATUS

UTC: {NOW}

- p=7：`PRINTED_LITERAL_CONTROL`（保留）
- p≈2N/3：**`MAX_ORDER_BOUND_SENSITIVITY`**（不再解释为论文主实现）
- Liang 实际 simulation order：**`PAPER_AR_ORDER_NOT_RECOVERABLE`**

本轮 **禁止 p sweep**。
""",
        encoding="utf-8",
    )

    # ---- 8) decision ----
    if material:
        decision = "C2_4B2_LOOK_ANGLE_AMBIGUITY_MATERIAL"
        why = f"修正 SNR 后 θ 仍材料性；{look_label}"
    elif few_poles and (r1["N_5PCT_IN_BAND"] <= 1).all():
        decision = "C2_4B2_P7_LITERAL_ORDER_INSUFFICIENT_AT_20DB"
        why = "p=7 在 20 dB 物理带内极点/峰均≈1，字面阶数不足"
    elif many_local_few_5:
        decision = "C2_4B2_BLOCKED_BY_PEAK_PICKING_AMBIGUITY"
        why = "存在多个 local maxima 但 5% 后仅 1–2 峰 → PAPER_PEAK_PICKING_RULE_NOT_EXPLICIT"
    else:
        decision = "C2_4B2_BLOCKED_BY_AR_ORDER_AMBIGUITY"
        why = "PAPER_AR_ORDER_NOT_RECOVERABLE；20 dB 多峰缺失与 p 字面/上界解释纠缠"

    (OUT / "R3_C2_4B_2H_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-C2.4B-2H",
                "decision": decision,
                "why": why,
                "look_angle_label": look_label,
                "snr_bug_fixed": True,
                "paper_ar_order": "PAPER_AR_ORDER_NOT_RECOVERABLE",
                "p7": "PRINTED_LITERAL_CONTROL",
                "p_2nr3": "MAX_ORDER_BOUND_SENSITIVITY",
                "many_local_few_5pct": many_local_few_5,
                "few_physical_poles_20db": few_poles,
                "frozen": [
                    "C2_4B2_PAPER_SPECTRUM_PARTIAL",
                    "MMAR_R0_MODE_ORDER_ADVANTAGE_OBSERVED",
                    "PAPER_TRUE_PEAK_STRUCTURE_CONSISTENT_AT_SPECTRUM_LEVEL",
                ],
                "not_done": ["D(z)", "Fig5-7", "5/-5dB", "p sweep", "E-STD", "MC", "P5"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_C2_4B_2H_REPORT.md").write_text(
        f"""# R3-C2.4B-2H AR order / 20dB diagnostic

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 1. SNR bug 已修

`Pn=Ps/10^(20/10)`；`SNR_AT_R0_DB` 见 look_angle_20db_fixed.csv。  
**{look_label}**（修正后）。

## 2–3. 极点 vs 选峰（p=7 未改）

见 `P7_POLE_STRUCTURE_R0_VS_20DB.csv`、`P7_LOCAL_MAXIMA_VS_PROMINENCE.csv`。

20 dB：N_POLES / N_LOCAL_MAXIMA / N_5PCT_PEAKS 分列。

## 4–5. AR order

- p=7 = `PRINTED_LITERAL_CONTROL`
- 2N/3 = **`MAX_ORDER_BOUND_SENSITIVITY`**（Ref.[18] 为 p≤2N/3 上限）
- **`PAPER_AR_ORDER_NOT_RECOVERABLE`**

## 停止

无 p sweep；不进 D(z)/Fig.5–7/E-STD/MC/P5。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(f"# R3-C2.4B-2H\\n\\n**{decision}**\\n\\n{why}\\n", encoding="utf-8")
    print("DECISION", decision)
    print("look", look_label, "many_local_few5", many_local_few_5, "few_poles", few_poles)
    print(pk_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
