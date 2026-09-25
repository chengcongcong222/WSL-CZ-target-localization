#!/usr/bin/env python3
"""R3-C2.4B-3A: Hankel amplitude bridge. No D(z), no Fig.5/6."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4B_3A_HANKEL_AMPLITUDE"
MOD_PATH = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2C" / "_kraken_paper" / "yang2014_f350.mod"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 350.0
R0 = 5010.0
DR = 2.5
ZR = 70.0
L_SIDE = 5
N_EL = 2 * L_SIDE + 1
P7 = 7
N_OMEGA = 16384
DW = 2 * np.pi / N_OMEGA
K_MIN, K_MAX = 1.20, 1.55
TWO_PI = 2 * np.pi
PROM = 0.05
SPANS = {"sufficient_4990m": 4990.0, "insufficient_1990m": 1990.0}
ZS_LIST = [4.0, 50.0]


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
    return a


def ar_spectrum(a, s2, omega):
    a = np.asarray(a, dtype=complex)
    acc = np.ones_like(omega, dtype=complex)
    for kk in range(1, a.size + 1):
        acc = acc + a[kk - 1] * np.exp(-1j * omega * kk)
    return s2 / np.maximum(np.abs(acc), 1e-30)


def detect_peaks_5pct(P, x):
    n = P.size
    prom = PROM * float(P.max())
    cand = [i for i in range(1, n - 1) if P[i] >= P[i - 1] and P[i] >= P[i + 1] and P[i] > 0]
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
    return keep, [float(x[i]) for i in keep], [float(P[i]) for i in keep]


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
    contrib_r0 = np.zeros(M)
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
        contrib_r0[m] = float(abs(A_m))  # proxy theoretical amplitude scale at r0-ish
        # store better: |A_m| at r0
        A0 = np.sqrt(2 * np.pi) * abs(ps * pr) / np.sqrt(km * R0)
        contrib_r0[m] = float(A0)
    return B, contrib_r0, k_re


def hankel_at_k(r, B, S, k_list):
    out = []
    for kr in k_list:
        g = np.trapezoid(B * np.exp(1j * kr * r) * S, r) * np.exp(1j * np.pi / 4) / np.sqrt(2 * np.pi * kr)
        out.append(complex(g))
    return out


def rank_corr(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.size < 3:
        return np.nan
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = np.sqrt((rx**2).sum() * (ry**2).sum())
    return float((rx * ry).sum() / den) if den > 0 else np.nan


def main() -> int:
    mod = parse_mod(MOD_PATH)
    k_model = mod["k"].real.copy()
    M = mod["M"]
    c_zr = ssp_c(ZR)
    d = (c_zr / FREQ) / 2
    omega = np.linspace(-np.pi, np.pi, N_OMEGA, endpoint=False)

    hankel_rows = []
    sanity_rows = []
    mode_count_rows = []

    for span_name, span in SPANS.items():
        n_r = int(round(span / DR)) + 1
        r = R0 + np.arange(n_r) * DR
        for zs in ZS_LIST:
            B, contrib_r0, k_re = gen_hla(mod, zs, r, d, 0.0, c_zr)
            for s_tag, S in [
                ("CURRENT_DATA_SHADING", np.full(n_r, np.mean(np.abs(B) ** 2) ** (-0.5))),
                ("ORACLE_SQRT_R_SPREADING_CONTROL", np.sqrt(r / R0)),
            ]:
                y = B * S
                a = complex_mcov_ar(y, P7)
                u = np.array([y[n] + np.dot(a, [y[n - kk] for kk in range(1, P7 + 1)]) for n in range(P7, n_r)])
                s2 = float(np.mean(np.abs(u) ** 2))
                P = ar_spectrum(a, s2, omega)
                _, xs, lvP = detect_peaks_5pct(P, omega)
                k_hats = []
                ar_levels = []
                for w, lv in zip(xs, lvP):
                    ib = map_omega_to_k(float(w), DR)
                    if ib:
                        k_hats.append(ib[0])
                        ar_levels.append(lv)

                # Eq24-style nearest unique theoretical mode (for bridge labels only)
                k_sorted = np.sort(k_re)
                assigned = []
                for kh in k_hats:
                    dist = np.abs(k_sorted - kh)
                    j = int(np.argmin(dist))
                    assigned.append(int(np.argsort(k_re)[j]))  # original mode index

                g_list = hankel_at_k(r, B, S, k_hats)
                for i, (kh, gm, midx) in enumerate(zip(k_hats, g_list, assigned)):
                    hankel_rows.append(
                        {
                            "span": span_name,
                            "zs": zs,
                            "shading": s_tag,
                            "k_hat": kh,
                            "ar_peak_level_NOT_amplitude": ar_levels[i] if i < len(ar_levels) else np.nan,
                            "g_re": gm.real,
                            "g_im": gm.imag,
                            "g_abs": abs(gm),
                            "g_phase": float(np.angle(gm)),
                            "matched_mode_id": midx + 1,
                            "k_model": float(k_re[midx]),
                            "abs_k_err": abs(kh - k_re[midx]),
                            "theory_contrib_r0": float(contrib_r0[midx]),
                            "note": "AR_PEAK_LEVEL_NOT_MODAL_AMPLITUDE",
                        }
                    )

                # sanity: rank correlation |g| vs theory among selected
                if len(k_hats) >= 3 and s_tag == "CURRENT_DATA_SHADING":
                    abs_g = [abs(g) for g in g_list]
                    th = [contrib_r0[m] for m in assigned]
                    rc = rank_corr(abs_g, th)
                    # dominant overlap: argmax
                    dom_g = int(np.argmax(abs_g))
                    dom_t = int(np.argmax(th))
                    overlap = int(assigned[dom_g] == assigned[dom_t])
                    # dynamic range
                    dr_g = float(max(abs_g) / max(min(abs_g), 1e-30))
                    # phase anomaly: none expected for noiseless matched case; report circular variance of g phases vs theory sign of phi_s*phi_r
                    phases = [np.angle(g) for g in g_list]
                    sanity_rows.append(
                        {
                            "span": span_name,
                            "zs": zs,
                            "n_modes_M0": len(k_hats),
                            "rank_corr_abs_g_vs_theory": rc,
                            "dominant_mode_overlap": overlap,
                            "amp_dynamic_range": dr_g,
                            "phase_note": "no systematic unwrap check; noiseless should be stable",
                        }
                    )
                mode_count_rows.append(
                    {
                        "span": span_name,
                        "zs": zs,
                        "shading": s_tag,
                        "M0_modes_into_future_DZ": len(k_hats),
                        "mode_ids": ";".join(str(m + 1) for m in assigned),
                    }
                )

    pd.DataFrame(hankel_rows).to_csv(OUT / "HANKEL_AT_AR_PEAKS.csv", index=False)
    san = pd.DataFrame(sanity_rows)
    san.to_csv(OUT / "AMPLITUDE_SANITY_CHECK.csv", index=False)
    pd.DataFrame(mode_count_rows).to_csv(OUT / "MODE_COUNT_FOR_DZ.csv", index=False)

    (OUT / "AR_VS_HANKEL_LEVEL_AUDIT.md").write_text(
        f"""# AR_VS_HANKEL_LEVEL_AUDIT

UTC: {NOW}

`AR_PEAK_LEVEL_NOT_MODAL_AMPLITUDE`

- AR 峰高只作诊断列 `ar_peak_level_NOT_amplitude`
- 模态幅度 **只** 用 Liang Eq.(7) generalized Hankel 在 k̂_m 处：
  g = e^{{jπ/4}}/√(2π k̂) ∫ B(r) e^{{j k̂ r}} S(r) dr

两条 shading：
- CURRENT_DATA_SHADING（全局 mean 常数）
- ORACLE_SQRT_R_SPREADING_CONTROL（机制上限）

本轮 **R0 noiseless only**。
""",
        encoding="utf-8",
    )

    (OUT / "DZ_INPUT_CONTRACT.md").write_text(
        f"""# DZ_INPUT_CONTRACT

UTC: {NOW}

未来 D(z) 输入（本轮只锁契约，不计算 b、D(z)）：

| 字段 | 来源 |
| --- | --- |
| matched theoretical mode IDs | Eq.(24) / nearest unique |
| k̂_m | AR p=7 printed branch |
| g(k̂_m,z_r) | **Eq.(7) Hankel**（复数） |
| φ_m(z_r), φ_m(z) | KRAKEN |
| δ | **0 = ORACLE_OFFSET_ALIGNMENT** |
| Δ regularizer | Liang ~½ max\\|φ\\| |

**禁止** AR peak height → amplitude。
**禁止**本轮计算 b=(Φ+U)^{-1}g 或 D(z)。

状态：机制验证，**非** Liang Fig.5/6 strict reproduction。
""",
        encoding="utf-8",
    )

    # decision
    if san.empty:
        decision = "HANKEL_AMPLITUDE_BRIDGE_FAILED"
        why = "no sanity rows"
    else:
        rc = float(san["rank_corr_abs_g_vs_theory"].mean())
        ov = float(san["dominant_mode_overlap"].mean())
        m0 = float(pd.DataFrame(mode_count_rows)["M0_modes_into_future_DZ"].mean())
        if rc == rc and rc >= 0.6 and ov >= 0.5:
            decision = "HANKEL_AMPLITUDE_BRIDGE_VALIDATED"
        elif rc == rc and rc >= 0.3:
            decision = "HANKEL_AMPLITUDE_BRIDGE_PARTIAL"
        else:
            decision = "HANKEL_AMPLITUDE_BRIDGE_FAILED"
        why = f"rank_corr={rc:.3f}, dominant_overlap={ov:.2f}, mean M0={m0:.1f}"

    (OUT / "R3_C2_4B_3A_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-C2.4B-3A",
                "decision": decision,
                "why": why,
                "frozen": [
                    "LIANG_PAPER_REPRO_PARTIAL_CLOSED",
                    "C2_4B2_PAPER_SPECTRUM_PARTIAL",
                    "PAPER_IMPLEMENTATION_AMBIGUITY",
                    "MMAR_R0_MODE_ORDER_ADVANTAGE_OBSERVED",
                    "PAPER_AR_ORDER_NOT_RECOVERABLE",
                    "PAPER_PEAK_PICKING_RULE_NOT_EXPLICIT",
                    "PAPER_LOOK_ANGLE_AMBIGUITY_NONMATERIAL",
                ],
                "task_nature": "mechanism_validation_not_Fig5_6_strict_repro",
                "ar": "AR_PEAK_LEVEL_NOT_MODAL_AMPLITUDE",
                "not_done": ["D(z)", "b=(Phi+U)^-1 g", "Fig5-6", "20/5/-5dB", "E-STD", "MC", "P5"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_C2_4B_3A_REPORT.md").write_text(
        f"""# R3-C2.4B-3A Hankel amplitude bridge

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 范围

- **机制验证**，非 Liang Fig.5/6 严格复现
- R0 无噪；p=7；span 4990/1990；z_s=4/50
- AR 只给 k̂；幅度 **Eq.(7) Hankel**
- 两种 shading 控制；不调 window

## 产物

`HANKEL_AT_AR_PEAKS.csv`（复数 g）· `AMPLITUDE_SANITY_CHECK.csv` · `MODE_COUNT_FOR_DZ.csv` · `DZ_INPUT_CONTRACT.md`

## 停止

不计算 b/D(z)/Fig.5–6；不进 E-STD/MC/P5。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(f"# R3-C2.4B-3A\\n\\n**{decision}**\\n\\n{why}\\n", encoding="utf-8")
    print("DECISION", decision, why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
