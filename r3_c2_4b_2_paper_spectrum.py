#!/usr/bin/env python3
"""R3-C2.4B-2: Liang paper-condition wavenumber spectrum reproduction (Fig.3/4 level).

No D(z), no Fig.5-7, no E-STD. Primary branch p=7.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4B_2_PAPER_SPECTRUM"
MOD_PATH = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2C" / "_kraken_paper" / "yang2014_f350.mod"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 350.0
R0 = 5010.0
SPANS = {"sufficient_4990m": 4990.0, "insufficient_1990m": 1990.0}
DR = 2.5  # REF7_CONSISTENT_SAMPLING_ASSUMPTION
ZS_LIST = [4.0, 50.0]
ZR = 70.0  # HLA depth
L_SIDE = 5  # 2L+1=11
N_EL = 2 * L_SIDE + 1
THETA = 0.0  # look=beam direction; ORACLE_LOOK_DIRECTION
P7 = 7
# p_moving = floor(2*N_r/3) — computed per sequence length (OUR_MOVING_SAMPLE_INTERPRETATION)
N_OMEGA = 16384
DW = 2 * np.pi / N_OMEGA
K_MIN, K_MAX = 1.20, 1.55  # slightly wider theory band for matching
TWO_PI = 2 * np.pi
PROM = 0.05


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
        return 1533.0 + (1478.0 - 1533.0) * (z - 10.0) / 30.0
    return 1478.0


def complex_mcov_ar(y, p):
    y = np.asarray(y, dtype=complex).ravel()
    if y.size <= p + 1:
        raise ValueError(f"y.size={y.size} too small for p={p}")
    N = y.size
    Af, bf, Ab, bb = [], [], [], []
    for n in range(p, N):
        Af.append(np.array([y[n - kk] for kk in range(1, p + 1)], dtype=complex))
        bf.append(complex(-y[n]))
        Ab.append(np.array([np.conj(y[n - p + kk]) for kk in range(1, p + 1)], dtype=complex))
        bb.append(complex(-np.conj(y[n - p])))
    A = np.vstack([np.asarray(Af, dtype=complex), np.asarray(Ab, dtype=complex)])
    b = np.concatenate([np.asarray(bf, dtype=complex), np.asarray(bb, dtype=complex)])
    a, *_ = np.linalg.lstsq(A, b, rcond=None)
    u = np.array([y[n] + np.dot(a, [y[n - kk] for kk in range(1, p + 1)]) for n in range(p, N)])
    return a, float(np.mean(np.abs(u) ** 2))


def ar_spectrum(a, s2, omega, power=1):
    a = np.asarray(a, dtype=complex)
    acc = np.ones_like(omega, dtype=complex)
    for kk in range(1, a.size + 1):
        acc = acc + a[kk - 1] * np.exp(-1j * omega * kk)
    return s2 / np.maximum(np.abs(acc) ** power, 1e-30)


def detect_peaks(P, omega, rel_prom=PROM):
    n = P.size
    prom = rel_prom * float(P.max())
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
        if all(abs(i - j) > 3 for j, _ in keep):
            keep.append((i, float(P[i])))
    keep.sort()
    return [i for i, _ in keep], np.array([float(omega[i]) for i, _ in keep])


def wrap(x):
    return (x + np.pi) % TWO_PI - np.pi


def map_omega_to_k(w, dr, kmin=K_MIN, kmax=K_MAX):
    ek = DW / (2 * dr) + 1e-12
    cands = [-(w + TWO_PI * m) / dr for m in range(-10, 11)]
    return sorted({round(c, 8) for c in cands if (kmin - ek) <= c <= (kmax + ek)})


def eq24_dp(k_est, k_model):
    k_est = np.asarray(k_est, float)
    k_model = np.asarray(k_model, float)
    M0, M = k_est.size, k_model.size
    if M0 == 0:
        return [], 0.0
    if M0 > M:
        return None, np.inf
    dp = np.full((M0, M), np.inf)
    prev = np.full((M0, M), -1, int)
    for j in range(M):
        dp[0, j] = (k_est[0] - k_model[j]) ** 2
    for i in range(1, M0):
        best, bj = np.inf, -1
        for j in range(M):
            if bj >= 0:
                dp[i, j] = (k_est[i] - k_model[j]) ** 2 + best
                prev[i, j] = bj
            if dp[i - 1, j] < best:
                best, bj = dp[i - 1, j], j
    j = int(np.argmin(dp[M0 - 1]))
    idx = [j]
    for i in range(M0 - 1, 0, -1):
        j = int(prev[i, j])
        idx.append(j)
    idx = idx[::-1]
    return idx, float(dp[M0 - 1, idx[-1]])


def beam_factor(L, d, X):
    """Eq.4 array factor including 1/(2L+1)."""
    den = np.sin(0.5 * d * X)
    num = np.sin((L + 0.5) * d * X)
    if abs(den) < 1e-12:
        return float(L + 0.5) / 0.5 if False else float(2 * L + 1) / (2 * L + 1)
    return float(num / den) / (2 * L + 1)


def gen_hla_modes(mod, zs, zr, r_samples, d, theta, c_ref_phase=1.0):
    """Eq.(1)-(2) multimode HLA pressure at each r_i and element l=0 (center).

    Returns B(r) via Eq.(3) with theta_hat=theta (ORACLE_LOOK_DIRECTION),
    plus per-mode contribution table (first range).
    """
    depths, phi, k = mod["depths"], mod["phi"], mod["k"]
    M = mod["M"]
    izs = int(np.argmin(np.abs(depths - zs)))
    izr = int(np.argmin(np.abs(depths - zr)))
    alpha = -k.imag
    k_re = k.real
    # HLA weights at theta_hat
    lvals = np.arange(-L_SIDE, L_SIDE + 1)
    w = np.exp(1j * (2 * np.pi * FREQ / c_ref_phase) * lvals * d * np.sin(theta))  # k=ω/c for steering

    r = np.asarray(r_samples, float)
    B = np.zeros(r.size, dtype=complex)
    mode_rows = []
    for m in range(M):
        km = k_re[m]
        am = alpha[m]
        X = -(km - 1j * am) * np.sin(theta) + (2 * np.pi * FREQ / c_ref_phase) * np.sin(theta)
        BF = beam_factor(L_SIDE, d, X)
        ps = phi[izs, m]
        pr = phi[izr, m]
        # A_m ~ sqrt(2pi) e^{-j pi/4} phi_s phi_r / sqrt(k_m r_i) — r dependence in amplitude
        # Eq(1) far field uses r_i in denom; Eq(2) freezes A_m at r_i
        # full element field summed then beamformed ≈ A_m e^{-j(km-j am)r_i} * array
        contrib = np.zeros(r.size, dtype=complex)
        for i, ri in enumerate(r):
            A_m = np.sqrt(2 * np.pi) * np.exp(-1j * np.pi / 4) * ps * pr / np.sqrt(km * ri)
            # beamforming of plane-wave-like mode across array at theta
            # sum_l e^{j k0 l d sin thetahat} e^{-j(km-j am)(ri + l d sin theta)}
            s = 0j
            for l in lvals:
                s += np.exp(1j * (2 * np.pi * FREQ / c_ref_phase) * l * d * np.sin(theta)) * np.exp(
                    -1j * (km - 1j * am) * (ri + l * d * np.sin(theta))
                )
            s = s / N_EL
            contrib[i] = A_m * s
        B += contrib
        if m < 25 or m in (0, 1, 2):
            mode_rows.append(
                {
                    "mode_id": m + 1,
                    "k_true": km,
                    "alpha": am,
                    "phi_zs": complex(ps),
                    "phi_zr": complex(pr),
                    "BF_real": float(np.real(BF)),
                    "BF_imag": float(np.imag(BF)),
                    "contrib_r0_abs": float(abs(contrib[0])),
                    "contrib_r0_phase": float(np.angle(contrib[0])),
                }
            )
    # also store all theoretical modes
    all_modes = []
    for m in range(M):
        all_modes.append(
            {
                "mode_id": m + 1,
                "k_re": float(k_re[m]),
                "k_im": float(k.imag[m]),
                "alpha": float(alpha[m]),
                "phi_zs_re": float(np.real(phi[izs, m])),
                "phi_zr_re": float(np.real(phi[izr, m])),
                "abs_phi_zs": float(abs(phi[izs, m])),
                "abs_phi_zr": float(abs(phi[izr, m])),
            }
        )
    return B, mode_rows, all_modes


def sab_spectrum(r, B, S, k_grid):
    """Generalized Hankel Eq.(7) on beam output."""
    out = np.zeros(len(k_grid), dtype=complex)
    for i, kr in enumerate(k_grid):
        integ = B * np.exp(1j * kr * r) * S
        out[i] = np.trapezoid(integ, r) * np.exp(1j * np.pi / 4) / np.sqrt(2 * np.pi * kr)
    return out


def score_mode_order(k_hats, k_model, tol_unique=5e-4):
    """MODE_ORDER_RECOVERY_RATE: Eq24 assign vs oracle unique nearest."""
    k_hats = np.asarray(k_hats, float)
    k_model = np.asarray(k_model, float)
    if k_hats.size == 0:
        return {
            "N_DETECTED_PEAKS": 0,
            "N_UNIQUE_MODE_IDS": 0,
            "N_ORACLE_AMBIGUOUS": 0,
            "N_CORRECT": 0,
            "MODE_ORDER_RECOVERY_RATE": np.nan,
            "median_k_err": np.nan,
            "max_k_err": np.nan,
            "eq24_cost": np.nan,
        }
    idx, cost = eq24_dp(k_hats, k_model)
    assigned = [int(i) for i in (idx or [])]
    correct = 0
    n_unique = 0
    n_amb = 0
    errs = []
    for i, kh in enumerate(k_hats):
        dist = np.abs(k_model - kh)
        jmin = int(np.argmin(dist))
        dmin = float(dist[jmin])
        # uniqueness: second nearest
        dist2 = dist.copy()
        dist2[jmin] = np.inf
        d2 = float(dist2.min())
        unique = (d2 - dmin) > tol_unique
        if not unique:
            n_amb += 1
            continue
        n_unique += 1
        errs.append(dmin)
        if assigned and i < len(assigned) and assigned[i] == jmin:
            correct += 1
    rate = correct / n_unique if n_unique else np.nan
    return {
        "N_DETECTED_PEAKS": int(k_hats.size),
        "N_UNIQUE_MODE_IDS": n_unique,
        "N_ORACLE_AMBIGUOUS": n_amb,
        "N_CORRECT": correct,
        "MODE_ORDER_RECOVERY_RATE": float(rate),
        "median_k_err": float(np.median(errs)) if errs else np.nan,
        "max_k_err": float(np.max(errs)) if errs else np.nan,
        "eq24_cost": float(cost),
    }


def main() -> int:
    mod = parse_mod(MOD_PATH)
    k_model = mod["k"].real
    M = mod["M"]
    c_zr = ssp_c(ZR)
    lam = c_zr / FREQ
    d = lam / 2
    k0_steer = 2 * np.pi * FREQ / c_zr  # OUR_LOCAL_SOUND_SPEED_LAMBDA_INTERPRETATION

    # provenance / audits
    (OUT / "PAPER_ENVIRONMENT_PROVENANCE.md").write_text(
        f"""# PAPER_ENVIRONMENT_PROVENANCE

UTC: {NOW}

| 参数 | 值 | 来源 |
| --- | --- | --- |
| f=350 Hz | 350 | LIANG_PRIMARY |
| HLA 2L+1=11, L=5 | 11 | LIANG_PRIMARY |
| d=λ/2 | λ=c(z_r)/f, c(70 m)={c_zr} | **OUR_LOCAL_SOUND_SPEED_LAMBDA_INTERPRETATION**（主文未给 c_ref） |
| HLA depth 70 m | 70 | LIANG_PRIMARY |
| z_s 4 / 50 m | 4,50 | LIANG_PRIMARY |
| r0=5010 m | 5010 | LIANG_PRIMARY |
| v=2.5 m/s | 2.5 | LIANG_PRIMARY |
| motion away along beam | yes | LIANG_PRIMARY |
| θ̂=θ | ORACLE_LOOK_DIRECTION | LIANG_PRIMARY |
| SSP / H=88 / bottom | Yang2014/Ref7 | YANG_REF7_REFERENCED |
| KRAKEN .mod | yang2014_f350.mod | YANG_REF7_REFERENCED (validated parser) |
| Δr=2.5 m, Δt=1 s | | **OUR REF7_CONSISTENT_SAMPLING_ASSUMPTION**（Liang 未明示） |
| R=4990 / 1990 m | | LIANG_PRIMARY（1990 来自正文 insufficient） |
""",
        encoding="utf-8",
    )
    (OUT / "RANGE_SAMPLING_AUDIT.md").write_text(
        f"""# RANGE_SAMPLING_AUDIT

UTC: {NOW}

`REF7_CONSISTENT_SAMPLING_ASSUMPTION`: Δt=1 s, Δr=2.5 m（非 Liang 明文）

| span | N_r = R/Δr+1 | r_end |
| --- | --- | --- |
| 4990 m | {int(4990/2.5)+1} | {R0+4990} |
| 1990 m | {int(1990/2.5)+1} | {R0+1990} |

含首尾点。禁止为改谱结果修改 Δr。
""",
        encoding="utf-8",
    )
    (OUT / "SHADING_IMPLEMENTATION_AUDIT.md").write_text(
        f"""# SHADING_IMPLEMENTATION_AUDIT

UTC: {NOW}

Liang Eq.(8): S(r)=⟨|B(r)|²⟩^{{-1/2}}

实现：`OUR_RANGE_AVERAGE_GLOBAL_MEAN` — 对 |B|² 在整个 range 序列上取均值后开方取负一次幂
（等价于尺度补偿的全局归一；主文未指定滑动窗）。

**不按 Fig.3/4 调 window。** 若后续需要局部平滑，单独标 interpretation。
来源：Liang Eq.(8) + Yang generalized-Hankel 实践。
""",
        encoding="utf-8",
    )

    # theoretical modes (from KRAKEN at zs=4 as representative dump — full table)
    _, _, all_modes = gen_hla_modes(mod, 4.0, ZR, np.array([R0 + 100.0]), d, THETA, c_zr)
    pd.DataFrame(all_modes).to_csv(OUT / "PAPER_THEORETICAL_MODES.csv", index=False)

    sab_rows = []
    mm_rows = []
    assign_rows = []
    metrics_rows = []
    contrib_dumped = False

    for span_name, span in SPANS.items():
        n_r = int(round(span / DR)) + 1
        r = R0 + np.arange(n_r) * DR
        for zs in ZS_LIST:
            for snr_tag, noise_amp in [("R0_NOISELESS", 0.0), ("R1_20dB", None)]:
                B, contrib, _ = gen_hla_modes(mod, zs, ZR, r, d, THETA, c_zr)
                if snr_tag.startswith("R1"):
                    # 20 dB SNR at r=r0 relative to signal power
                    sig_pow = np.mean(np.abs(B) ** 2)
                    rng = np.random.default_rng(0)
                    nse = np.sqrt(sig_pow / 100.0 / 2) * (
                        rng.standard_normal(r.size) + 1j * rng.standard_normal(r.size)
                    )
                    B = B + nse
                if not contrib_dumped:
                    pd.DataFrame(contrib).to_csv(OUT / "beam_mode_contributions_example.csv", index=False)
                    contrib_dumped = True

                # S(r) global mean shading
                S = np.full(r.size, np.mean(np.abs(B) ** 2) ** (-0.5))

                # --- SAB ---
                k_grid = np.arange(K_MIN, K_MAX, DW)
                g = sab_spectrum(r, B, S, k_grid)
                abs_g = np.abs(g)
                idx_s, ws_s = detect_peaks(abs_g, np.linspace(-np.pi, np.pi, abs_g.size, endpoint=False) * 0 + k_grid)
                # detect on |g| vs k_grid directly
                idx_s, _ = detect_peaks(abs_g, k_grid)
                k_sab = []
                for i in idx_s:
                    k_sab.append(float(k_grid[i]))
                    sab_rows.append(
                        {
                            "span": span_name,
                            "zs": zs,
                            "snr": snr_tag,
                            "k_peak": float(k_grid[i]),
                            "abs_g": float(abs_g[i]),
                        }
                    )
                sc_sab = score_mode_order(k_sab, k_model)

                # --- MMAR ---
                y = B * S
                p_moving = int(np.floor(2 * y.size / 3))
                for p, p_tag, role in [
                    (P7, "PRINTED_LITERAL_CONTROL", "PAPER_PRINTED_PRIMARY_BRANCH"),
                    (p_moving, "OUR_MOVING_SAMPLE_INTERPRETATION", "SENSITIVITY"),
                ]:
                    if p + 2 > y.size:
                        continue
                    a, s2 = complex_mcov_ar(y, p)
                    omega = np.linspace(-np.pi, np.pi, N_OMEGA, endpoint=False)
                    for power, pow_tag in [(1, "PRINTED_EQ19_POWER_1"), (2, "STANDARD_PSD_POWER_2")]:
                        P = ar_spectrum(a, s2, omega, power)
                        idx_a, ws_a = detect_peaks(P, omega)
                        k_mm = []
                        for w in ws_a:
                            ib = map_omega_to_k(float(w), DR)
                            if ib:
                                k_mm.append(ib[0])
                                mm_rows.append(
                                    {
                                        "span": span_name,
                                        "zs": zs,
                                        "snr": snr_tag,
                                        "p": p,
                                        "p_tag": p_tag,
                                        "role": role,
                                        "power": power,
                                        "power_tag": pow_tag,
                                        "omega": float(w),
                                        "k_peak": ib[0],
                                        "abs_P": float(P[int(np.argmax(P[idx_a]))]) if idx_a else np.nan,
                                    }
                                )
                        sc_mm = score_mode_order(k_mm, k_model)
                        # assignments
                        if k_mm:
                            idx24, cost24 = eq24_dp(np.array(k_mm), k_model)
                            for i, kh in enumerate(k_mm):
                                dist = np.abs(k_model - kh)
                                jmin = int(np.argmin(dist))
                                dist2 = dist.copy()
                                dist2[jmin] = np.inf
                                amb = not ((float(dist2.min()) - float(dist[jmin])) > 5e-4)
                                assign_rows.append(
                                    {
                                        "span": span_name,
                                        "zs": zs,
                                        "snr": snr_tag,
                                        "method": "MMAR",
                                        "p_tag": p_tag,
                                        "k_hat": kh,
                                        "assigned_mode_id": (idx24[i] + 1) if idx24 and i < len(idx24) else -1,
                                        "oracle_mode_id": jmin + 1,
                                        "oracle_status": "ORACLE_MODE_ID_AMBIGUOUS" if amb else "UNIQUE",
                                        "abs_k_err": float(dist[jmin]),
                                        "eq24_cost": cost24,
                                    }
                                )
                        metrics_rows.append(
                            {
                                "span": span_name,
                                "zs": zs,
                                "snr": snr_tag,
                                "method": f"MMAR_{p_tag}_{pow_tag}",
                                "role": role,
                                **sc_mm,
                            }
                        )
                metrics_rows.append(
                    {
                        "span": span_name,
                        "zs": zs,
                        "snr": snr_tag,
                        "method": "SAB_FOURIER",
                        "role": "BASELINE",
                        **sc_sab,
                    }
                )

    pd.DataFrame(sab_rows).to_csv(OUT / "SAB_WAVENUMBER_SPECTRUM.csv", index=False)
    pd.DataFrame(mm_rows).to_csv(OUT / "MMAR_WAVENUMBER_SPECTRUM.csv", index=False)
    pd.DataFrame(assign_rows).to_csv(OUT / "MODE_ORDER_ASSIGNMENTS.csv", index=False)
    met = pd.DataFrame(metrics_rows)
    met.to_csv(OUT / "SAB_VS_MMAR_METRICS.csv", index=False)

    # main comparison: p=7 power=1 vs SAB, R0
    main_mm = met[
        (met["method"] == "MMAR_PRINTED_LITERAL_CONTROL_PRINTED_EQ19_POWER_1")
        & (met["snr"] == "R0_NOISELESS")
    ]
    main_sab = met[(met["method"] == "SAB_FOURIER") & (met["snr"] == "R0_NOISELESS")]
    def avg(df, col):
        return float(df[col].mean()) if len(df) else np.nan

    sab_rate = avg(main_sab, "MODE_ORDER_RECOVERY_RATE")
    mm_rate = avg(main_mm, "MODE_ORDER_RECOVERY_RATE")
    sab_med = avg(main_sab, "median_k_err")
    mm_med = avg(main_mm, "median_k_err")

    # mechanism repro: MMAR mode-order >= SAB or peaks more complete
    if mm_rate == mm_rate and sab_rate == sab_rate and mm_rate >= sab_rate - 0.05 and mm_rate > 0.5:
        decision = "C2_4B2_PAPER_SPECTRUM_MECHANISM_REPRODUCED"
    elif mm_rate == mm_rate and mm_rate > 0.3:
        decision = "C2_4B2_PAPER_SPECTRUM_PARTIAL"
    else:
        decision = "C2_4B2_PAPER_SPECTRUM_NOT_REPRODUCED"

    (OUT / "FIG3_FIG4_FIDELITY_AUDIT.md").write_text(
        f"""# FIG3_FIG4_FIDELITY_AUDIT

UTC: {NOW}

`FIG3_FIG4_CAPTION_PANEL_LABEL_TYPO`：以正文+布局为准 SAB=(a,c,e), MMAR=(b,d,f)。

本阶段只要求**机制级**谱结构与 mode-order 一致性，不要求像素曲线重合。

- Fig.3 sufficient 4990 m SAB vs MMAR
- Fig.4 insufficient 1990 m SAB vs MMAR

R0 无噪主表见 `SAB_VS_MMAR_METRICS.csv`。

SAB MODE_ORDER_RECOVERY_RATE≈{sab_rate:.3f}；MMAR p=7≈{mm_rate:.3f}。
""",
        encoding="utf-8",
    )

    dec = {
        "stage": "R3-C2.4B-2",
        "decision": decision,
        "sab_mode_order_recovery_rate": sab_rate,
        "mmar_p7_mode_order_recovery_rate": mm_rate,
        "sab_median_k_err": sab_med,
        "mmar_p7_median_k_err": mm_med,
        "primary_branch": "P7_PRINTED_LITERAL_BRANCH / PRINTED_EQ19_MODULUS_POWER_1",
        "frozen": [
            "AR_CORE_IMPLEMENTATION_VALIDATED",
            "P1_COMPLEX_EXP_IDENTITY_VALIDATED",
            "OUR_EXACT_SOLVER_FOR_EQ24_VALIDATED",
            "BAND_EDGE_GRID_QUANTIZATION_CORRECTED",
        ],
        "not_claimed": ["D(z) reproduced", "Fig.5/6/7", "E-STD", "RC3 complete"],
        "not_done": ["D(z)", "Fig5-7", "5dB/-5dB", "E-STD", "MC", "P5"],
        "created_utc": NOW,
    }
    (OUT / "R3_C2_4B_2_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "R3_C2_4B_2_REPORT.md").write_text(
        f"""# R3-C2.4B-2 Liang paper-condition wavenumber spectrum

UTC: {NOW}

## 判定

### `{decision}`

## 指标（R0 无噪，p=7 主分支 vs SAB）

| | SAB | MMAR p=7 |
| --- | ---: | ---: |
| MODE_ORDER_RECOVERY_RATE | {sab_rate:.3f} | {mm_rate:.3f} |
| median \\|k̂−k'\\| | {sab_med:.3e} | {mm_med:.3e} |

几何：350 Hz，HLA 11×d=λ/2（c(70)={c_zr}），z_r=70，z_s=4/50，r0=5010，Δr=2.5（REF7 假设），span 4990/1990 m。

## 未做

D(z)、Fig.5–7、5/−5 dB、E-STD、MC、P5。**不宣称深度复现或 E-STD 可用。**
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# R3-C2.4B-2\\n\\n**{decision}**\\n\\nSAB rate={sab_rate:.3f} MMAR p7 rate={mm_rate:.3f}\\n",
        encoding="utf-8",
    )
    print("DECISION", decision, "SAB", sab_rate, "MMAR", mm_rate)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
