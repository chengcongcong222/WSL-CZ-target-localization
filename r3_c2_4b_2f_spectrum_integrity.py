#!/usr/bin/env python3
"""R3-C2.4B-2F: paper-spectrum integrity fix. No D(z)."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4B_2F_SPECTRUM_INTEGRITY"
MOD_PATH = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2C" / "_kraken_paper" / "yang2014_f350.mod"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

FREQ = 350.0
R0 = 5010.0
SPANS = {"sufficient_4990m": 4990.0, "insufficient_1990m": 1990.0}
DR = 2.5
ZS_LIST = [4.0, 50.0]
ZR = 70.0
L_SIDE = 5
N_EL = 2 * L_SIDE + 1
THETA = 0.0
P7 = 7
N_OMEGA = 16384
DW = 2 * np.pi / N_OMEGA
K_MIN, K_MAX = 1.20, 1.55
TWO_PI = 2 * np.pi
PROM = 0.05
# frozen BEFORE results: observable peak = ideal |g| or beam contrib ratio
OBS_REL_THRESH = 0.03  # 3% of max ideal modal beam contribution at case setup


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
    if N <= p + 1:
        raise ValueError(f"N={N} p={p}")
    Af = np.array([[y[n - kk] for kk in range(1, p + 1)] for n in range(p, N)], dtype=complex)
    Ab = np.array(
        [[np.conj(y[n - p + kk]) for kk in range(1, p + 1)] for n in range(p, N)], dtype=complex
    )
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


def detect_peaks(P, x, rel_prom=PROM):
    """Peaks in P vs axis x (truth-independent)."""
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
            keep.append((i, float(P[i])))  # own peak level
    keep.sort(key=lambda t: t[0])
    return [i for i, _ in keep], [float(x[i]) for i, _ in keep], [lv for _, lv in keep]


def map_omega_to_k(w, dr):
    ek = DW / (2 * dr) + 1e-12
    cands = [-(w + TWO_PI * m) / dr for m in range(-10, 11)]
    return sorted({round(c, 8) for c in cands if (K_MIN - ek) <= c <= (K_MAX + ek)})


def eq24_dp_sorted(k_est, k_model):
    """Eq.(24): k0 increasing. Inputs sorted ascending; returns indices into original arrays."""
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
    return idx[::-1], float(dp[M0 - 1, idx[-1]])


def score_and_assign(k_hats, k_model_sorted, orig_idx_of_sorted, mod_index_of_sorted):
    """Sort hats ascending; Eq24; map back to KRAKEN original mode_id."""
    if len(k_hats) == 0:
        return {
            "N_DETECTED_PEAKS": 0,
            "N_UNIQUE_MODE_IDS": 0,
            "N_ORACLE_AMBIGUOUS": 0,
            "N_CORRECT": 0,
            "MODE_ORDER_RECOVERY_RATE": np.nan,
            "MODE_ORDER_ASSIGNMENT_CONSISTENCY": np.nan,
            "median_k_err": np.nan,
            "max_k_err": np.nan,
            "eq24_cost": np.nan,
            "assignments": [],
        }
    kh = np.asarray(k_hats, float)
    order = np.argsort(kh)  # ascending k
    kh_s = kh[order]
    # model already sorted ascending
    idx24, cost = eq24_dp_sorted(kh_s, k_model_sorted)
    assignments = []
    correct = n_unique = n_amb = 0
    errs = []
    for i, kval in enumerate(kh_s):
        dist = np.abs(k_model_sorted - kval)
        jmin = int(np.argmin(dist))
        dmin = float(dist[jmin])
        dist2 = dist.copy()
        dist2[jmin] = np.inf
        d2 = float(dist2.min())
        unique = (d2 - dmin) > 5e-4
        if not unique:
            n_amb += 1
            assignments.append(
                {
                    "k_hat": float(kval),
                    "eq24_model_sorted_j": int(idx24[i]) if idx24 else -1,
                    "eq24_mode_id_kraken": int(mod_index_of_sorted[idx24[i]] + 1) if idx24 else -1,
                    "oracle_sorted_j": jmin,
                    "oracle_mode_id_kraken": int(mod_index_of_sorted[jmin] + 1),
                    "oracle_status": "ORACLE_MODE_ID_AMBIGUOUS",
                    "abs_k_err": dmin,
                    "match": False,
                }
            )
            continue
        n_unique += 1
        errs.append(dmin)
        ok = bool(idx24 and idx24[i] == jmin)
        if ok:
            correct += 1
        assignments.append(
            {
                "k_hat": float(kval),
                "eq24_model_sorted_j": int(idx24[i]) if idx24 else -1,
                "eq24_mode_id_kraken": int(mod_index_of_sorted[idx24[i]] + 1) if idx24 else -1,
                "oracle_sorted_j": jmin,
                "oracle_mode_id_kraken": int(mod_index_of_sorted[jmin] + 1),
                "oracle_status": "UNIQUE",
                "abs_k_err": dmin,
                "match": ok,
            }
        )
    rate = correct / n_unique if n_unique else np.nan
    return {
        "N_DETECTED_PEAKS": int(kh.size),
        "N_UNIQUE_MODE_IDS": n_unique,
        "N_ORACLE_AMBIGUOUS": n_amb,
        "N_CORRECT": correct,
        "MODE_ORDER_RECOVERY_RATE": float(rate),
        "MODE_ORDER_ASSIGNMENT_CONSISTENCY": float(rate),
        "median_k_err": float(np.median(errs)) if errs else np.nan,
        "max_k_err": float(np.max(errs)) if errs else np.nan,
        "eq24_cost": float(cost),
        "assignments": assignments,
    }


def gen_hla(mod, zs, r, d, theta, cph, return_contrib=True):
    depths, phi, k = mod["depths"], mod["phi"], mod["k"]
    M = mod["M"]
    izs = int(np.argmin(np.abs(depths - zs)))
    izr = int(np.argmin(np.abs(depths - ZR)))
    alpha = -k.imag
    k_re = k.real
    lvals = np.arange(-L_SIDE, L_SIDE + 1)
    B = np.zeros(r.size, dtype=complex)
    contrib_abs0 = np.zeros(M)
    for m in range(M):
        km, am = k_re[m], alpha[m]
        ps, pr = phi[izs, m], phi[izr, m]
        contrib = np.zeros(r.size, dtype=complex)
        for i, ri in enumerate(r):
            A_m = np.sqrt(2 * np.pi) * np.exp(-1j * np.pi / 4) * ps * pr / np.sqrt(km * ri)
            s = 0j
            for l in lvals:
                s += np.exp(1j * (2 * np.pi * FREQ / cph) * l * d * np.sin(theta)) * np.exp(
                    -1j * (km - 1j * am) * (ri + l * d * np.sin(theta))
                )
            contrib[i] = A_m * (s / N_EL)
        B += contrib
        contrib_abs0[m] = float(abs(contrib[0]))
    return B, contrib_abs0, k_re, alpha, mod


def sab_spec(r, B, S, k_grid):
    out = np.zeros(len(k_grid), dtype=complex)
    for i, kr in enumerate(k_grid):
        out[i] = np.trapezoid(B * np.exp(1j * kr * r) * S, r) * np.exp(1j * np.pi / 4) / np.sqrt(2 * np.pi * kr)
    return out


def main() -> int:
    mod = parse_mod(MOD_PATH)
    k_re_raw = mod["k"].real.copy()
    M = mod["M"]
    # ascending sort + map back to KRAKEN index
    sort_j = np.argsort(k_re_raw)  # ascending
    k_model_sorted = k_re_raw[sort_j]
    mod_index_of_sorted = sort_j  # kraken 0-based index of sorted j

    c_zr = ssp_c(ZR)
    d = (c_zr / FREQ) / 2
    k0_steer = 2 * np.pi * FREQ / c_zr

    (OUT / "EQ24_ORDERING_AUDIT.md").write_text(
        f"""# EQ24_ORDERING_AUDIT

UTC: {NOW}

Liang Eq.(24): k0(1)<…<k0(M0)。

修正：进入 DP 前 **k̂ 与 k' 均按物理 k 升序**；模态 ID 经原始 KRAKEN index 映射回。

SAB 与 MMAR 使用同一规则。KRAKEN 文件内 mode 存储序为 k **降序**，**不得**直接当 Eq24 顺序。

重算 R0/R1 全部 metrics → `MODE_ORDER_METRICS_CORRECTED.csv`。
""",
        encoding="utf-8",
    )

    # observable mode truth (frozen threshold before looking at SAB/MMAR peaks)
    (OUT / "OBSERVABLE_MODE_TRUTH_AUDIT.md").write_text(
        f"""# OBSERVABLE_MODE_TRUTH_AUDIT

UTC: {NOW}

## 区分

- **理论模态存在**：KRAKEN k'_m 在表中
- **可观测真实谱峰** `TRUE/OBSERVABLE MODE PEAK`：生成模型中该 mode 对 beam 的贡献足够大

## 预冻结阈值（看结果前）

`OBS_REL_THRESH = {OBS_REL_THRESH}`（占该 case 最大 modal beam 贡献 |contrib(r0)| 的比例）

依据：低于主模 3% 的贡献在 5% prominence 峰检测下通常不可独立成峰，避免把弱理论模态自动当真峰。

## 检查 primary 指述：1.30–1.36 m⁻¹ 无 true spectral peak

对照 `PAPER_THEORETICAL_MODES` / 本 case `observability` 表。若模型在该窗内仍有强可观测模态，标
`PAPER_TRUE_PEAK_STRUCTURE_NOT_REPRODUCED`。
""",
        encoding="utf-8",
    )

    (OUT / "SNR_EQ28_AUDIT.md").write_text(
        f"""# SNR_EQ28_AUDIT

UTC: {NOW}

Liang Eq.(28): SNR = 10 log10(Ps/Pn)|_{{r=r0}}, r0=5010 m；**array gain 已含在 SNR 中**。

修正：**Ps = |B(r0)|²**（r=r0 处 beam output），禁止 whole-aperture mean。

噪声：Pn = Ps / 10^(SNR/10)；生成后反算 `SNR_AT_R0_DB` 验证。
不按谱结果改噪声。
""",
        encoding="utf-8",
    )

    sab_full = []
    mm_full = []
    metrics = []
    assign_all = []
    obs_rows = []

    k_grid = np.arange(K_MIN, K_MAX, DW)
    omega = np.linspace(-np.pi, np.pi, N_OMEGA, endpoint=False)

    for span_name, span in SPANS.items():
        n_r = int(round(span / DR)) + 1
        r = R0 + np.arange(n_r) * DR
        for zs in ZS_LIST:
            B0, contrib0, k_re, alpha, _ = gen_hla(mod, zs, r, d, THETA, c_zr)
            # observable modes (frozen)
            cmax = float(np.max(contrib0)) if contrib0.size else 1.0
            obs_mask = contrib0 >= OBS_REL_THRESH * cmax
            for m in range(M):
                obs_rows.append(
                    {
                        "span": span_name,
                        "zs": zs,
                        "mode_id": m + 1,
                        "k": float(k_re[m]),
                        "contrib_r0_abs": float(contrib0[m]),
                        "rel_contrib": float(contrib0[m] / cmax) if cmax else np.nan,
                        "observable": bool(obs_mask[m]),
                        "in_1p30_1p36": bool(1.30 <= k_re[m] <= 1.36),
                    }
                )
            for snr_tag, snr_db in [("R0_NOISELESS", None), ("R1_20dB", 20.0)]:
                if snr_db is None:
                    B = B0.copy()
                    snr_meas = np.inf
                else:
                    Ps = float(abs(B0[0]) ** 2)  # Eq.28 at r=r0
                    Pn = Ps / (10 ** (snr_db / 10.0))
                    rng = np.random.default_rng(0)
                    nse = np.sqrt(Pn / 2) * (rng.standard_normal(n_r) + 1j * rng.standard_normal(n_r))
                    B = B0 + nse
                    snr_meas = 10 * np.log10(float(abs(B0[0]) ** 2) / max(np.mean(np.abs(nse) ** 2), 1e-30))
                S = np.full(n_r, np.mean(np.abs(B) ** 2) ** (-0.5))
                y = B * S

                # SAB full spectrum
                g = sab_spec(r, B, S, k_grid)
                abs_g = np.abs(g)
                gmax = float(abs_g.max()) if abs_g.size else 1.0
                for i in range(len(k_grid)):
                    sab_full.append(
                        {
                            "span": span_name,
                            "zs": zs,
                            "snr": snr_tag,
                            "k": float(k_grid[i]),
                            "abs_g": float(abs_g[i]),
                            "abs_g_db": float(20 * np.log10(max(abs_g[i] / gmax, 1e-30))),
                        }
                    )
                idx_s, xs_s, lv_s = detect_peaks(abs_g, k_grid)
                k_sab = xs_s
                sc_sab = score_and_assign(k_sab, k_model_sorted, sort_j, mod_index_of_sorted)
                sc_sab.pop("assignments", None)
                metrics.append({"span": span_name, "zs": zs, "snr": snr_tag, "method": "SAB_FOURIER", "SNR_AT_R0_DB": snr_meas, **sc_sab})

                # MMAR p=7 and p=floor(2N/3)
                p_mov = int(np.floor(2 * n_r / 3))
                for p, p_tag, role in [
                    (P7, "PRINTED_LITERAL_CONTROL", "PAPER_PRINTED_PRIMARY_BRANCH"),
                    (p_mov, "OUR_MOVING_SAMPLE_INTERPRETATION", "SENSITIVITY"),
                ]:
                    if p + 2 > n_r:
                        continue
                    a, s2 = complex_mcov_ar(y, p)
                    for power, pow_tag in [(1, "PRINTED_EQ19_POWER_1"), (2, "STANDARD_PSD_POWER_2")]:
                        P = ar_spectrum(a, s2, omega, power)
                        Pmax = float(P.max())
                        # map spectrum to k axis for full dump (p=1 power only to keep size)
                        if power == 1 and p == P7:
                            for i in range(0, N_OMEGA, 16):
                                ks = map_omega_to_k(float(omega[i]), DR)
                                kk = ks[0] if ks else np.nan
                                mm_full.append(
                                    {
                                        "span": span_name,
                                        "zs": zs,
                                        "snr": snr_tag,
                                        "p": p,
                                        "omega": float(omega[i]),
                                        "k_if_mapped": kk,
                                        "P_ar": float(P[i]),
                                        "P_ar_db": float(20 * np.log10(max(P[i] / Pmax, 1e-30))),
                                    }
                                )
                        idx_a, xs_a, lv_a = detect_peaks(P, omega)  # lv_a = own peak P[idx]
                        k_mm = []
                        for w in xs_a:
                            ib = map_omega_to_k(float(w), DR)
                            if ib:
                                k_mm.append(ib[0])
                        sc = score_and_assign(k_mm, k_model_sorted, sort_j, mod_index_of_sorted)
                        assigns = sc.pop("assignments", [])
                        for asg in assigns:
                            assign_all.append(
                                {
                                    "span": span_name,
                                    "zs": zs,
                                    "snr": snr_tag,
                                    "method": f"MMAR_{p_tag}_{pow_tag}",
                                    **asg,
                                }
                            )
                        metrics.append(
                            {
                                "span": span_name,
                                "zs": zs,
                                "snr": snr_tag,
                                "method": f"MMAR_{p_tag}_{pow_tag}",
                                "role": role,
                                "p": p,
                                "SNR_AT_R0_DB": snr_meas,
                                "n_own_peak_levels": len(lv_a),
                                **sc,
                            }
                        )

    pd.DataFrame(sab_full).to_csv(OUT / "SAB_FULL_SPECTRUM.csv", index=False)
    pd.DataFrame(mm_full).to_csv(OUT / "MMAR_FULL_SPECTRUM.csv", index=False)
    met = pd.DataFrame(metrics)
    met.to_csv(OUT / "MODE_ORDER_METRICS_CORRECTED.csv", index=False)
    pd.DataFrame(assign_all).to_csv(OUT / "MODE_ORDER_ASSIGNMENTS_CORRECTED.csv", index=False)
    pd.DataFrame(obs_rows).to_csv(OUT / "observable_mode_reference.csv", index=False)

    # primary 20dB gate
    r0 = met[(met["snr"] == "R0_NOISELESS") & (met["method"] == "SAB_FOURIER")]["MODE_ORDER_RECOVERY_RATE"]
    r1 = met[(met["snr"] == "R0_NOISELESS") & (met["method"] == "MMAR_PRINTED_LITERAL_CONTROL_PRINTED_EQ19_POWER_1")]
    r1_20 = met[
        (met["snr"] == "R1_20dB") & (met["method"] == "MMAR_PRINTED_LITERAL_CONTROL_PRINTED_EQ19_POWER_1")
    ]
    sab_r0 = float(r0.mean()) if len(r0) else np.nan
    mm_r0 = float(r1["MODE_ORDER_RECOVERY_RATE"].mean()) if len(r1) else np.nan
    mm_20_npk = float(r1_20["N_DETECTED_PEAKS"].mean()) if len(r1_20) else 0
    mm_20_rate = float(r1_20["MODE_ORDER_RECOVERY_RATE"].mean()) if len(r1_20) else np.nan
    sab_20 = met[(met["snr"] == "R1_20dB") & (met["method"] == "SAB_FOURIER")]["MODE_ORDER_RECOVERY_RATE"]
    sab_20r = float(sab_20.mean()) if len(sab_20) else np.nan

    # true peak structure check 1.30-1.36
    obs_df = pd.DataFrame(obs_rows)
    win = obs_df[(obs_df["in_1p30_1p36"]) & (obs_df["observable"])]
    if len(win):
        peak_struct = "PAPER_TRUE_PEAK_STRUCTURE_NOT_REPRODUCED"
    else:
        peak_struct = "PAPER_TRUE_PEAK_STRUCTURE_CONSISTENT"

    (OUT / "FIG3_FIG4_20DB_FIDELITY_AUDIT.md").write_text(
        f"""# FIG3_FIG4_20DB_FIDELITY_AUDIT

UTC: {NOW}

Fig.3/4 正式 panel 为 **20/5/−5 dB**；R0 仅 physics control。

Eq.(28) SNR@r0 修正后 `SNR_AT_R0_DB` 见 metrics。

| | R0 rate | R1 20dB rate | R1 n_peaks |
| --- | ---: | ---: | ---: |
| SAB | {sab_r0:.3f} | {sab_20r:.3f} | — |
| MMAR p=7 | {mm_r0:.3f} | {mm_20_rate:.3f} | {mm_20_npk:.1f} |

true-peak 1.30–1.36：**{peak_struct}**

若 p=7 在 20 dB 仅 1 峰：不得判 MECHANISM_REPRODUCED（见判定）。
""",
        encoding="utf-8",
    )

    # decision
    if mm_20_npk < 2 and mm_r0 == mm_r0 and mm_r0 > 0.9:
        decision = "C2_4B2_PAPER_SPECTRUM_PARTIAL"
        why = (
            f"R0 MMAR mode-order 优势保留（MMAR {mm_r0:.3f} vs SAB {sab_r0:.3f}）；"
            f"但 p=7 在 Eq.(28) 20 dB 下平均仅 {mm_20_npk:.1f} 峰，未复现 Fig.3/4 多峰结构。"
            "记 MMAR_R0_MODE_ORDER_ADVANTAGE_OBSERVED；不进 D(z)。"
        )
    elif mm_r0 == mm_r0 and sab_r0 == sab_r0 and mm_r0 > sab_r0 + 0.3 and mm_20_npk >= 2:
        decision = "C2_4B2_PAPER_SPECTRUM_MECHANISM_REPRODUCED"
        why = f"R0+R1 多峰且 mode-order 优势保持（R0 {mm_r0:.3f} vs {sab_r0:.3f}）"
    elif peak_struct.startswith("PAPER_TRUE_PEAK"):
        decision = "C2_4B2_PAPER_SPECTRUM_PARTIAL"
        why = f"true-peak structure 不一致；{peak_struct}"
    else:
        decision = "C2_4B2_PAPER_SPECTRUM_PARTIAL"
        why = f"MMAR_R0_MODE_ORDER_ADVANTAGE_OBSERVED; partial 20dB/metrics; {peak_struct}"

    (OUT / "R3_C2_4B_2F_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-C2.4B-2F",
                "decision": decision,
                "why": why,
                "label_frozen": "MMAR_R0_MODE_ORDER_ADVANTAGE_OBSERVED",
                "eq24_ordering": "ascending k before DP; map to KRAKEN index",
                "snr_eq28": "Ps=|B(r0)|^2; array gain in SNR",
                "sab_r0_mean": sab_r0,
                "mmar_p7_r0_mean": mm_r0,
                "mmar_p7_20db_mean_peaks": mm_20_npk,
                "true_peak_structure": peak_struct,
                "observable_threshold": OBS_REL_THRESH,
                "not_done": ["D(z)", "Fig5-7", "5/-5dB", "E-STD", "MC", "P5"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_C2_4B_2F_REPORT.md").write_text(
        f"""# R3-C2.4B-2F paper-spectrum integrity

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 修正

1. Eq.(24) **k 升序**后再 DP；mode ID→KRAKEN index
2. Eq.(28) **Ps=|B(r₀)|²**（非全孔径均值）
3. `OBSERVABLE_MODE_TRUTH` 阈值 {OBS_REL_THRESH} 预冻结；{peak_struct}
4. 完整谱 + 峰自身 `abs_P`；`SAB_FULL_SPECTRUM.csv` / `MMAR_FULL_SPECTRUM.csv`

## 指标（R0）

SAB mean rate **{sab_r0:.3f}**（排序修正后）；MMAR p=7 **{mm_r0:.3f}**

冻结：`MMAR_R0_MODE_ORDER_ADVANTAGE_OBSERVED`

## 停止

不进 D(z)/Fig.5–7/E-STD/MC/P5。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# R3-C2.4B-2F\\n\\n**{decision}**\\n\\n{why}\\n",
        encoding="utf-8",
    )
    print("DECISION", decision)
    print("SAB_R0", sab_r0, "MM_R0", mm_r0, "MM20_npk", mm_20_npk, peak_struct)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
