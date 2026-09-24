#!/usr/bin/env python3
"""R3-C2.4B-2G: paper implementation ambiguity audit (S(r), look angle). No D(z)."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4B_2G_IMPL_AMBIGUITY"
MOD_PATH = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2C" / "_kraken_paper" / "yang2014_f350.mod"
PREV = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4B_2F_SPECTRUM_INTEGRITY"
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
THETAS = [0.0, 30.0, 60.0, 90.0]  # deg, pre-fixed materiality scan — NOT tuning


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


def detect_peaks(P, x, rel_prom=PROM):
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
        if all(abs(i - j) > 3 for j in keep):
            keep.append(i)
    keep.sort()
    return keep, [float(x[i]) for i in keep]


def map_omega_to_k(w, dr):
    ek = DW / (2 * dr) + 1e-12
    cands = [-(w + TWO_PI * m) / dr for m in range(-10, 11)]
    return sorted({round(c, 8) for c in cands if (K_MIN - ek) <= c <= (K_MAX + ek)})


def eq24_dp_sorted(k_est, k_model):
    k_est = np.asarray(k_est, float)
    k_model = np.asarray(k_model, float)
    M0, M = k_est.size, k_model.size
    if M0 == 0:
        return [], 0.0
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


def score_order(k_hats, k_sorted):
    if not k_hats:
        return {"N_DETECTED_PEAKS": 0, "MODE_ORDER_RECOVERY_RATE": np.nan, "median_k_err": np.nan}
    kh = np.sort(np.asarray(k_hats, float))
    idx24, _ = eq24_dp_sorted(kh, k_sorted)
    correct = n_u = 0
    errs = []
    for i, kv in enumerate(kh):
        dist = np.abs(k_sorted - kv)
        jmin = int(np.argmin(dist))
        d2 = dist.copy()
        d2[jmin] = np.inf
        if (float(d2.min()) - float(dist[jmin])) <= 5e-4:
            continue
        n_u += 1
        errs.append(float(dist[jmin]))
        if idx24 and idx24[i] == jmin:
            correct += 1
    return {
        "N_DETECTED_PEAKS": int(kh.size),
        "N_UNIQUE": n_u,
        "MODE_ORDER_RECOVERY_RATE": float(correct / n_u) if n_u else np.nan,
        "median_k_err": float(np.median(errs)) if errs else np.nan,
    }


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
    bf_abs = np.zeros(M)
    for m in range(M):
        km, am = k_re[m], alpha[m]
        ps, pr = phi[izs, m], phi[izr, m]
        X = -(km - 1j * am) * np.sin(th) + (2 * np.pi * FREQ / cph) * np.sin(th)
        den = np.sin(0.5 * d * X)
        num = np.sin((L_SIDE + 0.5) * d * X)
        BF = (num / den / N_EL) if abs(den) > 1e-12 else 1.0
        bf_abs[m] = abs(BF)
        for i, ri in enumerate(r):
            A_m = np.sqrt(2 * np.pi) * np.exp(-1j * np.pi / 4) * ps * pr / np.sqrt(km * ri)
            s = 0j
            for l in lvals:
                s += np.exp(1j * (2 * np.pi * FREQ / cph) * l * d * np.sin(th)) * np.exp(
                    -1j * (km - 1j * am) * (ri + l * d * np.sin(th))
                )
            B[i] += A_m * (s / N_EL)
    return B, bf_abs, k_re


def sab_spec(r, B, S, k_grid):
    out = np.zeros(len(k_grid), dtype=complex)
    for i, kr in enumerate(k_grid):
        out[i] = np.trapezoid(B * np.exp(1j * kr * r) * S, r) * np.exp(1j * np.pi / 4) / np.sqrt(2 * np.pi * kr)
    return out


def count_peaks_in_window(k_peaks, lo, hi):
    return sum(1 for k in k_peaks if lo <= k <= hi)


def main() -> int:
    # ---- A) true-peak from saved R0 full SAB spectrum ----
    prev_spec = PREV / "SAB_FULL_SPECTRUM.csv"
    peak_audit = []
    if prev_spec.exists():
        df = pd.read_csv(prev_spec)
        r0df = df[df["snr"] == "R0_NOISELESS"]
        for (span, zs), g in r0df.groupby(["span", "zs"]):
            abs_g = g.sort_values("k")["abs_g"].to_numpy()
            kk = g.sort_values("k")["k"].to_numpy()
            _, xs = detect_peaks(abs_g, kk)
            n_win = count_peaks_in_window(xs, 1.30, 1.36)
            peak_audit.append(
                {
                    "span": span,
                    "zs": zs,
                    "n_peaks_total": len(xs),
                    "n_peaks_1p30_1p36": n_win,
                    "peaks_in_window": ";".join(f"{k:.4f}" for k in xs if 1.30 <= k <= 1.36),
                }
            )
    peak_df = pd.DataFrame(peak_audit)
    peak_df.to_csv(OUT / "true_peak_window_recount.csv", index=False)
    all_zero = bool(len(peak_df) and (peak_df["n_peaks_1p30_1p36"] == 0).all())
    true_peak_label = (
        "PAPER_TRUE_PEAK_STRUCTURE_CONSISTENT_AT_SPECTRUM_LEVEL"
        if all_zero
        else "PAPER_TRUE_PEAK_STRUCTURE_NOT_REPRODUCED"
    )

    (OUT / "TRUE_PEAK_AUDIT_CORRECTION.md").write_text(
        f"""# TRUE_PEAK_AUDIT_CORRECTION

UTC: {NOW}

撤销 `PAPER_TRUE_PEAK_STRUCTURE_NOT_REPRODUCED`（旧依据=3% modal contribution，**不是**谱峰）。

R0 全谱 + 冻结 5% prominence 重算 1.30–1.36：

{peak_df.to_string(index=False) if len(peak_df) else 'n/a'}

→ **{true_peak_label}**

`OBS_REL_THRESH=3%` 降级为 **`MODE_CONTRIBUTION_DIAGNOSTIC`**，不得称 true peak。
""",
        encoding="utf-8",
    )

    mod = parse_mod(MOD_PATH)
    k_sorted = np.sort(mod["k"].real.copy())
    c_zr = ssp_c(ZR)
    d = (c_zr / FREQ) / 2
    omega = np.linspace(-np.pi, np.pi, N_OMEGA, endpoint=False)
    k_grid = np.arange(K_MIN, K_MAX, DW)

    # ---- B) S(r) audit + two controls on R0/20dB, p=7, theta=0 ----
    (OUT / "SHADING_SOURCE_FIDELITY_AUDIT.md").write_text(
        f"""# SHADING_SOURCE_FIDELITY_AUDIT

UTC: {NOW}

## 当前

`OUR_RANGE_AVERAGE_GLOBAL_MEAN`：S(r)=const=[mean_aperture |B|²]^{{-1/2}}

**不能**补偿 1/√r 扩展（Liang Eq.8 / Yang 广义 Hankel：⟨·⟩=range averaging，使 S≈√r）。

## 控制（非 Liang 实现）

`ORACLE_SQRT_R_SPREADING_CONTROL`：S(r)=√(r/r0)

只用于判断“20 dB 多峰缺失是否主因在 shading”；**不得**冒充数据驱动 Eq.(8)。

`PAPER_RANGE_SMOOTHING_WINDOW_NOT_EXPLICIT`（primary/Ref7 未给窗口；禁止找最佳窗）。
""",
        encoding="utf-8",
    )

    shade_rows = []
    for span_name, span in SPANS.items():
        n_r = int(round(span / DR)) + 1
        r = R0 + np.arange(n_r) * DR
        for zs in ZS_LIST:
            B0, _, _ = gen_hla(mod, zs, r, d, 0.0, c_zr)
            for snr_tag, snr_db in [("R0_NOISELESS", None), ("R1_20dB", 20.0)]:
                if snr_db is None:
                    B = B0.copy()
                else:
                    Ps = float(abs(B0[0]) ** 2)
                    Pn = Ps / (10 ** (snr_db / 10.0))
                    rng = np.random.default_rng(0)
                    nse = np.sqrt(Pn / 2) * (rng.standard_normal(n_r) + 1j * rng.standard_normal(n_r))
                    B = B0 + nse
                for s_tag, S in [
                    ("GLOBAL_MEAN_CONTROL", np.full(n_r, np.mean(np.abs(B) ** 2) ** (-0.5))),
                    ("ORACLE_SQRT_R_SPREADING_CONTROL", np.sqrt(r / R0)),
                ]:
                    y = B * S
                    g = sab_spec(r, B, S, k_grid)
                    idx_s, xs_s = detect_peaks(np.abs(g), k_grid)
                    sc_s = score_order(xs_s, k_sorted)
                    a, s2 = complex_mcov_ar(y, P7)
                    P = ar_spectrum(a, s2, omega, 1)
                    idx_a, xs_a = detect_peaks(P, omega)
                    k_mm = [map_omega_to_k(float(w), DR)[0] for w in xs_a if map_omega_to_k(float(w), DR)]
                    sc_m = score_order(k_mm, k_sorted)
                    shade_rows.append(
                        {
                            "span": span_name,
                            "zs": zs,
                            "snr": snr_tag,
                            "shading": s_tag,
                            "sab_n_peaks": sc_s.get("N_DETECTED_PEAKS", 0),
                            "sab_order_rate": sc_s.get("MODE_ORDER_RECOVERY_RATE", np.nan),
                            "mmar_n_peaks": sc_m.get("N_DETECTED_PEAKS", 0),
                            "mmar_order_rate": sc_m.get("MODE_ORDER_RECOVERY_RATE", np.nan),
                            "mmar_n_20db_peaks": sc_m.get("N_DETECTED_PEAKS", 0) if snr_tag != "R0_NOISELESS" else np.nan,
                        }
                    )
    shade_df = pd.DataFrame(shade_rows)
    shade_df.to_csv(OUT / "shading_control_comparison.csv", index=False)

    mm20_g = shade_df[(shade_df["snr"] == "R1_20dB") & (shade_df["shading"] == "GLOBAL_MEAN_CONTROL")]["mmar_n_peaks"].mean()
    mm20_o = shade_df[(shade_df["snr"] == "R1_20dB") & (shade_df["shading"] == "ORACLE_SQRT_R_SPREADING_CONTROL")]["mmar_n_peaks"].mean()
    shading_primary = bool(mm20_o == mm20_o and mm20_g == mm20_g and mm20_o >= mm20_g + 1.5)

    # ---- D) look angle materiality (R0 + 20dB, p=7, global mean S) ----
    th_rows = []
    bf_rows = []
    for th in THETAS:
        for span_name, span in list(SPANS.items())[:1]:  # 4990 sufficient enough for materiality
            n_r = int(round(span / DR)) + 1
            r = R0 + np.arange(n_r) * DR
            for zs in ZS_LIST:
                B0, bf_abs, k_re = gen_hla(mod, zs, r, d, th, c_zr)
                for m in range(len(bf_abs)):
                    bf_rows.append({"theta_deg": th, "zs": zs, "mode_id": m + 1, "k": float(k_re[m]), "abs_BF": float(bf_abs[m])})
                for snr_tag, snr_db in [("R0_NOISELESS", None), ("R1_20dB", 20.0)]:
                    if snr_db is None:
                        B = B0.copy()
                    else:
                        Ps = float(abs(B0[0]) ** 2)
                        Pn = Ps / (10 ** 20.0)
                        rng = np.random.default_rng(0)
                        B = B0 + np.sqrt(Pn / 2) * (
                            rng.standard_normal(n_r) + 1j * rng.standard_normal(n_r)
                        )
                    S = np.full(n_r, np.mean(np.abs(B) ** 2) ** (-0.5))
                    y = B * S
                    a, s2 = complex_mcov_ar(y, P7)
                    P = ar_spectrum(a, s2, omega, 1)
                    _, xs_a = detect_peaks(P, omega)
                    k_mm = [map_omega_to_k(float(w), DR)[0] for w in xs_a if map_omega_to_k(float(w), DR)]
                    sc = score_order(k_mm, k_sorted)
                    th_rows.append(
                        {
                            "theta_deg": th,
                            "zs": zs,
                            "snr": snr_tag,
                            "mmar_n_peaks": sc.get("N_DETECTED_PEAKS", 0),
                            "mmar_order_rate": sc.get("MODE_ORDER_RECOVERY_RATE", np.nan),
                            "bf_spread": float(np.std(bf_abs) / max(np.mean(bf_abs), 1e-30)),
                        }
                    )
    th_df = pd.DataFrame(th_rows)
    th_df.to_csv(OUT / "look_angle_peak_sensitivity.csv", index=False)
    pd.DataFrame(bf_rows).to_csv(OUT / "look_angle_BF_distribution.csv", index=False)
    # material if peak counts differ by >=2 across theta
    spread = th_df.groupby(["zs", "snr"])["mmar_n_peaks"].agg(lambda s: s.max() - s.min())
    material = bool((spread >= 2).any())
    look_label = (
        "PAPER_LOOK_ANGLE_CONFIG_AMBIGUITY_MATERIAL"
        if material
        else "PAPER_LOOK_ANGLE_AMBIGUITY_NONMATERIAL"
    )

    (OUT / "LOOK_ANGLE_MATERIALITY_AUDIT.md").write_text(
        f"""# LOOK_ANGLE_MATERIALITY_AUDIT

UTC: {NOW}

`PAPER_LOOK_ANGLE_NOT_EXPLICIT`（主文仅 θ̂=θ，无数值）。

预固定扫描 θ=0/30/60/90°（**非调参**）。比较 BF 分布与 R0/20 dB p=7 峰数。

→ **{look_label}**（峰数跨 θ 极差≥2 为 material）

不选择“最佳 θ”作为论文真值。
""",
        encoding="utf-8",
    )

    # F) rename decimated if present
    dec = PREV / "MMAR_FULL_SPECTRUM.csv"
    if dec.exists():
        (PREV / "MMAR_SPECTRUM_DECIMATED_FOR_AUDIT.md").write_text(
            f"UTC: {NOW}\n\n旧文件名 FULL 有误：保存时 range(0,N_OMEGA,16) 抽稀。\n重命名解释：`MMAR_SPECTRUM_DECIMATED_FOR_AUDIT`（数据仍在 MMAR_FULL_SPECTRUM.csv 历史路径，已标注）。\n",
            encoding="utf-8",
        )

    # ---- G) decision ----
    if shading_primary:
        decision = "C2_4B2_SHADING_IMPLEMENTATION_IS_PRIMARY_LIMIT"
        why = f"oracle √r 下 20 dB MMAR 峰数 {mm20_o:.1f} vs global-mean {mm20_g:.1f}；shading 主导"
    elif material:
        decision = "C2_4B2_LOOK_ANGLE_AMBIGUITY_MATERIAL"
        why = f"θ 扫描峰数材料性敏感；{look_label}"
    elif shading_primary and material:
        decision = "C2_4B2_BLOCKED_BY_PAPER_IMPLEMENTATION_AMBIGUITY"
        why = "shading + look angle 双歧义"
    else:
        # if 20dB still ~1 peak under oracle sqrt and angle nonmaterial → still not ready
        if mm20_o == mm20_o and mm20_o < 2:
            decision = "C2_4B2_BLOCKED_BY_PAPER_IMPLEMENTATION_AMBIGUITY"
            why = f"oracle √r 后 20 dB 仍约 {mm20_o:.1f} 峰；shading 非主因，存在其它实现/配置歧义（PAPER_RANGE_SMOOTHING_WINDOW_NOT_EXPLICIT 等）"
        else:
            decision = "C2_4B2_SPECTRUM_LAYER_READY_FOR_DZ_MECHANISM"
            why = f"20 dB 多峰恢复且 θ 非材料；可进 Hankel→D(z)"

    (OUT / "R3_C2_4B_2G_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-C2.4B-2G",
                "decision": decision,
                "why": why,
                "frozen": [
                    "C2_4B2_PAPER_SPECTRUM_PARTIAL",
                    "MMAR_R0_MODE_ORDER_ADVANTAGE_OBSERVED",
                    true_peak_label,
                    look_label,
                ],
                "true_peak_label": true_peak_label,
                "look_angle_label": look_label,
                "mmar_20db_peaks_global_mean": float(mm20_g) if mm20_g == mm20_g else None,
                "mmar_20db_peaks_oracle_sqrtr": float(mm20_o) if mm20_o == mm20_o else None,
                "shading_is_primary_limit": shading_primary,
                "not_done": ["D(z)", "Fig5-7", "5/-5dB", "E-STD", "MC", "P5"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_C2_4B_2G_REPORT.md").write_text(
        f"""# R3-C2.4B-2G implementation ambiguity

UTC: {NOW}

## 判定

### `{decision}`

{why}

## A. true-peak 口径

→ **{true_peak_label}**（谱峰层 1.30–1.36 峰数见 recount）；3% 仅 MODE_CONTRIBUTION_DIAGNOSTIC

## B/C. S(r)

global mean（常数）vs oracle √(r/r0)：20 dB MMAR 峰数 **{mm20_g:.1f} → {mm20_o:.1f}**

`PAPER_RANGE_SMOOTHING_WINDOW_NOT_EXPLICIT`

## D. θ

→ **{look_label}**（0/30/60/90° 材料性，不选优）

## 冻结

`C2_4B2_PAPER_SPECTRUM_PARTIAL`；`MMAR_R0_MODE_ORDER_ADVANTAGE_OBSERVED`

## 停止

不进 D(z)/Fig.5–7/E-STD/MC/P5。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(f"# R3-C2.4B-2G\\n\\n**{decision}**\\n\\n{why}\\n", encoding="utf-8")
    print("DECISION", decision)
    print("true_peak", true_peak_label, "look", look_label)
    print("20dB peaks global", mm20_g, "oracle_sqrtr", mm20_o)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
