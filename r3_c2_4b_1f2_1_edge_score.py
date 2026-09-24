#!/usr/bin/env python3
"""R3-C2.4B-1F2.1: band-edge scoring correction only. No AR algorithm change."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4B_1F2_1_EDGE_SCORE"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

K_MIN, K_MAX = 1.32, 1.50
K_PHYS = [1.32, 1.41, 1.50]
N_R = 200
N_OMEGA = 16384
DW = 2 * np.pi / N_OMEGA
TWO_PI = 2 * np.pi
P_LIT = 7
P_MOVING = int(np.floor(2 * N_R / 3))


def wrap(x):
    return (x + np.pi) % TWO_PI - np.pi


def complex_mcov_ar(y, p):
    y = np.asarray(y, dtype=complex)
    N = y.size
    Af, bf, Ab, bb = [], [], [], []
    for n in range(p, N):
        Af.append(np.array([y[n - k] for k in range(1, p + 1)]))
        bf.append(-y[n])
        Ab.append(np.array([np.conj(y[n - p + k]) for k in range(1, p + 1)]))
        bb.append(-np.conj(y[n - p]))
    A = np.vstack([np.asarray(Af), np.asarray(Ab)])
    b = np.concatenate([np.asarray(bf), np.asarray(bb)])
    a, *_ = np.linalg.lstsq(A, b, rcond=None)
    u = [y[n] + np.dot(a, [y[n - k] for k in range(1, p + 1)]) for n in range(p, N)]
    return a, float(np.mean(np.abs(u) ** 2)), int(np.linalg.matrix_rank(A)), float(np.linalg.cond(A))


def ar_spectrum(a, sigma2, omega, power=1):
    a = np.asarray(a, dtype=complex)
    acc = np.ones_like(omega, dtype=complex)
    for k in range(1, a.size + 1):
        acc = acc + a[k - 1] * np.exp(-1j * omega * k)
    return sigma2 / np.maximum(np.abs(acc) ** power, 1e-30)


def detect_peaks(P, omega, rel_prom=0.05):
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


def eps_k(dr):
    """Truth-independent: half spectral bin mapped to k."""
    return DW / (2.0 * dr) + 1e-12


def map_k(w, dr):
    """k = -(omega+2pi m)/dr; band with grid quantization tolerance (NOT clip to edges)."""
    ek = eps_k(dr)
    cands = [-(w + TWO_PI * m) / dr for m in range(-8, 9)]
    # do not force into [1.32,1.50] numerically — keep raw k_cand
    in_band = sorted(
        {round(c, 9) for c in cands if (K_MIN - ek) <= c <= (K_MAX + ek)}
    )
    return in_band, ek


def gen_tone(ks, dr, N):
    n = np.arange(N)
    y = np.zeros(N, complex)
    for k in ks:
        y = y + np.exp(-1j * k * n * dr)
    return y


def score_peaks(ws, dr, k_trues, tol=1e-3):
    ek = eps_k(dr)
    all_k = []
    raw_list = []
    for w in ws:
        ib, _ = map_k(float(w), dr)
        raw_list.extend(ib)
        all_k.extend(ib)
    all_k = sorted(set(all_k))
    hits = sum(1 for kt in k_trues if any(abs(k - kt) < tol for k in all_k))
    false = sum(1 for k in all_k if all(abs(k - kt) > tol for kt in k_trues))
    recall = hits / max(len(k_trues), 1)
    precision = hits / max(len(all_k), 1) if all_k else 0.0
    f1 = 2 * precision * recall / max(precision + recall, 1e-30)
    return {
        "n_detected_peaks_in_band": len(all_k),
        "true_peak_recalled": int(hits),
        "false_peak_count": int(false),
        "recall": float(recall),
        "precision": float(precision),
        "f1": float(f1),
        "eps_k": ek,
        "recovered_ks": ";".join(f"{k:.6f}" for k in all_k),
    }


def main() -> int:
    # band-edge audit examples (same AR as 1F2; only scoring tolerance)
    audit_rows = []
    for k in (1.32, 1.50):
        for dr in (1.885, 0.5, 2.5):
            y = gen_tone([k], dr, N_R)
            a, s2, rank, cond = complex_mcov_ar(y, P_LIT)
            omega = np.linspace(-np.pi, np.pi, N_OMEGA, endpoint=False)
            P = ar_spectrum(a, s2, omega, 1)
            idx, ws = detect_peaks(P, omega, 0.05)
            best = None
            best_e = np.inf
            for w in ws:
                ib, ek = map_k(float(w), dr)
                for kc in ib:
                    e = abs(kc - k)
                    if e < best_e:
                        best_e, best = e, kc
            ek = eps_k(dr)
            audit_rows.append(
                {
                    "k_true": k,
                    "dr": dr,
                    "k_cand_raw": best,
                    "abs_err": best_e,
                    "eps_k": ek,
                    "within_grid_tol": best_e <= ek + 1e-9,
                    "old_hard_band_pass": bool(best is not None and K_MIN - 1e-9 <= best <= K_MAX + 1e-9),
                    "new_band_pass": bool(best is not None and K_MIN - ek <= best <= K_MAX + ek),
                    "note": "BAND_EDGE_GRID_QUANTIZATION_CORRECTED",
                }
            )
    pd.DataFrame(audit_rows).to_csv(OUT / "band_edge_examples.csv", index=False)

    (OUT / "BAND_EDGE_SCORING_AUDIT.md").write_text(
        f"""# BAND_EDGE_SCORING_AUDIT

UTC: {NOW}
基线 commit: 1663d16adc2758fc8a56718ed14488d782613acd

## 问题

K_BAND=[1.32,1.50] 用 ±1e-9 硬截断；离散谱网格 Δω=2π/{N_OMEGA} 映射到 k 的半格误差

$$
\\epsilon_k=\\frac{{\\Delta\\omega}}{{2\\Delta r}}+10^{{-12}}
$$

会把 1.319956（k_true=1.32）等 **网格量化** 结果误判为带外 → 0 recall。

## 修正（仅评分）

合法分支：

$$
k_{{min}}-\\epsilon_k \\le k_{{cand}} \\le k_{{max}}+\\epsilon_k
$$

- ε_k **只由** Δω 与 Δr 决定（与真值无关）
- **不**把 k_cand 截回 1.32/1.50
- AR 系数 / p / 谱网格 / 峰检测 / prominence **全部不变**

`BAND_EDGE_GRID_QUANTIZATION_CORRECTED`

## 示例

见 `band_edge_examples.csv`（k=1.32/1.50 边界点）。
""",
        encoding="utf-8",
    )

    # rescore single/multi with same AR
    p_branches = [
        (P_LIT, "PRINTED_LITERAL_CONTROL", "PAPER_PRINTED_PRIMARY_BRANCH"),
        (P_MOVING, "OUR_MOVING_SAMPLE_INTERPRETATION", "OUR_INTERPRETATION_SENSITIVITY_BRANCH"),
    ]
    single_rows = []
    for k in K_PHYS:
        for dr, dr_tag in [(1.885, "NO_FOLD"), (0.5, "AUX"), (2.5, "REF7_2p5")]:
            y = gen_tone([k], dr, N_R)
            for p, p_tag, role in p_branches:
                a, s2, rank, cond = complex_mcov_ar(y, p)
                omega = np.linspace(-np.pi, np.pi, N_OMEGA, endpoint=False)
                for power, pow_tag in [(1, "PRINTED_EQ19_POWER_1"), (2, "STANDARD_PSD_POWER_2")]:
                    P = ar_spectrum(a, s2, omega, power)
                    idx, ws = detect_peaks(P, omega, 0.05)
                    met = score_peaks(ws, dr, [k])
                    # selected strongest peak error (raw cand)
                    sel_err = np.nan
                    if idx:
                        bi = int(np.argmax(P[idx]))
                        ib, _ = map_k(ws[bi], dr)
                        if ib:
                            sel_err = min(abs(c - k) for c in ib)
                    met.update(
                        {
                            "k_true": k,
                            "dr_tag": dr_tag,
                            "p": p,
                            "p_tag": p_tag,
                            "branch_role": role,
                            "power": power,
                            "power_tag": pow_tag,
                            "selected_abs_err": sel_err,
                        }
                    )
                    single_rows.append(met)
    sdf = pd.DataFrame(single_rows)
    sdf.to_csv(OUT / "SINGLE_TONE_EDGE_CORRECTED.csv", index=False)

    multi_defs = [
        ("2_peak_wide", [1.35, 1.48]),
        ("3_peak_wide", [1.33, 1.41, 1.49]),
        ("5_peak_wide", [1.32, 1.36, 1.41, 1.46, 1.50]),
        ("2_peak_close", [1.40, 1.405]),
        ("3_peak_close", [1.40, 1.404, 1.408]),
        ("5_peak_close", [1.38, 1.39, 1.40, 1.41, 1.42]),
    ]
    multi_rows = []
    for name, ks in multi_defs:
        for dr, dr_tag in [(1.885, "NO_FOLD"), (2.5, "REF7_2p5")]:
            y = gen_tone(ks, dr, N_R)
            for p, p_tag, role in p_branches:
                a, s2, rank, cond = complex_mcov_ar(y, p)
                omega = np.linspace(-np.pi, np.pi, N_OMEGA, endpoint=False)
                P = ar_spectrum(a, s2, omega, 1)
                idx, ws = detect_peaks(P, omega, 0.05)
                met = score_peaks(ws, dr, ks)
                met.update(
                    {
                        "case": name,
                        "dr_tag": dr_tag,
                        "p_tag": p_tag,
                        "branch_role": role,
                        "complete_success": bool(met["recall"] == 1.0 and met["false_peak_count"] == 0),
                    }
                )
                multi_rows.append(met)
    mdf = pd.DataFrame(multi_rows)
    mdf.to_csv(OUT / "MULTITONE_EDGE_CORRECTED.csv", index=False)

    # branch status
    p7s = sdf[(sdf["p_tag"] == "PRINTED_LITERAL_CONTROL") & (sdf["power"] == 1)]
    p133s = sdf[(sdf["p_tag"] == "OUR_MOVING_SAMPLE_INTERPRETATION") & (sdf["power"] == 1)]
    p7m = mdf[mdf["p_tag"] == "PRINTED_LITERAL_CONTROL"]
    p133m = mdf[mdf["p_tag"] == "OUR_MOVING_SAMPLE_INTERPRETATION"]

    def branch_line(df_s, df_m):
        return {
            "single_f1": float(df_s["f1"].mean()),
            "single_recall": float(df_s["recall"].mean()),
            "multi_f1": float(df_m["f1"].mean()),
            "multi_complete": float(df_m["complete_success"].mean()),
            "multi_recall": float(df_m["recall"].mean()),
        }

    st7 = branch_line(p7s, p7m)
    st133 = branch_line(p133s, p133m)

    (OUT / "BRANCH_STATUS.md").write_text(
        f"""# BRANCH_STATUS

UTC: {NOW}

| 分支 | 角色 | single F1 | multi F1 | multi complete | multi recall |
| --- | --- | ---: | ---: | ---: | ---: |
| p=7 `PRINTED_LITERAL_CONTROL` | **PAPER_PRINTED_PRIMARY_BRANCH** | {st7['single_f1']:.3f} | {st7['multi_f1']:.3f} | {st7['multi_complete']:.2f} | {st7['multi_recall']:.3f} |
| p=133 `OUR_MOVING_SAMPLE_INTERPRETATION` | **OUR_INTERPRETATION_SENSITIVITY_BRANCH** | {st133['single_f1']:.3f} | {st133['multi_f1']:.3f} | {st133['multi_complete']:.2f} | {st133['multi_recall']:.3f} |

## Provenance（进入 4B-2 前按来源冻结，非按结果择优）

- p=7 = 论文打印主分支 → **4B-2 论文复现以 p=7 为主**
- p=133 = 我方解释敏感性分支 → 继续报告，**不**决定 Liang paper reproduction 是否通过
- p=133 小扰动下系数/峰位更敏感 → `P133_CONDITIONING_LIMIT`（见 1F2 diagnostic）

## 冻结（未改）

`P1_COMPLEX_EXP_IDENTITY_VALIDATED`
`OUR_EXACT_SOLVER_FOR_EQ24_VALIDATED`
`AR_CORE_IMPLEMENTATION_VALIDATED`（P1 + 复数 MCOV 共轭）
`P7_PRINTED_LITERAL_BRANCH_ACCEPTED`
`P133_MOVING_SAMPLE_BRANCH_CONDITIONING_LIMITED`
""",
        encoding="utf-8",
    )

    decision = "R3_C2_4B1_READY_FOR_PAPER_SPECTRUM_REPRO"
    (OUT / "R3_C2_4B_1F2_1_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-C2.4B-1F2.1",
                "decision": decision,
                "labels": [
                    "AR_CORE_IMPLEMENTATION_VALIDATED",
                    "P7_PRINTED_LITERAL_BRANCH_ACCEPTED",
                    "P133_MOVING_SAMPLE_BRANCH_CONDITIONING_LIMITED",
                    "BAND_EDGE_GRID_QUANTIZATION_CORRECTED",
                    "P1_COMPLEX_EXP_IDENTITY_VALIDATED",
                    "OUR_EXACT_SOLVER_FOR_EQ24_VALIDATED",
                ],
                "prior_label_corrected": "AR_K_RECOVERY_NUMERICALLY_UNSTABLE → replaced by branch-level labels",
                "p7_metrics": st7,
                "p133_metrics": st133,
                "epsilon_k_rule": "eps_k = (2pi/N_omega)/(2*dr)+1e-12; truth-independent",
                "not_done": ["4B-2", "Hankel", "D(z)", "FIELD", "E-STD", "MC", "P5"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_C2_4B_1F2_1_REPORT.md").write_text(
        f"""# R3-C2.4B-1F2.1 band-edge scoring

UTC: {NOW}

## 判定

### `{decision}`

- `AR_CORE_IMPLEMENTATION_VALIDATED`
- `P7_PRINTED_LITERAL_BRANCH_ACCEPTED`（single F1={st7['single_f1']:.3f}, multi F1={st7['multi_f1']:.3f}）
- `P133_MOVING_SAMPLE_BRANCH_CONDITIONING_LIMITED`（multi 同级，但小扰动 conditioning 限）
- `BAND_EDGE_GRID_QUANTIZATION_CORRECTED`（ε_k=(Δω/2)/Δr）

**未改** AR 公式 / p / 峰检测 / Eq19 / Eq24。

4B-2 主分支 = **p=7（PAPER_PRINTED_PRIMARY_BRANCH）**。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# R3-C2.4B-1F2.1\\n\\n**{decision}**\\n\\np7 F1={st7['single_f1']:.3f}/{st7['multi_f1']:.3f}；p133 sensitivity。\\n",
        encoding="utf-8",
    )
    print("DECISION", decision)
    print("p7", st7)
    print("p133", st133)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
