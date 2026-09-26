#!/usr/bin/env python3
"""R3-RC3-REANCHOR-5A: trackable-line availability boundary (15 frequency subsets)."""
from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_LINE_AVAILABILITY_BOUNDARY"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
PAIRS = ROOT / "results" / "P3_RC3_increment" / "selected_hard_pairs.csv"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0
ZR = 200.0
Z_TRUE = [180.0, 200.0, 220.0]
Z_GRID = np.arange(150.0, 250.0 + 1e-9, 5.0)
FREQS = [201.0, 235.0, 283.0, 338.0]


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


def platform_xy(t, turn_deg=0.0, u=U_PLAT):
    t = np.asarray(t, float)
    xp = u * t.copy()
    yp = np.zeros_like(t)
    if abs(turn_deg) > 1e-12:
        t_turn = float(np.max(t)) * 0.5
        d = math.radians(turn_deg)
        c, s = math.cos(d), math.sin(d)
        m = t > t_turn
        dt = t[m] - t_turn
        xp[m] = u * t_turn + u * dt * c
        yp[m] = u * dt * s
    return xp, yp


def range_traj(t, r0_km, theta0_deg, v, psi_deg, turn_deg):
    r0 = float(r0_km) * 1000.0
    th0 = math.radians(float(theta0_deg))
    psi = math.radians(float(psi_deg))
    xp, yp = platform_xy(t, turn_deg)
    xt = r0 * np.cos(th0) + v * t * np.cos(psi)
    yt = r0 * np.sin(th0) + v * t * np.sin(psi)
    return np.hypot(xt - xp, yt - yp)


def pressure_series(mod, r, zs, zr=ZR):
    depths, phi, k = mod["depths"], mod["phi"], mod["k"]
    izs = int(np.argmin(np.abs(depths - zs)))
    izr = int(np.argmin(np.abs(depths - zr)))
    kre, alpha = k.real, -k.imag
    ps, pr = phi[izs], phi[izr]
    r = np.asarray(r, float)
    p = np.zeros(r.size, dtype=complex)
    for m in range(len(kre)):
        p += (
            np.sqrt(2 * np.pi / (kre[m] * r))
            * ps[m]
            * pr[m]
            * np.exp(-1j * kre[m] * r - alpha[m] * r - 1j * np.pi / 4)
        )
    return p


def tl_shape(p):
    L = 20.0 * np.log10(np.maximum(np.abs(p), 1e-30))
    return L - float(np.mean(L))


def rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def all_subsets(freqs):
    out = []
    for k in range(1, len(freqs) + 1):
        for comb in itertools.combinations(freqs, k):
            out.append(comb)
    return out


def main() -> int:
    mods = {f: parse_mod(ZGRID / f"zgrid_f{int(f)}.mod") for f in FREQS}
    subsets = all_subsets(FREQS)

    (OUT / "LINE_AVAILABILITY_AXIS_LOCK.md").write_text(
        f"""# LINE_AVAILABILITY_AXIS_LOCK

UTC: {NOW}

`S1_BRIDGE_LINE_AVAILABILITY_ONLY` / `TRACKABLE_LINE_AVAILABILITY_AXIS`

- 环境：E0_NOMINAL（无 SSP mismatch / 无 IID TL error / 无频漂）
- F4 = {{201,235,283,338}} Hz；全部 **15** 个非空子集冻结运行
- 每频独立去均值 TL 形状；多频拼接
- 双侧 profile z∈150:5:250
- 主量：ΔJ = J_alt* − J_true*
- 子集 PASS：fraction_tested(ΔJ>0)≥0.8 且 median ΔJ>0（路线 Gate，非统计检验）
- 仍假定：所选线谱持续存在、中心频率稳定、幅度变化不污染 TL

**不是** 完整 S1；不含 168/204/232/279/320 Hz 迁移。
""",
        encoding="utf-8",
    )

    pairs = pd.read_csv(PAIRS)
    if "selected" in pairs.columns:
        pairs = pairs[pairs["selected"] == True]

    # precompute per-freq TL for each pair / traj / z
    # then evaluate all subsets
    subset_rows = []
    pair_rows = []
    single_corr_rows = []

    for _, row in pairs.iterrows():
        pid = row["pair_id"]
        mech = row.get("mechanism", "")
        T = float(row["T_s"])
        turn = float(row["turn_deg"])
        t = np.arange(0.0, T + 1e-9, 1.0)
        r_ref = range_traj(
            t, row["ref_r0_km"], row["ref_theta0_deg"], row["ref_v"], row["ref_psi_deg"], turn
        )
        r_alt = range_traj(
            t, row["alt_r0_km"], row["alt_theta0_deg"], row["alt_v"], row["alt_psi_deg"], turn
        )
        static_blind = pid == "A05"

        # TL_ref[f][z], TL_alt[f][z]
        TL_ref = {f: {} for f in FREQS}
        TL_alt = {f: {} for f in FREQS}
        for f in FREQS:
            for z in Z_GRID:
                TL_ref[f][float(z)] = tl_shape(pressure_series(mods[f], r_ref, float(z)))
                TL_alt[f][float(z)] = tl_shape(pressure_series(mods[f], r_alt, float(z)))
        for z_true in Z_TRUE:
            y_single = {f: TL_ref[f][float(z_true)] for f in FREQS}
            # single-freq ΔJ for correlation diagnostics
            dJ_single = {}
            for f in FREQS:
                Jt = min(rms(y_single[f], TL_ref[f][float(z)]) for z in Z_GRID)
                Ja = min(rms(y_single[f], TL_alt[f][float(z)]) for z in Z_GRID)
                dJ_single[f] = Ja - Jt
            if not static_blind:
                vals = [dJ_single[f] for f in FREQS]
                single_corr_rows.append(
                    {
                        "pair_id": pid,
                        "mechanism": mech,
                        "z_true": z_true,
                        "dJ_201": dJ_single[201.0],
                        "dJ_235": dJ_single[235.0],
                        "dJ_283": dJ_single[283.0],
                        "dJ_338": dJ_single[338.0],
                        "dJ_min": float(np.min(vals)),
                        "dJ_max": float(np.max(vals)),
                        "dJ_range": float(np.max(vals) - np.min(vals)),
                    }
                )

            for comb in subsets:
                y = np.concatenate([y_single[f] for f in comb])
                Jt = min(
                    rms(y, np.concatenate([TL_ref[f][float(z)] for f in comb]))
                    for z in Z_GRID
                )
                Ja = min(
                    rms(y, np.concatenate([TL_alt[f][float(z)] for f in comb]))
                    for z in Z_GRID
                )
                dJ = Ja - Jt
                pair_rows.append(
                    {
                        "pair_id": pid,
                        "mechanism": mech,
                        "z_true": z_true,
                        "n_lines": len(comb),
                        "frequencies": ";".join(str(int(f)) for f in comb),
                        "J_true_star": Jt,
                        "J_alt_star": Ja,
                        "delta_J": dJ,
                        "delta_J_gt_0": bool(dJ > 0),
                        "delta_J_gt_0p5": bool(dJ > 0.5),
                        "static_blind": static_blind,
                    }
                )

    df = pd.DataFrame(pair_rows)
    df.to_csv(OUT / "SUBSET_PROFILED_MARGIN_RESULTS.csv", index=False)
    pd.DataFrame(single_corr_rows).to_csv(
        OUT / "SINGLE_FREQUENCY_MARGIN_CORRELATION.csv", index=False
    )

    # subset summary (exclude A05)
    sum_rows = []
    for comb in subsets:
        key = ";".join(str(int(f)) for f in comb)
        g = df[(df["frequencies"] == key) & (~df["static_blind"])]
        dJ = g["delta_J"]
        frac0 = float((dJ > 0).mean())
        med = float(dJ.median())
        ok = frac0 >= 0.8 and med > 0
        sum_rows.append(
            {
                "n_lines": len(comb),
                "frequencies": key,
                "median_delta_J": med,
                "p10_delta_J": float(dJ.quantile(0.10)),
                "min_delta_J": float(dJ.min()),
                "fraction_tested_dJ_gt_0": frac0,
                "fraction_tested_dJ_gt_0p5": float((dJ > 0.5).mean()),
                "subset_pass": ok,
                "status": "SUBSET_TRAJECTORY_INFORMATION_SURVIVES" if ok else "SUBSET_FAIL",
            }
        )
    sdf = pd.DataFrame(sum_rows)
    sdf.to_csv(OUT / "SUBSET_SUMMARY.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(OUT / "ALL_15_FREQUENCY_SUBSETS.csv", index=False)

    # line count summary
    lc_rows = []
    for n in (1, 2, 3, 4):
        gs = sdf[sdf["n_lines"] == n]
        if not len(gs):
            continue
        worst = gs.loc[gs["fraction_tested_dJ_gt_0"].idxmin()]
        best = gs.loc[gs["fraction_tested_dJ_gt_0"].idxmax()]
        lc_rows.append(
            {
                "n_lines": n,
                "n_subsets": len(gs),
                "worst_subset": worst["frequencies"],
                "worst_frac_gt0": worst["fraction_tested_dJ_gt_0"],
                "worst_median_dJ": worst["median_delta_J"],
                "best_subset": best["frequencies"],
                "best_frac_gt0": best["fraction_tested_dJ_gt_0"],
                "median_of_subset_medians": float(gs["median_delta_J"].median()),
                "n_pass": int(gs["subset_pass"].sum()),
                "all_pass": bool(gs["subset_pass"].all()),
            }
        )
    ldf = pd.DataFrame(lc_rows)
    ldf.to_csv(OUT / "LINE_COUNT_SUMMARY.csv", index=False)

    mech_rows = []
    for n in (1, 2, 3, 4):
        for mech in ["A", "B", "C"]:
            g = df[(df["n_lines"] == n) & (df["mechanism"] == mech) & (~df["static_blind"])]
            if not len(g):
                continue
            dJ = g["delta_J"]
            mech_rows.append(
                {
                    "n_lines": n,
                    "mechanism": mech,
                    "median_delta_J": float(dJ.median()),
                    "fraction_tested_dJ_gt_0": float((dJ > 0).mean()),
                    "n": len(g),
                }
            )
    pd.DataFrame(mech_rows).to_csv(OUT / "MECHANISM_BY_LINE_COUNT.csv", index=False)

    (OUT / "STATIC_RANGE_BLIND_CASE.md").write_text(
        f"""# STATIC_RANGE_BLIND_CASE

UTC: {NOW}

`RELATIVE_TL_STATIC_RANGE_BLIND_CASE`：A05 参考距离恒定，去均值 TL≡0。
不进入主统计。结构边界，非“100%”。

B/C 的 180/200/220 签名非重复（4G）；A01 因距离恒定，深度签名重复属预期。
""",
        encoding="utf-8",
    )

    # decision
    def all_pass(n):
        g = sdf[sdf["n_lines"] == n]
        return len(g) > 0 and bool(g["subset_pass"].all())

    if all_pass(1):
        decision = "ANY_SINGLE_OF_FOUR_TRACKABLE_LINES_SUFFICIENT_IDEAL_CONTROL"
        why = "4 个单频子集全部 PASS（理想稳定线谱条件）"
    elif all_pass(2):
        decision = "ANY_TWO_OF_FOUR_TRACKABLE_LINES_SUFFICIENT_IDEAL_CONTROL"
        why = "6 个双频全部 PASS；单频未全过"
    elif all_pass(3):
        decision = "THREE_TRACKABLE_LINES_REQUIRED_IDEAL_CONTROL"
        why = "4 个三频全部 PASS"
    elif all_pass(4):
        # check if 1/2/3 mixed
        n_pass_low = int(sdf[sdf["n_lines"] < 4]["subset_pass"].sum())
        if n_pass_low > 0 and not all_pass(2):
            decision = "LINE_COMBINATION_DEPENDENT_IDEAL_CONTROL"
            why = "部分低线数子集 PASS、部分 FAIL"
        else:
            decision = "FOUR_TRACKABLE_LINES_REQUIRED_IDEAL_CONTROL"
            why = "仅四频组合稳定通过"
    else:
        n1 = int(sdf[sdf["n_lines"] == 1]["subset_pass"].sum())
        n2 = int(sdf[sdf["n_lines"] == 2]["subset_pass"].sum())
        if (n1 > 0 and n1 < 4) or (n2 > 0 and n2 < 6):
            decision = "LINE_COMBINATION_DEPENDENT_IDEAL_CONTROL"
            why = "同线数内组合分裂"
        else:
            decision = "LINE_AVAILABILITY_BOUNDARY_UNRESOLVED"

    (OUT / "R3_RC3_REANCHOR_5A_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-RC3-REANCHOR-5A",
                "decision": decision,
                "why": why,
                "axis": "S1_BRIDGE_LINE_AVAILABILITY_ONLY",
                "n_subsets": len(subsets),
                "subset_pass_by_n": {
                    str(n): int(sdf[sdf["n_lines"] == n]["subset_pass"].sum()) for n in (1, 2, 3, 4)
                },
                "frozen": [
                    "DEPTH_ROLE_MIXED_OR_UNRESOLVED",
                    "NO_STRONG_DEPTH_COMPENSATION_OBSERVED",
                    "DEPTH_PROFILED_NUISANCE_VARIABLE",
                    "PROFILED_RC3_MARGIN_ROBUST_TO_TESTED_ERROR",
                    "PROFILED_RC3_MARGIN_SURVIVES_TESTED_SSP_STRESS",
                ],
                "not_done": [
                    "freq drift",
                    "line amplitude drift",
                    "S1 full",
                    "SSP mismatch",
                    "IID TL",
                    "Liang",
                    "P5",
                    "new pairs",
                ],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_RC3_REANCHOR_5A_REPORT.md").write_text(
        f"""# R3-RC3-REANCHOR-5A trackable-line availability

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 轴

`S1_BRIDGE_LINE_AVAILABILITY_ONLY` · E0 · 15 非空子集 · 双侧 depth profile · A05 剔除

子集 PASS：fraction(ΔJ>0)≥0.8 且 median ΔJ>0

## 按线数

{ldf.to_string(index=False) if len(ldf) else 'n/a'}

## 含义边界

理想稳定线谱下的**结构性**下限；**不是** S1 已验证。
下一阶段：5B 线谱幅度/频率不稳定性。

## 文档修正

- index：4F→**4G**，目录 `..._FIX2`
- 4G 文案：B/C 签名非重复；A01 静态重复属预期（不改 CSV）
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# R3-RC3-REANCHOR-5A\\n\\n**{decision}**\\n\\n{why}\\n",
        encoding="utf-8",
    )
    print("DECISION", decision)
    print(ldf.to_string(index=False) if len(ldf) else "no lc")
    print(sdf[["n_lines", "frequencies", "median_delta_J", "fraction_tested_dJ_gt_0", "subset_pass"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
