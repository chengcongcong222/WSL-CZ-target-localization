#!/usr/bin/env python3
"""R3-RC3-REANCHOR-3: observation-error boundary for depth-profiled RC3 margin."""
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_OBSERVATION_ERROR_BOUNDARY"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
PAIRS = ROOT / "results" / "P3_RC3_increment" / "selected_hard_pairs.csv"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0
ZR = 200.0
Z_TRUE = [180.0, 200.0, 220.0]
Z_GRID = np.arange(150.0, 250.0 + 1e-9, 5.0)
FREQS = [201.0, 235.0, 283.0, 338.0]
EPS_AXIS = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0]
SEEDS = [0, 1, 2, 3, 4]


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
    kre = k.real
    alpha = -k.imag
    ps = phi[izs]
    pr = phi[izr]
    r = np.asarray(r, float)
    p = np.zeros(r.size, dtype=complex)
    for m in range(len(kre)):
        amp = np.sqrt(2 * np.pi / (kre[m] * r))
        p += amp * ps[m] * pr[m] * np.exp(-1j * kre[m] * r - alpha[m] * r - 1j * np.pi / 4)
    return p


def tl_shape(p):
    L = 20.0 * np.log10(np.maximum(np.abs(p), 1e-30))
    return L - float(np.mean(L))


def rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def main() -> int:
    pairs = pd.read_csv(PAIRS)
    if "selected" in pairs.columns:
        pairs = pairs[pairs["selected"] == True]
    mods = {f: parse_mod(ZGRID / f"zgrid_f{int(f)}.mod") for f in FREQS}

    (OUT / "OBSERVATION_ERROR_AXIS_LOCK.md").write_text(
        f"""# OBSERVATION_ERROR_AXIS_LOCK

UTC: {NOW}

`TL_OBSERVATION_ERROR_SENSITIVITY_AXIS`：

ε = {EPS_AXIS} dB（每频零均值 TL 形状扰动，扰动后再去均值）

seeds = {SEEDS}（有限敏感性，**非** Monte Carlo 性能概率）

信号条件：`S2_LIKE_TRACKABLE_FOUR_TONE_UPPER_BOUND`（201/235/283/338 Hz 持续可跟踪；每频独立去均值；不要求源级；不声称 S0/S1）

指标：ΔJ = J_alt* − J_true*（两者均对 z∈150:5:250 profile）

角色不变：`DEPTH_PROFILED_NUISANCE_VARIABLE`；标签 `IDEAL_PROFILED_PROPAGATION_MARGIN_OBSERVED`（非 96.7% 性能）。
""",
        encoding="utf-8",
    )

    rows = []
    for _, row in pairs.iterrows():
        pid = row["pair_id"]
        mech = row.get("mechanism", "")
        T = float(row["T_s"])
        turn = float(row["turn_deg"])
        t = np.arange(0.0, T + 1e-9, 1.0)
        r_ref = range_traj(t, row["ref_r0_km"], row["ref_theta0_deg"], row["ref_v"], row["ref_psi_deg"], turn)
        r_alt = range_traj(t, row["alt_r0_km"], row["alt_theta0_deg"], row["alt_v"], row["alt_psi_deg"], turn)

        # precompute TL signatures for all z (ref: only z_true needed for true traj; still need all z for J_true*)
        sig_true = {}  # z -> concat TL
        sig_alt = {}
        for z in Z_GRID:
            st = []
            sa = []
            for f in FREQS:
                st.append(tl_shape(pressure_series(mods[f], r_ref, z)))
                sa.append(tl_shape(pressure_series(mods[f], r_alt, z)))
            sig_true[float(z)] = np.concatenate(st)
            sig_alt[float(z)] = np.concatenate(sa)

        for z_true in Z_TRUE:
            y_true = sig_true[float(z_true)]
            static_blind = bool(pid == "A05" or (np.allclose(r_ref, r_ref[0]) and np.allclose(r_alt, r_alt[0])))
            for eps in EPS_AXIS:
                for seed in SEEDS:
                    rng = np.random.default_rng(1000 * seed + int(eps * 100) + int(z_true))
                    if eps == 0.0:
                        y_obs = y_true.copy()
                    else:
                        # per-frequency noise then re-demean each freq block
                        nseg = len(y_true) // len(FREQS)
                        y_obs = y_true.copy()
                        for fi in range(len(FREQS)):
                            sl = slice(fi * nseg, (fi + 1) * nseg)
                            e = rng.normal(0.0, eps, nseg)
                            e = e - e.mean()
                            y_obs[sl] = y_true[sl] + e
                    Jt = min(rms(y_obs, sig_true[float(z)]) for z in Z_GRID)
                    Ja = min(rms(y_obs, sig_alt[float(z)]) for z in Z_GRID)
                    dJ = Ja - Jt
                    rows.append(
                        {
                            "pair_id": pid,
                            "mechanism": mech,
                            "z_true": z_true,
                            "eps_db": eps,
                            "seed": seed,
                            "J_true_star": Jt,
                            "J_alt_star": Ja,
                            "delta_J": dJ,
                            "delta_J_gt_0": bool(dJ > 0),
                            "delta_J_gt_0p5": bool(dJ > 0.5),
                            "static_blind": static_blind,
                            "blind_label": "RELATIVE_TL_STATIC_RANGE_BLIND_CASE" if static_blind else "",
                        }
                    )

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "PROFILED_MARGIN_VS_ERROR.csv", index=False)

    # by mechanism and eps (exclude A05 blind from main mech stats? keep but flag)
    sum_rows = []
    for mech, g0 in df.groupby("mechanism"):
        g = g0[~g0["static_blind"]]
        for eps, ge in g.groupby("eps_db"):
            dJ = ge["delta_J"]
            sum_rows.append(
                {
                    "mechanism": mech,
                    "eps_db": eps,
                    "median_delta_J": float(dJ.median()),
                    "p10_delta_J": float(dJ.quantile(0.10)),
                    "frac_dJ_gt_0": float((dJ > 0).mean()),
                    "frac_dJ_gt_0p5": float((dJ > 0.5).mean()),
                    "median_J_true_star": float(ge["J_true_star"].median()),
                    "median_J_alt_star": float(ge["J_alt_star"].median()),
                    "n": len(ge),
                }
            )
    sdf = pd.DataFrame(sum_rows)
    sdf.to_csv(OUT / "MARGIN_BY_MECHANISM.csv", index=False)

    (OUT / "STATIC_RANGE_BLIND_CASE.md").write_text(
        f"""# STATIC_RANGE_BLIND_CASE

UTC: {NOW}

`RELATIVE_TL_STATIC_RANGE_BLIND_CASE`：A05（及任何 r(t) 恒定对）

相对 TL 形状去均值后 ≡ 0，**无法**用传播形状区分两个静态绝对距离。
不参与机制 margin 统计；支持“RC3=连续轨迹一致性约束，而非单帧绝对距离仪”。
""",
        encoding="utf-8",
    )

    # decision at each eps (all mech pooled, excluding blind)
    nd = df[~df["static_blind"]]
    eps_stats = []
    for eps in EPS_AXIS:
        ge = nd[nd["eps_db"] == eps]
        dJ = ge["delta_J"]
        med = float(dJ.median())
        frac0 = float((dJ > 0).mean())
        eps_stats.append({"eps": eps, "median_dJ": med, "frac_gt0": frac0, "frac_gt0p5": float((dJ > 0.5).mean())})
    pd.DataFrame(eps_stats).to_csv(OUT / "margin_vs_error_summary.csv", index=False)

    # robust if frac(dJ>0)>=0.8 at eps=2.0; conditional if good at 0.5-1.0 but not 2.0; collapse if <0.5 at 0.5
    f2 = next(s for s in eps_stats if s["eps"] == 2.0)
    f1 = next(s for s in eps_stats if s["eps"] == 1.0)
    f05 = next(s for s in eps_stats if s["eps"] == 0.5)
    if f2["frac_gt0"] >= 0.8 and f2["median_dJ"] > 0:
        decision = "PROFILED_RC3_MARGIN_ROBUST_TO_TESTED_ERROR"
    elif f1["frac_gt0"] >= 0.7:
        decision = "PROFILED_RC3_MARGIN_CONDITIONAL_ON_LOW_ERROR"
    elif f05["frac_gt0"] < 0.5:
        decision = "PROFILED_RC3_MARGIN_COLLAPSES_WITH_OBSERVATION_ERROR"
    else:
        decision = "PROFILED_RC3_OBSERVATION_ERROR_BOUNDARY_UNRESOLVED"

    why = (
        f"ΔJ>0 比例: ε=0.5→{f05['frac_gt0']:.2f}, 1.0→{f1['frac_gt0']:.2f}, 2.0→{f2['frac_gt0']:.2f}; "
        f"median ΔJ: {f05['median_dJ']:.2f}/{f1['median_dJ']:.2f}/{f2['median_dJ']:.2f} dB "
        f"(S2_LIKE_TRACKABLE_FOUR_TONE_UPPER_BOUND; 非实测噪声)"
    )

    (OUT / "R3_RC3_REANCHOR_3_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-RC3-REANCHOR-3",
                "decision": decision,
                "why": why,
                "signal_condition": "S2_LIKE_TRACKABLE_FOUR_TONE_UPPER_BOUND",
                "error_axis_label": "TL_OBSERVATION_ERROR_SENSITIVITY_AXIS",
                "margin_label": "IDEAL_PROFILED_PROPAGATION_MARGIN_OBSERVED",
                "roles_frozen": {
                    "science": "DEPTH_ROLE_MIXED_OR_UNRESOLVED",
                    "system": "DEPTH_PROFILED_NUISANCE_VARIABLE",
                },
                "eps_stats": eps_stats,
                "not_done": ["env mismatch", "S0/S1", "Liang", "P5", "depth RMSE"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_RC3_REANCHOR_3_REPORT.md").write_text(
        f"""# R3-RC3-REANCHOR-3 observation-error boundary

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 轴与样本

ε∈{EPS_AXIS} dB · seeds {SEEDS} · 30 对 × z_true 180/200/220 · z profile 150:5:250  
标签：`TL_OBSERVATION_ERROR_SENSITIVITY_AXIS`（非海试标定）  
信号：`S2_LIKE_TRACKABLE_FOUR_TONE_UPPER_BOUND`

## ΔJ = J_alt* − J_true*（均 profile）

| ε (dB) | median ΔJ | P(ΔJ>0) | P(ΔJ>0.5) |
| ---: | ---: | ---: | ---: |
| 0.5 | {f05['median_dJ']:.3f} | {f05['frac_gt0']:.3f} | {f05['frac_gt0p5']:.3f} |
| 1.0 | {f1['median_dJ']:.3f} | {f1['frac_gt0']:.3f} | {f1['frac_gt0p5']:.3f} |
| 2.0 | {f2['median_dJ']:.3f} | {f2['frac_gt0']:.3f} | {f2['frac_gt0p5']:.3f} |

A05：`RELATIVE_TL_STATIC_RANGE_BLIND_CASE`（已从机制统计剔除）

## 边界

角色不变；REANCHOR-2 的 96.7% 仅 `DECISION_THRESHOLD_SENSITIVITY` / 理想模型，非性能。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(f"# R3-RC3-REANCHOR-3\\n\\n**{decision}**\\n\\n{why}\\n", encoding="utf-8")
    print("DECISION", decision)
    print(pd.DataFrame(eps_stats).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
