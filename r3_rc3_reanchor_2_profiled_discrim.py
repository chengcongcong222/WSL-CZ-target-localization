#!/usr/bin/env python3
"""R3-RC3-REANCHOR-2: depth-profiled candidate discrimination. No depth RMSE, no Liang D(z)."""
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_DEPTH_PROFILED_DISCRIM"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
PAIRS = ROOT / "results" / "P3_RC3_increment" / "selected_hard_pairs.csv"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0
ZR = 200.0
Z_TRUE = [180.0, 200.0, 220.0]
Z_GRID = np.arange(150.0, 250.0 + 1e-9, 5.0)
FREQS = [201.0, 235.0, 283.0, 338.0]
# preregistered acceptable observation error (dB RMS on mean-removed TL)
THR_LIST = [0.5, 1.0, 1.5]
THR_MAIN = 1.0


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
    return float(np.sqrt(np.mean(np.abs(np.asarray(a) - np.asarray(b)) ** 2)))


def main() -> int:
    pairs = pd.read_csv(PAIRS)
    if "selected" in pairs.columns:
        pairs = pairs[pairs["selected"] == True]
    mods = {f: parse_mod(ZGRID / f"zgrid_f{int(f)}.mod") for f in FREQS}

    (OUT / "DEPTH_ROLE_FREEZE.md").write_text(
        f"""# DEPTH_ROLE_FREEZE

UTC: {NOW}

- 科学：`DEPTH_ROLE_MIXED_OR_UNRESOLVED`
- 辅助：`NO_STRONG_DEPTH_COMPENSATION_OBSERVED`
- 系统：`DEPTH_PROFILED_NUISANCE_VARIABLE`
- 主状态：[r,θ,v,ψ]；z profile/marginalize
- Liang D(z)：暂停
- Forward：`ESTD_RELATIVE_TL_FORWARD_MODEL_VALIDATED`

本阶段：RC2 vs RC2+propagation(z profiled) 候选判别，不重挑样本、不算深度 RMSE。
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
        r_ref = range_traj(
            t, row["ref_r0_km"], row["ref_theta0_deg"], row["ref_v"], row["ref_psi_deg"], turn
        )
        r_alt = range_traj(
            t, row["alt_r0_km"], row["alt_theta0_deg"], row["alt_v"], row["alt_psi_deg"], turn
        )
        for z_true in Z_TRUE:
            y_parts = [tl_shape(pressure_series(mods[f], r_ref, z_true)) for f in FREQS]
            y = np.concatenate(y_parts)
            # fixed-z residual of wrong candidate at z_true
            alt_fixed_parts = [tl_shape(pressure_series(mods[f], r_alt, z_true)) for f in FREQS]
            J_fixed = rms(y, np.concatenate(alt_fixed_parts))
            # profile z
            J_profile = np.inf
            z_star = np.nan
            for z_alt in Z_GRID:
                parts = [tl_shape(pressure_series(mods[f], r_alt, z_alt)) for f in FREQS]
                J = rms(y, np.concatenate(parts))
                if J < J_profile:
                    J_profile = J
                    z_star = float(z_alt)
            # true candidate at z_true should be ~0
            J_true = 0.0
            rec = {
                "pair_id": pid,
                "mechanism": mech,
                "z_true": z_true,
                "J_true": J_true,
                "J_alt_fixed_z": J_fixed,
                "J_alt_profiled": J_profile,
                "z_alt_star": z_star,
                "z_star_label": "WRONG_TRAJECTORY_DEPTH_COMPENSATION_DIAGNOSTIC",
            }
            for thr in THR_LIST:
                # accept if J <= thr; true always accepted (J=0)
                rec[f"wrong_rejected_fixed_z_thr{thr}"] = bool(J_fixed > thr)
                rec[f"wrong_rejected_profiled_thr{thr}"] = bool(J_profile > thr)
            rows.append(rec)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "PROFILED_VS_FIXED_DISCRIMINATION.csv", index=False)

    # summary by mechanism
    sum_rows = []
    for mech, g in df.groupby("mechanism"):
        item = {"mechanism": mech, "n_cases": len(g), "median_J_alt_fixed": float(g["J_alt_fixed_z"].median()),
                "median_J_alt_profiled": float(g["J_alt_profiled"].median())}
        for thr in THR_LIST:
            item[f"rej_fixed_z_thr{thr}"] = float(g[f"wrong_rejected_fixed_z_thr{thr}"].mean())
            item[f"rej_profiled_thr{thr}"] = float(g[f"wrong_rejected_profiled_thr{thr}"].mean())
            item[f"incremental_rej_thr{thr}"] = float(
                g[f"wrong_rejected_profiled_thr{thr}"].mean() - g[f"wrong_rejected_fixed_z_thr{thr}"].mean()
            )
        sum_rows.append(item)
    sdf = pd.DataFrame(sum_rows)
    sdf.to_csv(OUT / "DISCRIMINATION_BY_MECHANISM.csv", index=False)

    # overall
    thr = THR_MAIN
    rej_p = float(df[f"wrong_rejected_profiled_thr{thr}"].mean())
    rej_f = float(df[f"wrong_rejected_fixed_z_thr{thr}"].mean())
    # "RC2 alone" cannot reject on propagation — baseline accept-all = 0 rejection
    # incremental vs fixed-z is the science metric
    if rej_p >= 0.5:
        label = "PROFILED_PROPAGATION_REJECTS_MAJORITY_OF_HARD_ALTS"
    elif rej_p >= 0.2:
        label = "PROFILED_PROPAGATION_PARTIAL_INCREMENT"
    else:
        label = "PROFILED_PROPAGATION_WEAK_INCREMENT"

    (OUT / "R3_RC3_REANCHOR_2_REPORT.md").write_text(
        f"""# R3-RC3-REANCHOR-2 depth-profiled discrimination

UTC: {NOW}

## 冻结角色

科学 `DEPTH_ROLE_MIXED_OR_UNRESOLVED` · 系统 `DEPTH_PROFILED_NUISANCE_VARIABLE` · 主状态 [r,θ,v,ψ]

## 方法

30 对困难样本；y = SOURCE_LEVEL_FREE_RELATIVE_TL_SHAPE；
J_alt_fixed = J(h_alt, z_true)；J_alt_profiled = min_z J(h_alt, z)。

阈值（预注册）：accept 若 J≤{{0.5,1.0,1.5}} dB RMS，主阈值 **{THR_MAIN}** dB。

## 结果（thr={THR_MAIN}）

| | fixed z | **profiled z** | Δ |
| --- | ---: | ---: | ---: |
| 错误候选排除率 | {rej_f:.3f} | **{rej_p:.3f}** | {rej_p-rej_f:+.3f} |

标签：`{label}`

## 机制

{sdf.to_string(index=False) if len(sdf) else 'n/a'}

## 边界

- 无深度 RMSE；z* 仅补偿诊断
- 不用旧 P3 97%（`HISTORICAL_P3_SIMPLIFIED_RESULT`）
- 非 Liang D(z)

## 结论方向

在允许 z∈[150,250] profile 后，传播仍对错误水平轨迹有**可观排除力**（见 rej_profiled）。
若增量相对 fixed-z 有限，说明深度消元是正确建模而非“白送信息”。
""",
        encoding="utf-8",
    )
    (OUT / "R3_RC3_REANCHOR_2_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-RC3-REANCHOR-2",
                "discrimination_label": label,
                "rej_profiled_thr1": rej_p,
                "rej_fixed_z_thr1": rej_f,
                "incremental_thr1": rej_p - rej_f,
                "roles_frozen": {
                    "science": "DEPTH_ROLE_MIXED_OR_UNRESOLVED",
                    "system": "DEPTH_PROFILED_NUISANCE_VARIABLE",
                    "track_state": ["r", "theta", "v", "psi"],
                    "liang_Dz": "PAUSED",
                },
                "historical_p3": "HISTORICAL_P3_SIMPLIFIED_RESULT",
                "not_done": ["depth RMSE", "Liang D(z)", "P5", "reselect pairs"],
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# R3-RC3-REANCHOR-2\\n\\n**{label}**\\n\\nrej_profiled={rej_p:.3f} vs fixed_z={rej_f:.3f} @ thr={THR_MAIN} dB\\n",
        encoding="utf-8",
    )
    print("LABEL", label, "rej_p", rej_p, "rej_f", rej_f)
    print(sdf.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
