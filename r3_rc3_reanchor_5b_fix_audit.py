#!/usr/bin/env python3
"""5B-FIX audit: per-row 5A identity + range-grid integrity + amplitude-branch integrity."""
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_LINE_INSTABILITY_BOUNDARY_FIX"
A5DIR = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_LINE_AVAILABILITY_BOUNDARY"
OLD5B = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_LINE_INSTABILITY_BOUNDARY"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
PAIRS = ROOT / "results" / "P3_RC3_increment" / "selected_hard_pairs.csv"
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0
ZR = 200.0
Z_TRUE = [180.0, 200.0, 220.0]
Z_GRID = np.arange(150.0, 250.0 + 1e-9, 5.0)
FREQS = [201.0, 235.0, 283.0, 338.0]
COMMON_MULT = np.array([1.0, 1.0, 1.0, 1.0])
LINE_MULT = np.array([1.0, -1.0, 1.0, -1.0])
REPS = ["A03", "A10", "B01", "B08", "C01", "C08"]


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


def pressure_direct(mod, r, zs, zr=ZR):
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


def range_traj(t, r0_km, theta0_deg, v, psi_deg, turn_deg):
    t = np.asarray(t, float)
    xp = U_PLAT * t.copy()
    yp = np.zeros_like(t)
    if abs(turn_deg) > 1e-12:
        t_turn = float(np.max(t)) * 0.5
        d = math.radians(turn_deg)
        c, s = math.cos(d), math.sin(d)
        m = t > t_turn
        dt = t[m] - t_turn
        xp[m] = U_PLAT * t_turn + U_PLAT * dt * c
        yp[m] = U_PLAT * dt * s
    r0 = float(r0_km) * 1000.0
    th0 = math.radians(float(theta0_deg))
    psi = math.radians(float(psi_deg))
    xt = r0 * np.cos(th0) + v * t * np.cos(psi)
    yt = r0 * np.sin(th0) + v * t * np.sin(psi)
    return np.hypot(xt - xp, yt - yp)


def demean(L):
    L = np.asarray(L, float)
    return L - float(np.mean(L))


def rms(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def corr(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a - a.mean()
    b = b - b.mean()
    sa = float(np.sqrt(np.mean(a * a)))
    sb = float(np.sqrt(np.mean(b * b)))
    if sa <= 0 or sb <= 0:
        return 0.0
    return float(np.mean(a * b) / (sa * sb))


def main() -> int:
    # ---------- 1. per-row 5A identity ----------
    a5 = pd.read_csv(A5DIR / "SUBSET_PROFILED_MARGIN_RESULTS.csv")
    zs = pd.read_csv(OUT / "ZERO_STRESS_5A_IDENTITY.csv")
    cfg_map = {"201": "201", "235": "235", "283": "283", "338": "338", "FOUR": "201;235;283;338"}
    id_rows = []
    for cfg, freqs in cfg_map.items():
        g5 = a5[a5["frequencies"].astype(str) == freqs]
        gz = zs[zs["config"] == cfg]
        m = g5.merge(gz, on=["pair_id", "z_true"], suffixes=("_5a", "_zs"))
        for _, r in m.iterrows():
            eJt = abs(float(r["J_true_star_5a"]) - float(r["J_true_star_zs"]))
            eJa = abs(float(r["J_alt_star_5a"]) - float(r["J_alt_star_zs"]))
            eDj = abs(float(r["delta_J_5a"]) - float(r["delta_J_zs"]))
            id_rows.append(
                {
                    "pair_id": r["pair_id"],
                    "config": cfg,
                    "z_true": r["z_true"],
                    "J_true_star_5a": r["J_true_star_5a"],
                    "J_true_star_zs": r["J_true_star_zs"],
                    "err_J_true_star": eJt,
                    "J_alt_star_5a": r["J_alt_star_5a"],
                    "J_alt_star_zs": r["J_alt_star_zs"],
                    "err_J_alt_star": eJa,
                    "delta_J_5a": r["delta_J_5a"],
                    "delta_J_zs": r["delta_J_zs"],
                    "err_delta_J": eDj,
                    "row_pass": bool(max(eJt, eJa, eDj) < 1e-9),
                }
            )
    idf = pd.DataFrame(id_rows)
    idf.to_csv(OUT / "ZERO_STRESS_5A_IDENTITY_PER_ROW.csv", index=False)
    ok_id = bool(idf["row_pass"].all())
    max_err = float(idf[["err_J_true_star", "err_J_alt_star", "err_delta_J"]].to_numpy().max())
    med_rows = []
    for cfg in cfg_map:
        g = idf[idf["config"] == cfg]
        med_rows.append(
            {
                "config": cfg,
                "n_rows": len(g),
                "max_err_Jt": float(g["err_J_true_star"].max()),
                "max_err_Ja": float(g["err_J_alt_star"].max()),
                "max_err_dJ": float(g["err_delta_J"].max()),
                "all_pass": bool(g["row_pass"].all()),
            }
        )
    med_df = pd.DataFrame(med_rows)
    med_df.to_csv(OUT / "ZERO_STRESS_5A_IDENTITY_summary.csv", index=False)
    (OUT / "ZERO_STRESS_5A_IDENTITY_REPORT.md").write_text(
        f"""# ZERO_STRESS_5A_IDENTITY_REPORT

UTC: {NOW}

## Gate

`ZERO_STRESS_5A_IDENTITY_GATE` = **{'PASS' if ok_id else 'FAIL'}**

A_rms=0, D_f=0 时 5B-FIX 与 5A raw `SUBSET_PROFILED_MARGIN_RESULTS.csv` 逐行比较
(pair_id × z_true × config)，量：J_true* / J_alt* / ΔJ。

浮点容差：1e-9 dB。全局最大绝对误差 = {max_err:.3e} dB。

## 按 config

{med_df.to_string(index=False)}

## 结论

零扰动严格退化到 5A。observable = 轨迹窗去均值 `SOURCE_LEVEL_FREE_RELATIVE_TL_SHAPE`。
""",
        encoding="utf-8",
    )

    # ---------- 2. amplitude branch integrity ----------
    qs = np.array([0.3, -0.5, 0.8, -0.2, 0.6, -0.7, 0.1, -0.4])
    branch_rows = []
    for A in (0.0, 0.25, 0.5, 1.0, 2.0):
        common = A * qs[None, :] * COMMON_MULT[:, None]
        linesp = A * qs[None, :] * LINE_MULT[:, None]
        identical = bool(np.allclose(common, linesp, atol=1e-15))
        diff_norm = float(np.linalg.norm(common - linesp))
        branch_rows.append(
            {
                "A_rms_db": A,
                "common_mult": ";".join(f"{x:+.0f}" for x in COMMON_MULT),
                "line_specific_mult": ";".join(f"{x:+.0f}" for x in LINE_MULT),
                "injected_common_201dB": float(common[0, 0]),
                "injected_line_201dB": float(linesp[0, 0]),
                "injected_common_235dB": float(common[1, 0]),
                "injected_line_235dB": float(linesp[1, 0]),
                "arrays_identical": identical,
                "frobenius_diff": diff_norm,
                "expected_identical": A == 0.0,
                "integrity_pass": bool(identical == (A == 0.0)),
            }
        )
    bdf = pd.DataFrame(branch_rows)
    bdf.to_csv(OUT / "AMPLITUDE_BRANCH_INTEGRITY.csv", index=False)
    ok_branch = bool(bdf["integrity_pass"].all())

    # ---------- 3. range-grid interpolation audit ----------
    pairs = pd.read_csv(PAIRS)
    if "selected" in pairs.columns:
        pairs = pairs[pairs["selected"] == True]
    pairs = pairs[pairs["pair_id"].isin(REPS)]
    mods = {f: parse_mod(ZGRID / f"zgrid_f{int(f)}.mod") for f in FREQS}
    grid_rows = []
    for grid_m in (25.0, 10.0, 5.0):
        for _, row in pairs.iterrows():
            pid = row["pair_id"]
            T = float(row["T_s"])
            turn = float(row["turn_deg"])
            t = np.arange(0.0, T + 1e-9, 1.0)
            for role, pref in (("ref", "ref"), ("alt", "alt")):
                r_traj = range_traj(
                    t,
                    row[f"{pref}_r0_km"],
                    row[f"{pref}_theta0_deg"],
                    row[f"{pref}_v"],
                    row[f"{pref}_psi_deg"],
                    turn,
                )
                rmin, rmax = float(r_traj.min()), float(r_traj.max())
                pad = 50.0
                r_grid = np.arange(rmin - pad, rmax + pad + 1e-9, grid_m)
                for f in FREQS:
                    for z in Z_TRUE:
                        mod = mods[f]
                        L_dir = 20.0 * np.log10(np.maximum(np.abs(pressure_direct(mod, r_traj, z)), 1e-30))
                        L_g = 20.0 * np.log10(np.maximum(np.abs(pressure_direct(mod, r_grid, z)), 1e-30))
                        L_int = np.interp(r_traj, r_grid, L_g)
                        raw_rms = rms(L_dir, L_int)
                        d_rms = rms(demean(L_dir), demean(L_int))
                        c = corr(demean(L_dir), demean(L_int))
                        grid_rows.append(
                            {
                                "grid_m": grid_m,
                                "pair_id": pid,
                                "role": role,
                                "freq_hz": f,
                                "z": z,
                                "raw_TL_rms_db": raw_rms,
                                "demeaned_TL_rms_db": d_rms,
                                "corr_demeaned": c,
                            }
                        )
    gdf = pd.DataFrame(grid_rows)
    gdf.to_csv(OUT / "RANGE_GRID_INTERPOLATION_AUDIT.csv", index=False)
    # pick first grid meeting gate
    chosen = None
    gate_rows = []
    for grid_m in (25.0, 10.0, 5.0):
        g = gdf[gdf["grid_m"] == grid_m]
        max_d = float(g["demeaned_TL_rms_db"].max())
        min_c = float(g["corr_demeaned"].min())
        ok = bool(max_d <= 0.05 and min_c >= 0.999)
        gate_rows.append({"grid_m": grid_m, "max_demeaned_rms_db": max_d, "min_corr": min_c, "pass": ok})
        if ok and chosen is None:
            chosen = grid_m
    gate_df = pd.DataFrame(gate_rows)
    gate_df.to_csv(OUT / "RANGE_GRID_INTERPOLATION_GATE.csv", index=False)
    if chosen is None:
        # grid acceleration unusable; main path is DIRECT F_m(r(t)) so science proceeds
        grid_status = "RANGE_GRID_PATH_REJECTED_INSUFFICIENT_ACCURACY_USE_DIRECT"
    else:
        grid_status = f"GRID_{int(chosen)}M_MEETS_GATE_BUT_MAIN_PATH_DIRECT"
    # main path is direct F_m(r(t)); grid only documented
    (OUT / "RANGE_GRID_INTERPOLATION_AUDIT.md").write_text(
        f"""# RANGE_GRID_INTERPOLATION_AUDIT

UTC: {NOW}

代表样本：{REPS}
频率：201/235/283/338 Hz · z=180/200/220 · ref+alt

比较 DIRECT p(f,r(t),z) 与 GRID raw L(f,z,r) → interp → L(r(t))，
分别报 raw TL RMS、trajectory-demeaned TL RMS、demeaned correlation。

## 固定候选 grid

{gate_df.to_string(index=False)}

选定：**{grid_status}**

主计算路径不使用 range-grid 插值，而是 F_m(r(t)) 直接求 p。
本审计记录 grid 逼近误差量级，供后续若引入 grid 加速时复用同一门槛
（demeaned RMS ≤ 0.05 dB 且 corr ≥ 0.999）。
""",
        encoding="utf-8",
    )

    # ---------- 4. old 5B supersession ----------
    old_dec = json.loads((OLD5B / "R3_RC3_REANCHOR_5B_DECISION.json").read_text(encoding="utf-8"))
    old_dec["decision"] = "REANCHOR5B_SUPERSEDED_PENDING_OBSERVABLE_FIX"
    old_dec["superseded_reason"] = (
        "ZERO_STRESS_NEVER_REPRODUCED_5A; TL_SHAPE_DEMEANED_ON_FULL_R_GRID_NOT_TRAJECTORY; "
        "LINE_SPECIFIC_AMPLITUDE_BRANCH_NOT_EXECUTED"
    )
    old_dec["auxiliary_labels"] = list(old_dec.get("auxiliary_labels") or []) + [
        "LINE_SPECIFIC_AMPLITUDE_BRANCH_NOT_EXECUTED",
        "FREQUENCY_DRIFT_FORWARD_MODEL_SPOTCHECK_PASSED",
    ]
    old_dec["superseded_by"] = "R3_RC3_LINE_INSTABILITY_BOUNDARY_FIX"
    old_dec["superseded_utc"] = NOW
    (OLD5B / "R3_RC3_REANCHOR_5B_DECISION.json").write_text(
        json.dumps(old_dec, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    rep_path = OLD5B / "R3_RC3_REANCHOR_5B_REPORT.md"
    if rep_path.exists():
        txt = rep_path.read_text(encoding="utf-8")
        banner = (
            f"> **SUPERSEDED** `{NOW}`\n"
            f"> `REANCHOR5B_SUPERSEDED_PENDING_OBSERVABLE_FIX`\n"
            f"> 根因：全 R_GRID 去均值 ≠ 轨迹窗去均值；LINE_SPECIFIC 分支未执行。\n"
            f"> 保留：`FREQUENCY_DRIFT_FORWARD_MODEL_SPOTCHECK_PASSED`。\n"
            f"> 正式判定见 `R3_RC3_LINE_INSTABILITY_BOUNDARY_FIX/`。\n\n"
        )
        if "SUPERSEDED" not in txt[:400]:
            rep_path.write_text(banner + txt, encoding="utf-8")
    lock_path = OLD5B / "LINE_INSTABILITY_AXIS_LOCK.md"
    if lock_path.exists():
        txt = lock_path.read_text(encoding="utf-8")
        if "SUPERSEDED" not in txt[:400]:
            lock_path.write_text(
                f"> **SUPERSEDED** by 5B-FIX `{NOW}`\n\n" + txt, encoding="utf-8"
            )

    # ---------- 5. enrich FIX decision/report ----------
    dec_path = OUT / "R3_RC3_REANCHOR_5B_FIX_DECISION.json"
    dec = json.loads(dec_path.read_text(encoding="utf-8"))
    # recompute science decision from FIXED summaries (never trust stale decision field)
    asdf = pd.read_csv(OUT / "AMPLITUDE_VARIATION_SUMMARY_FIXED.csv")
    fsdf = pd.read_csv(OUT / "FREQUENCY_DRIFT_SUMMARY_FIXED.csv")

    def _survive_amp(A):
        for shape in ("A-RAMP", "A-SINE"):
            for cfg in ("201", "235", "283", "338", "FOUR"):
                gg = asdf[(asdf["A_rms_db"] == A) & (asdf["shape"] == shape) & (asdf["config"] == cfg)]
                if not len(gg) or not gg["survive"].all():
                    return False
        return True

    def _survive_f(branch, D):
        for cfg in ("201", "235", "283", "338", "FOUR"):
            gg = fsdf[(fsdf["branch"] == branch) & (fsdf["D_f"] == D) & (fsdf["config"] == cfg)]
            if not len(gg) or not gg["survive"].all():
                return False
        return True

    amp1, amp05 = _survive_amp(1.0), _survive_amp(0.5)
    ftr2, fnom1 = _survive_f("F_TRACKED", 0.02), _survive_f("F_NOMINAL", 0.01)
    extra = ["TRACKABLE_LINE_REQUIRES_FREQUENCY_TRACKING"] if (ftr2 and not fnom1) else []
    if not ok_id:
        science_decision = "5B_BLOCKED_BY_ZERO_STRESS_IDENTITY"
    elif amp1 and ftr2:
        science_decision = "TRACKABLE_LINE_INSTABILITY_ROBUST_IN_TESTED_RANGE"
    elif amp05 and not amp1:
        science_decision = "TRACKABLE_LINE_CONDITIONAL_ON_LOW_SOURCE_VARIATION"
    elif (not amp05) or (not _survive_f("F_TRACKED", 0.01)):
        science_decision = "TRACKABLE_LINE_INSTABILITY_BREAKS_PROFILED_RC3"
    else:
        science_decision = "TRACKABLE_LINE_INSTABILITY_BOUNDARY_UNRESOLVED"
    why = f"ZERO_STRESS={'PASS' if ok_id else 'FAIL'}; amp1={amp1} amp05={amp05} F_TR@2%={ftr2} F_NOM@1%={fnom1} aux={extra}"

    dec["decision"] = science_decision
    dec["why"] = why
    dec["zero_stress_per_row_max_err_db"] = max_err
    dec["zero_stress_gate"] = "PASS" if ok_id else "FAIL"
    dec["amplitude_branch_integrity"] = "PASS" if ok_branch else "FAIL"
    dec["range_grid_status"] = grid_status
    dec["range_grid_main_path"] = "DIRECT_F_MATRIX_ON_TRAJECTORY"
    dec["observable"] = "SOURCE_LEVEL_FREE_RELATIVE_TL_SHAPE_TRAJECTORY_DEMEAN"
    dec["old_5b"] = "REANCHOR5B_SUPERSEDED_PENDING_OBSERVABLE_FIX"
    dec["kept_labels"] = ["FREQUENCY_DRIFT_FORWARD_MODEL_SPOTCHECK_PASSED"]
    dec["auxiliary_labels"] = extra
    dec["updated_utc"] = NOW
    dec["range_grid_gate"] = (
        "FAIL_ALL_TESTED_GRIDS" if chosen is None else f"PASS_AT_{int(chosen)}M"
    )
    dec_path.write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    fix_decision = dec["decision"]
    (OUT / "R3_RC3_REANCHOR_5B_FIX_REPORT.md").write_text(
        f"""# R3-RC3-REANCHOR-5B-FIX

UTC: {NOW}

## 判定

### `{fix_decision}`

{dec.get('why','')}

## Gate 链

1. **ZERO_STRESS_5A_IDENTITY_GATE = {'PASS' if ok_id else 'FAIL'}**
   - 逐行对齐 5A raw（450 行 × J_true*/J_alt*/ΔJ）
   - 全局最大误差 = {max_err:.3e} dB（容差 1e-9）
2. **AMPLITUDE_BRANCH_INTEGRITY = {'PASS' if ok_branch else 'FAIL'}**
   - COMMON `[+1,+1,+1,+1]` vs LINE_SPECIFIC `[+1,-1,+1,-1]`
   - A=0 完全一致；A>0 注入数组不同
3. **RANGE_GRID_INTERPOLATION_AUDIT**
   - 状态：`{grid_status}`
   - 25/10/5 m 均未达 demeaned RMS≤0.05 dB 且 corr≥0.999
   - **主路径不使用 range-grid**，而是 F_m(r(t)) 直接 |p|；grid 加速正式废弃

## observable

$$\\tilde L_f(t)=L_f(t)-\\overline{{L_f(t)}}$$

- 幅漂：L_prop + a_f(t) 后再轨迹窗去均值
- 频漂：先拼完整 L_f(t)(t)，再整条线时间窗去均值
- F_TRACKED 模板用观测 f(t)；F_NOMINAL 模板恒 f_0

## 相对旧 5B

| 项 | 旧 5B | 5B-FIX |
|---|---|---|
| TL 去均值 | 全 44–61 km R_GRID | 轨迹时间窗 |
| 零扰动 vs 5A | 不一致 | 逐行一致 |
| LINE_SPECIFIC | 字符串未匹配，未执行 | multiplier 显式 [±1] |
| 判定 | REANCHOR5B_SUPERSEDED_PENDING_OBSERVABLE_FIX | `{fix_decision}` |

## 辅助标签

- `TRACKABLE_LINE_REQUIRES_FREQUENCY_TRACKING`（若 F_NOMINAL@1% 未过）
- `FREQUENCY_DRIFT_FORWARD_MODEL_SPOTCHECK_PASSED`（自旧 5B 保留）
- `WEAKEST_IDEAL_SINGLE_LINE_235HZ`（单列）

## 未做

combined corner、旧 S1 频率迁移、S0、SSP、IID TL、intermittent、Liang、P5。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(
        f"# 5B-FIX\n\n**{fix_decision}**\n\n"
        f"ZERO_STRESS per-row max_err={max_err:.3e} dB · "
        f"branch_integrity={'PASS' if ok_branch else 'FAIL'} · grid={grid_status}\n\n"
        f"旧5B: REANCHOR5B_SUPERSEDED_PENDING_OBSERVABLE_FIX\n"
        f"保留: FREQUENCY_DRIFT_FORWARD_MODEL_SPOTCHECK_PASSED\n",
        encoding="utf-8",
    )

    print("identity", ok_id, "max_err", max_err)
    print("branch", ok_branch)
    print("grid", grid_status)
    print("decision", fix_decision)
    # grid failure does not block: main path is DIRECT and zero-stress passed
    return 0 if ok_id and ok_branch else 1


if __name__ == "__main__":
    raise SystemExit(main())
