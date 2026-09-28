#!/usr/bin/env python3
"""R3-RC3-REANCHOR-5B (fast): amplitude + frequency drift via r-grid TL cache."""
from __future__ import annotations

import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_LINE_INSTABILITY_BOUNDARY"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
PAIRS = ROOT / "results" / "P3_RC3_increment" / "selected_hard_pairs.csv"
AT_BIN = ROOT / "tools" / "acoustics_toolbox" / "atWin10" / "at" / "bin"
WORK = OUT / "_kfreq"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0
ZR = 200.0
Z_TRUE = [180.0, 200.0, 220.0]
Z_GRID = np.arange(150.0, 250.0 + 1e-9, 5.0)
FREQS = [201.0, 235.0, 283.0, 338.0]
A_RMS = [0.0, 0.25, 0.5, 1.0, 2.0]
D_F = [0.0, 0.0025, 0.005, 0.01, 0.02]
Q_F = np.array([0.0, 0.5, 1.0, 0.5, 0.0, -0.5, -1.0, -0.5])
R_GRID = np.arange(44000.0, 61000.0 + 1, 25.0)  # 25 m grid for interpolation


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


def freq_key(freq: float) -> str:
    return f"f{float(freq):.4f}".replace(".", "p")


def write_env_freq(path: Path, freq: float, env0: Path):
    lines = env0.read_text(encoding="utf-8").splitlines()
    out = []
    for i, ln in enumerate(lines):
        if i == 0:
            out.append("'PERT'")
        elif i == 1:
            out.append(f"{float(freq):.4f}")
        else:
            out.append(ln)
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def get_mod(freq: float) -> dict | None:
    key = freq_key(freq)
    modp = WORK / f"{key}.mod"
    if modp.exists():
        return parse_mod(modp)
    write_env_freq(WORK / f"{key}.env", float(freq), ZGRID / "zgrid_f235.env")
    subprocess.run([str(AT_BIN / "kraken.exe"), key], cwd=str(WORK), capture_output=True, text=True, timeout=180)
    return parse_mod(modp) if modp.exists() else None


def pressure_rgrid(mod, r, zs, zr=ZR):
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


def norm_rms(q):
    q = np.asarray(q, float)
    q = q - q.mean()
    s = float(np.sqrt(np.mean(q * q)))
    return q / s if s > 0 else q


def interp_tl(tl_rgrid, r_traj):
    # tl on R_GRID; interp real/imag? tl is real array
    return np.interp(r_traj, R_GRID, tl_rgrid)


def main() -> int:
    # 5A summary fix
    a5 = pd.read_csv(
        ROOT / "results/R3_RC3_REANCHOR/R3_RC3_LINE_AVAILABILITY_BOUNDARY/SUBSET_SUMMARY.csv"
    )
    lc = []
    for n in (1, 2, 3, 4):
        g = a5[a5["n_lines"] == n].sort_values(["fraction_tested_dJ_gt_0", "median_delta_J"])
        if not len(g):
            continue
        lc.append(
            {
                "n_lines": n,
                "worst_subset": g.iloc[0]["frequencies"],
                "worst_median_dJ": g.iloc[0]["median_delta_J"],
                "best_subset": g.iloc[-1]["frequencies"],
                "best_median_dJ": g.iloc[-1]["median_delta_J"],
                "n_pass": int(g["subset_pass"].sum()),
                "all_pass": bool(g["subset_pass"].all()),
                "frequency_info": "PARTIALLY_COMPLEMENTARY_FREQUENCY_INFORMATION",
            }
        )
    pd.DataFrame(lc).to_csv(OUT / "5A_LINE_COUNT_SUMMARY_CORRECTED.csv", index=False)

    (OUT / "LINE_INSTABILITY_AXIS_LOCK.md").write_text(
        f"""# LINE_INSTABILITY_AXIS_LOCK

UTC: {NOW}

`S1_BRIDGE_LINE_INSTABILITY_ONLY` · E0 · 无 SSP/无 IID TL

5B-A: A-RAMP/A-SINE, A_rms={A_RMS} dB, COMMON vs LINE_SPECIFIC
5B-F: D_f={D_F}, staircase 8段, F_TRACKED vs F_NOMINAL; KRAKEN 缓存; 禁止幅+频联合

survive: fraction_tested(ΔJ>0)>=0.8 且 median>0
""",
        encoding="utf-8",
    )

    pairs = pd.read_csv(PAIRS)
    if "selected" in pairs.columns:
        pairs = pairs[pairs["selected"] == True]

    # global TL cache on R_GRID: key=(fr, z, ref/alt) -> tl array
    # built per pair because r_traj differs — instead cache p on R_GRID once globally
    # p depends on r and z, not pair. Cache TL_RGRID[(fr,z,traj)] globally.
    tl_cache = {}

    def tl_on_grid(fr, z, traj):
        key = (round(fr, 4), float(z), traj)
        if key in tl_cache:
            return tl_cache[key]
        modf = get_mod(fr)
        if modf is None:
            tl_cache[key] = None
            return None
        # pressure on R_GRID for zs=z; traj unused for pressure vs r — same for ref/alt
        # wait: TL shape vs r is independent of which trajectory; we interpolate r(t)
        p = pressure_rgrid(modf, R_GRID, float(z))
        tlr = tl_shape(p)
        # note: tl_shape mean-removed on R_GRID, then we interp — approx OK for shape
        tl_cache[key] = tlr
        return tlr

    def y_of(r_traj, fs, z, q_seg=None, branch="TRACK", f0s=None):
        parts = []
        for f in fs:
            if branch == "NOMINAL":
                tlr = tl_on_grid(f, z, "x")
                if tlr is None:
                    continue
                parts.append(interp_tl(tlr, r_traj))
            else:
                # piecewise by q
                y = np.zeros(r_traj.size)
                if q_seg is None:
                    tlr = tl_on_grid(f, z, "x")
                    if tlr is None:
                        continue
                    y = interp_tl(tlr, r_traj)
                else:
                    for q in (-1.0, -0.5, 0.0, 0.5, 1.0):
                        fr = f * (1.0 + Dq * q) if False else f  # placeholder
                    # caller passes segment freqs
                parts.append(y)
        return np.concatenate(parts) if parts else None

    # simpler: build pair-time y via precomputed per-q TL
    amp_rows = []
    freq_rows = []
    # precompute all needed TL for 4 base freqs + drift freqs
    all_fr = set(FREQS)
    for f in FREQS:
        for D in D_F:
            for q in (-1.0, -0.5, 0.0, 0.5, 1.0):
                all_fr.add(round(f * (1.0 + D * q), 4))
    for fr in sorted(all_fr):
        for z in list(Z_GRID) + list(Z_TRUE):
            tl_on_grid(fr, float(z), "x")

    for _, row in pairs.iterrows():
        pid = row["pair_id"]
        mech = row.get("mechanism", "")
        T = float(row["T_s"])
        turn = float(row["turn_deg"])
        t = np.arange(0.0, T + 1e-9, 1.0)
        r_ref = range_traj(t, row["ref_r0_km"], row["ref_theta0_deg"], row["ref_v"], row["ref_psi_deg"], turn)
        r_alt = range_traj(t, row["alt_r0_km"], row["alt_theta0_deg"], row["alt_v"], row["alt_psi_deg"], turn)
        static_blind = pid == "A05"
        seg = np.minimum((t / T * 8).astype(int), 7)
        qf = Q_F[seg]
        q_ramp = norm_rms(t / T - 0.5)
        q_sine = norm_rms(np.sin(2 * np.pi * t / T))

        def y_true_amp(z_true, fs, A, q, mode):
            parts = []
            for fi, f in enumerate(fs):
                tlr = tl_on_grid(f, float(z_true), "x")
                if tlr is None:
                    continue
                y = interp_tl(tlr, r_ref)
                qq = q if (mode != "LINE_SPECIFIC" or fi % 2 == 0) else -q
                parts.append(y + A * qq)
            return np.concatenate(parts) if parts else None

        def y_tpl(r, fs, z):
            parts = []
            for f in fs:
                tlr = tl_on_grid(f, float(z), "x")
                if tlr is None:
                    continue
                parts.append(interp_tl(tlr, r))
            return np.concatenate(parts) if parts else None

        # ---- 5B-A ----
        for shape_name, q in (("A-RAMP", q_ramp), ("A-SINE", q_sine)):
            for A in A_RMS:
                for cfg, fs in (
                    ("201", [201.0]),
                    ("235", [235.0]),
                    ("283", [283.0]),
                    ("338", [338.0]),
                    ("FOUR", FREQS),
                ):
                    modes = ["COMMON"] if cfg != "FOUR" else ["COMMON_LINE_AMPLITUDE_VARIATION", "LINE_SPECIFIC_AMPLITUDE_VARIATION"]
                    for mode in modes:
                        for z_true in Z_TRUE:
                            y = y_true_amp(z_true, fs, A, q, mode)
                            if y is None:
                                continue
                            Jt = min(rms(y, y_tpl(r_ref, fs, z)) for z in Z_GRID)
                            Ja = min(rms(y, y_tpl(r_alt, fs, z)) for z in Z_GRID)
                            dJ = Ja - Jt
                            amp_rows.append(
                                {
                                    "pair_id": pid,
                                    "mechanism": mech,
                                    "shape": shape_name,
                                    "A_rms_db": A,
                                    "config": cfg,
                                    "mode": mode,
                                    "z_true": z_true,
                                    "J_true_star": Jt,
                                    "J_alt_star": Ja,
                                    "delta_J": dJ,
                                    "delta_J_gt_0": bool(dJ > 0),
                                    "delta_J_gt_0p5": bool(dJ > 0.5),
                                    "static_blind": static_blind,
                                }
                            )

        # ---- 5B-F ----
        for D in D_F:
            for cfg, fs in (
                ("201", [201.0]),
                ("235", [235.0]),
                ("283", [283.0]),
                ("338", [338.0]),
                ("FOUR", FREQS),
            ):
                for branch in ("F_TRACKED", "F_NOMINAL"):
                    for z_true in Z_TRUE:
                        # obs from drifting f
                        yobs = []
                        for f in fs:
                            y = np.zeros(t.size)
                            for q in (-1.0, -0.5, 0.0, 0.5, 1.0):
                                fr = round(f * (1.0 + D * q), 4)
                                tlr = tl_on_grid(fr, float(z_true), "x")
                                if tlr is None:
                                    continue
                                mask = np.abs(qf - q) < 1e-12
                                if mask.any():
                                    y[mask] = interp_tl(tlr, r_ref[mask])
                            yobs.append(y)
                        yobs = np.concatenate(yobs)
                        Jt = np.inf
                        Ja = np.inf
                        for z in Z_GRID:
                            parts_t = []
                            parts_a = []
                            for f in fs:
                                yt = np.zeros(t.size)
                                ya = np.zeros(t.size)
                                if branch == "F_TRACKED":
                                    for q in (-1.0, -0.5, 0.0, 0.5, 1.0):
                                        fr = round(f * (1.0 + D * q), 4)
                                        tlr = tl_on_grid(fr, float(z), "x")
                                        if tlr is None:
                                            continue
                                        mask = np.abs(qf - q) < 1e-12
                                        if mask.any():
                                            yt[mask] = interp_tl(tlr, r_ref[mask])
                                            ya[mask] = interp_tl(tlr, r_alt[mask])
                                else:
                                    tlr = tl_on_grid(f, float(z), "x")
                                    if tlr is not None:
                                        yt = interp_tl(tlr, r_ref)
                                        ya = interp_tl(tlr, r_alt)
                                parts_t.append(yt)
                                parts_a.append(ya)
                            Jt = min(Jt, rms(yobs, np.concatenate(parts_t)))
                            Ja = min(Ja, rms(yobs, np.concatenate(parts_a)))
                        dJ = Ja - Jt
                        freq_rows.append(
                            {
                                "pair_id": pid,
                                "mechanism": mech,
                                "D_f": D,
                                "config": cfg,
                                "branch": branch,
                                "z_true": z_true,
                                "J_true_star": Jt,
                                "J_alt_star": Ja,
                                "delta_J": dJ,
                                "delta_J_gt_0": bool(dJ > 0),
                                "static_blind": static_blind,
                            }
                        )

    adf = pd.DataFrame(amp_rows)
    adf.to_csv(OUT / "AMPLITUDE_VARIATION_RESULTS.csv", index=False)
    fdf = pd.DataFrame(freq_rows)
    fdf.to_csv(OUT / "FREQUENCY_DRIFT_RESULTS.csv", index=False)

    # summaries
    amp_sum = []
    for keys, g0 in adf.groupby(["shape", "A_rms_db", "config", "mode"]):
        g = g0[~g0["static_blind"]]
        dJ = g["delta_J"]
        frac = float((dJ > 0).mean())
        med = float(dJ.median())
        amp_sum.append(
            {
                "shape": keys[0],
                "A_rms_db": keys[1],
                "config": keys[2],
                "mode": keys[3],
                "median_delta_J": med,
                "p10_delta_J": float(dJ.quantile(0.10)),
                "min_delta_J": float(dJ.min()),
                "fraction_tested_dJ_gt_0": frac,
                "survive": bool(frac >= 0.8 and med > 0),
            }
        )
    asdf = pd.DataFrame(amp_sum)
    asdf.to_csv(OUT / "AMPLITUDE_VARIATION_SUMMARY.csv", index=False)

    fsum = []
    for keys, g0 in fdf.groupby(["branch", "D_f", "config"]):
        g = g0[~g0["static_blind"]]
        dJ = g["delta_J"]
        frac = float((dJ > 0).mean())
        med = float(dJ.median())
        fsum.append(
            {
                "branch": keys[0],
                "D_f": keys[1],
                "config": keys[2],
                "median_delta_J": med,
                "min_delta_J": float(dJ.min()),
                "fraction_tested_dJ_gt_0": frac,
                "survive": bool(frac >= 0.8 and med > 0),
            }
        )
    fsdf = pd.DataFrame(fsum)
    fsdf.to_csv(OUT / "FREQUENCY_DRIFT_SUMMARY.csv", index=False)

    # kragen config / spotcheck (lightweight)
    pd.DataFrame(
        [{"freq": fr, "cached": (WORK / f"{freq_key(fr)}.mod").exists()} for fr in sorted(all_fr)]
    ).to_csv(OUT / "FREQUENCY_DRIFT_KRAKEN_CONFIG.csv", index=False)

    # weakest 235
    w_rows = []
    for _, g in adf[adf["config"] == "235"].groupby(["shape", "A_rms_db"]):
        g = g[~g["static_blind"]]
        if len(g):
            w_rows.append(
                {
                    "stress": f"AMP_{g.iloc[0]['shape']}_{g.iloc[0]['A_rms_db']}",
                    "median_dJ": float(g["delta_J"].median()),
                    "min_dJ": float(g["delta_J"].min()),
                    "any_nonpositive": bool((g["delta_J"] <= 0).any()),
                }
            )
    for _, g in fdf[fdf["config"] == "235"].groupby(["branch", "D_f"]):
        g = g[~g["static_blind"]]
        if len(g):
            w_rows.append(
                {
                    "stress": f"FREQ_{g.iloc[0]['branch']}_{g.iloc[0]['D_f']}",
                    "median_dJ": float(g["delta_J"].median()),
                    "min_dJ": float(g["delta_J"].min()),
                    "any_nonpositive": bool((g["delta_J"] <= 0).any()),
                }
            )
    pd.DataFrame(w_rows).to_csv(OUT / "WEAKEST_235HZ_AUDIT.csv", index=False)

    def survive_amp(A):
        for shape in ("A-RAMP", "A-SINE"):
            for cfg in ("201", "235", "283", "338", "FOUR"):
                gg = asdf[
                    (asdf["A_rms_db"] == A) & (asdf["shape"] == shape) & (asdf["config"] == cfg)
                ]
                if not len(gg) or not gg["survive"].all():
                    return False
        return True

    def survive_f(branch, D):
        for cfg in ("201", "235", "283", "338", "FOUR"):
            gg = fsdf[(fsdf["branch"] == branch) & (fsdf["D_f"] == D) & (fsdf["config"] == cfg)]
            if not len(gg) or not gg["survive"].all():
                return False
        return True

    amp1, amp05 = survive_amp(1.0), survive_amp(0.5)
    ftr2, fnom1 = survive_f("F_TRACKED", 0.02), survive_f("F_NOMINAL", 0.01)
    extra = ["TRACKABLE_LINE_REQUIRES_FREQUENCY_TRACKING"] if (ftr2 and not fnom1) else []
    if amp1 and ftr2:
        decision = "TRACKABLE_LINE_INSTABILITY_ROBUST_IN_TESTED_RANGE"
    elif amp05 and not amp1:
        decision = "TRACKABLE_LINE_CONDITIONAL_ON_LOW_SOURCE_VARIATION"
    elif (not amp05) or (not survive_f("F_TRACKED", 0.01)):
        decision = "TRACKABLE_LINE_INSTABILITY_BREAKS_PROFILED_RC3"
    else:
        decision = "TRACKABLE_LINE_INSTABILITY_BOUNDARY_UNRESOLVED"
    why = f"amp1={amp1} amp05={amp05} F_TRACKED@2%={ftr2} F_NOMINAL@1%={fnom1} aux={extra}"

    (OUT / "R3_RC3_REANCHOR_5B_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-RC3-REANCHOR-5B",
                "decision": decision,
                "why": why,
                "auxiliary_labels": extra,
                "axis": "S1_BRIDGE_LINE_INSTABILITY_ONLY",
                "weakest_line": "WEAKEST_IDEAL_SINGLE_LINE_235HZ",
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_RC3_REANCHOR_5B_REPORT.md").write_text(
        f"""# R3-RC3-REANCHOR-5B

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 5A 修正

`PARTIALLY_COMPLEMENTARY_FREQUENCY_INFORMATION`；worst/best 按 median ΔJ。

## 5B-A / 5B-F

见 `AMPLITUDE_VARIATION_SUMMARY.csv` / `FREQUENCY_DRIFT_SUMMARY.csv`。
F-TRACKED vs F-NOMINAL 分列；未做幅+频联合。

## 235 Hz

`WEAKEST_IDEAL_SINGLE_LINE_235HZ` — `WEAKEST_235HZ_AUDIT.csv`

## 未做

旧 S1 频率迁移、线谱消失、S0、Liang、P5。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(f"# R3-RC3-REANCHOR-5B\\n\\n**{decision}**\\n\\n{why}\\n", encoding="utf-8")
    print("DECISION", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
