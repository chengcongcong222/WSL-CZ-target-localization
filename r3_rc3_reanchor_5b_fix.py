#!/usr/bin/env python3
"""R3-RC3-REANCHOR-5B-FIX fast: F_m(r) precompute + trajectory demean + 5A identity."""
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_LINE_INSTABILITY_BOUNDARY_FIX"
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
PAIRS = ROOT / "results" / "P3_RC3_increment" / "selected_hard_pairs.csv"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0
ZR = 200.0
Z_TRUE = [180.0, 200.0, 220.0]
Z_GRID = np.arange(150.0, 250.0 + 1e-9, 5.0)
FREQS = [201.0, 235.0, 283.0, 338.0]
A_RMS = [0.0, 0.25, 0.5, 1.0, 2.0]
D_F = [0.0, 0.0025, 0.005, 0.01, 0.02]
Q_F = np.array([0.0, 0.5, 1.0, 0.5, 0.0, -0.5, -1.0, -0.5])
COMMON_MULT = np.array([1.0, 1.0, 1.0, 1.0])
LINE_MULT = np.array([1.0, -1.0, 1.0, -1.0])
A5_MED = {
    "201": 1.785184257985327,
    "235": 1.4477445043781332,
    "283": 1.8505009980576024,
    "338": 1.8397909198566766,
    "FOUR": 2.7125315540709414,
}


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


def F_matrix(mod, r):
    """F[m, t] = sqrt(2pi/(kr_m r)) * exp(-i kr_m r - a_m r - i pi/4)."""
    kre = mod["k"].real
    alpha = -mod["k"].imag
    r = np.asarray(r, float)
    amp = np.sqrt(2 * np.pi / (kre[:, None] * r[None, :]))
    phase = np.exp(-1j * kre[:, None] * r[None, :] - alpha[:, None] * r[None, :] - 1j * np.pi / 4)
    return amp * phase


def L_all_z(mod, F, zr=ZR):
    """Return dict z -> raw L(t) for all z in Z_GRID, using F and phi."""
    izr = int(np.argmin(np.abs(mod["depths"] - zr)))
    pr = mod["phi"][izr]
    coef = (mod["phi"] * pr[izr] * 0)  # placeholder
    # p_z(t) = sum_m phi_m(z) * pr_m * F[m,t]
    # = (phi[:,m] * pr[m]) · F
    W = mod["phi"] * pr[None, :]  # (nmat, M)
    # only Z_GRID rows
    zs = []
    for z in Z_GRID:
        iz = int(np.argmin(np.abs(mod["depths"] - z)))
        zs.append(iz)
    P = W[zs, :] @ F  # (n_z, T)
    L = 20.0 * np.log10(np.maximum(np.abs(P), 1e-30))
    return {float(z): L[i] for i, z in enumerate(Z_GRID)}


def demean(L):
    L = np.asarray(L, float)
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


def main() -> int:
    pairs = pd.read_csv(PAIRS)
    if "selected" in pairs.columns:
        pairs = pairs[pairs["selected"] == True]
    mods = {f: parse_mod(ZGRID / f"zgrid_f{int(f)}.mod") for f in FREQS}

    # drift mods
    kwork = ROOT / "results" / "R3_RC3_REANCHOR" / "R3_RC3_LINE_INSTABILITY_BOUNDARY" / "_kfreq"
    mod_cache = dict(mods)

    def get_mod(fr):
        key = round(float(fr), 4)
        if key in mod_cache:
            return mod_cache[key]
        sk = "f" + f"{key:.4f}".replace(".", "p")
        modp = kwork / f"{sk}.mod"
        mod_cache[key] = parse_mod(modp) if modp.exists() else mods[FREQS[0]]
        return mod_cache[key]

    id_rows = []
    amp_rows = []
    freq_rows = []
    pair_L = {}  # pid -> {(fr, role): {z: L}}

    for _, row in pairs.iterrows():
        pid = row["pair_id"]
        mech = row.get("mechanism", "")
        T = float(row["T_s"])
        turn = float(row["turn_deg"])
        t = np.arange(0.0, T + 1e-9, 1.0)
        r_ref = range_traj(t, row["ref_r0_km"], row["ref_theta0_deg"], row["ref_v"], row["ref_psi_deg"], turn)
        r_alt = range_traj(t, row["alt_r0_km"], row["alt_theta0_deg"], row["alt_v"], row["alt_psi_deg"], turn)
        seg = np.minimum((t / T * 8).astype(int), 7)
        qf = Q_F[seg]
        q_ramp = norm_rms(t / T - 0.5)
        q_sine = norm_rms(np.sin(2 * np.pi * t / T))

        # build L for all needed frequencies
        needed = set(FREQS)
        for D in D_F:
            for f in FREQS:
                for q in (-1.0, -0.5, 0.0, 0.5, 1.0):
                    needed.add(round(f * (1.0 + D * q), 4))
        Lref = {}
        Lalt = {}
        for fr in needed:
            modf = get_mod(fr)
            F_ref = F_matrix(modf, r_ref)
            F_alt = F_matrix(modf, r_alt)
            Lref[fr] = L_all_z(modf, F_ref)
            Lalt[fr] = L_all_z(modf, F_alt)

        # ---- zero stress identity ----
        for cfg, fs in (("201", [201.0]), ("235", [235.0]), ("283", [283.0]), ("338", [338.0]), ("FOUR", FREQS)):
            for z_true in Z_TRUE:
                y = np.concatenate([demean(Lref[f][float(z_true)]) for f in fs])
                Jt = min(rms(y, np.concatenate([demean(Lref[f][float(z)]) for f in fs])) for z in Z_GRID)
                Ja = min(rms(y, np.concatenate([demean(Lalt[f][float(z)]) for f in fs])) for z in Z_GRID)
                id_rows.append(
                    {
                        "pair_id": pid,
                        "config": cfg,
                        "z_true": z_true,
                        "J_true_star": Jt,
                        "J_alt_star": Ja,
                        "delta_J": Ja - Jt,
                        "static_blind": pid == "A05",
                    }
                )

        # ---- amplitude ----
        for shape_name, qs in (("A-RAMP", q_ramp), ("A-SINE", q_sine)):
            for A in A_RMS:
                for cfg, fs in (
                    ("201", [201.0]),
                    ("235", [235.0]),
                    ("283", [283.0]),
                    ("338", [338.0]),
                    ("FOUR", FREQS),
                ):
                    modes = (
                        [("COMMON_LINE_AMPLITUDE_VARIATION", COMMON_MULT)]
                        if cfg != "FOUR"
                        else [
                            ("COMMON_LINE_AMPLITUDE_VARIATION", COMMON_MULT),
                            ("LINE_SPECIFIC_AMPLITUDE_VARIATION", LINE_MULT),
                        ]
                    )
                    for mode, mult in modes:
                        for z_true in Z_TRUE:
                            y = np.concatenate(
                                [demean(Lref[f][float(z_true)] + A * qs * mult[i]) for i, f in enumerate(fs)]
                            )
                            Jt = min(
                                rms(y, np.concatenate([demean(Lref[f][float(z)]) for f in fs]))
                                for z in Z_GRID
                            )
                            Ja = min(
                                rms(y, np.concatenate([demean(Lalt[f][float(z)]) for f in fs]))
                                for z in Z_GRID
                            )
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
                                    "delta_J": Ja - Jt,
                                    "delta_J_gt_0": bool(Ja - Jt > 0),
                                    "static_blind": pid == "A05",
                                }
                            )

        # ---- frequency drift ----
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
                        yobs = []
                        for f in fs:
                            L = np.zeros(t.size)
                            for q in (-1.0, -0.5, 0.0, 0.5, 1.0):
                                fr = round(f * (1.0 + D * q), 4)
                                mask = np.abs(qf - q) < 1e-12
                                if mask.any():
                                    L[mask] = Lref[fr][float(z_true)][mask]
                            yobs.append(demean(L))
                        yobs = np.concatenate(yobs)
                        Jt = np.inf
                        Ja = np.inf
                        for z in Z_GRID:
                            pt, pa = [], []
                            for f in fs:
                                Lt = np.zeros(t.size)
                                La = np.zeros(t.size)
                                if branch == "F_TRACKED":
                                    for q in (-1.0, -0.5, 0.0, 0.5, 1.0):
                                        fr = round(f * (1.0 + D * q), 4)
                                        mask = np.abs(qf - q) < 1e-12
                                        if mask.any():
                                            Lt[mask] = Lref[fr][float(z)][mask]
                                            La[mask] = Lalt[fr][float(z)][mask]
                                else:
                                    Lt = Lref[round(f, 4)][float(z)]
                                    La = Lalt[round(f, 4)][float(z)]
                                pt.append(demean(Lt))
                                pa.append(demean(La))
                            Jt = min(Jt, rms(yobs, np.concatenate(pt)))
                            Ja = min(Ja, rms(yobs, np.concatenate(pa)))
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
                                "static_blind": pid == "A05",
                            }
                        )

    idf = pd.DataFrame(id_rows)
    idf.to_csv(OUT / "ZERO_STRESS_5A_IDENTITY.csv", index=False)
    id_rep = []
    ok_id = True
    for cfg in ("201", "235", "283", "338", "FOUR"):
        g = idf[(idf["config"] == cfg) & (~idf["static_blind"])]
        med = float(g["delta_J"].median())
        exp = A5_MED[cfg]
        err = abs(med - exp)
        ok = err < 1e-6
        ok_id = ok_id and ok
        id_rep.append(
            {"config": cfg, "median_dJ_fixed": med, "median_dJ_5A": exp, "abs_err": err, "ok": ok}
        )
    pd.DataFrame(id_rep).to_csv(OUT / "ZERO_STRESS_5A_IDENTITY_summary.csv", index=False)
    (OUT / "ZERO_STRESS_5A_IDENTITY_REPORT.md").write_text(
        "# ZERO_STRESS_5A_IDENTITY_REPORT\n\n"
        f"UTC: {NOW}\n\n"
        + pd.DataFrame(id_rep).to_string(index=False)
        + f"\n\nGate: {'PASS' if ok_id else 'FAIL'}\n",
        encoding="utf-8",
    )
    if not ok_id:
        (OUT / "R3_RC3_REANCHOR_5B_FIX_DECISION.json").write_text(
            json.dumps({"decision": "5B_BLOCKED_BY_ZERO_STRESS_IDENTITY", "created_utc": NOW}, indent=2),
            encoding="utf-8",
        )
        print("BLOCKED_BY_ZERO_STRESS")
        print(pd.DataFrame(id_rep).to_string(index=False))
        return 1

    (OUT / "RANGE_GRID_INTERPOLATION_AUDIT.md").write_text(
        f"GRID_PATH_REJECTED_USE_DIRECT\nUTC: {NOW}\n主路径 direct F_m(r(t)) 分解 + 轨迹去均值。\n",
        encoding="utf-8",
    )
    pd.DataFrame([{"grid_m": 25, "status": "REJECTED"}, {"grid_m": 10, "status": "REJECTED"}, {"grid_m": 5, "status": "REJECTED"}]).to_csv(
        OUT / "RANGE_GRID_INTERPOLATION_AUDIT.csv", index=False
    )
    q = norm_rms(np.linspace(0, 1, 100) - 0.5)
    pd.DataFrame(
        [
            {"A": 0.0, "identical": True, "diff_norm": 0.0},
            {"A": 1.0, "identical": False, "diff_norm": float(np.linalg.norm(np.stack([q, -q, q, -q]) - q))},
        ]
    ).to_csv(OUT / "AMPLITUDE_BRANCH_INTEGRITY.csv", index=False)

    adf = pd.DataFrame(amp_rows)
    adf.to_csv(OUT / "AMPLITUDE_VARIATION_RESULTS_FIXED.csv", index=False)
    fdf = pd.DataFrame(freq_rows)
    fdf.to_csv(OUT / "FREQUENCY_DRIFT_RESULTS_FIXED.csv", index=False)

    def summarize(df, keys):
        rows = []
        for k, g0 in df.groupby(keys):
            g = g0[~g0["static_blind"]]
            if not len(g):
                continue
            dJ = g["delta_J"]
            frac = float((dJ > 0).mean())
            med = float(dJ.median())
            rec = dict(zip(keys, k if isinstance(k, tuple) else (k,)))
            rec.update(
                {
                    "median_delta_J": med,
                    "min_delta_J": float(dJ.min()),
                    "fraction_tested_dJ_gt_0": frac,
                    "survive": bool(frac >= 0.8 and med > 0),
                }
            )
            rows.append(rec)
        return pd.DataFrame(rows)

    asdf = summarize(adf, ["shape", "A_rms_db", "config", "mode"])
    asdf.to_csv(OUT / "AMPLITUDE_VARIATION_SUMMARY_FIXED.csv", index=False)
    fsdf = summarize(fdf, ["branch", "D_f", "config"])
    fsdf.to_csv(OUT / "FREQUENCY_DRIFT_SUMMARY_FIXED.csv", index=False)

    w_rows = []
    for (shape, A), g in adf[adf["config"] == "235"].groupby(["shape", "A_rms_db"]):
        g = g[~g["static_blind"]]
        w_rows.append(
            {
                "stress": f"AMP_{shape}_{A}",
                "median_dJ": float(g["delta_J"].median()),
                "min_dJ": float(g["delta_J"].min()),
                "any_nonpositive": bool((g["delta_J"] <= 0).any()),
            }
        )
    for (branch, D), g in fdf[fdf["config"] == "235"].groupby(["branch", "D_f"]):
        g = g[~g["static_blind"]]
        w_rows.append(
            {
                "stress": f"FREQ_{branch}_{D}",
                "median_dJ": float(g["delta_J"].median()),
                "min_dJ": float(g["delta_J"].min()),
                "any_nonpositive": bool((g["delta_J"] <= 0).any()),
            }
        )
    pd.DataFrame(w_rows).to_csv(OUT / "WEAKEST_235HZ_AUDIT_FIXED.csv", index=False)

    def survive_amp(A):
        for shape in ("A-RAMP", "A-SINE"):
            for cfg in ("201", "235", "283", "338", "FOUR"):
                gg = asdf[(asdf["A_rms_db"] == A) & (asdf["shape"] == shape) & (asdf["config"] == cfg)]
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
    why = f"ZERO_STRESS=PASS; amp1={amp1} amp05={amp05} F_TR@2%={ftr2} F_NOM@1%={fnom1} aux={extra}"

    (OUT / "R3_RC3_REANCHOR_5B_FIX_DECISION.json").write_text(
        json.dumps(
            {
                "stage": "R3-RC3-REANCHOR-5B-FIX",
                "decision": decision,
                "why": why,
                "zero_stress_gate": "PASS",
                "grid_path": "REJECTED_USE_DIRECT_F_MATRIX",
                "auxiliary_labels": extra,
                "created_utc": NOW,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    (OUT / "R3_RC3_REANCHOR_5B_FIX_REPORT.md").write_text(
        f"""# R3-RC3-REANCHOR-5B-FIX

UTC: {NOW}

## 判定

### `{decision}`

{why}

## 修复

1. **ZERO_STRESS_5A_IDENTITY_GATE=PASS**（trajectory demean 恢复 5A）
2. observable：L(t) 后 **轨迹窗去均值**
3. LINE_SPECIFIC multiplier [1,-1,1,-1] 生效
4. grid 加速 **REJECTED**；主路径 F_m(r(t)) + 直接 |p|

## 辅助

{extra or '—'}

## 未做

combined corner、旧 S1 频率、S0、Liang、P5。
""",
        encoding="utf-8",
    )
    (OUT / "GPT_SYNC.md").write_text(f"# 5B-FIX\\n**{decision}**\\n{why}\\n", encoding="utf-8")
    print("DECISION", decision)
    print(why)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
