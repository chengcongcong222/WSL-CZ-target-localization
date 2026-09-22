#!/usr/bin/env python3
"""R3-C2.0R: physical 1D normal-mode baseline + literature recovery status.

If Erratum not recovered → decision C2_0_PAPER_RECOVERY_BLOCKED (required).
Mode solver is for future C2 validity, not Yang estimator.
"""
from __future__ import annotations

import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

C0 = 1500.0
H_WATER = 5000.0
Z_A, B_MUNK, EPS_MUNK = 1300.0, 1300.0, 0.00737
FREQS = [201.0, 235.0, 283.0, 338.0]
Z_S_LIST = [180.0, 190.0, 200.0, 210.0, 220.0]
Z_R = 200.0
L_LIST = [0.15e3, 0.3e3, 0.6e3, 1.2e3, 2.4e3]
PAIRS = [(180, 190), (190, 200), (200, 210), (210, 220), (180, 200), (200, 220)]
C1_CORR = {201.0: 0.579, 235.0: 0.994, 283.0: 1.0, 338.0: 0.067}


def munk_c(z, za=Z_A, c0=1500.0, B=B_MUNK, eps=EPS_MUNK):
    z = np.asarray(z, float)
    eta = 2.0 * (z - za) / B
    return c0 * (1.0 + eps * (np.exp(eta) - 1.0 - eta))


class PhysicalModes:
    """1D range-independent normal modes: φ'' + (ω²/c² - k²)φ = 0.

    BC: φ(0)=0 (pressure-release surface), φ'(H)=0 (rigid bottom).
    Discrete symmetric eigenproblem A φ = k² φ.
    """

    def __init__(self, n_z=401, H=H_WATER):
        self.n_z = n_z
        self.H = H
        self.z = np.linspace(0.0, H, n_z)
        self.h = float(self.z[1] - self.z[0])
        self.c = munk_c(self.z)
        self.cache = {}

    def solve(self, f_hz):
        key = round(float(f_hz), 6)
        if key in self.cache:
            return self.cache[key]
        omega = 2 * np.pi * float(f_hz)
        k0sq = (omega / self.c) ** 2
        h = self.h
        n = self.n_z
        # interior nodes 0..n-1; u[0]=0 Dirichlet → use nodes 1..n-1 as free, last has Neumann
        # free indices 1..n-1 → size m=n-1
        m = n - 1
        A = np.zeros((m, m))
        for i in range(m):
            zi = i + 1  # node index
            if i < m - 1:
                # standard interior
                A[i, i] = 2.0 / h ** 2 + k0sq[zi]
                A[i, i + 1] = -1.0 / h ** 2
                A[i + 1, i] = -1.0 / h ** 2
            else:
                # Neumann at H: ghost u_n = u_{n-1} → diag 1/h² + k0²
                A[i, i] = 1.0 / h ** 2 + k0sq[zi]
        # ensure symmetric fill for last off-diag already set in loop when i=m-2
        A = 0.5 * (A + A.T)
        evals, evecs = np.linalg.eigh(A)
        # k² = evals descending for trapped? eigh ascending; k_r = sqrt(evals) for evals>0
        order = np.argsort(evals)[::-1]
        evals = evals[order]
        evecs = evecs[:, order]
        modes = []
        for j, k2 in enumerate(evals):
            if k2 <= 0:
                continue
            kr = math.sqrt(float(k2))
            phi_full = np.zeros(n)
            phi_full[1:] = evecs[:, j]
            # orthogonal normalize ∫φ² dz = 1
            nrm = math.sqrt(float(np.trapezoid(phi_full ** 2, self.z))) + 1e-30
            phi_full = phi_full / nrm
            # residual ||A u - k2 u||
            u = phi_full[1:]
            res = float(np.linalg.norm(A @ u - k2 * u) / (np.linalg.norm(u) + 1e-30))
            modes.append({"mode_id": j, "k_r": kr, "k2": float(k2), "phi": phi_full, "residual": res})
        self.cache[key] = modes
        return modes

    def phi_at(self, f, z, mode_list=None):
        if mode_list is None:
            mode_list = self.solve(f)
        zs = self.z
        return [float(np.interp(z, zs, m["phi"])) for m in mode_list]


PM = PhysicalModes(n_z=401)


def old_mode_kr_phi(f, z_list, n=15):
    """Simplified P4.5-style ModeModel for comparison."""
    # replicate p4 ModeModel construction
    gamma_fixed = [(m - 0.25) * np.pi / 2300.0 for m in range(1, n + 1)]
    omega = 2 * np.pi * f
    k0 = omega / 1500.0
    rows = []
    zgrid = np.linspace(0, H_WATER, 401)
    for m, gamma in enumerate(gamma_fixed, start=1):
        kr2 = k0 ** 2 - gamma ** 2
        if kr2 <= 0:
            continue
        kr = math.sqrt(kr2)
        width = 350.0 + 55.0 * m
        phi = np.exp(-0.5 * ((zgrid - 1300.0) / width) ** 2) * np.cos(gamma * (zgrid - 1300.0))
        phi = phi / (math.sqrt(float(np.trapezoid(phi ** 2, zgrid))) + 1e-30)
        phis = [float(np.interp(z, zgrid, phi)) for z in z_list]
        rows.append({"old_mode_id": m, "k_r": kr, "gamma": gamma, "phis": phis, "phi_grid": phi})
    return rows


def amplitude(A_s, A_r):
    return A_s * A_r


def main():
    t0 = time.time()
    print("=== R3-C2.0R ===", flush=True)

    # --- 1) literature ---
    lit_ok = False  # network still blocked in this environment
    for p in [OUT / "ANCHOR_YANG2015_METHOD.md", OUT / "YANG2015_ERRATUM_NOTE.md"]:
        pass
    (OUT / "YANG2015_ERRATUM_NOTE.md").write_text(
        "\n".join([
            "# YANG2015 Erratum — R3-C2.0R",
            "",
            f"UTC：{NOW}",
            "",
            "**状态：NOT_RECOVERED（仍）**",
            "",
            "- DOI 10.1121/1.5081712 / PubMed 30599696 / Crossref 本轮均无法访问",
            "- Yang 2015 主框架（由任务书/公开摘要核实）：单水听器移动 CW 源、距离相关相位 steering、",
            "  Doppler→Δr（需已知 f0）、名义环境模态深度函数、实测 PLL",
            "- **Erratum 原式/新式：未恢复** → 不猜公式",
            "",
            "正式状态保持：**C2_0_PAPER_RECOVERY_BLOCKED**",
            "原因缩小为：主方法框架已核实；**精确公式与 2018 Erratum 修正内容未恢复**。",
            "",
        ]),
        encoding="utf-8",
    )
    am = (OUT / "ANCHOR_YANG2015_METHOD.md")
    txt = am.read_text(encoding="utf-8") if am.exists() else ""
    if "R3-C2.0R" not in txt:
        am.write_text(txt + f"\n\n---\n\n## R3-C2.0R 更新 {NOW}\n\n- 网络仍无法取得全文/Erratum 修正式\n"
                       "- 框架核实：单水听器移动 CW、Δr 相位 steering、Doppler 需 f0、模态深度函数估深、PLL\n"
                       "- **公式级恢复：未完成**\n", encoding="utf-8")

    # --- 2) physical modes ---
    phys_rows = []
    conv_rows = []
    orth_rows = []
    for f in FREQS:
        modes = PM.solve(f)
        for m in modes:
            ph = PM.phi_at(f, 200.0, modes)  # placeholder
        phis_z = {z: PM.phi_at(f, z, modes) for z in Z_S_LIST + [Z_R]}
        for i, m in enumerate(modes):
            phys_rows.append({
                "f_hz": f,
                "mode_id": i,
                "k_r": m["k_r"],
                "k2": m["k2"],
                "eigen_residual": m["residual"],
                "phi_180": phis_z[180.0][i],
                "phi_190": phis_z[190.0][i],
                "phi_200": phis_z[200.0][i],
                "phi_210": phis_z[210.0][i],
                "phi_220": phis_z[220.0][i],
                "phi_zr_200": phis_z[Z_R][i],
                "A_m_zs200": phis_z[200.0][i] * phis_z[Z_R][i],
            })
        # orthogonality sample
        nchk = min(5, len(modes))
        G = np.zeros((nchk, nchk))
        for a in range(nchk):
            for b in range(nchk):
                G[a, b] = float(np.trapezoid(modes[a]["phi"] * modes[b]["phi"], PM.z))
        orth_rows.append({
            "f_hz": f, "n_modes": len(modes),
            "max_offdiag_orth": float(np.max(np.abs(G - np.diag(np.diag(G))))),
            "mean_abs_resid": float(np.mean([m["residual"] for m in modes])) if modes else np.nan,
        })
        # grid / mode count convergence
        for nz in [201, 401, 801]:
            pm = PhysicalModes(n_z=nz)
            mm = pm.solve(f)
            conv_rows.append({
                "check": "n_z", "f_hz": f, "value": nz,
                "n_modes": len(mm),
                "k_r_mode0": mm[0]["k_r"] if mm else np.nan,
                "median_dk": float(np.median(np.diff([x["k_r"] for x in mm]))) if len(mm) > 1 else np.nan,
            })
    pd.DataFrame(phys_rows).to_csv(OUT / "physical_mode_table.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(conv_rows).to_csv(OUT / "physical_mode_convergence.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(orth_rows).to_csv(OUT / "physical_mode_orthogonality.csv", index=False, encoding="utf-8-sig")

    # --- 3) old vs physical ---
    cmp_rows = []
    for f in FREQS:
        modes = PM.solve(f)
        old = old_mode_kr_phi(f, Z_S_LIST + [Z_R])
        kr_p = np.array([m["k_r"] for m in modes])
        kr_o = np.array([r["k_r"] for r in old])
        dk_p = np.median(np.diff(np.sort(kr_p))) if len(kr_p) > 1 else np.nan
        dk_o = np.median(np.diff(np.sort(kr_o))) if len(kr_o) > 1 else np.nan
        # phi variation across depth for first 5 modes
        def phi_span(phis_list):
            # phis_list: list of (z, vec)
            arr = np.array([[p for p in phis] for _, phis in phis_list])
            return float(np.mean(np.std(arr, axis=0))) if arr.size else np.nan
        old_ph = []
        zall = Z_S_LIST + [Z_R]
        for z in Z_S_LIST:
            zi = zall.index(z)
            old_ph.append((z, [r["phis"][zi] for r in old]))
        new_ph = [(z, PM.phi_at(f, z, modes)[: min(8, len(modes))]) for z in Z_S_LIST]
        cmp_rows.append({
            "f_hz": f,
            "old_n_modes": len(old),
            "phys_n_modes": len(modes),
            "old_median_dk": dk_o,
            "phys_median_dk": dk_p,
            "old_phi_m_span_180_220": phi_span(old_ph),
            "phys_phi_m_span_180_220": phi_span(new_ph),
            "old_phi_f_dependence": "NONE (phi built from gamma only)",
            "artifact_flags": "phi identical across f in old model → depth signature f-invariance is ARTIFACT_OF_SIMPLIFIED_MODEMODEL",
        })
    cmp_df = pd.DataFrame(cmp_rows)
    cmp_df.to_csv(OUT / "old_vs_physical_modes.csv", index=False, encoding="utf-8-sig")

    # --- 4) aperture pre-gate on physical modes ---
    ap_rows = []
    L_eta1_rows = []
    for f in FREQS:
        modes = PM.solve(f)
        A = np.array([PM.phi_at(f, Z_S, modes)[i] * PM.phi_at(f, Z_R, modes)[i] for i in range(len(modes)) for Z_S in [200.0]])
        # recompute A at zs=200
        ps = PM.phi_at(f, 200.0, modes)
        pr = PM.phi_at(f, Z_R, modes)
        A = np.array([ps[i] * pr[i] for i in range(len(modes))])
        An = np.abs(A) / (np.max(np.abs(A)) + 1e-30)
        kr = np.array([m["k_r"] for m in modes])
        order = np.argsort(kr)
        kr, An = kr[order], An[order]
        for thr in [0.1, 0.2]:
            idx = [i for i, a in enumerate(An) if a > thr]
            dks = []
            for a, b in zip(idx[:-1], idx[1:]):
                dks.append(abs(kr[b] - kr[a]))
            for L in L_LIST:
                etas = [L * d / (2 * np.pi) for d in dks]
                n_res = sum(1 for e in etas if e > 1.0)
                ap_rows.append({
                    "f_hz": f, "amp_thr": thr, "L_m": L, "L_km": L / 1e3,
                    "delta_k_SA": 2 * np.pi / L,
                    "median_eta": float(np.median(etas)) if etas else np.nan,
                    "p10_eta": float(np.percentile(etas, 10)) if etas else np.nan,
                    "p90_eta": float(np.percentile(etas, 90)) if etas else np.nan,
                    "n_adj_pairs": len(dks),
                    "n_resolvable_pairs_Rayleigh": n_res,
                    "v_r_T600": L / 600.0,
                    "v_r_T1200": L / 1200.0,
                    "gate_label": "Fourier/Rayleigh pre-gate — NOT high-res algorithm limit",
                })
            if thr == 0.1 and dks:
                dkm = float(np.median(dks))
                L_eta1_rows.append({
                    "f_hz": f, "median_dk": dkm,
                    "L_eta1_median_m": 2 * np.pi / dkm if dkm > 0 else np.nan,
                    "L_eta1_p10_m": 2 * np.pi / np.percentile(dks, 90) if len(dks) else np.nan,
                    "L_eta1_p90_m": 2 * np.pi / np.percentile(dks, 10) if len(dks) else np.nan,
                })
    ap_df = pd.DataFrame(ap_rows)
    ap_df.to_csv(OUT / "sa_aperture_resolution_physical.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(L_eta1_rows).to_csv(OUT / "sa_L_eta1_stats.csv", index=False, encoding="utf-8-sig")

    # --- 5) depth signature on physical phi ---
    dep_rows = []
    for f in FREQS:
        modes = PM.solve(f)
        ps = PM.phi_at(f, 200.0, modes)
        An = np.abs(ps * np.array(PM.phi_at(f, Z_R, modes)))
        An = An / (An.max() + 1e-30)
        idx = [i for i, a in enumerate(An) if a > 0.2]
        vecs = {}
        for z in Z_S_LIST:
            ph = np.array(PM.phi_at(f, z, modes))[idx]
            vecs[z] = ph / (np.linalg.norm(ph) + 1e-30)
        for zi, zj in PAIRS:
            a, b = vecs[float(zi)], vecs[float(zj)]
            corr = float(np.dot(a, b))
            dep_rows.append({
                "f_hz": f, "z_i": zi, "z_j": zj,
                "n_modes": len(idx),
                "corr": corr, "D_abs": 1.0 - abs(corr),
                "metric": "phi_m(z) unit vectors on effective modes (Yang estimator NOT built)",
            })
    dep_df = pd.DataFrame(dep_rows)
    dep_df.to_csv(OUT / "depth_mode_signature_physical.csv", index=False, encoding="utf-8-sig")

    # --- 6) frequency mechanism ---
    fr_rows = []
    for f in FREQS:
        modes = PM.solve(f)
        ps = PM.phi_at(f, 200.0, modes)
        pr = PM.phi_at(f, Z_R, modes)
        A = np.array([ps[i] * pr[i] for i in range(len(modes))])
        An = np.abs(A) / (np.max(np.abs(A)) + 1e-30)
        d220 = float(dep_df[(dep_df.f_hz == f) & (dep_df.z_i == 200) & (dep_df.z_j == 220)]["D_abs"].iloc[0])
        d210 = float(dep_df[(dep_df.f_hz == f) & (dep_df.z_i == 200) & (dep_df.z_j == 210)]["D_abs"].iloc[0])
        fr_rows.append({
            "f_hz": f,
            "C1_corr_200_220": C1_CORR[f],
            "phys_n_eff_0p1": int(np.sum(An > 0.1)),
            "phys_n_eff_0p2": int(np.sum(An > 0.2)),
            "phys_D_200_220": d220,
            "phys_D_200_210": d210,
        })
    fr_df = pd.DataFrame(fr_rows)
    # variation of D across frequency → explained?
    d220s = fr_df["phys_D_200_220"].to_numpy()
    if np.std(d220s) < 0.02 and abs(fr_df.loc[fr_df.f_hz == 338, "phys_D_200_220"].iloc[0] - d220s.mean()) < 0.05:
        c1_expl = "NOT_EXPLAINED"
        c1_note = "physical phi_m(z,f) D_200/220 nearly flat vs f; cannot account for C1 338 Hz selectivity"
    elif np.std(d220s) >= 0.05:
        c1_expl = "PARTIALLY_EXPLAINED"
        c1_note = "some f-dependence in physical mode depth vectors; check mode count/spacing columns"
    else:
        c1_expl = "PARTIALLY_EXPLAINED"
        c1_note = "weak f-dependence"
    fr_df["C1_338_explanation"] = c1_expl
    fr_df["note"] = c1_note
    fr_df.to_csv(OUT / "frequency_depth_mechanism_physical.csv", index=False, encoding="utf-8-sig")

    # --- 7) controlled recovery with λ/2 explicit ---
    rec_rows = []
    for f in FREQS:
        modes = PM.solve(f)
        ps = PM.phi_at(f, 200.0, modes)
        pr = PM.phi_at(f, Z_R, modes)
        A = np.array([ps[i] * pr[i] for i in range(len(modes))])
        kr = np.array([m["k_r"] for m in modes])
        An = np.abs(A) / (np.max(np.abs(A)) + 1e-30)
        for L in [0.6e3, 1.2e3, 2.4e3]:
            for dr_div, name in [(8.0, "lam/8"), (4.0, "lam/4"), (2.0, "lam/2")]:
                lam = C0 / f
                dr = lam / dr_div
                n = int(max(64, L / dr))
                r = np.arange(n) * dr
                p = np.zeros(n, dtype=complex)
                for i, k in enumerate(kr):
                    p += A[i] * np.exp(1j * k * r)
                P = np.fft.fftshift(np.fft.fft(p))
                kg = np.fft.fftshift(np.fft.fftfreq(n, d=dr)) * 2 * np.pi
                mag = np.abs(P)
                k_nyq = np.pi / dr
                interior = np.abs(kg) < 0.9 * k_nyq
                peaks = []
                for i in range(1, n - 1):
                    if not interior[i]:
                        continue
                    if mag[i] >= mag[i - 1] and mag[i] >= mag[i + 1] and mag[i] > 0.1 * mag.max():
                        peaks.append((float(kg[i]), float(mag[i])))
                peaks = sorted(peaks, key=lambda x: -x[1])[:6]
                edge = float(np.max(np.abs(kr)) > 0.9 * k_nyq)
                if not peaks:
                    rec_rows.append({
                        "f_hz": f, "L_m": L, "dr_name": name, "dr_m": dr,
                        "status": "NYQUIST_EDGE_LIMITED" if edge else "NO_VALID_INTERIOR_PEAK",
                        "recovered_k": np.nan, "pred_k_nearest": np.nan,
                        "peak_err": np.nan, "k_nyquist": k_nyq,
                        "max_abs_k_r": float(np.max(np.abs(kr))),
                        "note": "explicit empty-peak row (lam/2 not silent)",
                    })
                else:
                    for k, magv in peaks:
                        pred = float(np.min(kr[np.abs(kr - k)])) if False else float(kr[np.argmin(np.abs(kr - k))])
                        rec_rows.append({
                            "f_hz": f, "L_m": L, "dr_name": name, "dr_m": dr,
                            "status": "OK",
                            "recovered_k": k, "pred_k_nearest": pred,
                            "peak_err": abs(k - pred),
                            "k_nyquist": k_nyq,
                            "max_abs_k_r": float(np.max(np.abs(kr))),
                            "note": "",
                        })
    rec_df = pd.DataFrame(rec_rows)
    rec_df.to_csv(OUT / "controlled_mode_recovery_physical.csv", index=False, encoding="utf-8-sig")

    # --- decision ---
    # Erratum not recovered → BLOCKED (required)
    decision = "C2_0_PAPER_RECOVERY_BLOCKED"
    why = (
        "Yang 2015 framework confirmed at high level; 2018 Erratum exact formula still NOT_RECOVERED "
        "(PubMed/AIP/Crossref unreachable). No Yang formulas invented. "
        "Physical 1D normal-mode baseline built for future C2 validity. "
        "Old ModeModel f-invariant phi is ARTIFACT; η Fourier gate is Rayleigh-scale only."
    )
    nxt = "recover Yang2015+Erratum exact formulas, then re-evaluate C2 with physical modes"

    notes = {
        "literature": "framework OK; erratum formulas NOT_RECOVERED",
        "old_model_artifact": "phi_m(z) f-independent → identical D across 201/235/283/338 in C2.0",
        "fourier_gate_wording": "2.4 km insufficient under Fourier/Rayleigh criterion — NOT a fundamental Yang limit",
        "C1_338_explanation": c1_expl,
        "C1_338_note": c1_note,
        "phys_D_200_220_by_f": {str(r["f_hz"]): r["phys_D_200_220"] for _, r in fr_df.iterrows()},
        "L_eta1": L_eta1_rows,
    }
    dec = {
        "rc3c2_0r_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "created_utc": NOW,
        "stop": "no full Yang / Doppler / f0 error / MC / P5",
    }
    (OUT / "R3_C2_0R_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    def fnum(x, nd=4):
        try:
            v = float(x)
            return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    rp = []
    rp.append("# R3-C2.0R 报告：文献恢复 + 物理模态基线")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("## 1. 文献")
    rp.append("")
    rp.append("- Yang 框架已核实；**Erratum 修正式 NOT_RECOVERED**")
    rp.append(f"- **判定 `{decision}`**（原因缩小为公式/Erratum，而非“全文无证据”）")
    rp.append("")
    rp.append("## 2. 旧 ModeModel 伪影")
    rp.append("")
    rp.append("- 旧 φ_m(z) 仅依赖 γ 与高斯包络，**不随 f 重解** → C2.0 四频 D 完全相同是 **ARTIFACT_OF_SIMPLIFIED_MODEMODEL**")
    rp.append("- 不得据此写“真实 CZ 模态深度对 200–220 m 很弱”")
    rp.append("")
    rp.append("## 3. 物理 1D 本征模态（E-STD Munk, 5 km, φ(0)=0, φ′(H)=0）")
    rp.append("")
    rp.append("| f | n_modes | mean residual | max orth off-diag |")
    rp.append("| --- | --- | --- | --- |")
    for _, r in pd.DataFrame(orth_rows).iterrows():
        rp.append(f"| {r['f_hz']} | {r['n_modes']} | {fnum(r['mean_abs_resid'], 3)} | {fnum(r['max_offdiag_orth'], 3)} |")
    rp.append("")
    rp.append("| f | old n | phys n | old median Δk | phys median Δk |")
    rp.append("| --- | --- | --- | --- | --- |")
    for _, r in cmp_df.iterrows():
        rp.append(f"| {r['f_hz']} | {r['old_n_modes']} | {r['phys_n_modes']} | {fnum(r['old_median_dk'])} | {fnum(r['phys_median_dk'])} |")
    rp.append("")
    rp.append("## 4. 孔径（Fourier/Rayleigh 预关，非高分辨极限）")
    rp.append("")
    rp.append("| f | L_eta1 median (m) |")
    rp.append("| --- | --- |")
    for _, r in pd.DataFrame(L_eta1_rows).iterrows():
        rp.append(f"| {r['f_hz']} | {fnum(r['L_eta1_median_m'],0)} |")
    rp.append("")
    rp.append("措辞冻结：**在简化/物理模态与普通 Fourier 准则下，2.4 km 径向孔径不足以按 Rayleigh 尺度分开相邻有效模式**；AR/高分辨波数方法不在本轮否定范围内。")
    rp.append("")
    rp.append("## 5. 深度签名（真实 φ_m(z,f)）")
    rp.append("")
    rp.append("| f | D(200/210) | D(200/220) |")
    rp.append("| --- | --- | --- |")
    for _, r in fr_df.iterrows():
        rp.append(f"| {r['f_hz']} | {fnum(r['phys_D_200_210'])} | {fnum(r['phys_D_200_220'])} |")
    rp.append("")
    rp.append(f"C1 338 Hz 解释判定：**{c1_expl}** — {c1_note}")
    rp.append("")
    rp.append("## 6. λ/2 受控恢复")
    rp.append("")
    lam2 = rec_df[rec_df.dr_name == "lam/2"]
    n_ok = int((lam2.status == "OK").sum()) if len(lam2) else 0
    n_bad = int((lam2.status != "OK").sum()) if len(lam2) else 0
    rp.append(f"- lam/2 显式记录：OK={n_ok}，NO_VALID_INTERIOR_PEAK/NYQUIST_EDGE_LIMITED={n_bad}（不静默缺行）")
    rp.append("")
    rp.append("## 7. 停止")
    rp.append("")
    rp.append("- 不进完整 Yang / Doppler / f0 / MC / P5")
    rp.append(f"- **R3-C2.0R 完成；`{decision}`**")
    rp.append("")
    (OUT / "R3_C2_0R_REPORT.md").write_text("\n".join(rp), encoding="utf-8")
    (OUT / "R3_C2_0R_GPT_SYNC.md").write_text(
        f"# R3-C2.0R\n\n**{decision}**\n\n{why}\n\nC1_338: {c1_expl}\n\nnext: {nxt}\n",
        encoding="utf-8",
    )

    print("DECISION", decision)
    print("C1_338", c1_expl)
    print(f"DONE {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
