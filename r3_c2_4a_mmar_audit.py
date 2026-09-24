#!/usr/bin/env python3
"""R3-C2.4A: C2.3A wording fix + E-STD mode-set robustness. NO MMAR code, NO invented formulas."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
ZGRID = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_1" / "_kraken_zgrid"
A3 = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_3A"
OUT = ROOT / "results" / "R3_C2_AR_MMAR" / "R3_C2_4A"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()
FREQS = [201.0, 235.0, 283.0, 338.0]
L_RAD = 2400.0
DK_R = 2 * np.pi / L_RAD
OMEGA = lambda f: 2 * np.pi * f


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


def main() -> int:
    # ---- 0) C2.3A wording freeze ----
    dec3 = json.loads((A3 / "R3_C2_3A_DECISION.json").read_text(encoding="utf-8"))
    dec3["baseline_strict_success"] = "0/20"
    dec3["baseline_breakdown"] = {
        "19/20": "EMPTY_UNIQUE_MODE",
        "1/20": "NONINFORMATIVE_SINGLE_UNIQUE_MODE (283Hz, ztrue=220m, zhat=250m, margin<1)",
    }
    dec3["unique_group_counts"] = {
        "UNIQUE": 2,
        "UNRESOLVED_GROUP": 1777,
        "scope": "AGGREGATED_OVER_ALL_SCANS (not baseline-20 only)",
    }
    dec3["output_artifact_missing"] = "estd_wavenumber_peaks.csv empty; peaks live in estd_peak_mode_mapping.csv"
    dec3["eq5_note"] = "ORACLE_EQ5_DEPTH_SIGNATURE_SURVIVES (uses true zs for b_m ideal peaks; not recovered performance)"
    dec3["frozen_sentence"] = (
        "E-STD 中模态深度签名在理想上限下仍然存在，但在 L<=2.4 km 的普通 Fourier 合成孔径下，"
        "谱峰通常覆盖多个理论模态，无法可靠完成 Yang Eq.(6) 所要求的谱峰→模态编号对应。"
    )
    (A3 / "R3_C2_3A_DECISION.json").write_text(json.dumps(dec3, ensure_ascii=False, indent=2), encoding="utf-8")

    # ---- 1) PRIMARY PDF attempt log ----
    pdf_note = f"""# MMAR_PRIMARY_SOURCE

UTC: {NOW}

## 目标文献（身份锁定）

Liang Guo-Long, Zhang Yi-Feng, Zou Nan, Wang Jin-Jin,
*Match-Mode Autoregressive Method for Moving Source Depth Estimation in Shallow Water Waveguides*,
Mathematical Problems in Engineering 2018, Article ID 7824671,
DOI **10.1155/2018/7824671**

## PRIMARY_SOURCE_AVAILABLE = FALSE

本轮尝试（全部失败）：

- downloads.hindawi.com/journals/mpe/2018/7824671.pdf → 403
- onlinelibrary.wiley.com/doi/pdf/10.1155/2018/7824671 → 403
- onlinelibrary.wiley.com/doi/pdfdirect/... → 403
- doi.org/10.1155/2018/7824671 → 403

## 状态

**C2_4A_PRIMARY_PDF_REQUIRED**

请用户提供完整 PDF（作者公开稿亦可）。

## 禁止

- 不根据 abstract / 网页片段补 Eq.(1)–(27)
- 不编码 AR / MMAR
- 不把 MMAR 写成 “Yang2015 的 AR 版”

## 方法身份（仅任务书已确认的定性描述，非公式）

MMAR ≠ 单水听器 Yang Eq.(6) 换谱估计器；链路为：

HLA beamforming → 改进 AR 估 k̂ → generalized Hankel/FT 取模态幅度 → 模态匹配定深。

**AR 峰高不得当作模态幅度。**
"""
    (OUT / "MMAR_PRIMARY_SOURCE.md").write_text(pdf_note, encoding="utf-8")
    (OUT / "MMAR_METHOD_IDENTITY.md").write_text(
        f"""# MMAR_METHOD_IDENTITY

UTC: {NOW}

（等待 PRIMARY PDF 后细锁；以下仅任务书已给定性边界）

| 模块 | 来源倾向 | 状态 |
| --- | --- | --- |
| HLA spatial beamforming | Liang2018 新增于 Yang 单听器路线 | WAIT_PDF |
| moving-source range sequence | 与 Yang SAB 同类 | WAIT_PDF |
| AR high-res wavenumber | Liang2018 核心新增 | WAIT_PDF |
| generalized Hankel amplitudes | 与 Yang/Ref.7 同源 | WAIT_PDF |
| matched-mode depth | Liang2018 Eq.(20)+ | WAIT_PDF |
| mode-order matching Eq.(24)–(27) | Liang2018 关键 | WAIT_PDF |

**不得称为 Yang2015 已验证的 AR 版本。**
""",
        encoding="utf-8",
    )
    for name, body in {
        "MMAR_EQUATION_LOCK.md": "# MMAR_EQUATION_LOCK\n\nNOT_RECOVERED — PRIMARY_PDF_REQUIRED\n禁止猜公式。\n",
        "MMAR_EQUATION_TABLE.csv": "equation_id,exact_formula,paper_page,symbols,input,output,purpose,implementation_required,confidence\n"
        + "".join([f"Eq.({i}),NOT_RECOVERED,,,,,,,NOT_RECOVERED\n" for i in range(1, 28)]),
        "MMAR_SYMBOL_TABLE.csv": "symbol,definition,status\nALL,NOT_RECOVERED,WAIT_PDF\n",
        "MMAR_OBSERVABLE_CHAIN.md": "# MMAR_OBSERVABLE_CHAIN\n\nWAIT_PDF：AR 只定 k̂；幅度来自 Hankel/FT（任务书约束，细节待原文）。\n禁止 AR peak height 作为 amplitude。\n",
        "MMAR_MODE_ORDER_MATCHING.md": "# MMAR_MODE_ORDER_MATCHING\n\nWAIT_PDF：Eq.(24)–(27) 未恢复前禁止发明 Hungarian/nearest matching。\n",
        "MMAR_PAPER_CONFIG.md": "# MMAR_PAPER_CONFIG\n\nWAIT_PDF 部分；任务书已给 simulation 骨架（f=350, zs=4/50, HLA 11 el, d=λ/2, HLA z=70m, v=2.5 m/s, r0=5010m, SNR 20/5/-5）待原文逐项核验。\n",
        "MMAR_PAPER_FIDELITY_TARGETS.md": "# MMAR_PAPER_FIDELITY_TARGETS\n\nWAIT_PDF：Fig3–7 目标待原文核验后锁定。\n",
        "MMAR_ESTD_TRANSFER_AUDIT.md": "# MMAR_ESTD_TRANSFER_AUDIT\n\nWAIT_PDF 时先列任务书已知假设边界（range-independent / fixed-depth / known f0,Δr / HLA / far-field / shallow vs CZ modal density）。细节待原文。\n",
    }.items():
        (OUT / name).write_text(body, encoding="utf-8")
        if name.endswith(".csv") and name != "MMAR_EQUATION_TABLE.csv":
            pass

    # ---- 2) E-STD mode family classification (zs-independent) ----
    fam_rows = []
    set_rows = []
    res_rows = []
    for f in FREQS:
        M = parse_mod(ZGRID / f"zgrid_f{int(f)}.mod")
        depths, phi, k = M["depths"], M["phi"], M["k"]
        nmode = M["M"]
        k_re, k_im = k.real, k.imag
        alpha = -k_im
        omega = OMEGA(f)
        izr = int(np.argmin(np.abs(depths - 200)))
        # search band 150-250 on zgrid nodes
        trapped_idx = []
        iband = np.where((depths >= 150) & (depths <= 250))[0]
        max_band = np.max(np.abs(phi[iband]), axis=0) if iband.size else np.zeros(nmode)
        phi_r = np.abs(phi[izr])
        max_phi_all = float(np.max(np.abs(phi)))
        # attenuation to 50 km
        att50 = np.exp(-alpha * 50000.0)
        for m in range(nmode):
            cp = omega / k_re[m] if k_re[m] != 0 else np.nan
            # turning classification by depth of max |phi| and near-boundary energy
            izmax = depths[int(np.argmax(np.abs(phi[:, m])))]
            e_surf = float(np.mean(np.abs(phi[depths <= 40, m]) ** 2)) if np.any(depths <= 40) else 0
            e_bot = float(np.mean(np.abs(phi[depths >= 800, m]) ** 2)) if np.any(depths >= 800) else 0
            e_mid = float(np.mean(np.abs(phi[(depths > 200) & (depths < 1300), m]) ** 2))
            if e_mid >= max(e_surf, e_bot):
                cls = "internally_refracted_trapped"
            elif e_surf >= e_bot and e_surf > 0.5 * e_mid:
                cls = "surface_interacting"
            elif e_bot > e_surf and e_bot > 0.5 * e_mid:
                cls = "bottom_interacting"
            else:
                cls = "unresolved_classification"
            fam_rows.append(
                {
                    "f_hz": f,
                    "mode_id": m + 1,
                    "k_r": k_re[m],
                    "phase_speed": cp,
                    "alpha": alpha[m],
                    "phi_zr200": phi[izr, m],
                    "abs_phi_zr200": phi_r[m],
                    "max_abs_phi_search_band": max_band[m],
                    "z_at_max_phi": izmax,
                    "atten_50km": att50[m],
                    "turning_class": cls,
                }
            )
            if cls == "internally_refracted_trapped":
                trapped_idx.append(m)

        # candidate sets (zs-independent)
        rel_r = phi_r / max(phi_r.max(), 1e-30)
        rel_b = max_band / max(max_band.max(), 1e-30)
        for thr in (0.05, 0.10, 0.20):
            set_rows.append({"f_hz": f, "set": "SET-B_RECEIVER_OBSERVABLE", "threshold": thr, "n_modes": int((rel_r >= thr).sum())})
            set_rows.append({"f_hz": f, "set": "SET-C_SEARCH_BAND_SUPPORTED", "threshold": thr, "n_modes": int((rel_b >= thr).sum())})
        trapped = list(trapped_idx)
        set_rows.append({"f_hz": f, "set": "SET-D_REFRACTED_TRAPPED", "threshold": np.nan, "n_modes": len(trapped)})
        set_rows.append({"f_hz": f, "set": "SET-A_ALL_KRAKEN", "threshold": np.nan, "n_modes": nmode})

        # Fourier robustness on each set
        def rayleigh_unique(idx):
            if not idx:
                return 0
            kk = np.sort(k_re[idx])
            uniq = 0
            for i in range(len(kk)):
                dl = kk[i] - kk[i - 1] if i > 0 else np.inf
                dr = kk[i + 1] - kk[i] if i < len(kk) - 1 else np.inf
                if dl > DK_R and dr > DK_R:
                    uniq += 1
            return uniq

        sets = {
            "SET-A": list(range(nmode)),
            "SET-B_0.05": list(np.where(rel_r >= 0.05)[0]),
            "SET-B_0.10": list(np.where(rel_r >= 0.10)[0]),
            "SET-B_0.20": list(np.where(rel_r >= 0.20)[0]),
            "SET-C_0.05": list(np.where(rel_b >= 0.05)[0]),
            "SET-C_0.10": list(np.where(rel_b >= 0.10)[0]),
            "SET-C_0.20": list(np.where(rel_b >= 0.20)[0]),
            "SET-D_trapped": trapped,
        }
        for sn, idx in sets.items():
            kk = np.sort(k_re[idx]) if idx else np.array([])
            if len(kk) >= 2:
                dk = np.diff(kk)
                p10, p50, p90 = np.percentile(dk, [10, 50, 90])
                dmin = float(dk.min())
            else:
                p10 = p50 = p90 = dmin = np.nan
            set_rows.append(
                {
                    "f_hz": f,
                    "set": sn + "_FOURIER",
                    "threshold": np.nan,
                    "n_modes": len(idx),
                    "n_rayleigh_unique": rayleigh_unique(idx),
                    "dk_min": dmin,
                    "dk_p10": p10,
                    "dk_p50": p50,
                    "dk_p90": p90,
                }
            )
            # required L = 2pi/dk
            for stat, val in (("min", dmin), ("p10", p10), ("p50", p50), ("p90", p90)):
                if val and np.isfinite(val) and val > 0:
                    res_rows.append({"f_hz": f, "set": sn, "dk_stat": stat, "dk": float(val), "L_required_m": float(2 * np.pi / val)})

    pd.DataFrame(fam_rows).to_csv(OUT / "estd_mode_family_classification.csv", index=False)
    pd.DataFrame(set_rows).to_csv(OUT / "estd_mode_set_fourier_robustness.csv", index=False)
    pd.DataFrame(res_rows).to_csv(OUT / "ESTD_REQUIRED_WAVENUMBER_RESOLUTION.csv", index=False)

    # robustness label
    fs = pd.DataFrame(set_rows)
    fs_uniq = fs[fs["set"].str.endswith("_FOURIER")]
    max_uniq = int(fs_uniq["n_rayleigh_unique"].fillna(0).max()) if not fs_uniq.empty else 0
    fourier_label = "FOURIER_MODE_IDENTITY_LIMIT_ROBUST_TO_MODE_SET" if max_uniq <= 1 else "FOURIER_LIMIT_DEPENDS_ON_MODE_SET"

    # ---- 3) HLA migration table (geometry only; no AR) ----
    # paper: 11 elements (L=5 => 2L+1), d=lambda/2, HLA depth 70 m (from task brief — WAIT_PDF for formula-level)
    hla_rows = []
    for f in FREQS:
        lam = 1500.0 / f  # c_min based; paper uses lambda/2 spacing
        for name, n, d_m, note in [
            ("A_PAPER_EQUIV", 11, 0.5 * lam, "paper: 2*5+1 elements, d=lambda/2"),
            ("B_PROJECT_14M", 14, 1.0, "14 m physical HLA representative (array pitch TBD from project freeze)"),
            ("C_4EL_ENHANCE", 4, 1.0, "4-element condition-enhancement if frozen config requires"),
        ]:
            aperture = (n - 1) * d_m
            # array gain upper bound ~ N for incoherent noise (ideal)
            ag = float(n)
            # far-field: 2 D^2 / lambda
            ff = 2.0 * aperture**2 / lam if lam > 0 else np.nan
            hla_rows.append(
                {
                    "f_hz": f,
                    "config": name,
                    "N_elem": n,
                    "d_m": d_m,
                    "physical_aperture_m": aperture,
                    "array_gain_upper_bound": ag,
                    "far_field_min_range_m": ff,
                    "note": note,
                }
            )
    pd.DataFrame(hla_rows).to_csv(OUT / "MMAR_HLA_MIGRATION_TABLE.csv", index=False)

    # ---- 4) decision ----
    decision = "C2_4A_PRIMARY_PDF_REQUIRED"
    why = (
        "Liang2018 MMAR 主文 PDF 未取得（Hindawi/Wiley/DOI 均 403）。"
        "按规则停止公式恢复与编码。C2.3A 口径已修正；E-STD 模式族/Fourier 稳健性与 AR 分辨率需求已审计（不依赖 MMAR 公式）。"
    )
    dec = {
        "stage": "R3-C2.4A",
        "rc3_c2_4a_decision": decision,
        "why": why,
        "fourier_mode_set_label": fourier_label,
        "max_rayleigh_unique_across_sets": max_uniq,
        "c2_3a_baseline": "0/20 effective (19 EMPTY_UNIQUE_MODE + 1 NONINFORMATIVE_SINGLE_UNIQUE_MODE)",
        "yang_route": "YANG_ROUTE_UNDECIDED",
        "mmar_route": "RC3-C2-AR/MMAR independent candidate; not Yang-AR",
        "not_done": ["AR code", "MMAR run", "new FIELD", "E-STD depth", "MC", "RC2-delta", "P5"],
        "created_utc": NOW,
    }
    (OUT / "R3_C2_4A_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8")

    (OUT / "R3_C2_4A_REPORT.md").write_text(
        f"""# R3-C2.4A MMAR primary-source recovery + E-STD 准入前审计

UTC: {NOW}

## 判定

### `{decision}`

{why}

**请用户提供 Liang et al. 2018 PDF（DOI 10.1155/2018/7824671）。**

## C2.3A 口径冻结（不重跑）

- baseline strict：**0/20 有效**（19 EMPTY_UNIQUE_MODE；1 例 283 Hz/z=220 仅 1 unique 且 ẑ=250、margin<1 → NONINFORMATIVE）
- UNIQUE=2 / GROUP=1777：**AGGREGATED_OVER_ALL_SCANS**
- `estd_wavenumber_peaks.csv` 空：`OUTPUT_ARTIFACT_MISSING`
- EQ5 20/20：**ORACLE_EQ5_DEPTH_SIGNATURE_SURVIVES**（非可恢复性能）

冻结句：E-STD 理想深度签名存在，但 L≤2.4 km 普通 Fourier 无法可靠完成谱峰→唯一模态编号。

## 模式族稳健性（无 true zs）

见 `estd_mode_family_classification.csv`、`estd_mode_set_fourier_robustness.csv`。

**{fourier_label}**（max n_rayleigh_unique={max_uniq}）

## AR 需要分开的 Δk

见 `ESTD_REQUIRED_WAVENUMBER_RESOLUTION.csv`（由相邻模态 Δk 分位数换算 L_required=2π/Δk）。

## HLA 迁移（仅几何）

见 `MMAR_HLA_MIGRATION_TABLE.csv`。**HLA 不是免费**；后续须区分 MMAR-PAPER vs SINGLE_SENSOR_ABLATION。

## 停止

不编码 AR、不跑 MMAR、不进 P5。等待用户 PDF + GPT 审计。
""",
        encoding="utf-8",
    )
    (OUT / "R3_C2_4A_GPT_SYNC.md").write_text(
        f"# R3-C2.4A\\n\\n**{decision}**\\n\\n请提供 Liang2018 PDF。\\n{fourier_label}\\nC2.3A 0/20 effective.\\n",
        encoding="utf-8",
    )
    print("DECISION", decision, fourier_label, "max_uniq", max_uniq)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
