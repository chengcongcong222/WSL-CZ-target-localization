#!/usr/bin/env python3
"""R3-0 / RC3 Literature Method Recovery + R3-A start (distance baseline).

Outputs under results/R3_RC3_method_recovery/.
No P5, no Bellhop audit, no auto R3-B/C.
"""
from __future__ import annotations

import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import p4_g4_performance_boundary as p4

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_RC3_method_recovery"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

# E-STD
C0 = 1500.0
Z_S_TRUE = 200.0
Z_R_TRUE = 200.0
R_TRUE = 50e3
Z_TRUE = 200.0
R_GRID = np.linspace(45e3, 60e3, 61)       # 0.25 km
Z_GRID = np.linspace(150.0, 250.0, 41)     # 2.5 m
F_MFP = np.linspace(150.0, 375.0, 24)     # broadband MFP
F_FRINGE = np.linspace(150.0, 375.0, 46)  # fringe spectrogram
S2_LINES = [166.0, 201.0, 235.0, 283.0, 338.0]
S1_LINES = [168.0, 204.0, 232.0, 279.0, 320.0]

CONFIG = {
    "package": "R3_RC3_method_recovery",
    "created_utc": NOW,
    "reanchor": "RC3 formal literature-method validation; P4.5 is simplified baseline only, NOT RC3 upper bound",
    "project_progress": {
        "RC1": "scene frozen, basis formed",
        "RC2": "kinematic baseline + hard candidates largely done",
        "RC3": "formal method validation just entering main phase",
        "five_state_perf": "prior results = simplified baseline only",
    },
    "E_STD": {
        "r_true_km": 50.0,
        "z_true_m": 200.0,
        "zs_zr_m": 200.0,
        "r_search_km": [45.0, 60.0],
        "z_search_m": [150.0, 250.0],
        "band_hz": [150.0, 375.0],
        "propagation": "normal-mode theory model consistent with P3 (no Bellhop)",
    },
    "frozen_mainlines": {
        "upper": "Broadband MFP (range-depth ambiguity) — physics ceiling, not final engineering method",
        "RC3-A": "CZ broadband interference fringes / waveguide-invariant distance (r, r-change)",
        "RC3-B": "MMAC multipath angle+delay joint match (r,z)",
        "RC3-C": "HLA modal/wavenumber spectrum + synthetic-aperture depth discrimination (z)",
    },
    "protocol": [
        "1 paper-like reproduction (original observables/formulas)",
        "2 migrate to E-STD",
        "3 RC2 hard-pair comparison",
        "4 parameter RMSE/coverage (P4.5-style)",
    ],
    "allowed_conclusions": ["可用", "有条件可用", "与RC2冗余", "当前场景不适用"],
    "stop_this_round": [
        "literature matrix", "freeze three RC3 mainlines",
        "MFP upper bound", "RC3-A first paper-like reproduction",
    ],
    "not_this_round": ["R3-B full", "R3-C full", "P5", "Bellhop audit"],
}


# ---------------------------------------------------------------------------
# Literature method matrix (project-recovered + standard literature mapping)
# ---------------------------------------------------------------------------
METHODS = [
    dict(
        method_id="MFP",
        method_name="宽带匹配场处理 Matched-Field Processing",
        literature="Bucker 1976; Tolstoy/Clay 1987 ocean acoustics matched field; deep-water range-depth MFP classics",
        observables="阵元复声压/宽带频谱向量",
        target_params="r, z (可扩展 θ)",
        formula_flow="p_rep(f,x;r,z)=Σ_m A_m φ_m(z_s)φ_m(z_r)e^{ik_mr}/√k_m; J_Bartlett=|w^H p_rep|²/(norms); unknown S(f) via spectral normalization/whitening",
        array_req="垂直阵最经典；水平阵也可用但深度模糊更重；单水听器多频可做降维MFP",
        band_signal_req="宽带有利；S0连续谱信息最多；S2多线谱次之",
        env_prior_req="需要环境/SSP 做前向场（匹配环境）",
        paper_perf_metric="模糊面主瓣宽、旁瓣比、定位RMSE",
        e_std_transfer_gap="深海第一会聚区45–60km、zs=zr=200m；水平拖曳阵 vs 文献常用VLA；环境匹配假设强",
        project_impl_degree="未用完整复声压MFP；P3/P4.5仅用CZ幅度轮廓简化特征",
        mainline_role="上限基准",
        selected=True,
        selection_reason="定义E-STD完整传播场理论上r/z信息上限，便于判断物理是否本身就难",
    ),
    dict(
        method_id="WI_FRINGE",
        method_name="波导不变量/距离–频率干涉条纹 (Waveguide invariant / range-frequency fringes)",
        literature="S.D. Chuprov 1980s waveguide invariant; Ratilal/Tindle; deep-water CZ fringe literature; LOFARgram fringe slope methods",
        observables="距离–频率强度图 I(r,f) 或运动平台时–频谱 I(t,f)→r(t)",
        target_params="r, Δr, β(波导不变量)",
        formula_flow="干涉条纹 dr/df≈β r/f; 2D FFT/Radon/条纹斜率提取→β; 由斜率或条纹间距反演距离/距离变化",
        array_req="单水听器/单阵即可（能量图）；水平阵可增强",
        band_signal_req="宽带连续谱最佳；稳定线谱可在离散频点采样条纹",
        env_prior_req="弱–中：β理想波导≈1；深海CZ需数值/经验修正",
        paper_perf_metric="距离RMSE、β估计误差、条纹可见度",
        e_std_transfer_gap="文献多浅海/中距离；深海CZ条纹间距与β不同；平台运动导致r(t)时变",
        project_impl_degree="未按原文条纹提取/β公式复现；仅有简化CZ轮廓",
        mainline_role="RC3-A距离主线",
        selected=True,
        selection_reason="与甲方最初关注的会聚区干涉条纹直接对应；RC2最确定缺口是距离",
    ),
    dict(
        method_id="MMAC",
        method_name="多途到达特征匹配 MMAC (multipath angle + delay)",
        literature="deep-water multipath arrival-angle/delay localization; HLA multipath spatial structure; angle more range-sensitive, delay more depth-sensitive",
        observables="各路径到达角(AOA)+相对时延τ；可加幅度",
        target_params="r, z",
        formula_flow="前向模型给出路径i: (θ_i(r,z), τ_i(r,z)); J=Σ_i w_i[(θ_i-θ̂_i)²/σ_θ²+(τ_i-τ̂_i)²/σ_τ²]; 网格/优化搜索(r,z)",
        array_req="需要阵列分辨多途角度→水平拖曳阵/HLA关键",
        band_signal_req="宽带利于时延分辨；稳定谱线可辅助",
        env_prior_req="需要传播模型给出多途结构",
        paper_perf_metric="(r,z)模糊面、联合RMSE",
        e_std_transfer_gap="深海CZ多途与浅海不同；200m/200m深度；需先RC2增强方位再用传播结构",
        project_impl_degree="未正式实现角+时延联合；历史M4包有多途相对延迟候选但环境关联门失败",
        mainline_role="RC3-B距离–深度主线",
        selected=True,
        selection_reason="文献指出AOA偏距、时延偏差；最符合“方位增强→传播结构→r/z”原始技术链",
    ),
    dict(
        method_id="HLA_MODAL_SA",
        method_name="水平阵模态/波数谱 + 合成孔径深度判别",
        literature="T.C. Yang et al. horizontal line array wavenumber spectrum / synthetic aperture depth discrimination; modal beamforming on HLA",
        observables="HLA阵元序列→波数谱(k_r或k_z结构)；运动合成孔径扩展孔径",
        target_params="z (主), 辅助r",
        formula_flow="阵列FFT/波束→k谱; 模态横向波数与深度结构相关; 合成孔径增大有效孔径提升深度模态可分性",
        array_req="长水平拖曳阵 + 平台运动（合成孔径）",
        band_signal_req="宽带/多频有利",
        env_prior_req="需要模态/波导模型解释k谱",
        paper_perf_metric="深度判别正确率、z RMSE、模糊面z切片宽度",
        e_std_transfer_gap="文献场景/阵长/频带需核对；E-STD 45–60km CZ；200m深度是否可分待验证",
        project_impl_degree="仅讨论过、未正式实现；P4.5深度UNRESOLVED不能替代本方法验证",
        mainline_role="RC3-C深度主线",
        selected=True,
        selection_reason="专门补当前最弱深度维；Yang水平阵路线此前未进入统一验证",
    ),
    dict(
        method_id="MATCHED_MODE",
        method_name="多线谱匹配模态 Matched-mode / modal matching",
        literature="matched-mode processing in deep water; normal-mode filtering",
        observables="阵列模态滤波输出或单点多频复谱",
        target_params="r, z",
        formula_flow="估计模态幅度→与候选(r,z)预测模态幅度匹配",
        array_req="有利于VLA/能分辨模态的阵；HLA模态处理变体",
        band_signal_req="多频/多线谱",
        env_prior_req="强：需要模态模型",
        paper_perf_metric="r,z RMSE",
        e_std_transfer_gap="与MFP信息同源；工程上模态提取更难",
        project_impl_degree="未实现",
        mainline_role="非主线（与MFP信息高度重叠）",
        selected=False,
        selection_reason="本轮用MFP作上限即可覆盖模态匹配的信息论角色",
    ),
    dict(
        method_id="ARRAY_INVARIANT",
        method_name="阵不变量 Array invariant (ω–k / slope)",
        literature="array invariant χ from array element wavenumber-frequency plot",
        observables="阵元–频率–相位/波数",
        target_params="等效距离/斜距相关量",
        formula_flow="从ω–k图提取不变量χ→距离类估计",
        array_req="阵列",
        band_signal_req="宽带",
        env_prior_req="相对弱",
        paper_perf_metric="距离RMSE",
        e_std_transfer_gap="深海CZ适用性待核；与条纹类信息部分重叠",
        project_impl_degree="未实现",
        mainline_role="后备（控制范围）",
        selected=False,
        selection_reason="本轮优先WI/条纹作距离主线，阵不变量暂不扩展",
    ),
    dict(
        method_id="DOPPLER_TMA",
        method_name="谱线+Doppler TMA",
        literature="bearing-only/Doppler TMA; frequency-rate tracking",
        observables="方位 + 线谱频率轨迹",
        target_params="r, v, ψ",
        formula_flow="f_obs=f0(1+v_rad/c); 与方位联合TMA",
        array_req="方位需阵列；频谱单点即可",
        band_signal_req="稳定谱线 S1/S2；S0不可用",
        env_prior_req="弱（运动学）",
        paper_perf_metric="TMA RMSE",
        e_std_transfer_gap="未知f0时恒定v_rad可被吸收（P3已示）",
        project_impl_degree="P3对照支路已做，增量近零（恒定v_rad对）",
        mainline_role="非主线（运动学对照）",
        selected=False,
        selection_reason="已作为对照；非传播RC3主线",
    ),
    dict(
        method_id="PBPR_SHALLOW",
        method_name="浅海被动宽带距离/深度 PBPR类",
        literature="shallow-water passive broadband range/depth (thermocline, active-like features in literature)",
        observables="宽带能量/干涉结构",
        target_params="r, z",
        formula_flow="浅海波导特征匹配",
        array_req="阵或单点",
        band_signal_req="宽带",
        env_prior_req="强浅海环境",
        paper_perf_metric="浅海r,z RMSE",
        e_std_transfer_gap="浅海温跃层/主动目标 vs 深海CZ 200m/200m，物理跨度过大",
        project_impl_degree="不迁移",
        mainline_role="后备",
        selected=False,
        selection_reason="场景不匹配，暂不进入主线",
    ),
]

FROZEN_MAINLINES = [
    dict(role="上限基准", method_id="MFP", name="宽带匹配场MFP", targets="r,z", why="判断理想声场信息上限"),
    dict(role="RC3-A距离主线", method_id="WI_FRINGE", name="会聚区宽带干涉条纹/波导不变量", targets="r,Δr", why="对应甲方会聚区条纹关注；补RC2距离缺口"),
    dict(role="RC3-B距离–深度主线", method_id="MMAC", name="多途角+时延联合MMAC", targets="r,z", why="AOA偏距/τ偏差；HLA方位增强后用传播结构"),
    dict(role="RC3-C深度主线", method_id="HLA_MODAL_SA", name="HLA模态/波数谱+合成孔径深度", targets="z", why="专门补深度维；Yang路线需正式验证"),
]


# ---------------------------------------------------------------------------
# Modal forward (reuse E-STD theory model)
# ---------------------------------------------------------------------------
def ensure_mode():
    if not p4.AMP_TABS:
        # build amp tables for MFP frequencies separately via ModeModel
        pass
    return p4.MODE_ENV["E0"]


def pressure_spectrum(mode, freqs, r_m, z_s, z_r):
    """Complex modal pressure vector (unknown absolute source)."""
    p = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        ms = mode.modes(float(f))
        acc = 0j
        rr = max(float(r_m), 1.0)
        for kr, phi in ms:
            a_s = float(np.interp(z_s, mode.z, phi))
            a_r = float(np.interp(z_r, mode.z, phi))
            att = np.exp(-2e-5 * (f / 200.0) * rr / 1000.0)
            acc += (a_s * a_r / math.sqrt(kr * rr)) * np.exp(1j * kr * rr) * att
        p[i] = acc
    return p


def unknown_source_match(p_obs, p_rep):
    """Eliminate complex S(f): residual after LS scaling per frequency — use broadband joint gain.
    For single snapshot broadband: maximize |<p_obs,p_rep>| / (||p_obs|| ||p_rep||).
    Also try per-freq amplitude only (incoherent Bartlett).
    """
    # coherent broadband matched field with single global complex S
    num = np.vdot(p_rep, p_obs)  # sum conj(p_rep)*p_obs
    den = (np.linalg.norm(p_rep) * np.linalg.norm(p_obs)) + 1e-30
    coherent = abs(num) / den
    # incoherent: per-frequency |p| match after removing S(f) amplitude
    a_obs = np.abs(p_obs)
    a_rep = np.abs(p_rep)
    # optimal S_i = a_obs_i / a_rep_i → residual on amplitudes
    s = a_obs / (a_rep + 1e-30)
    resid = a_obs - s * a_rep
    # better: compare normalized amplitude spectra
    no = a_obs / (np.linalg.norm(a_obs) + 1e-30)
    nr = a_rep / (np.linalg.norm(a_rep) + 1e-30)
    amp_match = float(np.dot(no, nr))
    # J: 1 - similarity (lower better)
    J_coh = 1.0 - coherent
    J_amp = 1.0 - amp_match
    return J_coh, J_amp, coherent, amp_match


def run_mfp_upper_bound():
    mode = ensure_mode()
    # Truth observations: S0 / S1 / S2 frequency sets
    source_sets = {
        "S0": F_MFP.copy(),
        "S1": np.array(S1_LINES, float),
        "S2": np.array(S2_LINES, float),
    }
    truths = [
        dict(tag="ON_GRID", r_m=50e3, z=200.0),
        dict(tag="OFF_GRID", r_m=52.3e3, z=213.0),
        dict(tag="OFF_GRID2", r_m=47.6e3, z=187.0),
    ]
    rows = []
    amb_saves = {}
    for src, freqs in source_sets.items():
        for tr in truths:
            p_obs = pressure_spectrum(mode, freqs, tr["r_m"], Z_S_TRUE, tr["z"])
            Jc = np.zeros((len(Z_GRID), len(R_GRID)))
            Ja = np.zeros_like(Jc)
            for iz, z in enumerate(Z_GRID):
                for ir, r in enumerate(R_GRID):
                    p_rep = pressure_spectrum(mode, freqs, r, Z_S_TRUE, z)
                    jc, ja, _, _ = unknown_source_match(p_obs, p_rep)
                    Jc[iz, ir] = jc
                    Ja[iz, ir] = ja
            # use incoherent amplitude match as primary (more robust to phase/unknown S)
            J = Ja
            # best
            iz, ir = np.unravel_index(np.argmin(J), J.shape)
            r_hat = float(R_GRID[ir])
            z_hat = float(Z_GRID[iz])
            # mainlobe widths at +0.05 from min (relative ambiguity threshold)
            jmin = float(J.min())
            def width_1d(vec, grid, thr):
                idx = np.where(vec <= thr)[0]
                if len(idx) < 2:
                    return 0.0
                return float(grid[idx.max()] - grid[idx.min()])
            thr = jmin + 0.05
            # r slice at best z, z slice at best r
            wr = width_1d(J[iz, :], R_GRID, thr)
            wz = width_1d(J[:, ir], Z_GRID, thr)
            # multimodal: local minima count in r at best z
            slice_r = J[iz, :]
            n_loc = 0
            for i in range(1, len(slice_r) - 1):
                if slice_r[i] <= slice_r[i - 1] and slice_r[i] <= slice_r[i + 1] and slice_r[i] <= jmin + 0.08:
                    n_loc += 1
            # truth global?
            it_r = int(np.argmin(np.abs(R_GRID - tr["r_m"])))
            it_z = int(np.argmin(np.abs(Z_GRID - tr["z"])))
            truth_is_best = bool(iz == it_z and ir == it_r)
            # r-z correlation of ambiguity basin
            mask = J <= jmin + 0.05
            corr = np.nan
            if mask.sum() >= 5:
                rr = R_GRID[np.where(mask)[1]] / 1e3
                zz = Z_GRID[np.where(mask)[0]]
                # need same length — use all mask cells
                idx = np.where(mask)
                rr = R_GRID[idx[1]] / 1e3
                zz = Z_GRID[idx[0]]
                if np.std(rr) > 1e-9 and np.std(zz) > 1e-9:
                    corr = float(np.corrcoef(rr, zz)[0, 1])
            rows.append({
                "source": src,
                "truth_tag": tr["tag"],
                "r_true_km": tr["r_m"] / 1e3,
                "z_true_m": tr["z"],
                "r_hat_km": r_hat / 1e3,
                "z_hat_m": z_hat,
                "err_r_km": r_hat / 1e3 - tr["r_m"] / 1e3,
                "err_z_m": z_hat - tr["z"],
                "J_min": jmin,
                "r_mainlobe_km": wr / 1e3,
                "z_mainlobe_m": wz,
                "n_local_minima_r": n_loc,
                "truth_is_global_opt": truth_is_best,
                "corr_rz_basin": corr,
                "n_freqs": int(len(freqs)),
            })
            # save ambiguity for S0 on-grid
            if src == "S0" and tr["tag"] == "ON_GRID":
                amb_saves["S0"] = J
                # write full grid CSV
                rec = []
                for iz, z in enumerate(Z_GRID):
                    for ir, r in enumerate(R_GRID):
                        rec.append({"r_km": r / 1e3, "z_m": z, "J_amp": J[iz, ir], "J_coh": Jc[iz, ir]})
                pd.DataFrame(rec).to_csv(OUT / "MFP_UPPER_BOUND_rz.csv", index=False, encoding="utf-8-sig")
            print(f"  MFP {src} {tr['tag']}: r_hat={r_hat/1e3:.2f} wr={wr/1e3:.2f}km "
                  f"z_hat={z_hat:.1f} wz={wz:.1f}m best={truth_is_best} nmin={n_loc}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "MFP_UPPER_BOUND_summary.csv", index=False, encoding="utf-8-sig")

    # figure ambiguity
    if "S0" in amb_saves:
        write_ambiguity_svg(FIG / "MFP_ambiguity_S0.svg", amb_saves["S0"], R_GRID, Z_GRID)
    return df


def write_ambiguity_svg(path, J, R, Z):
    w, h = 760, 420
    ml, mr, mt, mb = 70, 90, 48, 56
    pw, ph = w - ml - mr, h - mt - mb
    # show similarity-like: 1-J
    S = 1.0 - J
    smin, smax = float(S.min()), float(S.max())
    nx, ny = len(R), len(Z)
    cw, ch = pw / nx, ph / ny
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="28" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600">MFP 上限模糊面 J(r,z)（S0，匹配环境，无噪）</text>']
    for iz in range(ny):
        for ir in range(nx):
            u = (S[iz, ir] - smin) / (smax - smin + 1e-15)
            r = int(240 * (1 - u) + 31 * u)
            g = int(230 * (1 - u) + 111 * u)
            b = int(220 * (1 - u) + 139 * u)
            x = ml + ir * cw
            y = mt + (ny - 1 - iz) * ch
            p.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(cw,0.5):.1f}" height="{max(ch,0.5):.1f}" fill="rgb({r},{g},{b})"/>')
    # truth mark
    ir = int(np.argmin(np.abs(R - 50e3)))
    iz = int(np.argmin(np.abs(Z - 200)))
    x = ml + ir * cw + cw / 2
    y = mt + (ny - 1 - iz) * ch + ch / 2
    p.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="8" fill="none" stroke="#fff" stroke-width="2"/>')
    p.append(f'<text x="{x:.1f}" y="{y-12:.1f}" text-anchor="middle" font-size="10" fill="#fff" font-family="sans-serif">真值</text>')
    p.append(f'<text x="{ml+pw/2}" y="{h-14}" text-anchor="middle" font-size="12" font-family="sans-serif">r (km) 45→60</text>')
    p.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90 16 {mt+ph/2})">z (m) 150→250</text>')
    p.append('</svg>')
    path.write_text("\n".join(p), encoding="utf-8")


# ---------------------------------------------------------------------------
# RC3-A: waveguide invariant / range-frequency fringe — paper-like reproduction
# ---------------------------------------------------------------------------
def intensity_map(mode, freqs, r_grid, z_s, z_r):
    """I(r,f) modal intensity."""
    I = np.zeros((len(freqs), len(r_grid)))
    for i, f in enumerate(freqs):
        for j, r in enumerate(r_grid):
            p = pressure_spectrum(mode, [f], r, z_s, z_r)[0]
            I[i, j] = abs(p) ** 2
    return I


def fringe_beta_estimate(I, freqs, r_grid, beta_grid=None):
    """Paper-like: remove smooth range trend, 2D FFT/Radon-like slope search for β.

    Waveguide-invariant striation model (idealized):
      I(r,f) ~ trend(r,f) * cos(2π/Λ * (r - β * r0 * ln(f/f0)) )  approx
    or slope in (r,f): dr/df ≈ β * r / f
    We search β maximizing fringe coherence after range-frequency de-trend.
    """
    if beta_grid is None:
        beta_grid = np.linspace(0.3, 2.5, 45)
    # log intensity
    L = np.log(I + 1e-20)
    # remove per-frequency mean over range (broadband trend)
    L0 = L - L.mean(axis=1, keepdims=True)
    # remove per-range smooth trend with poly fit in f
    L1 = L0.copy()
    for j in range(L0.shape[1]):
        x = freqs
        y = L0[:, j]
        c = np.polyfit(x, y, 2)
        L1[:, j] = y - np.polyval(c, x)
    # For each beta, warp range axis: r_eff = r * (f/f_ref)^{beta-1} or use
    # phase-alignment: stripe phase φ = 2π * r / Λ(f) with Λ ∝ f^{1} / something
    # Practical paper-like: Radon on (r,f) after noting stripe slope s = dr/df = β r/f
    # Integrate along lines r = r0 + β r0 ln(f/f0)  (integral of β r/f df)
    f0 = float(np.mean(freqs))
    best = dict(beta=np.nan, score=-np.inf)
    scores = []
    for beta in beta_grid:
        # warp: for each range bin r_ref, integrate L1 along r(f)=r_ref*(f/f0)^beta
        # discrete: interpolate
        score = 0.0
        n = 0
        for jref in range(0, len(r_grid), 4):
            r_ref = r_grid[jref]
            acc = np.zeros(len(freqs))
            for i, f in enumerate(freqs):
                r_f = r_ref * (f / f0) ** beta
                acc[i] = np.interp(r_f, r_grid, L1[i, :])
            # fringe coherence: variance of aligned profile (higher = better alignment)
            score += float(np.var(acc))
            n += 1
        score /= max(n, 1)
        scores.append((beta, score))
        if score > best["score"]:
            best = dict(beta=float(beta), score=float(score))
    return best, pd.DataFrame(scores, columns=["beta", "score"])


def range_from_fringe(I, freqs, r_grid, beta, r_search=None):
    """Given β, estimate range by matched fringe alignment vs candidate r_ref."""
    if r_search is None:
        r_search = r_grid
    L = np.log(I + 1e-20)
    L1 = L.copy()
    for j in range(L.shape[1]):
        c = np.polyfit(freqs, L[:, j], 2)
        L1[:, j] = L[:, j] - np.polyval(c, freqs)
    f0 = float(np.mean(freqs))
    costs = []
    for r_ref in r_search:
        acc = np.zeros(len(freqs))
        for i, f in enumerate(freqs):
            r_f = r_ref * (f / f0) ** beta
            acc[i] = np.interp(r_f, r_grid, L1[i, :])
        # want high variance (stripes aligned) → minimize negative variance
        costs.append(-float(np.var(acc)))
    costs = np.asarray(costs)
    k = int(np.argmin(costs))
    return float(r_search[k]), costs


def run_rc3a_paper_like():
    """Paper-like WI/fringe reproduction on E-STD modal field, then metric summary."""
    mode = ensure_mode()
    freqs = F_FRINGE
    r_grid = R_GRID
    rows = []
    beta_rows = []

    cases = [
        dict(tag="PAPER_LIKE_IDEAL", r_km=50.0, z=200.0, note="E-STD modal field, matched env, noiseless — algorithm fidelity"),
        dict(tag="OFFGRID", r_km=52.3, z=213.0, note="off-grid truth"),
        dict(tag="S2_LINES", r_km=50.0, z=200.0, note="fringes sampled on S2 line freqs only", freqs=np.array(S2_LINES, float)),
    ]
    for case in cases:
        f_use = case.get("freqs", freqs)
        r_true = case["r_km"] * 1e3
        I = intensity_map(mode, f_use, r_grid, Z_S_TRUE, case["z"])
        best, beta_df = fringe_beta_estimate(I, f_use, r_grid)
        beta_df.insert(0, "case", case["tag"])
        beta_rows.append(beta_df)
        r_hat, costs = range_from_fringe(I, f_use, r_grid, best["beta"], r_grid)
        # mainlobe width of -cost
        c = -costs
        cmin = c.min()
        idx = np.where(c <= cmin + 0.05 * (c.max() - cmin + 1e-15))[0]
        wr = float(r_grid[idx.max()] - r_grid[idx.min()]) / 1e3 if len(idx) >= 2 else 0.0
        rows.append({
            "case": case["tag"],
            "note": case["note"],
            "r_true_km": case["r_km"],
            "z_true_m": case["z"],
            "beta_hat": best["beta"],
            "beta_score": best["score"],
            "r_hat_km": r_hat / 1e3,
            "err_r_km": r_hat / 1e3 - case["r_km"],
            "fringe_align_width_km": wr,
            "n_freqs": int(len(f_use)),
            "method": "WI-like fringe alignment after range-frequency de-trend (paper-like beta search + range warp)",
        })
        print(f"  RC3A {case['tag']}: beta={best['beta']:.2f} r_hat={r_hat/1e3:.2f} "
              f"err={r_hat/1e3-case['r_km']:.2f}km width={wr:.2f}km", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "RC3A_REPRO_RESULTS.csv", index=False, encoding="utf-8-sig")
    pd.concat(beta_rows, ignore_index=True).to_csv(OUT / "RC3A_beta_scores.csv", index=False, encoding="utf-8-sig")

    # save one I(r,f) and fringe figure for paper-like case
    I = intensity_map(mode, F_FRINGE, R_GRID, Z_S_TRUE, 200.0)
    write_fringe_svg(FIG / "RC3A_range_frequency_intensity.svg", I, R_GRID, F_FRINGE, df.iloc[0].to_dict())
    return df


def write_fringe_svg(path, I, R, F, meta):
    w, h = 780, 420
    ml, mr, mt, mb = 70, 100, 48, 56
    pw, ph = w - ml - mr, h - mt - mb
    # log I normalized
    L = np.log(I + 1e-20)
    L = (L - L.min()) / (L.max() - L.min() + 1e-15)
    nf, nr = L.shape
    cw, ch = pw / nr, ph / nf
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="28" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600">RC3-A 干涉条纹：I(r,f)（E-STD模态场）</text>']
    for i in range(nf):
        for j in range(nr):
            u = L[i, j]
            r = int(240 * (1 - u) + 31 * u)
            g = int(230 * (1 - u) + 111 * u)
            b = int(220 * (1 - u) + 139 * u)
            x = ml + j * cw
            y = mt + (nf - 1 - i) * ch
            p.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(cw,0.4):.1f}" height="{max(ch,0.4):.1f}" fill="rgb({r},{g},{b})"/>')
    p.append(f'<text x="{ml+pw/2}" y="{h-20}" text-anchor="middle" font-size="12" font-family="sans-serif">r (km)</text>')
    p.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90 16 {mt+ph/2})">f (Hz)</text>')
    p.append(f'<text x="{w-mr+10}" y="{mt+20}" font-size="11" font-family="sans-serif">β̂={meta.get("beta_hat", float("nan")):.2f}</text>')
    p.append(f'<text x="{w-mr+10}" y="{mt+40}" font-size="11" font-family="sans-serif">r̂={meta.get("r_hat_km", float("nan")):.2f} km</text>')
    p.append(f'<text x="{w-mr+10}" y="{mt+60}" font-size="11" font-family="sans-serif">err={meta.get("err_r_km", float("nan")):.2f} km</text>')
    p.append('</svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def write_reports(mfp_df, rc3a_df):
    mat = pd.DataFrame(METHODS)
    mat.to_csv(OUT / "RC3_LITERATURE_METHOD_MATRIX.csv", index=False, encoding="utf-8-sig")

    gap_rows = []
    for m in METHODS:
        gap_rows.append({
            "method_id": m["method_id"],
            "selected_mainline": m["selected"],
            "project_impl_degree": m["project_impl_degree"],
            "gap_vs_paper": "P4.5 simplified CZ profile ≠ paper algorithm" if m["method_id"] in ("WI_FRINGE", "MFP", "MMAC", "HLA_MODAL_SA") else "not started / out of scope",
            "next_validation_step": {
                "MFP": "R3-0 this round: paper-like MFP ambiguity on E-STD",
                "WI_FRINGE": "R3-0 this round: paper-like fringe/β reproduction; then E-STD+RC2 hard pairs",
                "MMAC": "R3-C later: angle+delay joint after HLA AOA",
                "HLA_MODAL_SA": "R3-B later: k-spectrum/synthetic aperture depth",
            }.get(m["method_id"], "defer"),
        })
    pd.DataFrame(gap_rows).to_csv(OUT / "RC3_IMPLEMENTATION_GAP.csv", index=False, encoding="utf-8-sig")

    # Selection markdown
    sel = []
    sel.append("# RC3 方法选择冻结（R3-0）")
    sel.append("")
    sel.append(f"UTC：{NOW}")
    sel.append("")
    sel.append("## 项目重锚定")
    sel.append("")
    sel.append("$$RC1\\ 环境/会聚区先验 \\rightarrow RC2\\ 运动学候选 \\rightarrow RC3\\ 传播观测约束r/z \\rightarrow 连续时间五维更新$$")
    sel.append("")
    sel.append("- **P4.5 仅作 RC3 简化基线估计器，不是 RC3 方法能力上限。**")
    sel.append("- 项目进度：RC1 基础形成；RC2 基线与困难候选基本完成；**RC3 正式方法验证刚进入主体**；五维综合性能此前结果不作最终结论。")
    sel.append("- 整体理论验证大约 **三分之一量级**（按科学问题而非实验次数）。")
    sel.append("")
    sel.append("## 统一验证协议（每条 RC3 方法）")
    sel.append("")
    sel.append("1. **原方法复现**（观测量/公式/流程/论文条件）")
    sel.append("2. **场景迁移** E-STD（允许迁移失败，须记录原因）")
    sel.append("3. **RC2 困难候选** RC2 vs RC2+方法")
    sel.append("4. **参数 RMSE/覆盖率**（P4.5 口径）")
    sel.append("")
    sel.append("结论仅四类：**可用 / 有条件可用 / 与RC2冗余 / 当前场景不适用**")
    sel.append("")
    sel.append("## 正式冻结的主线")
    sel.append("")
    sel.append("| 角色 | 方法 | 目标 | 选择依据 |")
    sel.append("| --- | --- | --- | --- |")
    for f in FROZEN_MAINLINES:
        sel.append(f"| {f['role']} | **{f['name']}** (`{f['method_id']}`) | {f['targets']} | {f['why']} |")
    sel.append("")
    sel.append("后备（不进主线）：匹配模态（与MFP信息重叠）、阵不变量、Doppler TMA（已对照）、浅海PBPR（场景跨度过大）。")
    sel.append("")
    sel.append("## 完整方法矩阵")
    sel.append("")
    sel.append("见 `RC3_LITERATURE_METHOD_MATRIX.csv`（含文献、观测量、公式、阵型、频带、环境先验、原性能、E-STD迁移差异、实现程度）。")
    sel.append("")
    sel.append("## 本轮范围")
    sel.append("")
    sel.append("- 完成：方法矩阵、主线冻结、**MFP 上限**、**RC3-A 条纹/波导不变量首轮复现**")
    sel.append("- **不**自动进入 RC3-B / RC3-C / P5")
    (OUT / "RC3_METHOD_SELECTION.md").write_text("\n".join(sel), encoding="utf-8")

    # MFP report
    def fnum(x, nd=3):
        try:
            if x is None:
                return "n/a"
            v = float(x)
            if not np.isfinite(v):
                return "n/a"
            return f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    mr = []
    mr.append("# MFP 理论上限基准（R3-0）")
    mr.append("")
    mr.append(f"UTC：{NOW}")
    mr.append("")
    mr.append("使用与 P3 一致的 **E-STD 正常模态理论前向模型**（不使用旧 CZ 时间轮廓评分，不进 Bellhop）。")
    mr.append("未知源：对候选谱做归一化/幅度匹配，**不**假设绝对源幅度已知。")
    mr.append("首轮：**无噪、环境匹配**。")
    mr.append("")
    mr.append("真值：r=50 km, z=200 m（及两组离网对照）。搜索：r∈[45,60] km, z∈[150,250] m。")
    mr.append("")
    mr.append("## 结果摘要")
    mr.append("")
    mr.append("| 源 | 真值 | r̂ km | ẑ m | r主瓣宽 km | z主瓣宽 m | 全局最优? | r局部峰 | corr(r,z) |")
    mr.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for _, r in mfp_df.iterrows():
        mr.append(
            f"| {r['source']} | {r['truth_tag']} r={r['r_true_km']} z={r['z_true_m']} | "
            f"{fnum(r['r_hat_km'])} | {fnum(r['z_hat_m'],1)} | {fnum(r['r_mainlobe_km'])} | "
            f"{fnum(r['z_mainlobe_m'],1)} | {r['truth_is_global_opt']} | {r['n_local_minima_r']} | {fnum(r['corr_rz_basin'])} |"
        )
    mr.append("")
    mr.append("## 物理解读（上限含义）")
    mr.append("")
    s0 = mfp_df[(mfp_df["source"] == "S0") & (mfp_df["truth_tag"] == "ON_GRID")]
    if len(s0):
        s0 = s0.iloc[0]
        mr.append(f"- **S0 宽带**：r 主瓣宽约 **{fnum(s0['r_mainlobe_km'])} km**，z 主瓣宽约 **{fnum(s0['z_mainlobe_m'],1)} m**。")
        mr.append(f"- 真值是否全局最优：**{s0['truth_is_global_opt']}**；r 方向局部峰数：**{int(s0['n_local_minima_r'])}**。")
        if s0["r_mainlobe_km"] >= 4:
            mr.append("- **r 主瓣仍宽** → 即便理想匹配场，在 E-STD 该孔径/频带/深度几何下，**传播场本身对 50 km 附近距离的分辨有限**；RC3 工程方法只能逼近、不能超过该物理上限。")
        else:
            mr.append("- **r 主瓣较窄** → 完整声场理论上含较丰富距离信息；若条纹类方法差，说明算法空间仍在。")
        if s0["z_mainlobe_m"] >= 40:
            mr.append(f"- **z 主瓣约 {fnum(s0['z_mainlobe_m'],1)} m** → 深度在本场景 MFP 上限也偏宽，与 P4.5 深度 UNRESOLVED 一致，**不能**用简化 CZ 弱深度否定全部深度方法，但需 R3-B/C 专门验证。")
        else:
            mr.append(f"- z 主瓣约 {fnum(s0['z_mainlobe_m'],1)} m → 深度信息在完整场中存在，简化特征未利用充分。")
    s2 = mfp_df[(mfp_df["source"] == "S2") & (mfp_df["truth_tag"] == "ON_GRID")]
    if len(s2):
        s2 = s2.iloc[0]
        mr.append(f"- **S2 多线谱**：r 主瓣宽 **{fnum(s2['r_mainlobe_km'])} km**，z 主瓣宽 **{fnum(s2['z_mainlobe_m'],1)} m**（频点少于 S0）。")
    mr.append("")
    mr.append("数据：`MFP_UPPER_BOUND_rz.csv`（S0 模糊面）、`MFP_UPPER_BOUND_summary.csv`；图：`figures/MFP_ambiguity_S0.svg`。")
    mr.append("")
    mr.append("**MFP 不是最终工程算法**，只定义当前标准场景完整传播场的 r/z 信息上限。")
    (OUT / "MFP_UPPER_BOUND_REPORT.md").write_text("\n".join(mr), encoding="utf-8")

    # RC3-A paper reproduction
    ra = []
    ra.append("# RC3-A 首轮复现：会聚区干涉条纹 / 波导不变量类距离估计")
    ra.append("")
    ra.append(f"UTC：{NOW}")
    ra.append("")
    ra.append("## 1. 选定文献方法线")
    ra.append("")
    ra.append("- 方法：**波导不变量 / 距离–频率干涉条纹**（Chuprov 波导不变量；深海/波导干涉条纹距离估计文献谱系）")
    ra.append("- 观测量：距离–频率强度 \(I(r,f)\)（论文级处理对象）")
    ra.append("- 目标参数：\(r,\\Delta r,\\beta\)")
    ra.append("- 核心关系：条纹斜率 \(dr/df\\approx \\beta\\, r/f\)；由斜率/条纹间距反演距离")
    ra.append("- **禁止**：替换成自造“单个干涉标量”")
    ra.append("")
    ra.append("## 2. 复现算法（paper-like）")
    ra.append("")
    ra.append("1. 由 E-STD 模态前向模型生成 \(I(r,f)\)（未知绝对源幅度：对强度取 log 后去趋势）")
    ra.append("2. 去除每频点距离均值与每距离的二次频率趋势（保留条纹）")
    ra.append("3. 沿 \(r(f)=r_{ref}(f/f_0)^{\\beta}\) 做条纹对齐搜索 \(\\beta\\in[0.3,2.5]\)（对应 \(\\int \\beta r/f\\,df\) 型条纹）")
    ra.append("4. 以对齐后剖面方差为条纹可见度评分，得 \(\\hat\\beta\)")
    ra.append("5. 固定 \(\\hat\\beta\)，搜索 \(r_{ref}\) 使条纹最对齐 → \(\\hat r\)")
    ra.append("")
    ra.append("## 3. 复现结果（`RC3A_REPRO_RESULTS.csv`）")
    ra.append("")
    ra.append("| case | r真值 km | β̂ | r̂ km | 误差 km | 条纹对齐宽 km |")
    ra.append("| --- | --- | --- | --- | --- | --- |")
    for _, r in rc3a_df.iterrows():
        ra.append(
            f"| {r['case']} | {r['r_true_km']} | {fnum(r['beta_hat'])} | {fnum(r['r_hat_km'])} | "
            f"{fnum(r['err_r_km'])} | {fnum(r['fringe_align_width_km'])} |"
        )
    ra.append("")
    ra.append("## 4. 与 MFP 上限对照（距离）")
    ra.append("")
    if len(mfp_df):
        m0 = mfp_df[(mfp_df["source"] == "S0") & (mfp_df["truth_tag"] == "ON_GRID")]
        if len(m0):
            ra.append(f"- MFP(S0) r 主瓣宽 ≈ **{fnum(m0.iloc[0]['r_mainlobe_km'])} km**")
        a0 = rc3a_df[rc3a_df["case"] == "PAPER_LIKE_IDEAL"]
        if len(a0):
            ra.append(f"- RC3-A 条纹对齐宽 ≈ **{fnum(a0.iloc[0]['fringe_align_width_km'])} km**，误差 ≈ **{fnum(a0.iloc[0]['err_r_km'])} km**")
    ra.append("")
    ra.append("## 5. 协议位置与下一步（本轮不执行）")
    ra.append("")
    ra.append("| 步骤 | 状态 |")
    ra.append("| --- | --- |")
    ra.append("| 1 原方法复现 | **本轮完成（首轮）** |")
    ra.append("| 2 E-STD 迁移诊断 | 部分（已在 E-STD 模态场上复现算法流程） |")
    ra.append("| 3 RC2 困难候选比较 | **未做（下轮 R3-A 深化）** |")
    ra.append("| 4 参数 RMSE/覆盖率 | **未做** |")
    ra.append("")
    ra.append("结论类型待 3–4 步完成后按：可用 / 有条件可用 / 与RC2冗余 / 当前场景不适用。")
    ra.append("")
    ra.append("**本轮停止**：不进入 RC3-B/C，不进入 P5。")
    (OUT / "RC3A_PAPER_REPRODUCTION.md").write_text("\n".join(ra), encoding="utf-8")

    # GPT sync
    gs = []
    gs.append("# R3-0 — GPT 同步稿")
    gs.append("")
    gs.append(f"- UTC: {NOW}")
    gs.append("- 重锚定：RC3 文献方法正式验证；P4.5≠RC3上限；不进入P5")
    gs.append("")
    gs.append("## 冻结主线")
    gs.append("")
    for f in FROZEN_MAINLINES:
        gs.append(f"- **{f['role']}**：{f['name']} → {f['targets']}")
    gs.append("")
    gs.append("## MFP 上限（摘录）")
    gs.append("")
    gs.append(mfp_df.to_string(index=False))
    gs.append("")
    gs.append("## RC3-A 条纹复现（摘录）")
    gs.append("")
    gs.append(rc3a_df.to_string(index=False))
    gs.append("")
    gs.append("## 停止")
    gs.append("")
    gs.append("本轮仅 R3-0 + RC3-A 首轮复现启动；**不**自动 R3-B/C，**不** P5。")
    gs.append("")
    (OUT / "R3_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    # decision json
    decision = {
        "package": "R3-0",
        "created_utc": NOW,
        "literature_matrix_rows": len(METHODS),
        "mainlines_frozen": [f["method_id"] for f in FROZEN_MAINLINES],
        "mfp_done": True,
        "rc3a_first_repro_done": True,
        "stop": "after R3-0; no auto R3-B/C; no P5",
        "p45_role": "simplified RC3 baseline estimator only, not method ceiling",
        "mfp_summary": mfp_df.to_dict(orient="records"),
        "rc3a_summary": rc3a_df.to_dict(orient="records"),
    }
    (OUT / "R3_0_decision.json").write_text(json.dumps(decision, ensure_ascii=False, indent=2, default=str), encoding="utf-8")


def main():
    t0 = time.time()
    print("=== R3-0 RC3 method recovery ===")
    print(f"OUT={OUT}")
    (OUT / "R3_0_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    print("[1] literature matrix + freeze ...", flush=True)
    print("[2] MFP upper bound ...", flush=True)
    mfp_df = run_mfp_upper_bound()
    print("[3] RC3-A paper-like fringe repro ...", flush=True)
    rc3a_df = run_rc3a_paper_like()
    print("[4] reports ...", flush=True)
    write_reports(mfp_df, rc3a_df)
    print(f"DONE in {time.time()-t0:.1f}s")
    print("files", [p.name for p in sorted(OUT.iterdir())])


if __name__ == "__main__":
    main()
