# YANG2015 Erratum — 已恢复（R3-C2.0R-FIX）

UTC：以 R3_C2_0R_FIX_REPORT.md 为准

## 状态：RECOVERED（2018 Erratum 关键修正）

来源：JASA/AIP Erratum 页面（DOI 10.1121/1.5081712，JASA 144(6), 3075, 2018）。

## 核心修正内容

引入初始距离偏移量：

δ = r₁ − r₁_data = r₂ − r₂_data

模态波数处谱峰应含额外相位：

g(k_m, z_r) = b_m φ_m(z_r) exp(i k_m δ)

（相对原论文 Eq.(5) 缺失 exp(i k_m δ) 的形式。）

## 冻结解释

1. exp(i k_m δ) **不改变** wavenumber spectral density；
2. **会改变** source-depth ambiguity function；
3. 原论文 simulation / data processing 相当于 δ = 0；
4. 实际数据必须知道或搜索 initial source range / offset δ；
5. 因此 Yang 方法**并非**严格“无需源距离”——与 RC2 距离候选 + RC3 收缩架构一致。

## 项目含义（正确链路）

RC2 运动/距离候选 → Δr(t) → k_m → (z, δ) 联合或条件搜索

## 未恢复部分

Yang 2015 主文若干正文公式细节可能仍记为 **YANG_MAIN_EQUATIONS_PARTIAL**（单独标注，**不得**再用 Erratum 卡住 normal-mode 基线验证）。

本轮**不**自行补齐未恢复的 2015 公式。
