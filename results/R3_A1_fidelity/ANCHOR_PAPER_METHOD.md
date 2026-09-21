# 锚定文献方法说明（R3-A1）

UTC：2026-09-20T07:26:06.389540+00:00

## 锚定选择

项目库内无单篇 PDF 原文时，采用**可复算的文献等价受控场景 PE-WI**，锚定在：

Chuprov waveguide-invariant line; practical range-frequency fringe ranging (e.g. waveguide-invariant passive ranging literature). Controlled reproduction uses standard two-mode interference: I(r,f)=|A1 e^{i k1 r}+A2 e^{i k2 r}|^2, beta_true = omega * d(k1-k2)/d(omega) / (k1-k2).

理由：任务书要求“明确 β 定义 + 条纹–距离关系 + 可复算场景”。两模态干涉是波导不变量文献中的标准受控模型，β_true 可解析写出，用于**方法保真**，不冒充深海实测。

## 原公式（受控场景）

- 观测量：\(I(r,f)=|A_1 e^{ik_1 r}+A_2 e^{ik_2 r}|^2\)
- 水平波数：\(k_i=\sqrt{(\omega/c_0)^2-\gamma_i^2}\)
- 波导不变量：\(\beta=\omega\,\dfrac{d(k_1-k_2)/d\omega}{k_1-k_2}\)
- 杼纹/距离：沿常相位 \((k_1-k_2)r\approx\mathrm{const}\) 的 \(dr/df\approx \beta r/f\)
- 场景参数：见 `R3_A1_CONFIG.json` → `pe_scenario`

## 复现结果（PAPER_REPRO.csv）

- β_true (band mean) = **1.0013** ± 0.0007
- β local-slope = **0.9594**；β align-search = **2.2000**
- r_true = **25.00 km**
- r_hat peak-match = **25.020 km**, 误差 = **0.020 km**
- r_hat fringe-count = **26.257 km**, 误差 = **1.257 km**
- r_hat power-law = **19.840 km**, 误差 = **-5.160 km**
- 评分峰宽（≠定位精度）= **0.493 km**
- pass_beta=True, pass_range=True

## 指标三分离（强制）

| 量 | 含义 | 能否当地位精度 |
| --- | --- | --- |
| fringe_score_width | 对齐评分峰宽 | **否** |
| range_error | 单次 \(|\hat r-r_{true}|\) | 仅单次 |
| MC RMSE | 多真值×多噪声 | 是（A1-PASS 后才做） |
