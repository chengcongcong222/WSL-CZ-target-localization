# MMAR_PRIMARY_SOURCE

UTC: 2026-09-24T06:44:02.063185+00:00

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
