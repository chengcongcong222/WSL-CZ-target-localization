# YANG2015_OBSERVABLE_CHAIN

UTC: 2026-09-23T08:46:33.049990+00:00

## 必须区分的 observable（主文未锁定前均不得实现）

| 候选对象 | 本轮状态 |
| --- | --- |
| complex beam output | **NOT_RECOVERED**（摘要出现 beam output 一词） |
| spectral amplitude | **NOT_RECOVERED** |
| spectral energy | **NOT_RECOVERED** |
| normalized modal coefficient | **NOT_RECOVERED** |
| modal shading coefficient | **NOT_RECOVERED** |

摘要仅说明 depth 从 **beam output** 估计，并用 nominal 环境的 **mode depth functions**。
深度匹配到底比较哪一个量：**NOT_RECOVERED**。

## b_m 分解（决定 KRAKEN mode table 如何变成 Yang observable）

必须明确 b_m 是否包含：

- φ_m(z_s)：**NOT_RECOVERED**
- 1/sqrt(k_m)：**NOT_RECOVERED**
- source level：**NOT_RECOVERED**
- attenuation：**NOT_RECOVERED**
- range spreading：**NOT_RECOVERED**
- other mode factors：**NOT_RECOVERED**

**禁止**继续把 `A_m = φ_m(z_s) φ_m(z_r)` 写成 Yang 正式公式。

## source-level nuisance

- 是否归一化：**NOT_RECOVERED**
- 是否只使用模态相对幅度：**NOT_RECOVERED**
- 是否存在比例因子解析消除：**NOT_RECOVERED**
- 是否通过相关/内积形成 ambiguity：**NOT_RECOVERED**

## Doppler / PLL（只恢复用途，不仿真）

| 项目 | 摘要级 | 公式级 |
| --- | --- | --- |
| Doppler → range increment Δr | 由数据估计 Doppler，已知原始频率 f0 | **NOT_RECOVERED** |
| PLL 输入 | 含随机+确定性相位的实数据相位 | **NOT_RECOVERED** |
| PLL 输出 | 去掉快变随机分量后的确定性 range-dependent 相位 | **NOT_RECOVERED** |

## 次级文献

- Liang et al. 2018, Match-Mode Autoregressive Method..., DOI 10.1155/2018/7824671：仅可交叉核验对 Yang 的描述与 SA 孔径限制；**本轮 PDF 403，未做公式对照**
- 不得用 AR/Hankel 公式冒充 Yang 2015 原公式
