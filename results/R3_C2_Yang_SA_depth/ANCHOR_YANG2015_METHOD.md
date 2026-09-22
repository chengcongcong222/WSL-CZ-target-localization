# ANCHOR — Yang 2015 / Erratum 2018

UTC：2026-09-22T05:17:48.914782+00:00

## 文献

- T.C. Yang, *Source depth estimation based on synthetic aperture beamforming for a moving source*, JASA 138(3), 1678–1686 (2015). DOI 10.1121/1.4929748
- Erratum: JASA 144(6), 3075 (2018). DOI 10.1121/1.5081712

## 检索状态

- **PubMed / AIP / DOI 本轮网络访问失败**（cookie/transport）
- **Erratum 具体修正公式：NOT_RECOVERED**
- 因此 **不得编写/猜测 Yang 正式波束形成公式细节**

## 任务书已给出的物理量（本轮唯一使用）

- 模态场：\(p(r,f)=\sum_m A_m e^{jk_{rm}r}\)，\(A_m\propto\phi_m(z_s)\phi_m(z_r)\)
- 合成孔径：\(L=|r(t_{end})-r(t_{start})|\)（**径向相对距离变化**）
- 分辨尺度：\(\Delta k_{SA}\approx 2\pi/L\)，\(\eta_m=L|k_{r,m+1}-k_{r,m}|/(2\pi)\)
- 可控恢复：对 \(p(r)\) 做空间谱（FFT）比对预测 \(k_r\) 峰

## 任务书未覆盖、本轮禁止臆造

- Yang 正文 steering 相位约定、SA beam 具体核、Doppler–Δr 精确式、Erratum 符号修正
- PLL 实测细节

## f0 条件

- **KNOWN_F0 / KNOWN_TRUE_RANGE_INCREMENT**（理论上限）
