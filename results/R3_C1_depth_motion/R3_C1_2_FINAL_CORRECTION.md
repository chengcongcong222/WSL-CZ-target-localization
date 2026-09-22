# R3-C1.2 FINAL CORRECTION 报告

UTC：2026-09-22T03:47:21.743561+00:00

## 1. Steering 符号修正

- \(p_m\sim e^{-jk_rx_m}\)（\(r_m=r-x_m\)）
- 正确：\(w_m=e^{-jk_0x_m}\)，则 \(w_m^*p_m\sim e^{+jk_0x_m}e^{-jk_rx_m}\) 相干
- 旧错误：\(w_m=e^{+jk_0x_m}\) 相位梯度加倍，B_steered 塌缩

## 2. 平面波单元测试

- pass = **True**
- out_correct=1, out_wrong=3.224e-05, out_bs=0.02251

## 3. B/D 重跑（A/C 冻结保留）

| 组 | T | D200/210 | D200/220 | mean FWHM |
| --- | --- | --- | --- | --- |
| A | 600 | 0.0031 | 0.0067 | 47.8 |
| A | 1200 | 0.0023 | 0.0062 | 64.2 |
| B | 600 | 0.0031 | 0.0067 | 47.8 |
| B | 1200 | 0.0023 | 0.0062 | 64.2 |
| C | 600 | 0.0012 | 0.0093 | 80.0 |
| C | 1200 | 0.0090 | 0.3401 | 76.6 |
| D | 600 | 0.0012 | 0.0093 | 80.0 |
| D | 1200 | 0.0090 | 0.3401 | 76.6 |

- HLA 相对 single 增量 T=600：ΔD220=**0.0000**，ΔFWHM=**0.0 m**
- T=1200 多频 20 m 信号：D220=**0.340**，采样稳定=**True**

## 4. 终判

### `C1_QM_MULTIFREQ_ONLY_CONDITIONAL`

HLA increment still weak (ΔD220_600=0.0000). Multifreq long-track 20 m signal holds under sampling (D220_T1200=0.340, conv_stable=True); 10 m still not separable (D210=0.0090). Frequency-selective depth cue (e.g. 338 Hz) noted — aligns with Yang modal-depth idea.

**下一步**：Zhu route PERMANENTLY_CLOSED; Yang 2015 modal synthetic aperture next (not this round)

Zhu **PERMANENTLY_CLOSED**（本轮后无 C1.x）。