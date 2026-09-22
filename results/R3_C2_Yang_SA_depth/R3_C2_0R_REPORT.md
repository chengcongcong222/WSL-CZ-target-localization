# R3-C2.0R 报告：文献恢复 + 物理模态基线

UTC：2026-09-22T06:33:52.618776+00:00

## 1. 文献

- Yang 框架已核实；**Erratum 修正式 NOT_RECOVERED**
- **判定 `C2_0_PAPER_RECOVERY_BLOCKED`**（原因缩小为公式/Erratum，而非“全文无证据”）

## 2. 旧 ModeModel 伪影

- 旧 φ_m(z) 仅依赖 γ 与高斯包络，**不随 f 重解** → C2.0 四频 D 完全相同是 **ARTIFACT_OF_SIMPLIFIED_MODEMODEL**
- 不得据此写“真实 CZ 模态深度对 200–220 m 很弱”

## 3. 物理 1D 本征模态（E-STD Munk, 5 km, φ(0)=0, φ′(H)=0）

| f | n_modes | mean residual | max orth off-diag |
| --- | --- | --- | --- |
| 201.0 | 400.0 | 0.000 | 0.000 |
| 235.0 | 400.0 | 0.000 | 0.000 |
| 283.0 | 400.0 | 0.000 | 0.000 |
| 338.0 | 400.0 | 0.000 | 0.000 |

| f | old n | phys n | old median Δk | phys median Δk |
| --- | --- | --- | --- | --- |
| 201.0 | 15 | 400 | 0.0000 | 0.0006 |
| 235.0 | 15 | 400 | 0.0000 | 0.0007 |
| 283.0 | 15 | 400 | 0.0000 | 0.0009 |
| 338.0 | 15 | 400 | 0.0000 | 0.0010 |

## 4. 孔径（Fourier/Rayleigh 预关，非高分辨极限）

| f | L_eta1 median (m) |
| --- | --- |
| 201.0 | 27995 |
| 235.0 | 35195 |
| 283.0 | 23108 |
| 338.0 | 21451 |

措辞冻结：**在简化/物理模态与普通 Fourier 准则下，2.4 km 径向孔径不足以按 Rayleigh 尺度分开相邻有效模式**；AR/高分辨波数方法不在本轮否定范围内。

## 5. 深度签名（真实 φ_m(z,f)）

| f | D(200/210) | D(200/220) |
| --- | --- | --- |
| 201.0 | 0.3814 | 0.2221 |
| 235.0 | 0.4169 | 0.4237 |
| 283.0 | 0.2858 | 0.4778 |
| 338.0 | 0.3618 | 0.6711 |

C1 338 Hz 解释判定：**PARTIALLY_EXPLAINED** — some f-dependence in physical mode depth vectors; check mode count/spacing columns

## 6. λ/2 受控恢复

- lam/2 显式记录：OK=0，NO_VALID_INTERIOR_PEAK/NYQUIST_EDGE_LIMITED=12（不静默缺行）

## 7. 停止

- 不进完整 Yang / Doppler / f0 / MC / P5
- **R3-C2.0R 完成；`C2_0_PAPER_RECOVERY_BLOCKED`**
