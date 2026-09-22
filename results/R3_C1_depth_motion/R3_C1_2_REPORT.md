# R3-C1.2 报告：Zhu-QM 第一 CZ 有限适配终判

UTC：2026-09-22T02:32:48.390339+00:00

- **QF 直接迁移：CLOSED**（CZ 中 \(s=\sin\theta\) 几何变量失效，不改造成新方法）
- **QM matched intensity**：保留做本轮有限适配（HLA 波束 + 多频联合）
- 已知真实轨迹 = **理论上限**（RC2 误差只会更难）

## 四组（T=600 s 摘要）

| 组 | mean FWHM m | mean PSL | D 200/210 | D 200/220 | pair D>0.05 |
| --- | --- | --- | --- | --- | --- |
| A single + 235 Hz | 47.8 | 0.919 | 0.003 | 0.007 | 0.500 |
| B HLA 14m endfire beamformed + 235 Hz | 47.8 | 0.919 | 0.003 | 0.007 | 0.500 |
| C single + multi-line {201,235,283,338} Hz | 80.0 | 0.000 | 0.007 | 0.031 | 0.000 |
| D HLA 14m + multi-line | 80.0 | 0.000 | 0.000 | 0.002 | 0.000 |

## 深度对 D=1−|corr|（T=600, 组 A vs D）

| pair | Δz | A | D(HLA+多频) |
| --- | --- | --- | --- |
| 180/190 | 10 | 0.501 | 0.000 |
| 190/200 | 10 | 0.050 | 0.000 |
| 200/210 | 10 | 0.003 | 0.000 |
| 210/220 | 10 | 0.001 | 0.001 |
| 180/200 | 20 | 0.255 | 0.000 |
| 200/220 | 20 | 0.007 | 0.002 |

- 深度分辨下限：**None m**（10 m 对不过则记 >10 m；20 m 过则记 ~20 m）

## 终判（Zhu 路线此后永久关闭）

### `C1_QM_CZ_DEPTH_WEAK`

Four QM groups fail to separate 200±10/20 m clearly on first CZ. A: FWHM=47.8m D200/210=0.0031 D200/220=0.0068; D joint: FWHM=80.0m D200/210=0.0002 D200/220=0.0018; HLA help=False, multifreq help=False. Known-track ceiling already optimistic.

**下一步**：permanently close Zhu-QM route; proceed to Yang 2015 modal synthetic aperture next (not this round)

## 停止

- 不再 C1.x；不进 P5；本轮不进 Yang 2015
- Zhu（QF+QM）**PERMANENTLY_CLOSED**
