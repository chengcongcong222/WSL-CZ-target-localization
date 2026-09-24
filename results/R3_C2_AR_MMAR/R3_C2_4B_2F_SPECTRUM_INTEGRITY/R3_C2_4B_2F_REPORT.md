# R3-C2.4B-2F paper-spectrum integrity

UTC: 2026-09-24T16:02:27.087144+00:00

## 判定

### `C2_4B2_PAPER_SPECTRUM_PARTIAL`

R0 MMAR mode-order 优势保留（MMAR 1.000 vs SAB 0.071）；但 p=7 在 Eq.(28) 20 dB 下平均仅 1.2 峰，未复现 Fig.3/4 多峰结构。记 MMAR_R0_MODE_ORDER_ADVANTAGE_OBSERVED；不进 D(z)。

## 修正

1. Eq.(24) **k 升序**后再 DP；mode ID→KRAKEN index
2. Eq.(28) **Ps=|B(r₀)|²**（非全孔径均值）
3. `OBSERVABLE_MODE_TRUTH` 阈值 0.03 预冻结；PAPER_TRUE_PEAK_STRUCTURE_NOT_REPRODUCED
4. 完整谱 + 峰自身 `abs_P`；`SAB_FULL_SPECTRUM.csv` / `MMAR_FULL_SPECTRUM.csv`

## 指标（R0）

SAB mean rate **0.071**（排序修正后）；MMAR p=7 **1.000**

冻结：`MMAR_R0_MODE_ORDER_ADVANTAGE_OBSERVED`

## 停止

不进 D(z)/Fig.5–7/E-STD/MC/P5。
