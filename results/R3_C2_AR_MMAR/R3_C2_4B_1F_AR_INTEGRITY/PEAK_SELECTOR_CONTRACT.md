# PEAK_SELECTOR_CONTRACT

UTC: 2026-09-24T10:18:32.049426+00:00

## 禁止

`choose peak closest to k_true` — **oracle peak selection 已删除**。

## Truth-independent 流程

1. 用固定 prominence=5% global max 检测 **全部** AR 谱峰（不知 k_true）
2. **冻结** 峰集合
3. 用冻结映射 k=sign·ω/Δr + 物理带分支 → k̂ 集合
4. `selected_peak` = **谱峰高度最大**的映射峰（与真值无关）
5. 真值**只做事后评分**：recall / precision / F1 / abs err / false_peak_count / peak_count_err

## p=1 SIGN_MAPPING_DIAGNOSTIC_ONLY

只用于冻结 ω↔k 符号；**不进**性能统计、**不是**新性能分支。

性能分支仍仅：

- `PRINTED_LITERAL_CONTROL` p=7
- `OUR_MOVING_SAMPLE_INTERPRETATION` p=133

## Eq.(19)

幂次 1 与 2 同时跑；报告峰位置是否一致（单调变换，不重复计成功样本）。
