# R3-C2.4B-1F AR implementation integrity

UTC: 2026-09-24T10:18:32.049426+00:00

## 判定

### `AR_K_RECOVERY_NUMERICALLY_UNSTABLE`

p=133 rank/cond unstable; single F1=0.444

## 修复

1. **真 modified covariance**：backward 为 y[n−p]=−Σa[k]y[n−p+k]（与 forward 独立）；`backward_dup_forward=False`
2. **oracle 选峰删除**：峰集先冻结；selected=最强峰；真值只评分
3. **符号冻结**（p=1 diagnostic）：**k=-omega/dr**
4. alias 类按 **kΔr 与物理带** 分类，不按误差改标签

## 指标（power=1，不与 power=2 重复计数）

| 项 | 值 |
| --- | --- |
| single F1 / recall / precision | 0.444 / 0.444 / 0.444 |
| single mean FP | 0.00 |
| multi complete / F1 | 0.67 / 0.843 |
| p=133 unstable | True |
| Eq24 solver | True |

N_r=200,p=133：2(N_r−p)=134 约束 / 133 参数。

## 停止

不进 4B-2 / Hankel / D(z) / FIELD / E-STD / MC / P5。
