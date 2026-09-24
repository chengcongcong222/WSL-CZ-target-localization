# AR_ORDER_PREREGISTRATION

UTC: 2026-09-24T09:59:41.308845+00:00
基线 commit：4c3dbc9ae539bbb2e732c107097baaa6489f2a62

在任何 AR 数值结果出现前冻结。**禁止**根据 D(z) 或 MODE_ORDER_RECOVERY_RATE 调整 p。

## 打印（PRIMARY）

Liang Eq.(18)：p often set to (2/3)(2L+1)；仿真 L=5 → 2L+1=11。
记号冲突：`AR_SAMPLE_COUNT_NOTATION_CONFLICT` / `PRIMARY_NOTATION_CONFLICT`（2L+1 同时为阵元数与 y[i] 长度）。

## 预注册两个解释分支（4B 都跑）

### A. `PRINTED_LITERAL_CONTROL`

```latex
p=\left\lfloor\frac{2}{3}(2L+1)\right\rfloor=\left\lfloor\frac{2\times 11}{3}\right\rfloor=7
```

仅“原文字面控制”，**不得**声称正确物理解释。

### B. `OUR_MOVING_SAMPLE_INTERPRETATION`

将 Eq.(17)/(18) 碰撞的 2L+1 解释为移动距离样本数 **N_r**：

```latex
p=\left\lfloor\frac{2}{3}N_r\right\rfloor
```

明确：**我们的实现解释**，不得反写成论文原文。

## 规则

- 两分支均可测试
- **禁止**按最终深度 / mode-order 结果再选 p
- 结果报告必须分列 A/B
