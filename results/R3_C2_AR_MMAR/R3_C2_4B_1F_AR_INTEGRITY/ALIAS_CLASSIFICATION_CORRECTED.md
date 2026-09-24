# ALIAS_CLASSIFICATION_CORRECTED

UTC: 2026-09-24T10:18:32.049426+00:00

符号映射（p=1 diagnostic 冻结）：**k=-omega/dr**（sign=-1）

分类仅由几何+物理带决定（**不**看恢复误差）：

| 类 | 条件 |
| --- | --- |
| `NO_FOLD_DIRECT_RECOVERY` | kΔr < π |
| `FOLDED_BUT_BRANCH_RECOVERABLE` | kΔr > π 且物理带 [1.32,1.50] 内**唯一**分支 |
| `FOLDED_BRANCH_AMBIGUOUS` | 物理带内多分支 |

修复前把 no-fold 标成 folded 的错误已废除。

sign_rows / single 表见同目录 csv。
