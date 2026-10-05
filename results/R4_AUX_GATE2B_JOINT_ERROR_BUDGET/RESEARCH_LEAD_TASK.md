这轮可以正式接受：

\[
\boxed{\texttt{AUX1\_NONIDEAL\_REQUIREMENTS\_ESTABLISHED\_CLIENT\_CONFIRMATION\_REQUIRED}}
\]

而且这一步比前两轮更接近真正的申报指标。

目前对两个候选架构，已经能说出**单因素、理想其他条件固定时**的 5% 距离 P95 要求：

| 条件 | Anchor A | Anchor B |
|---|---:|---:|
| 横向基线 | 5 km | 7 km |
| 随机定向方位 RMS | 0.05° | 0.075° |
| 共模系统偏差 | ≤±0.10° | ≤±0.10° |
| 差模半差 | ≤±0.05° | ≤±0.05° |
| 单节点每轴导航 \(1\sigma\) | ≤50 m | ≤50 m |
| 布放横向角偏差 | 至少测试到 ±20°仍通过 | 同左 |

但最后一句尤其重要：

\[
\boxed{\text{这些数现在绝对不能直接拼成一个联合工程指标。}}
\]

例如不能写：

> “0.05°随机方位 + 0.10°共模 + 0.05°差模 + 50 m导航 + 20°布放误差下仍≤5%”。

因为这一组合**从未验证过**。

另外，`±20°` 只是测试网格上限仍通过，所以正确表述是：

> 在已测 ±20° 范围内未观察到 5% Gate 失效。

不能说“容限就是20°”。

当前最大的现实缺口也已经非常明确：**真实辅助阵、DOA、姿态、导航、同步等能力全部 UNKNOWN。**

因此我建议再做最后一轮 geometry-only 数值 Gate：**联合误差预算 Gate-2B**。做完这一轮后，无论结果好坏，都应该停止纯仿真，向甲方/716团队确认硬件参数；不能继续无限算下去。

继续项目：

`chengcongcong222/WSL-CZ-target-localization`

当前正式基线：

`5d8b9c9f2d8efdb74a95da4ea546aae11457a24f`

研究负责人独立审计接受：

```text
AUX1_NONIDEAL_REQUIREMENTS_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED
```

Gate-2A 的准确证据为 **SINGLE_FACTOR ONLY**。

不得将其改写成联合误差能力。

当前：

```text
actual auxiliary bearing hardware capability = UNKNOWN
actual navigation / attitude capability = UNKNOWN

R4-A1-NEW = NOT_OPENED
R4 = 0%
```

本轮执行：

# `R4_AUX_GATE2B_JOINT_STATIC_ERROR_BUDGET`

这是在联系甲方确认硬件之前的**最后一个 geometry-only numerical Gate**。

完成后不得继续自动增加 geometry sweep。

---

# 1. 本轮科学问题

Gate-2A 分别证明：

```text
bearing systematic bias
navigation position error
deployment azimuth error
```

在单因素作用时都有一定 T5/T10 余量。

本轮只回答：

> 当这些非理想因素同时存在时，两个候选辅助节点架构还能保留多少联合误差预算？

不能使用：

```text
single-factor maximums
```

直接拼成工程 requirement。

必须真正联合验证。

---

# 2. 两个 architecture anchors 保持不变

## Anchor A

```text
baseline B = 5 km
random directed-bearing RMS = 0.05°
```

## Anchor B

```text
baseline B = 7 km
random directed-bearing RMS = 0.075°
```

继续：

```text
R = 50:1:60 km
both mirror deployments
```

不得新增第三个 anchor。

不得因为结果不好改 baseline/random RMS。

---

# 3. Gate-2A 的 T5 单因素预算

两个 anchor 均为：

```text
|b_common| <= 0.10°
|b_diff|   <= 0.05°
sigma_pos  <= 50 m per node per Cartesian axis
|beta|     <= 20°
```

其中：

```text
b_common = (b1+b2)/2
b_diff   = (b2-b1)/2
```

注意：

```text
b_diff=0.05°
```

表示两节点 systematic bias 的总差值可以达到：

```text
b2-b1 = 0.10°
```

不得把 half-difference误写成总差。

---

# 4. 不直接测试“全额单因素上限拼接”一个点

本轮采用固定 scaling：

```text
lambda =
0
0.25
0.50
0.75
1.00
```

对 T5 单因素预算同时缩放。

即：

```text
|b_common| = lambda × 0.10°
|b_diff|   = lambda × 0.05°
sigma_pos  = lambda × 50 m
|beta|     = lambda × 20°
```

例如：

## lambda = 0.50

```text
|b_common| = 0.05°
|b_diff|   = 0.025°
sigma_pos  = 25 m
|beta|     = 10°
```

这定义为：

```text
HALF_T5_JOINT_ENGINEERING_POINT
```

---

# 5. 所有 deterministic signs 必须全组合

对于每个：

```text
lambda > 0
```

完整测试：

```text
sign(b_common) = ±
sign(b_diff)   = ±
sign(beta)     = ±
```

即：

```text
2 × 2 × 2 = 8
```

个 sign combinations。

不得只用同号。

不得只汇报最好组合。

primary aggregate使用：

```text
worst sign combination
×
worst mirror
×
worst R
```

---

# 6. Position uncertainty

继续 Gate-2A 定义：

每个节点 position estimate：

```text
delta p_i ~ N(0, sigma_pos² I)
```

每 Cartesian axis 的 \(1\sigma\) 为：

```text
lambda × 50 m
```

MAIN 与 AUX 独立。

使用实际真实 AUX deployment：

```text
beta
```

生成 bearing。

Estimator 使用 noisy node positions。

不得把 deployment error 与 navigation error混为一项。

---

# 7. Bearing observation

保持：

```text
theta1_meas = theta1_true + eps1 + b1
theta2_meas = theta2_true + eps2 + b2
```

其中 random：

Anchor A：

```text
eps_i ~ N(0,0.05°²)
```

Anchor B：

```text
eps_i ~ N(0,0.075°²)
```

systematic：

```text
b1 = b_common - b_diff
b2 = b_common + b_diff
```

---

# 8. Monte Carlo

每：

```text
anchor
× lambda
× sign combination
× range
× mirror
```

使用至少：

```text
30,000 trials
```

与 Gate-2A 相同数量即可。

Commit A 前冻结一个新 seed。

共享 frozen standard-normal draws 可以用于 geometry comparisons。

必须保存：

- bearing draws；
- MAIN position draws；
- AUX position draws。

---

# 9. lambda=0

只需一个 deterministic sign identity。

它必须重现：

```text
Gate-2A zero-error control
```

至 MC sampling tolerance / exact shared-draw equivalence。

若不能重现：

```text
EXECUTION_INVALID
```

---

# 10. Primary joint metric

对每：

```text
anchor × lambda
```

定义：

```text
worst_joint_P95
```

为全部：

```text
11 ranges
× mirrors
× sign combinations
```

中最大：

```text
P95 relative range error
```

所有：

- parallel；
- nonfinite；
- behind-sensor；

继续计 infinite error。

不得删除 outliers。

---

# 11. 两个 Gate

## T5

```text
worst_joint_P95 <= 5%
```

## T10

```text
worst_joint_P95 <= 10%
```

同时保存：

```text
P99
failure rate
2D position P95
```

---

# 12. Joint-budget frontier

输出：

```text
AUX_JOINT_LAMBDA_FRONTIER.csv
```

至少：

```text
anchor
lambda
worst_joint_P95
worst_joint_P99
worst_range
worst_sign_common
worst_sign_diff
worst_sign_beta
worst_mirror
full_range_T5_pass
full_range_T10_pass
```

---

# 13. Primary frozen joint target

正式要求：

```text
lambda = 0.50
```

即 HALF_T5_JOINT_ENGINEERING_POINT。

### Gate-2B primary PASS

若至少一个 anchor 的：

```text
lambda=0.50
```

在全部 sign / range / mirror 下满足：

```text
P95 <=5%
```

则：

```text
AUX1_HALF_T5_JOINT_BUDGET_ESTABLISHED
```

这仍然不是 hardware admission。

---

# 14. 如果 lambda=0.50 失败

不得改 primary标准。

继续完整报告已预注册：

```text
lambda=0.25
lambda=0.75
lambda=1.0
```

并给：

```text
maximum tested lambda satisfying T5
maximum tested lambda satisfying T10
```

但科学判定：

```text
AUX1_HALF_T5_JOINT_BUDGET_NOT_ESTABLISHED
```

不得将0.25重新声明成 primary PASS。

---

# 15. T10 secondary joint budget

Gate-2A T10 单因素 allowance不同。

本轮额外做一个独立的 **T10 ladder**，不是把 T5 λ 结果冒充 T10 engineering budget。

## Anchor A T10 single-factor vector

```text
|b_common| = 0.10°
|b_diff|   = 0.10°
sigma_pos  = 100 m
|beta|     = 20°
```

## Anchor B

```text
|b_common| = 0.10°
|b_diff|   = 0.10°
sigma_pos  = 200 m
|beta|     = 20°
```

同样：

```text
lambda =
0
0.25
0.50
0.75
1.00
```

以及所有 signed combinations。

该 ladder 只判：

```text
P95 <=10%
```

输出：

```text
AUX_JOINT_T10_LAMBDA_FRONTIER.csv
```

---

# 16. 不要把 T5/T10 ladders 混起来

命名：

```text
T5_BUDGET_FAMILY
T10_BUDGET_FAMILY
```

一个 T5 lambda=1 point 和 T10 lambda=0.5 point 即使数值某些维度相同，也必须保存清晰 provenance。

---

# 17. 必须形成“联合工程 requirement”

如果 T5 half-budget通过，例如：

```text
Anchor A:
random RMS 0.05°
common bias <= ±0.05°
differential half-bias <= ±0.025°
navigation <=25 m/axis 1sigma
deployment <=±10°
```

正确表述：

> 在冻结联合误差模型、所有 signed bias/deployment combinations、50–60 km 11个距离点及双镜像下，P95 距离相对误差不超过5%。

不能写：

> 实际系统已经达到这些误差。

---

# 18. Hardware evidence仍然 UNKNOWN

不得修改 Gate-2A：

```text
actual auxiliary array/bearing capability = UNKNOWN
actual navigation/attitude capability = UNKNOWN
```

除非本轮新提供了真正的设备资料。

Repository 内已有资料不足，不能通过更多代码把 UNKNOWN变成 CONFIRMED。

---

# 19. 更新时间同步/target association状态

继续：

```text
NOT_NUMERICALLY_VALIDATED
```

本轮不做 target motion/time skew。

原因：

当前几何是 instantaneous snapshot。

这些是后续动态 tracking architecture 才进入的问题。

---

# 20. 本轮 route decision

## Outcome A

如果：

```text
T5 lambda=0.50
```

至少一个 anchor通过：

```text
AUX1_JOINT_STATIC_REQUIREMENT_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED
```

下一动作：

```text
STOP NUMERICAL ARCHITECTURE WORK
REQUEST CLIENT HARDWARE CONFIRMATION
```

**不得自动重开 A1。**

---

## Outcome B

如果 T5 half point失败，但 T10 half point至少一个 anchor通过：

```text
AUX1_JOINT_SUPPORTS_10PCT_NOT_5PCT
```

这非常重要，因为原项目目标本身是：

```text
多参数联合误差 <10%
```

由负责人决定是否接受 distance submetric 10%。

---

## Outcome C

如果连 T10 half point都失败：

```text
AUX1_STATIC_GEOMETRY_TOO_FRAGILE_FOR_PROJECT_TARGET
```

则应考虑：

```text
TDOA / external ranging / active-cooperative / target revision
```

而不是继续调 baseline。

---

# 21. 这轮之后必须止住纯几何仿真

无论 A/B/C：

```text
NO AUX Gate2C numerical run automatically
```

下一步必须是：

```text
CLIENT_REQUIREMENT_CONFIRMATION
```

或项目负责人改变架构。

这条 STOP 规则写进 decision JSON。

---

# 22. Client requirements更新

更新 Gate-2A 的：

```text
AUX_CLIENT_REQUIREMENTS.md
```

增加：

```text
joint validated requirement
```

与：

```text
single-factor maximum
```

两列。

必须明确二者不能互换。

建议形成：

| quantity | single-factor tested max | jointly validated requirement | client value | status |
|---|---:|---:|---:|---|

client value当前保持：

```text
UNKNOWN
```

---

# 23. 输出目录

```text
results/R4_AUX_GATE2B_JOINT_ERROR_BUDGET/
```

Commit A 至少：

```text
AUX_GATE2B_DESIGN_FREEZE.json
AUX_T5_JOINT_GRID.csv
AUX_T10_JOINT_GRID.csv
AUX_GATE2B_CRITERIA.md
```

Commit B 至少：

```text
AUX_T5_JOINT_RESULTS.csv
AUX_T10_JOINT_RESULTS.csv

AUX_JOINT_LAMBDA_FRONTIER.csv
AUX_JOINT_T10_LAMBDA_FRONTIER.csv

AUX_GATE2B_DECISION.json
AUX_GATE2B_REPORT.md
AUX_CLIENT_REQUIREMENTS_UPDATED.md

VALIDATION.json
GPT_SYNC.md
```

至少画：

```text
lambda vs worst P95 (T5 family)
lambda vs worst P95 (T10 family)
worst sign combination map
Anchor A vs B comparison
```

---

# 24. 两提交模式

## Commit A

任何 joint MC 前：

```text
R4 AUX Gate2B: freeze joint static error budget
```

push，确认 remote/main。

## Commit B

一次执行全部 frozen grids：

```text
R4 AUX Gate2B: evaluate joint static error budget
```

push后停止。

---

# 25. 独立验证

至少：

- frozen draw regeneration；
- full signed Cartesian combinations；
- 11 ranges；
- mirrors；
- no dropped failures；
- nearest-rank quantiles；
- Gate2A zero controls；
- independent slope-intersection implementation；
- correct common/differential bias algebra；
- correct position/noisy-estimator geometry；
- T5/T10 family provenance；
- policy pushed before MC；
- historical Gate0/1/2A artifacts unchanged。

---

# 26. 本轮禁止

禁止：

```text
acoustic propagation
depth
SSP
TDOA
tracker
moving-target time sequence
time-skew
target-association failures
new baseline or new random-bearing RMS
R4-A1 target-state runs
```

---

# 27. R4进度

无论结果：

```text
R4 = 0%
```

因为当前仍是 architecture-reset qualification。

R4-A1-NEW 只有：

1. joint static requirement成立；
2. client/hardware能力被确认；
3. time/association architecture可接受；

之后才由研究负责人正式打开。

---

# 28. 最终只返回

```text
parent AUX Gate-2A SHA:
5d8b9c9f2d8efdb74a95da4ea546aae11457a24f

Gate-2A independent audit:
ACCEPTED_SINGLE_FACTOR_REQUIREMENTS

AUX Gate-2B design SHA:
...

AUX Gate-2B execution SHA:
...

remote/main SHA:
...

T5 Anchor A:
lambda 0.25 -> ...
lambda 0.50 -> ...
lambda 0.75 -> ...
lambda 1.00 -> ...
maximum tested lambda passing 5% = ...

T5 Anchor B:
...

T10 Anchor A:
maximum tested lambda passing 10% = ...

T10 Anchor B:
...

HALF_T5 joint point:
Anchor A PASS / FAIL
Anchor B PASS / FAIL

HALF_T10 joint point:
Anchor A PASS / FAIL
Anchor B PASS / FAIL

worst joint failure mode:
...

Gate-2B decision:
...

actual hardware capability:
UNKNOWN / CONFIRMED

client confirmation required:
[...]

automatic further numerical work:
STOP

R4-A1-NEW:
NOT_OPENED

R4:
0%
```

提交后停止。

这轮如果通过，我们就不再缺“数学上的联合误差要求”了。之后真正缺的只剩一件事：**甲方能否提供满足这些要求的辅助阵、导航、姿态和同步能力。**这时继续仿真已经没有价值，应该拿 requirement sheet 去问设备条件。