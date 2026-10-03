# Controlling staged execution instruction

可以。现在让 Codex 正式执行，但要采用**阶段式指挥**，不能给它一句“按附件完成”后让它自由发挥。尤其是 SHGO `maxev` 已确认不是严格调用上限，因此第一关就是把**独立 objective 计数与硬封顶**做对。

把下面整段原样交给 Codex：

你现在接手 `chengcongcong222/WSL-CZ-target-localization` 的下一阶段研究执行。

# 0. 当前冻结状态

正式基线：

`7ab24845e6e1551b75287fefb1ab662e92b395b8`

必须首先核对：

```bash
git status
git rev-parse HEAD
git rev-parse origin/main
```

要求：

- `HEAD == origin/main == 7ab24845e6e1551b75287fefb1ab662e92b395b8`
- 工作区干净；
- 不允许改写历史；
- R3、原始 A1、FIX1、FIX2 数值证据保持冻结。

当前科学状态：

`A1_FIX2_BLOCKED_BY_UNCLOSED_ACOUSTIC_COVERAGE`

并继续保持：

- `R4 = 0%`
- A2 不开放；
- depth/B 不开放；
- SSP 不开放；
- P5 不开放。

本轮唯一允许执行的阶段：

`A1_SEARCH_TRACTABILITY_AUDIT`

---

# 1. 本阶段研究问题

FIX2 已确认：

1. matched acoustic basin 存在；
2. `J<0.001 dB` 的局部匹配截面可能极窄；
3. 主搜索与独立搜索命中了不同 basin；
4. 因此尚未闭合的是 **search coverage / tractability**，而不是简单的 forward-model 错误。

本轮只回答：

> 在不使用 truth、不降低 `J<0.001 dB` 门限、不事后增加预算的条件下，observation-only 的确定性全局搜索能否稳定恢复 matched basin，并表现出 raw budget convergence 和独立 solver agreement？

以下科学口径永久保持：

- likelihood section width ≠ optimizer capture basin；
- noiseless RC2 singleton 不算困难 noisy acoustic-search 强证据；
- finite candidate catalogs 未发现 alias ≠ global uniqueness；
- optimization failure ≠ physical non-identifiability。

---

# 2. 先创建阶段目录，但不要运行 noisy case

创建：

```text
results/R4_A1_SEARCH_TRACTABILITY_AUDIT/
```

先写：

```text
A1_SEARCH_TRACTABILITY_AUDIT_DESIGN.md
TRACTABILITY_METHOD.md
```

将现有任务定义中的科学 Gate、停止规则、禁止事项完整保存。

建议新代码隔离为：

```text
r4_a1_search_tractability.py
r4_a1_search_tractability_audit.py
tests/test_r4_a1_search_tractability.py
```

禁止修改：

- R3 scripts/results；
- A1/FIX1/FIX2 原始 evidence；
- FIX2 case results；
- 原 observations。

---

# 3. 搜索空间固定

继续使用已经接受的 FIX2 continuous RC2 bearing continuation。

每个 nuisance-depth branch 上，全局搜索坐标固定为：

\[
q=(r_0,u_r)
\]

其中：

- `r0` = initial range；
- `u_r` = radial target-motion coordinate。

其余水平状态必须通过 observation-only continuation 重建。

不能：

- 使用 truth coordinate；
- 使用 truth error；
- 使用 oracle basin center；
- 使用 FIX2 recovered state 做初始化；
- 使用 Hessian/SVD truth-local direction 做初始化；
- 返回粗四维 truth-neighborhood grid。

depth 仍是 inherited finite nuisance branch，不是新的 continuous-depth estimator。

---

# 4. 两个 solver 固定

Primary：

```python
scipy.optimize.shgo
```

Independent：

```python
scipy.optimize.direct
```

在 scientific run 前完成 API smoke test。

如果任一 solver 无法在冻结 objective 上正确工作：

```text
A1_TRACTABILITY_ENVIRONMENT_BLOCKED
```

提交证据并停止。

禁止替换为：

- DE；
- Sobol/Nelder-Mead；
- CMA-ES；
- Bayesian optimization；
- 其它随机搜索。

---

# 5. 最重要修正：不能使用 SHGO maxev 作为严格预算

本机已经实测：

> SHGO 的 `maxev` 不是严格 objective-call 上限。

因此 `maxev` 只能作为内部 solver option，**不能作为 T1/T2/T3 实际预算证据**。

必须实现 solver 外部独立预算控制器。

至少记录四个计数：

```text
n_objective_requests
n_exact_forward_evaluations
n_unique_states_evaluated
n_cache_hits
```

并记录：

```text
hard_budget_limit
termination_reason
hard_budget_triggered
```

## 硬预算定义

统一冻结：

> `hard_budget_limit` 作用于 objective 请求次数 `n_objective_requests`。

每次 solver 请求 objective 时：

1. 先增加 request counter；
2. 检查是否超过 hard cap；
3. 未超过才允许进入 cache / exact forward；
4. cache hit 仍然计为一次 solver objective request；
5. 另行记录真正执行 direct-modal forward 的次数。

这样：

- optimizer 行为预算可比较；
- computational reuse 又可以单独审计。

不得把 cache hit 从搜索预算中隐去。

---

# 6. Hard-cap 实现必须先单测

实现例如：

```python
class ObjectiveBudgetExceeded(RuntimeError):
    pass
```

objective wrapper 的原则：

```python
def __call__(self, x):
    if self.n_objective_requests >= self.limit:
        self.hard_budget_triggered = True
        raise ObjectiveBudgetExceeded(...)
    self.n_objective_requests += 1

    # then cache/exact evaluation
```

注意边界必须明确：

若 limit = N：

- 第 1...N 次请求允许；
- 第 N+1 次请求触发；
- 最终 `n_objective_requests == N`。

必须分别测试 SHGO 和 DIRECT 对自定义异常的传播行为。

如果某 solver：

- 吞掉异常；
- 无限继续调用；
- 无法保证硬 cap；

不要设计隐蔽 workaround。

直接判：

```text
HARD_BUDGET_ENFORCEMENT_FAILED
```

该 solver 不允许进入 scientific run。

在此问题解决前，不运行 21-case development。

---

# 7. Exact objective

Gate-bearing score 必须最终来自：

```python
direct modal forward model
```

不能使用 spline-only score 决策。

可以：

- cache exact objective；
- vectorize；
- spline 作 proposal/debug。

但每个 exported candidate 的：

- `J_exact`
- ranking
- recovered

必须由 direct modal calculation 重建。

要求测试：

1. truth-state exact self-match；
2. repeated same state exact determinism；
3. cached vs uncached exact identity；
4. geometry reconstruction；
5. bearing feasibility；
6. depth branch mapping；
7. state bounds；
8. truth-information exclusion。

---

# 8. 先做静态和运行时 truth-leak 审计

Gate-bearing search functions 的 signature/source 中不得出现：

```text
truth
errors
neighborhood
oracle
basin_width
basin_direction
```

如果某名称只是诊断输出，也必须确认不进入 solver path。

增加 runtime sentinel：

给 truth/evaluation helper monkeypatch 成一旦搜索阶段调用就抛异常。

至少覆盖：

- SHGO；
- DIRECT；
- candidate extraction；
- exact verification。

---

# 9. 预算只能根据运行成本预冻结

先运行：

- synthetic objective API smoke；
- exact objective unit test；
- 固定的 objective-cost calibration。

**不得查看任何 development recovery result。**

然后确定 SHGO 和 DIRECT 各自的：

```text
T1
T2
T3
```

允许两个 solver 的 numerical hard caps 不相同，因为算法结构不同。

但必须：

```text
T1 < T2 < T3
```

而且预算理由只能来自：

- objective 单次成本；
- solver 调用模式；
- 可接受执行规模。

不能来自：

- “T1 没找到，所以再加”；
- “P09 很难，所以给更多”；
- “到多少刚好 21/21”。

写入：

```text
TRACTABILITY_BUDGET_FREEZE.json
```

至少包括：

```json
{
  "baseline_commit": "...",
  "recovery_threshold_db": 0.001,
  "shgo": {
    "T1": "...",
    "T2": "...",
    "T3": "...",
    "solver_options": "..."
  },
  "direct": {
    "T1": "...",
    "T2": "...",
    "T3": "...",
    "solver_options": "..."
  },
  "hard_budget_definition": "n_objective_requests",
  "dedup_tolerances": "...",
  "state_agreement_tolerances": "...",
  "depth_agreement_rule": "...",
  "budget_choice_basis": "runtime/API calibration only"
}
```

---

# 10. 在第一个 noisy development case 之前做 PRE-RUN FREEZE

冻结：

- design；
- implementation；
- tests；
- budget JSON；
- solver configuration；
- threshold；
- tolerances。

生成 SHA256 manifest。

然后执行：

```bash
pytest ...
```

要求所有新测试通过。

建议此时做一个独立 git commit：

```text
R4 A1 tractability: freeze pre-run deterministic search design
```

如果有 GitHub 写权限，先 push。

这个 commit 必须发生在第一个 noisy development result 产生之前。

---

# 11. Development 固定为旧 21 个 nominal cases

只能使用 FIX2：

```text
SEARCH_BUDGET_CONVERGENCE.csv
```

中的 21 个：

```text
sigma_deg = 0.1
```

cases。

所有历史：

- P；
- Q；
- FIX2 fresh panel；

对于本方法都属于 development/regression。

无噪 case 可以跑 regression，但：

**不能计入 tractability Gate。**

---

# 12. Development 执行顺序

为了避免一个 solver 结果影响另一个：

## 先 SHGO

全部 21 cases：

```text
T1
T2
T3
```

每个 budget 必须 raw 独立执行。

禁止：

- T2 使用 T1 candidate 初始化；
- T3 使用 T2 candidate 初始化；
- cumulative retention 影响 raw result。

导出：

```text
DEVELOPMENT_SHGO_RESULTS.csv
SHGO_MINIMA_CATALOG.csv
SHGO_BUDGET_LOG.csv
```

## 再 DIRECT

同样：

```text
T1
T2
T3
```

独立执行。

不能读 SHGO minima/state。

导出：

```text
DEVELOPMENT_DIRECT_RESULTS.csv
DIRECT_MINIMA_CATALOG.csv
DIRECT_BUDGET_LOG.csv
```

---

# 13. 每次运行必须保存这些字段

至少：

```text
case_id
solver
budget
depth_branch
best_exact_J
recovered
r_km
theta_deg
v_mps
psi_deg
z_star_label_m
bearing_cost
bearing_cutoff
n_objective_requests
n_exact_forward_evaluations
n_unique_states_evaluated
n_cache_hits
hard_budget_limit
hard_budget_triggered
termination_reason
solver_success
solver_message
elapsed_seconds
```

所有失败项保留。

---

# 14. Raw budget convergence 的硬 Gate

SHGO 和 DIRECT **分别**判断。

必须同时：

### A

raw T2：

```text
21/21 recovered
```

### B

raw T3：

```text
21/21 recovered
```

### C

每个 case T2/T3 top state 在继承 FIX2 tolerance 内一致：

```text
[r, theta, v, psi]
=
[1e-5 km, 0.001 deg, 1e-4 m/s, 0.01 deg]
```

### D

nuisance-depth label 一致。

### E

T3 不得出现一个 T2 完全不存在的新 distinct：

```text
J < 0.001 dB
```

basin。

### F

结果不得依赖 cumulative lower-budget retention。

如果任一 solver 不满足：

```text
DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED = false
```

---

# 15. Development 失败后的强制 STOP

如果 development raw convergence 不通过：

最终 decision：

```text
A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED
```

然后：

1. 完整导出 failed cases；
2. 做 integrity reconstruction；
3. 写 report；
4. 更新 master ledger；
5. commit + push；
6. 停止。

**不得：**

- 增加 T4；
- 再加 maxiter；
- 扩大 budget；
- 改 tolerance；
- 换 solver；
- 生成 fresh confirmation；
- 创建 FIX3/FIX4。

---

# 16. Development convergence 通过后才比较 dual solver

T3 下要求：

```text
SHGO 21/21
DIRECT 21/21
```

然后独立比较 minima catalogs。

输出：

```text
DUAL_SOLVER_AGREEMENT.csv
```

检查：

- basin/state matching；
- exact J；
- depth branch；
- bearing feasibility。

如果：

- 一方恢复、一方不恢复；
- 两方均恢复但不同 basin；
- 存在无法匹配的 `J<0.001` distinct minima；

则：

```text
DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED = false
A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED
```

然后 STOP。

---

# 17. Alias 特殊 STOP

若两个独立 solver 找到 state-separated candidates，并且独立 direct-modal reconstruction 满足项目严格 exact joint match：

```text
bearing RMS < 1e-8 deg
acoustic RMS < 1e-6 dB
```

则进入：

```text
A1_CONTINUOUS_ALIAS_FINDING_STOP
```

这是科学结果。

不能再归因于 optimizer coverage。

如果没找到，只允许写：

> no unresolved exact alias found in tested solver catalogs

不得写：

> global uniqueness established

---

# 18. 只有 development 两个 Gate 都 PASS，才生成 fresh confirmation

顺序必须是：

```text
method/code/budgets/tolerances frozen
↓
freeze hash
↓
generate fresh panel
↓
freeze panel/design
↓
generate observations
```

新 panel：

```text
8 truths
```

每 truth：

```text
2 × sigma_theta = 0.1 deg
```

总计：

```text
16 noisy cases
```

必须使用此前从未使用的新 seed。

不能复用：

- 2026100203；
- 2026100204；
- 410001...；
- 420001...；
- 430001/430002。

选择一个新的固定 seed，并在 observation generation 前写入 freeze manifest。

---

# 19. Fresh Gate

仅运行 frozen raw T3。

要求：

```text
SHGO: 16/16
DIRECT: 16/16
```

且：

```text
dual-solver agreement = 16/16
```

不得：

- 看失败后增加预算；
- 删除 case；
- 改 panel；
- 改 threshold；
- 改 tolerance。

失败即：

```text
FRESH_CONFIRMATION_DUAL_SOLVER_RECOVERY_CONFIRMED = false
A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED
```

关闭 A1 search-repair chain。

---

# 20. 最终六个 Gate

必须严格输出：

```text
OBSERVATION_ONLY_SEARCH_VERIFIED
EXACT_OBJECTIVE_RECONSTRUCTED
DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED
DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED
FRESH_CONFIRMATION_DUAL_SOLVER_RECOVERY_CONFIRMED
NO_EXACT_ALIAS_FINDING_IN_TESTED_CATALOGS
```

注意：

如果 development 失败，fresh confirmation 没有执行，则第五个 Gate 不得伪造为 true。

应明确写：

```text
false / NOT_REACHED_DUE_TO_DEVELOPMENT_GATE
```

最终逻辑：

### 六项 PASS

```text
A1_SEARCH_TRACTABILITY_AUDIT_PASS
```

但：

```text
R4 = 0%
```

只表示后续可由研究负责人决定是否恢复 full A1 statistical experiment。

### 搜索/收敛/confirmation 任一失败

```text
A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED
```

关闭 A1 search-repair chain。

### Exact alias

```text
A1_CONTINUOUS_ALIAS_FINDING_STOP
```

### 环境无法实施

```text
A1_TRACTABILITY_ENVIRONMENT_BLOCKED
```

---

# 21. 必须产生的最终证据

目录：

```text
results/R4_A1_SEARCH_TRACTABILITY_AUDIT/
```

最终至少包含：

```text
A1_SEARCH_TRACTABILITY_AUDIT_DESIGN.md
TRACTABILITY_METHOD.md
TRACTABILITY_BUDGET_FREEZE.json
METHOD_FREEZE.json
DEVELOPMENT_SHGO_RESULTS.csv
DEVELOPMENT_DIRECT_RESULTS.csv
SHGO_MINIMA_CATALOG.csv
DIRECT_MINIMA_CATALOG.csv
DUAL_SOLVER_AGREEMENT.csv
ALIAS_DIAGNOSTIC.csv
INTEGRITY_AUDIT.csv
LOCAL_VALIDATION.md
A1_SEARCH_TRACTABILITY_REPORT.md
A1_SEARCH_TRACTABILITY_DECISION.json
GPT_SYNC.md
```

如果进入 fresh confirmation，再增加：

```text
FRESH_HOLDOUT_PANEL.csv
FRESH_HOLDOUT_DESIGN_FREEZE.json
FRESH_HOLDOUT_OBSERVATIONS.npz
FRESH_CONFIRMATION_RESULTS.csv
```

---

# 22. Integrity audit

最终必须独立重建：

- frozen hashes；
- budget freeze；
- objective count；
- exact scores；
- bearing cost；
- candidate ranking；
- recovery flag；
- T2/T3 state agreement；
- cross-solver agreement；
- alias flags；
- fresh panel generation seed/rule（若进入）。

**代码自己输出 PASS 不算审计。**

---

# 23. Git 提交要求

禁止覆盖 FIX2。

最终：

```bash
git status
git diff --check
pytest ...
git add ...
git commit -m "R4 A1: audit deterministic search tractability"
git push origin main
git rev-parse HEAD
git ls-remote origin refs/heads/main
```

核对：

```text
local HEAD == remote main
```

若没有写权限：

- 不伪称已 push；
- 输出本地 commit SHA；
- 明确 `REMOTE_PUSH_BLOCKED`。

---

# 24. 最终只向研究负责人汇报这些内容

```text
baseline SHA
pre-run freeze SHA
final SHA
remote/main SHA

decision

六个 Gate

SHGO:
T1 x/21
T2 x/21
T3 x/21

DIRECT:
T1 x/21
T2 x/21
T3 x/21

SHGO raw T2/T3 convergence x/21
DIRECT raw T2/T3 convergence x/21

dual solver T3 agreement x/21

fresh:
是否 reached
SHGO x/16
DIRECT x/16
agreement x/16

exact alias:
YES / NO_IN_TESTED_CATALOGS

hard budget:
是否严格执行
是否发生超限
objective-request counts

unit tests:
x passed

integrity/reconstruction:
x passed

failed cases:
[...]

R4:
0%

A2/depth/P5:
UNOPENED
```

不要继续执行任何下一阶段。

---

# 25. 第一个检查点

现在先只完成以下内容：

1. baseline/working-tree 核对；
2. 新阶段文件结构；
3. SHGO/DIRECT hard-budget wrapper；
4. hard-cap unit tests；
5. exact objective + truth-leak tests；
6. runtime calibration；
7. 冻结 T1/T2/T3；
8. 生成并提交 PRE-RUN FREEZE。

**到这里先停止并汇报。**

在我审计这个 PRE-RUN FREEZE 之前，禁止运行任何一个 21-case noisy development case。

这一点我建议比上一版更严格：**先让 Codex 只做到 PRE-RUN FREEZE，暂时不要让它一口气跑完。**原因是本轮最容易再次出现的问题已经不是物理模型，而是预算定义和 solver API 行为；尤其 SHGO 的硬封顶实现必须先独立审计。

Codex 第一轮完成后，把它返回的 **pre-run commit SHA、`TRACTABILITY_BUDGET_FREEZE.json`、hard-cap 测试结果和新增代码文件名**发给我。我会先审计这一关，再给它放行 21-case development 的下一条指令。

# Prior task definition (latest staged instruction takes precedence)

# R4-A1 Search Tractability Audit
## 冻结设计 + MIMO/Codex 执行任务

**基线 commit：** `7ab24845e6e1551b75287fefb1ab662e92b395b8`

**已接受状态：**
`A1_FIX2_BLOCKED_BY_UNCLOSED_ACOUSTIC_COVERAGE`

**管理状态：**
R4 = 0%；R3、原始 A1、FIX1、FIX2 证据全部冻结。
不得开放 A2、Depth/B、SSP、完整 bearing-sigma sweep 或 P5。

---

# 一、下一阶段唯一目标

阶段名：

`A1_SEARCH_TRACTABILITY_AUDIT`

FIX2 已经说明：

1. matched acoustic basin 确实存在；
2. 精确声学似然的局部盆地可以极窄；
3. 主搜索与独立搜索能命中不同的 basin；
4. 因而当前阻塞不是“是否有正确解”，而是：
   **observation-only 搜索是否能稳定、可审计地捕获正确 basin。**

本阶段不再通过 B4/B5 继续加密网格，也不以“把恢复率继续调高”为目标。

核心科学问题：

> 在不使用 truth、不降低 J<0.001 dB 门限、不事后调预算的条件下，
> 是否存在一套 observation-only、确定性、可复核的全局搜索结构，
> 能够对相关 matched acoustic basins 给出稳定的预算收敛和独立求解器一致性？

---

# 二、搜索空间

继续使用 FIX2 已接受的 RC2 observation-conditioned continuation。

全局搜索在每个 nuisance-depth branch 上只搜索：

`q = (r0, u_r)`

其中：

- `r0`：初始距离；
- `u_r`：FIX2 已使用的径向目标运动坐标。

其余水平状态由已接受的连续 bearing continuation 重构，并继续受到：

- 原连续 bearing likelihood；
- 物理状态边界；

约束。

深度仍然只是继承的有限 nuisance profile branch。

**禁止重新退回粗四维真值邻域网格。**

---

# 三、禁止信息流

所有 Gate-bearing estimator/search API 严禁读取或使用：

- truth coordinates；
- truth error；
- oracle basin width；
- Hessian/SVD truth-local direction；
- basin wall coordinate；
- FIX2 truth-near recovered state 作为初始化；
- 旧求解器 pass/fail label 作为初始化或分支选择依据。

oracle 信息只能在 observation-only 求解结束后：

- 计算 recovery；
- 做误差审计；
- 做 basin/alias 诊断。

必须同时做：

1. static source inspection；
2. runtime sentinel；

证明上述信息没有进入搜索器。

---

# 四、两套冻结的确定性全局求解器

在任何 noisy development case 运行前先检查环境：

## Solver A — primary

`scipy.optimize.shgo`

用途：

- deterministic simplicial global search；
- 在每个 nuisance-depth branch 上枚举 global/local minima；
- 最终候选必须再用 direct-modal exact objective 验证。

## Solver B — independent

`scipy.optimize.direct`

用途：

- DIRECT-style deterministic rectangle subdivision；
- 与 SHGO 相同物理域；
- 但不得读取 SHGO candidate/state；
- 不共享 optimizer state。

如果任一 API 不存在，或与当前 objective 无法兼容：

`A1_TRACTABILITY_ENVIRONMENT_BLOCKED`

然后停止。

**禁止静默替换为：**

- Differential Evolution；
- Sobol + Nelder-Mead；
- CMA-ES；
- Bayesian optimization；
- 其它随机优化器。

---

# 五、Gate-bearing objective

最终 Gate 使用：

**direct modal forward model**

不得使用 spline-only objective 作为最终判据。

允许：

- cache；
- batch/vectorization；
- spline 作为 proposal/debug；

前提是最终：

- exact J；
- exact ranking；
- exact recovery；

全部重新由 direct modal model 计算。

需要单测：

- exact truth self-match；
- duplicate call determinism；
- coordinate transformation；
- bearing feasibility；
- direct/saved candidate score reconstruction。

---

# 六、预算必须预冻结

允许先做：

- API smoke test；
- unit test；
- 与 recovery 无关的单次 objective runtime calibration。

但在第一个 noisy development case 之前，必须生成并冻结：

`TRACTABILITY_BUDGET_FREEZE.json`

至少包含：

- SHGO T1/T2/T3；
- DIRECT T1/T2/T3；
- n/evaluation caps；
- stopping tolerances；
- local solver options；
- basin/minimum dedup tolerance；
- state agreement tolerance；
- nuisance-depth agreement rule；
- library/internal seeds（如有）。

T1 < T2 < T3 必须单调增加。

预算大小只能根据：

- 算法结构；
- objective cost；
- 可运行性；

确定。

**不能根据某 case 是否恢复来调 T1/T2/T3。**

---

# 七、开发集

使用 FIX2 的同一组：

**21 个 nominal noisy development cases**

所有：

- P regression；
- Q；
- FIX2 fresh holdout；

从本方法开始都视为 development/regression material。

无噪 case 只用于 regression control，
**不计入 acoustic tractability Gate。**

每个：

- solver；
- budget；
- case；

至少导出：

- best exact J；
- recovered flag；
- exact state；
- nuisance depth；
- bearing cost/cutoff；
- n exact evaluations；
- termination status；
- deduplicated local/global minima catalog。

---

# 八、开发预算收敛 Gate

SHGO 与 DIRECT **分别**必须同时满足：

1. raw T2：21/21 recover；
2. raw T3：21/21 recover；
3. T2/T3 top recovered states 在 FIX2 继承 tolerance 内一致；
4. T2/T3 nuisance-depth 一致；
5. T3 不得出现 T2 完全未发现的新 distinct `J<0.001 dB` basin；
6. T2/T3 的成功不得依赖低预算候选累计保留。

这里判断的是：

**raw independent budget convergence**

而不是：

“之前某个预算曾经找到过，因此后面累计池里还在”。

如果 raw T2/T3 不闭合：

`DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED = false`

---

# 九、双求解器一致性 Gate

在 T3：

- SHGO 必须 21/21 recover；
- DIRECT 必须 21/21 recover。

然后比较两套独立 minima catalogs。

需要检查：

- matched basin 是否可一一匹配；
- state coordinate；
- nuisance depth；
- exact J；
- bearing feasibility。

任何 case 出现：

- 一方 recover、一方不 recover；
- 两方 recover 但落入不同 distinct basin；

均视为：

`DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED = false`

不能因为其中一个“更接近 truth”就判通过。

---

# 十、exact alias stop

如果两套 solver 独立发现：

- state-separated；
- bearing exact-match；
- acoustic exact-match；

的 distinct basin，并通过 direct modal independent reconstruction：

立即进入：

`A1_CONTINUOUS_ALIAS_FINDING_STOP`

这属于科学结果，不再是搜索器问题。

如果没发现，只能写：

> no unresolved exact alias found in tested solver catalogs

不能写：

> globally unique

---

# 十一、开发 Gate 失败时直接停止

如果开发阶段任一关键 Gate 失败：

`A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED`

立即停止。

特别是：

- 不再增加 T4；
- 不修改 solver；
- 不提高 maxiter；
- 不重新调 tolerance；
- 不新建 FIX3/FIX4；
- 不生成 fresh confirmation 去“救”开发失败。

这是本阶段最重要的停止规则。

---

# 十二、全新 confirmation

只有当 development 全部 PASS 后才能生成。

必须在：

- code；
- solver；
- T1/T2/T3；
- tolerance；
- objective；
- development workflow；

完全冻结之后，再生成：

**8 个全新 interior off-grid truths**

每个 truth：

- 两个 nominal `sigma_theta = 0.1 deg` realization；

共：

**16 noisy confirmation cases**

要求：

1. truth panel seed 冻结；
2. generation rule 冻结；
3. observation seeds 冻结；
4. panel hash 在 observation generation 之前冻结；
5. observation 生成后不得改 panel；
6. 不得根据结果再改 budget/tolerance。

noiseless 可以附带跑，
但只算 regression，不计 Gate。

---

# 十三、fresh Gate

在 frozen raw T3 下：

- SHGO：16/16；
- DIRECT：16/16；
- 两者 16/16 均需 state/minima-set agreement。

否则：

`FRESH_CONFIRMATION_DUAL_SOLVER_RECOVERY_CONFIRMED = false`

阶段结论：

`A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED`

并关闭 A1 search-repair chain。

---

# 十四、最终六个 Gate

必须输出：

1. `OBSERVATION_ONLY_SEARCH_VERIFIED`
2. `EXACT_OBJECTIVE_RECONSTRUCTED`
3. `DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED`
4. `DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED`
5. `FRESH_CONFIRMATION_DUAL_SOLVER_RECOVERY_CONFIRMED`
6. `NO_EXACT_ALIAS_FINDING_IN_TESTED_CATALOGS`

判定：

## 全部 true

`A1_SEARCH_TRACTABILITY_AUDIT_PASS`

注意：

PASS 只代表：

> 可以在研究负责人审计后重新进入完整 frozen A1 statistical experiment。

**R4 仍然是 0%。**

## search / convergence / confirmation 任一 false

`A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED`

然后：

> 关闭 A1 search-repair chain。
> 禁止靠加大 global-search budget 再建 FIX3/FIX4。

## exact alias 被证实

`A1_CONTINUOUS_ALIAS_FINDING_STOP`

## 必需 solver 环境不可用

`A1_TRACTABILITY_ENVIRONMENT_BLOCKED`

---

# 十五、要求的目录与文件

新建：

`results/R4_A1_SEARCH_TRACTABILITY_AUDIT/`

至少包含：

- `A1_SEARCH_TRACTABILITY_AUDIT_DESIGN.md`
- `TRACTABILITY_METHOD.md`
- `TRACTABILITY_BUDGET_FREEZE.json`
- `METHOD_FREEZE.json`
- `DEVELOPMENT_SHGO_RESULTS.csv`
- `DEVELOPMENT_DIRECT_RESULTS.csv`
- `SHGO_MINIMA_CATALOG.csv`
- `DIRECT_MINIMA_CATALOG.csv`
- `DUAL_SOLVER_AGREEMENT.csv`
- fresh panel/design freeze 文件
- `FRESH_CONFIRMATION_RESULTS.csv`
- `ALIAS_DIAGNOSTIC.csv`
- `INTEGRITY_AUDIT.csv`
- `LOCAL_VALIDATION.md`
- `A1_SEARCH_TRACTABILITY_REPORT.md`
- `A1_SEARCH_TRACTABILITY_DECISION.json`
- `GPT_SYNC.md`

所有失败 case、失败 solver run、termination failure 必须保留。

---

# 十六、MIMO / Codex 执行顺序

严格按以下顺序：

1. checkout / verify baseline `7ab24845e6e1551b75287fefb1ab662e92b395b8`；
2. 确认 R3、A1、FIX1、FIX2 frozen evidence 未改；
3. 创建上述阶段目录；
4. 保存本设计文件；
5. 检查 `scipy.optimize.shgo` 和 `scipy.optimize.direct`；
6. 若任一不可用：输出 environment-blocked，提交，停止；
7. 实现 observation-conditioned `(r0,u_r)` objective；
8. 加 unit tests；
9. 做与 recovery 无关的 runtime calibration；
10. **冻结 T1/T2/T3**；
11. 运行全部 21 development cases 的 SHGO T1/T2/T3；
12. 独立运行全部 21 development cases 的 DIRECT T1/T2/T3；
13. direct-modal 重建所有 top/minima candidate；
14. 判 raw T2/T3 convergence；
15. 判 dual-solver agreement；
16. 若 development FAIL：写最终 report/decision，提交，停止；
17. 若 development PASS：freeze method/code/config；
18. 生成全新的 8-truth confirmation panel；
19. freeze panel/design；
20. 生成 observations；
21. SHGO raw T3 跑 16 noisy cases；
22. DIRECT raw T3 跑 16 noisy cases；
23. 做 dual-solver confirmation + alias audit；
24. 生成 report / decision / GPT_SYNC；
25. 更新 `results/R4_MASTER/R4_PLAN.md`；
26. 更新 `results/R4_MASTER/R4_EVIDENCE_LEDGER.csv`；
27. commit + push `main`；
28. 核对 remote SHA；
29. 停止。

**不得在同一次执行里进入完整 A1、A2、B1、depth、SSP 或 P5。**

---

# 十七、最终回复给 GPT 的格式

完成后只需给出：

- commit SHA；
- remote/main SHA 核对；
- 最终 decision string；
- 六个 Gate boolean；
- SHGO T1/T2/T3 raw recovery；
- DIRECT T1/T2/T3 raw recovery；
- T2/T3 convergence count；
- SHGO vs DIRECT agreement count；
- fresh 16-case 结果（如实际进入）；
- exact alias 是否发现；
- unit test 数；
- reconstruction/integrity check 数；
- failed case 列表；
- 明确说明 R4 是否仍为 0%；
- 明确说明 A2/B/P5 是否仍未开放。

不要继续下一阶段。


# 2026-10-04 controlling execution-harness appendix

The research lead accepted the pre-run infrastructure at b4839b295777127ec0c8ade56b76c776db99148b, then authorized only DEVELOPMENT_EXECUTION_HARNESS_FREEZE, with zero noisy development execution. The complete new instruction is saved in DEVELOPMENT_EXECUTION_HARNESS_DESIGN.md. It supersedes earlier minimum/basin terminology for this audit: Gates are threshold-hit / discrete witness-cluster convergence, not local-stationarity certificates. WITNESS_SEMANTICS.md freezes the representative and comparison rules. The original task and earlier instructions above remain historical records. Budgets/core/thresholds/tolerances remain unchanged. The unreleased harness must be independently audited before any real case is run.
