# Controlling research-lead instruction: execution harness freeze

我已独立审计当前 PRE-RUN checkpoint。结论是：

\[
\boxed{\texttt{PRE\_RUN\_INFRASTRUCTURE\_FREEZE\_ACCEPTED}}
\]

当前远端 `main` 确为：

`b4839b295777127ec0c8ade56b76c776db99148b`

相对初始冻结 `5fa6309...` 只有一笔 CRLF/LF portability 修正；**核心搜索代码 `r4_a1_search_tractability.py` 未改变，预算文件 SHA 仍为 `40cd9c...`，noisy development=0。**

硬预算机制也可以接受：N 次允许、N+1 拒绝且拒绝前无 forward；cache hit 计 request；DIRECT 的 `SystemError -> ObjectiveBudgetExceeded` 因果链明确，并没有被伪装成 solver success。

不过目前还**不能直接放行 21-case**。原因不是 PRE-RUN 有错误，而是当前代码故意没有 development driver，而且 `candidate_catalog()` 明确标注为：

`EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM`

因此必须先冻结“**development 怎么定义 raw result、怎么比较 T2/T3、怎么比较 SHGO/DIRECT**”，否则后续很容易把低 J witness 静默升级成“局部极小值/盆地”。

下一条给 Codex：

你现在从已通过独立 PRE-RUN 基础设施审计的 commit 开始：

`b4839b295777127ec0c8ade56b76c776db99148b`

本轮仍然 **禁止运行任何真实 21-case noisy development case**。

目标只有一个：

# `DEVELOPMENT_EXECUTION_HARNESS_FREEZE`

完成 development driver、结果语义和 Gate 逻辑冻结，然后提交并停止，等待第二次独立审计。

---

## 1. 当前已接受且禁止修改的内容

以下内容已经通过独立审计，全部视为冻结：

- `r4_a1_search_tractability.py` 中的 `ExactObjective`
- `BudgetObjective`
- SHGO / DIRECT 两个 solver family
- `(r0,u_r)` 搜索坐标
- 21 个 inherited depth branches
- `J_exact < 0.001 dB`
- hard-cap 定义
- request/cache/exact-forward/unique-state counters
- T1/T2/T3 数值
- solver options
- state agreement tolerance
- dedup tolerance
- alias tolerance
- truth-information prohibition

尤其：

```text
SHGO   T1/T2/T3 = 16 / 64 / 256
DIRECT T1/T2/T3 = 16 / 64 / 256
```

单位继续是：

> admitted objective requests / case / depth branch

严禁修改 `TRACTABILITY_BUDGET_FREEZE.json`。

其 SHA 必须继续是：

`40cd9c0b2cbc7a6fced83faf97009c25d652223bf7cbb77e7ebc6a981089672d`

如果任何工作导致这个 SHA 改变：

**立即停止，不得运行 development。**

---

# 2. 必须先修正“minimum / basin”术语

当前：

```text
candidate_catalog()
```

产生的是：

```text
EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM
```

这个定义是正确的，必须保持。

因此在后续 development 中：

**禁止把 wrapper/cache 中的低 J evaluated point 写成：**

- local minimum
- global minimum
- basin minimum
- recovered basin

除非确有额外数学证据证明 stationarity/minimum。

本轮不新增额外 local-polish 算法，也不增加新的 evaluation budget。

统一术语：

### `raw witness`

solver 在冻结 hard budget 内实际 evaluation 过、且由 direct-modal exact objective 重建的可行状态。

### `raw top witness`

某 solver / case / budget 的全部 21 depth branches 中：

```text
J_exact 最小的 feasible evaluated witness
```

### `recovered`

严格定义：

```text
raw_top_witness.J_exact < 0.001 dB
```

这只表示：

> 该 capped observation-only search 实际命中了 threshold region。

不表示：

> 找到了局部极小值。

---

# 3. 新增 artifact 名称必须反映这一点

新增：

```text
DEVELOPMENT_SHGO_RESULTS.csv
DEVELOPMENT_DIRECT_RESULTS.csv

SHGO_WITNESS_CATALOG.csv
DIRECT_WITNESS_CATALOG.csv

SHGO_BUDGET_LOG.csv
DIRECT_BUDGET_LOG.csv
```

原设计中的：

```text
SHGO_MINIMA_CATALOG.csv
DIRECT_MINIMA_CATALOG.csv
```

本轮禁止伪造。

如果为了兼容目录规范必须保留这两个文件，则：

- 文件可以存在；
- 必须有 `candidate_kind`；
- 未获得 native certified minimum 时不得填入 witness；
- 可以为空；
- scientific Gate 不依赖它们。

同时在 `TRACTABILITY_METHOD.md` 和设计附录中声明：

> 本 tractability audit 的 Gate 是 threshold-hit / witness convergence，不是 local-stationarity certificate。

这是在 **0 个 noisy development result 之前**进行的方法语义澄清，不是结果后改 Gate。

---

# 4. Development case manifest 必须先隔离 truth 字段

不要让 scientific driver 直接拿完整：

- truth panel；
- truth error table；
- FIX2 recovery/error columns；

作为 runtime dataframe。

先从已有冻结 evidence 构造：

```text
DEVELOPMENT_CASE_MANIFEST.csv
```

只允许包含搜索实际需要的非 oracle 信息，例如：

```text
case_id
panel_id                  # 仅 identifier，可选
sigma_deg
seed
bearing_cutoff
observation_index
origin_group              # 仅 provenance，可选
```

不得包含：

```text
truth r/theta/v/psi
z_true
top1 error
survivor error
old recovered flag
old J result
old solver state
oracle basin information
```

必须通过 `CASE_OBSERVATIONS.npz` 的 `case_ids` 显式映射 observation，不得假定 CSV row order 静默一致。

最终 manifest 要 hash freeze。

---

# 5. Development 集固定为 exactly 21 noisy cases

只能选择：

```text
sigma_deg == 0.1
```

的 21 个 nominal development cases。

必须 assert：

```text
n_unique_case_id == 21
```

同时 assert：

- 没有 sigma=0 case；
- 没有 duplicate；
- observation 中 21 个 case 均唯一存在；
- 不读取 fresh future panel；
- 不读取 truth state。

无噪 regression 本轮 execution harness 不必运行。

---

# 6. Raw run 必须完全独立

每个：

```text
solver × case × budget × depth_branch
```

都重新创建：

- `ExactObjective`
- `BudgetObjective`
- solver instance/state
- cache

严格禁止：

- T2 复用 T1 cache；
- T3 复用 T2 cache；
- SHGO 读取 DIRECT witness；
- DIRECT 读取 SHGO witness；
- lower-budget cumulative candidate retention。

需要 unit test 用 object identity / sentinel 验证。

---

# 7. Hard-cap exhaustion 的科学解释固定

由于 hard cap 是本实验人为定义的 raw computational budget，因此：

```text
HARD_BUDGET_EXCEEDED
```

本身**不是 scientific failure**。

同样：

```text
solver_success = false
```

如果原因只是 frozen hard cap exhaustion，也不能自动令 recovery=false。

scientific recovery 只由：

```text
best feasible direct-modal raw witness J_exact < 0.001
```

决定。

但是以下属于 execution invalid：

- `HARD_BUDGET_ENFORCEMENT_FAILED`
- unrelated uncaught exception
- counter inconsistency
- objective reconstruction failure
- truth leak
- corrupted observation
- modified frozen budget/code

这种 case 必须：

```text
execution_valid = false
```

并使对应 scientific Gate 无法 PASS。

---

# 8. 每个 branch 保存完整 provenance

每个 solver/case/budget/depth branch 保存：

```text
case_id
solver
budget
depth_label_m

termination_reason
solver_success
execution_valid
solver_message
exception_chain

n_objective_requests
n_objective_attempts
n_blocked_requests
n_exact_forward_evaluations
n_unique_states_evaluated
n_cache_hits
hard_budget_limit
hard_budget_triggered

best_branch_J_exact
best_branch_r_km
best_branch_theta_deg
best_branch_v_mps
best_branch_psi_deg
best_branch_bearing_cost

elapsed_seconds
```

要求：

```text
n_objective_requests <= hard_budget_limit
```

且如果 wrapper hard-cap exhaustion：

```text
n_objective_requests == hard_budget_limit
n_blocked_requests == 1
```

---

# 9. Case-level raw result 固定

完成 21 个 depth branches 后，对同一：

```text
solver / case / budget
```

把所有 branch 的 feasible exact witnesses 合并。

选择：

```text
raw_top_witness = argmin J_exact
```

输出：

```text
case_id
solver
budget
best_exact_J
recovered
r_km
theta_deg
v_mps
psi_deg
z_star_label_m
bearing_cost
bearing_cutoff

n_depth_branches_completed
n_valid_branches
total_objective_requests
total_exact_forward_evaluations
total_unique_states
total_cache_hits

all_branches_execution_valid
```

必须：

```text
n_depth_branches_completed == 21
```

否则该 case execution invalid。

---

# 10. Witness catalog

保存所有可行 exact evaluated witness。

至少字段：

```text
case_id
solver
budget
depth_label_m
candidate_kind
J_exact
r_km
theta_deg
v_mps
psi_deg
bearing_cost
```

统一：

```text
candidate_kind =
EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM
```

对 catalog 做 frozen `DEDUP`：

```text
[0.005 km, 0.001 deg, 0.005 m/s, 0.05 deg]
```

depth label 必须相同才允许视为同一 cluster。

---

# 11. “basin” Gate 改为 subthreshold witness-cluster Gate

因为没有 stationarity certificate，因此旧句：

> T3 不得出现 T2 未发现的新 distinct J<0.001 basin

在代码和正式报告中必须改成：

> T3 不得出现 T2 未发现的新 distinct subthreshold witness cluster。

这里的 cluster 只是：

- `J_exact < 0.001`
- 同 depth label
- frozen DEDUP tolerance

形成的离散 witness cluster。

**不能称 basin。**

---

# 12. T2/T3 raw convergence 固定算法

对每个 solver 分别判断。

一个 case `raw_T2_T3_converged = true` 必须同时：

1. T2 `execution_valid=true`
2. T3 `execution_valid=true`
3. T2 `recovered=true`
4. T3 `recovered=true`
5. T2/T3 raw top witness 满足 STATE_AGREEMENT：

```text
r       <= 1e-5 km
theta   <= 0.001 deg
v       <= 1e-4 m/s
psi     <= 0.01 deg
```

6. depth label 完全一致
7. 每个 T3 subthreshold witness cluster 在 T2 中存在 matching cluster

第 7 项使用 frozen DEDUP tolerance，不使用 truth。

SHGO：

```text
21/21
```

DIRECT：

```text
21/21
```

才允许：

```text
DEVELOPMENT_RAW_BUDGET_CONVERGENCE_VALIDATED = true
```

---

# 13. Dual solver T3 agreement 固定算法

只有在两边 development raw convergence 均 21/21 后才评估。

每 case 必须：

1. SHGO T3 recovered
2. DIRECT T3 recovered
3. 两者 top witness STATE_AGREEMENT
4. depth label 一致
5. SHGO 每个 subthreshold witness cluster 可在 DIRECT catalog 中找到 matching cluster
6. DIRECT 每个 subthreshold witness cluster 可在 SHGO catalog 中找到 matching cluster

只有：

```text
21/21
```

才：

```text
DEVELOPMENT_DUAL_SOLVER_AGREEMENT_CONFIRMED=true
```

这里报告必须写：

```text
dual-solver subthreshold witness-cluster agreement
```

不得写：

```text
global-minimum set agreement
```

或：

```text
basin uniqueness established
```

---

# 14. Alias audit 仍保留原严格定义

Alias 不使用普通 witness-cluster threshold。

只有 state-separated candidates 独立重建满足：

```text
bearing RMS < 1e-8 deg
acoustic RMS < 1e-6 dB
```

才能进入：

```text
A1_CONTINUOUS_ALIAS_FINDING_STOP
```

有限 catalog 中没有发现，只能写：

```text
NO_EXACT_ALIAS_FINDING_IN_TESTED_CATALOGS
```

不能写 global uniqueness。

---

# 15. Development driver 本轮必须不可执行真实数据

新增例如：

```text
r4_a1_search_tractability_development.py
```

但本检查点中必须设：

```text
DEVELOPMENT_RELEASED = False
```

或采用等效硬 stop。

真实入口若检测未 release：

```text
raise SystemExit(
  "DEVELOPMENT_NOT_RELEASED: research-lead audit required"
)
```

测试必须验证真实 21-case driver 在当前 checkpoint 下无法启动。

---

# 16. 需要新增的单元测试

至少覆盖：

1. exactly 21-case manifest；
2. manifest 不含 truth/error/recovered/J-old/state-old；
3. observation case-id 显式映射；
4. 21 depth branches；
5. solver/budget/branch fresh cache；
6. T1/T2/T3 无状态继承；
7. SHGO/DIRECT 无交叉状态；
8. hard-cap exhaustion 不等于 recovery failure；
9. invalid enforcement 导致 execution invalid；
10. raw top witness selection；
11. state circular-angle agreement；
12. depth-label exact agreement；
13. witness dedup；
14. T3→T2 cluster containment；
15. bidirectional SHGO↔DIRECT cluster matching；
16. development fail 后 fresh generation function 不可调用；
17. runtime truth sentinels；
18. unreleased scientific driver 硬停止。

不得运行真实 noisy case。

---

# 17. 新的 RELEASE FREEZE

完成 harness 后重新运行：

```bash
python -m pytest tests/test_r4_a1_search_tractability.py -q
python r4_a1_search_tractability_audit.py audit
```

更新：

```text
METHOD_FREEZE.json
DEVELOPMENT_EXECUTION_FREEZE.json
LOCAL_VALIDATION.md
A1_SEARCH_TRACTABILITY_REPORT.md
GPT_SYNC.md
```

`DEVELOPMENT_EXECUTION_FREEZE.json` 必须记录：

```text
parent_pre_run_commit =
b4839b295777127ec0c8ade56b76c776db99148b

budget_sha256 =
40cd9c0b2cbc7a6fced83faf97009c25d652223bf7cbb77e7ebc6a981089672d

noisy_development_runs = 0

development_released = false
```

并 hash freeze：

- development driver；
- Gate logic；
- case manifest；
- witness semantics；
- tests；
- frozen budget；
- core objective code。

---

# 18. Git checkpoint

提交：

```text
R4 A1 tractability: freeze development execution harness
```

push `main`。

核对：

```bash
git rev-parse HEAD
git ls-remote origin refs/heads/main
```

必须一致。

---

# 19. 本轮停止位置

提交后立即停止。

**仍然不得运行 21 个 noisy development cases。**

向研究负责人只返回：

```text
parent PRE-RUN SHA
new execution-harness freeze SHA
remote/main SHA

budget SHA
core objective SHA
development driver SHA
Gate-logic SHA

test count
integrity check count

DEVELOPMENT_CASE_MANIFEST:
21/21
truth columns absent = true
observation explicit mapping = true

candidate semantics:
EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM

noisy development cases:
0

development_released:
false

R4:
0%

A2/depth/SSP/P5:
UNOPENED
```

不要继续运行 scientific development。

这个 checkpoint 通过后，我再放行真正的 **21-case SHGO → DIRECT development**。目前正式状态仍然是：

\[
\boxed{\text{R4}=0\%,\quad \text{scientific tractability NOT EVALUATED}}
\]

但 PRE-RUN 的底层预算、objective 与环境审计已经可以冻结，不需要返工。