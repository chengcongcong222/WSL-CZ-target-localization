# CAP06 Git、状态恢复与低频审阅协议

## 分支和保护

冻结继承父SHA：`dddad1f66e6d61668ecd5af38cad3e3af9fc9245`，研究分支：`research/env-cap-cap06`。

若研究分支已由上游预置设计文件，Codex须**fetch后checkout该分支**，验证REQ05父提交为其祖先，并记录初始HEAD。不得因为分支已存在而创建同名新分支或reset其他工作。若尚不存在才从父提交新建。若本地工作树含其他用户变更，使用独立`git worktree`或全新clone，不清理、reset、覆盖用户文件。

成果只写`results/ENV_CAP_CAP06/`与独立D盘目录（由本地preflight自行发现），旧历史与`main`只读；大二进制留D并存索引/字节数/SHA，Git≤120MiB。仅推科研分支，不force-push，不自动合并main。若远端HEAD有未经本任务创建的并发推进，先审计差异，停止写入并询问负责人。

## 推荐远端节点（最多三次科学检查点）

1. **Commit A**：P0物理合同、测量接口、所有使用的信号模型、O2/O6各自Gate与资源上限、留出抽样规则全部冻结。它可以跟在本任务包预置的设计种子commit之后。Commit A之前可反复运行解析与合成**开发**零控制并记录代码版本；不可窥看正式场景结果后调参。
2. **CP1**：P1/P2两条支线各已按冻结计划执行或合法`NOT_EVALUATED`关闭；将数字原始索引、冷审计、关键失败和是否具备P3资格推送。通常自动继续，不逐轮求GPT审批。
3. **CP2**：如果P3合法执行则推有条件联合/水平4D统计，否则推正式未评价/必要条件。自动继续P4。
4. **FINAL**：P4最终科研报告、证据矩阵、代码及审计、来源哈希、预算、RUN_STATE、STOP全部推送，检查本地HEAD与远端科研分支SHA一致。

注：上面是**Commit A + CP1 + CP2 + FINAL**；CP1/CP2/FINAL为主动审阅检查点最多三次，不代表只有三次Git commit。Codex可将一支线失败封存并继续独立支线，避免反复让用户参与。

## 需要人工审批的仅四类情况

- 变更硬背景（不再第一CZ、绕开单平台单阵基础资源）或核心研究对象；
- 需要**未登记**传感器/算法谱系、真实用户/第三方对外访问，或声学求解需要借用旧T2修补；
- 科学累计预算超限、或必须改变已经正式冻结的数值/统计Gate；
- 出现会导致历史记录被篡改、oracle混入估计器或数据来源造假的实质矛盾。

普通数值失败、缺实测数据、某支线未准入、仪器未知不能使整个任务反复问人；执行负结果允许路径和P4资源边界即可。

## 跨会话恢复

在每次阶段转移后更新`RUN_STATE.json`与`SESSION_HANDOFF.md`：`baseSHA`, `currentSHA`, `phase_status`, `SCIENCE_STARTED units`, `frozen_artifact_SHA`, `total_external_solver_calls`, `total_MC`, `wall_s`, `D_bytes`, `open_issues`, `next_exact_command`。恢复时先读取，再执行下一未完成状态；禁止再次跑已经SCIENCE_STARTED的实验。

本包`scripts/stage_state.py`可作轻量日志生成器，但不能用它的状态值替代科研证据。用户最终只需转发新阶段的Git最终SHA、`GPT_SYNC_CAP06_FINAL.md`及关键Gate即可进行整体审核。
