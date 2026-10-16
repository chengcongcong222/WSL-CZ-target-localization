# ENV-CAP 全阶段研究任务包 v1.0

## 一句话目标

**在“第一会聚区必须真实存在”的不可变背景下，构建可追溯环境族；辨别Yang 2015 / Liang 2018方法的信息上限、可提取性及跨环境性能边界；将精度需求反推为可检验的观测、环境与硬件资源条件；最后形成面向716合作项目申请的科研Demo与负结果报告。**

本包是交给具有长任务执行能力的 Codex 的**整体授权与科研执行合同**，不是一轮一轮等待GPT确认的短任务。用户许可Codex在预算内自主完成E0—E5，完成独立审计、记录负结果、按条件分支继续。**最多在必须改变硬约束、追加资源、无任何可用第一会聚区环境等实质性决策时请求研究负责人/ChatGPT审阅。** 中途允许向远端科研分支推送可审计检查点，但不自动合并main。

## 起动顺序

1. 打开 `01_CODEX_MASTER_PROMPT.md`，严格按照其中的执行规则办。
2. 阅读 `02_SCIENTIFIC_CONTRACT.md`、`03_PHASES_AND_GATES.md`、`04_ENVIRONMENT_PROTOCOL.md`、`05_METHOD_TRANSFER_PROTOCOL.md`、`06_PERFORMANCE_SENSITIVITY.md`、`07_RESOURCE_INVERSE_DESIGN.md`。
3. 读取 `08_GIT_CHECKPOINT_AUTONOMY.md`、`09_REPORT_ACCEPTANCE.md`、`10_REFERENCE_SOURCE_INDEX.csv`、`11_PUBLIC_DATA_SOURCES.md`、`12_STAGE_SPEC.json`、模板、脚本。
4. 先运行 `python scripts/verify_manifest.py --root .` 检查工作包自身；在项目仓库根运行 `python <包目录>/scripts/repo_preflight.py --repo . --expected-parent a9a80191ee4d0f0cecd9db296198606b63a68080`。
5. 在 `research/env-cap-20261016` 科研分支冻结研究设计，按门槛滚动执行，不碰历史结果；报告完成后推送分支，等待GPT最终整体审阅。

## 固定背景与最新有效证据

- **硬边界**：第一会聚区必须存在；是物理场景准入门槛，不是算法读取到的已知目标距离。
- **非硬边界**：5 km海深、50–60 km位置、200 m源/收深度仅为专家经验代表性场景；可合法变化，但必须标记环境与资源代价。
- 单平台单拖曳水平阵是核心基础资源；辅助节点、非共线横向孔径、垂向通道、环境资料都是可评估的**新增条件**，不得写成原硬件自有能力。
- 目标状态 `x=[r,theta,z,v,psi]`；不要求同一种观测独立确定五维。
- Yang论文浅海条件下的机理复现已有；E-STD中理想深度签名存在，但2.4km普通Fourier模态身份恢复受限；Liang/MMAR论文含噪和完整深度估计链未闭合。
- T2现有数值修补链**已STOP**（FUNC2解析准入失败，不是所有高精度方法均不可行）；T3-W1只有有符号方位前提下的条件水平Demo，B3偏置触界依赖明显；R4总进度0%，不改历史Gate。

## 项目科研索引与代码复用说明

`10_REFERENCE_SOURCE_INDEX.csv`列出的GitHub文件只是在当前冻结parent上**预期存在**的历史证据；Codex仍需按Git对象与SHA验明身份。旧环境接手文件《验证环境鲁棒性.txt》属于较早的阶段快照，可以辅助阅读，但其“当前进度”不得覆盖GitHub中较新的FUNC2和T3-W1正式判定。

## 最终输出形态

一份可交研究负责人/ChatGPT整体审核的 `GPT_SYNC_ENV_CAP_FINAL.md`、`ENV_CAP_FINAL_SCIENTIFIC_REPORT.md`、环境和方法证据库、性能与资源界图表、完整代码与结果哈希、独立审计、负结果、研究Go/No-Go、项目申报可用的图文片段。**没有取得的误差不许填0，没证明的资源下界不许称“最低需求”。**

## 本包完整性

`FILE_MANIFEST_SHA256.json`为**唯一**包内文件SHA登记；启动提示与README不复制其子文件哈希。验证脚本独立哈希所有文件，包含零失败才准入。压缩包SHA仅在GPT交付回复中提供，不能用其它已失效的外部记录替代包内校验。
