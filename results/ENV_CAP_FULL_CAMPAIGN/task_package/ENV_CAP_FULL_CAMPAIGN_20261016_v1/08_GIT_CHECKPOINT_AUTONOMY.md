# Codex滚动执行、Git检查点及异常升级协议

## 版本与隔离

父main固定SHA `a9a80191ee4d0f0cecd9db296198606b63a68080`；从该版本创建新分支 `research/env-cap-20261016`，持续把科学代码、配置、报告、图、哈希索引提交到**新结果目录** `results/ENV_CAP_FULL_CAMPAIGN/`。历史结果目录只读，`results/R4_MASTER/R4_EVIDENCE_LEDGER.csv`也只读；如需ENV-CAP汇总，创建自己的`ENV_CAP_EVIDENCE_LEDGER.csv`。不能改写别的任务的STOP/Gate，不能force-push或重置main。

D盘中间文件推荐`D:/ProjectStorage/WSL-CZ/ENV_CAP/`；原始大型场、声音、模型、可复现包索引保留于D盘；Git只保留小体积科学源、CSV汇总、图、必要中间数据索引。所有跨环境数据以路径+SHA256+字节数+来源ID登记。如果Git输出>120MiB，应转大文件到D并保持哈希。

## 开始

```
git fetch origin main
git rev-parse origin/main  # 应为父SHA；如不符，暂停仅做身份核查
git switch -c research/env-cap-20261016 <frozen-parent-sha>
```

如果本地`main`包含未推送内容、历史结果未提交或remote已经变化，不能把当前工作强行覆盖、也不能借他人新HEAD修改原合同；记录冲突、停下索要合并选择。

## 先冻后算

建立设计提交`ENV-CAP: freeze full-stage contracts and source registry`，包括本包快照及派生科学合同、数据索引、代码单测、阶段预算、预注册环境/轨迹/种子和Gate。源数据若尚未下载，提前登记获取规则、URL和身份审核方式；下载后freeze SHA作为单独准入步骤，不能事后因性能改来源。

审计每个科学任务：`SCIENCE_STARTED.json`包含阶段/单元ID/代码SHA/参数SHA/输入SHA/seed/purpose/预算；失败后该单元不能改阈值重试。若预先登记多个独立物理环境，可按预登记配置逐个执行，环境失败不作为改良其它环境的借口。

## 推荐三个主动远端检查点（其余本地提交）

**CP1：环境准入（E1）**：`ENV-CAP/CP1_CZ_REVIEW.md`一页说明源/环境族、CZ2数量、ray-TL数值审计、无CZ负例、是否能够继续E2；完成commit并`git push -u origin research/env-cap-20261016`。如至少一合法CZ2，直接按冻结合同继续；如果零CZ2且要重选未冻结海域/下载超许可/追加环境上限，发人工审阅请求。

**CP2：Yang/MMAR机制与性能（E3）**：`ENV-CAP/CP2_METHOD_REVIEW.md`列每环境L0–L6准入与否、负结果、提取链/条件性能、尚未打开的指标，commit/push。若方法Gate失败且已预设E4理论资源边界可继续，继续而不请求重新许可；只有需要变更核心模型和研究预算才人工审核。

**FINAL：整个阶段（E5）**：保证独立冷审计`PASS`或如实`FAIL`，提交全部结果、`GPT_SYNC_ENV_CAP_FINAL.md`和`REVIEW_ENTRYPOINT.md`，push分支后STOP；用户转发审阅入口给GPT，不自动合并main。

如果只触发早期科学NO-GO，允许CP1后直接做负结果界/资料需求与E5；不要求凑够三个推送。

## 什么时候必须向用户/GPT提问？

只要满足任一项，冻结当前工作状态、生成`CHECKPOINT_REQUEST_TO_GPT.md`并报告有确证的当前SHA和可选方案（不得伪造自动通过）：
- 第一CZ物理准入0个，且希望改变事前冻结的候选范围、源深或场景类型。
- 需要超出72h、8GiB、30GiB、240求解调用或打开未注册算法/硬件。
- 关键资料的原始来源、许可、原理矛盾阻碍有效研究，无法仅以`NOT_EVALUATED`继续下游必要条件分析。
- 用户要求修改固定硬背景，或者候选需不可追溯“环境调好看”的虚假条件。
- 科研伦理/证据冲突或可能会把oracle/前向错误导入申请书。 

**不需要**每运行一次仿真就询问；也**不需要**每个负面Gate都询问。负结果是成果，继续已有合同允许的、独立的分析。

## 错误恢复边界

正式SCIENCE_STARTED前：允许修复依赖、路径、Git行尾、图表排版、解析单元控制的实现bug，但每项改动留记录并重新冻结代码SHA。正式SCIENCE_STARTED后：仅允许只读恢复已封存输出的统计/绘图错误，不能修改科学方法或重采原始观测。科学程序失败即该试验单元STOP，不得跳过失败重跑；阶段可执行合同明确批准的其它独立单元。

## 跨会话状态

每次阶段结束/交互中断，都更新 `RUN_STATE.json` 与 `SESSION_HANDOFF.md`：当前分支SHA、完成/未做步骤、下一确定性命令、资源累计、关口状态、已推检查点、历史不可触碰资料和已知异常。**Codex重启时先读这两个文件**，不得重新运行完成的科学单元。

## 给GPT的低频审阅消息模板

`项目：ENV-CAP；分支/commit SHA: ...；本次必须决策：...；冻结条件：...；已完成/未执行：...；核心证据：...；负结果：...；推荐方案A/B和代价：...；请仅审该决策，不要求重复整个项目。`
