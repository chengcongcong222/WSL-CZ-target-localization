# 任务包移交与最终给ChatGPT的最少交互流程

## 第一次：立即交给本地高级Codex

把压缩包原样交给Codex，粘贴 `15_CODEX_START_复制.txt` 的完整正文。要求Codex读完整研究合同后自主执行到E5，并维护跨会话恢复状态。不要拆成几十个短任务或重做旧历史T2科研实验。

## 中途：最多按必要决策转给GPT

正常`CP1`和`CP2`进度提交不需要停下来向GPT逐项请示。只有触发需要改变首CZ物理场景准入、数据源、核心方法、硬成本或新增计算预算时，把`17_REVIEW_CHECKPOINT_TEMPLATE.md`填写后的内容连同Git branch SHA发给GPT。没有实质决策的checkpoint让Codex自行继续。

## 最后：一次性整体审核

Codex推送FINAL至 `research/env-cap-20261016` 后，转给GPT这些信息：
- Git科研分支和完整最终SHA（必须不是一个从未验证的假链接）；
- `REVIEW_ENTRYPOINT.md`和`GPT_SYNC_ENV_CAP_FINAL.md`的内容或可获取的Git引用；
- `ENV_CAP_FINAL_SCIENTIFIC_REPORT.md`；
- `ENV_CAP_FINAL_DECISION.json`；
- `INDEPENDENT_AUDIT_SUMMARY.json`；
- 如果大文件只在D盘，则提供`raw_local_archive_index.csv`与不可访问项清单。

GPT最后一次审核：方法身份/物理CZ/环境来源、敏感性反推的证据层级、数值可靠性、项目申请能写/不能写、下一方向是否值得继续。审阅通过之后由用户决定是否合并科研分支至main并开始申请书写作，**自动合并不在本阶段授权范围**。

## 失败仍然能有正确的交付

一旦没有合格CZ或模态数值仍无法认证，Codex不需要为了“完成项目”伪造精度曲线；可以交付严谨的物理/数值准入报告、必要的孔径/源频/SSP资源条件、主要不足和项目选型建议。最终报告要说明这类结果为理论或资料边界，不是端到端定位Demo。
