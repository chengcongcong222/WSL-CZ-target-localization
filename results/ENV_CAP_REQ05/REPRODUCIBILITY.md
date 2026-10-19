# 复核入口

首先阅读GPT_SYNC_REQ05.md和RESOURCE_AND_RESEARCH_DECISION.md。历史证据索引为EVIDENCE_INDEX.csv（检出SHA256及父提交blob SHA256分列，CRLF差异单列），方法原文使用封存PDF；不重新下载数据。

Commit A冻结代码与18个无量纲矩阵。analyze.py使用NumPy正交投影；cold_audit.py不导入主程序，使用Fraction精确消元/逆矩阵及独立标量差分。ALGEBRA_RESULTS.json和EXACT_NULLSPACES.json提供全部矩阵和零空间，INDEPENDENT_COLD_AUDIT.json逐条记录124项。

run_once.py使用D盘排他启动标记；正式主分析和冷审计已各执行一次，不得在当前阶段删除标记或重试。开发脚本仅用于追溯冻结前构建过程，不是继续科研入口。build_delivery.py只格式化已保存结果/编写报告，不调用科学分析。后续外部审核如需独立重算，应使用另行批准的审计副本及独立输出位置，禁止覆盖本轮正式结果。

运行环境：Windows，Python 3.10，NumPy 1.26.4，psutil 5.9.8。PYTHONUTF8=1、PYTHONDONTWRITEBYTECODE=1；环境与临时目录来自D:/ProjectStorage/WSL-CZ/runtime-environment.ps1；依赖只读路径D:/ProjectStorage/WSL-CZ/T2_PREC3_FUNC2/runtime_dependencies。使用该依赖目录不运行FUNC2。全部计算为小型确定性代数控制，无传播求解、MC或实际观测生成。

冻结前两次文档/身份工具错误已记录：CRLF/LF检出差异与Windows GBK默认解码。均在正式Commit A前处理；未修改任何历史输入和科学阈值。预检与正式研究分别记账。

结果文件不自引用其所在提交SHA；最终执行SHA以Git中包含报告的Commit B及远端核对为准。DOC_DELIVERY_VALIDATION是交付文件/论断一致性检查，不是第二次科学冷审计。
