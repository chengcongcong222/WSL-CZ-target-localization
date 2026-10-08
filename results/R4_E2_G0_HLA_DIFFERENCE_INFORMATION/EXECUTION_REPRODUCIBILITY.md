科学执行入口为r4_e2_g0_launch.py；原始冻结科学文件r4_e2_g0.py未改变。EXECUTION_STARTED.json禁止重复科学执行。

审计入口：在D盘项目目录应用runtime-environment.ps1后执行 `python -X utf8 r4_e2_g0_audit.py --checks-only`。仅重建已有模态/记录并更新当前阶段审计文件，不调用传播求解器、不重跑场景、不改阈值，也不重写报告/追加主账本。脚本没有调用原始科学函数。研究负责人可另外核对冻结文件与Git提交。
