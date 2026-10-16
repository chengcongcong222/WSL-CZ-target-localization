# 运行环境与工程化滚动执行说明（Windows / WSL）

## Windows主运行路径

历史工作盘通常在 `D:\XiaomiMiMoProjects\WSL系统汇聚区目标定位` 或 `D:\ProjectStorage\WSL-CZ`。本包**不假设能够由GPT访问用户D盘**。Codex在本机首先自动定位正确Git仓库并以远端SHA核对；不要让新任务覆盖旧源码目录。建议科研本地缓存 `D:\ProjectStorage\WSL-CZ\ENV_CAP`，分块保存原始海洋文件和solver输出，Git仅commit索引和轻量图表。

若在WSL中执行，可用对应挂载 `/mnt/d/ProjectStorage/WSL-CZ/ENV_CAP`；不要同时在Windows和WSL用两个工作树编辑同一个受控文件。新阶段自动记录Python解释器、NumPy/SciPy、GSW、Bellhop/其他exe编译指纹、操作系统、BLAS线程和可用CPU/内存/磁盘。

## 启动与预检命令（示意）

```
python scripts/verify_manifest.py --root <解压包目录>
python <解压包目录>/scripts/repo_preflight.py --repo <Git仓库根> --expected-parent a9a80191ee4d0f0cecd9db296198606b63a68080
python <解压包目录>/scripts/stage_state.py init <新Git研究结果目录>
```

如果 `origin/main` 远端SHA在正式E0前已变化，先比较是否可快进并保持原父提交历史身份；禁止自动改成新的父SHA后混用旧任务包。若其他人提交了无关修订，允许在独立worktree以冻结SHA建立分支，不重置他人的main。

## 前向环境求解器注意事项

- 在Windows/WSL加载旧`tools/acoustics_toolbox`先检查可执行架构、license、Windows/WSL指令差别、I/O路径编码、浮点和输出版本。是否可运行不等于计算可信。
- BELLHOP BOX的深度和水平距离单位要根据程序语法审查，保留过去`ZBOX=5001m,RBOX=70km`的正确案例作为单元测试，不把`70km`误写成`70000km`或`5000km`。
- WOA23标准深度与目标海深一致性、压强/温度/绝对盐度单位和经度经向取样必须事前单元验证。
- 旧`m452/m453/m454`代码中绝对WSL路径不可直接执行；先抽取可核验算法而非修改老文件。明确资料原样复用和代码迁移的差异。
- 若没有足够数值设备执行模态相干方法，切换合法报告`L0/L1/L2`与资料需求分支，不能通过绕开数值Gate做L6性能。

## 复现与孤儿任务防止

每次科学求解前写`SCIENCE_STARTED.json`，结束后写 `JOB_EXIT.json`，包含PID、子进程、status、wall/RSS/路径/输入hash。崩溃后只允许只读恢复已完成不可变原始数据，**不重新生成同名科学生成**；未完成的独立试验单元标`FAILED`并继续预注册可独立阶段。

长输出日志可以滚动写盘，但摘要必须可阅读。Codex上下文不足时写 `SESSION_HANDOFF.md`，刷新任务契约和哈希后接续已允许任务，不因为换会话重新抽种子。
