# ENV-CAP-CZROB-02：最终同步

**结论：CZROB02_NUMERICAL_ADMISSION_NOT_CLOSED_AND_ENERGY_ATTRIBUTION_UNRESOLVED。** 本轮负结果分支交付完成并STOP，不建议现在启动相干提供器开发。

- 父SHA：`4e53bf035625aeddf4d8ab6ee009242c06726a9c`。
- 设计SHA：`3d2a56e1821c3c94949346bb5d9e80170dfc221b`。
- 最终SHA：本报告所在Commit B；完整身份见Git科研分支HEAD。
- 分支：`research/env-cap-czrob02`；不合并main。

## 核心科研结果

同一WOA派生节点上的C-linear和C1 PCHIP，在9个单元中都出现首深回返候选；保存的首候选移动3.94–32.21 m，转向深度移动1.81–24.53 m。18个首候选点的独立参考比较通过。局部极小值计数由C-linear的4–8个变为PCHIP的1–3个，显示细节表示依赖。这些是**局部诊断，非稳健CZ通过**。

18个全角域积分阶次Gate全部失败，最大差6.53–48.12 m；因此0/9取得本轮数值准入，代表射线ODE按条件门取消。**0/9不是CZ消失或物理不可行。** 冷审计39,601项中6项失败，集中在三个近节点转向角，独立差5.29–19.90 m。冻结实现对跨节点近转向的积分精度不足，没有修补、重试或放宽标准。

18组存档射线前缀分类统计经独立解析核对一致：原焦区同时有零边界与边界接触回返。**RAY计数不能换算能量占比。** 全部9单元能量归因仍EXCLUSIVE_ENERGY_ATTRIBUTION_UNRESOLVED；PCHIP总TL与各支能量均NOT_EVALUATED。

## 研究去留

现有局部诊断没有显示C1对照令首候选整体消失，但全域精度和主要路径能量归因均未闭合。暂缓高精度相干提供器；不启动Yang/MMAR性能，不重开T2/FUNC2。原REC01已接受的条件证据保持，三套气候组合仍只算一个来源族。新研究没有证明真实海域CZ、可检测性或测深精度。

## 执行核对

922项历史输入身份检查通过，冻结后931项执行身份检查；25项合成控制通过。主分析一次，冷审计一次，科学壁钟106.64 s，采样RSS峰值486.0 MiB。新增外部声学求解0、MC0、新接收/SHD0，所有新增文件在D盘。main、历史STOP及R4=0%不变。

详细结论见[CZROB02_RESEARCH_DECISION.md](CZROB02_RESEARCH_DECISION.md)，所有失败见[CZROB02_INDEPENDENT_COLD_AUDIT.json](CZROB02_INDEPENDENT_COLD_AUDIT.json)，两种表示全曲线、原节点/系数、边界路径穿越和原始哈希索引均一并交付。交付完整性检查与科学Gate必须分别阅读。
