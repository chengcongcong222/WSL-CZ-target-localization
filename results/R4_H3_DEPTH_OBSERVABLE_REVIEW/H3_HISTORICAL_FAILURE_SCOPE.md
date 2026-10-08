# 历史失败作用域

所有旧报告、代码、Gate及失败样本保持不可变。以下是研究负责人已接受结论的继承，不重跑、不计新成绩。

| 历史证据 | 保留事实 | H3必须避免的重复 |
|---|---|---|
| [R3-C2.3A](../R3_C2_Yang_SA_depth/R3_C2_3A/R3_C2_3A_REPORT.md) | 0.6/1.2/2.4 km普通Fourier合成孔径在E-STD密集模态下有身份瓶颈；2.4 km Rayleigh独立模式为0，理想系数20/20有深度签名 | 不把FIELD输出的k、φ、模式编号或幅度直接交给提取器。理想签名不证明实际恢复 |
| [MMAR幅度桥](../R3_C2_AR_MMAR/R3_C2_4B_3A_HANKEL_AMPLITUDE/R3_C2_4B_3A_REPORT.md)及[D(z)合同](../R3_C2_AR_MMAR/R3_C2_4B_3A_HANKEL_AMPLITUDE/DZ_INPUT_CONTRACT.md) | 部分浅海机制和AR→Hankel桥已验证，完整深度函数和E-STD未知运动链未闭合 | AR用于位置，不用峰高作能量；不继承δ=0或真实mode-ID作为实际算法条件 |
| [B1](../R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY/B1_REPORT.md) | 理想已知水平状态下有深度签名；不是带水平误差的测深 | oracle只作上限，不能作主性能 |
| [B1A](../R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY/B1A_REPORT.md) | 注册0.25 km/0.05 m/s范围13/150实际通过；最小0.125 km/0.025 m/s也仅7/54；全部16矩形失败，369个非单调见证 | 不再缩小步长找通过范围，不把局部水平容差当可实现上游精度 |
| [B1B收口](../R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT/R4_B1_CLOSEOUT.md) | 保存格点MIN含中心6/6成功依赖真值支撑；去中心、Q25/MEDIAN不支持原联合工程路线 | 不把固定有限点网格MIN变成完整水平域证书；H3不得恢复同一exact-TL点值处理器 |
| [E2修复停止](../R4_E2_NINE_PAIR_FORMULA_REPAIR/GPT_SYNC.md) | E2_FORMULA_REPAIR_FAILED限定为FROZEN_PRE_DEVELOPMENT_PRESSURE_IDENTITY_GUARD_FAILED；P1解析秩48/48为4，零方向1.43e-17，完整信息回归0次 | 不将坐标一致性守卫失败解释成解析径向导数失败或E2物理No-Go |

E2拆分：九频对数值可恢复性已有支持；解析水平链及结构秩已有修复证据；完整公式Gate与未知校准/深度下科学增量未建立。当前PAUSED_WITH_FORMULA_REPAIR_FAILURE是资源分配状态，保留将来经重新预注册统一Cartesian实现复用的可能。本轮不修E2。

H3改变的是观测量/接收维度及深度输出合同，不能只给旧处理器换名字。深度区间或垂向响应即使成功，也不自动解决1200 s纯方位速度瓶颈。原E2-G0=FAIL_UNCHANGED，原pilot=IMPLEMENTATION_INVALID，R4=0%。
