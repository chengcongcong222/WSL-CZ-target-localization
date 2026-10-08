# 水平支持来源审计

## 判定与定义

FULL_RC2_HORIZONTAL_SUPPORT_NOT_ESTABLISHED。对H3三组主阵几何，没有已有输入能同时提供观测生成、无真值泄露的完整多分支水平支持、门槛的覆盖解释和全域候选保持证明。

| 支持层 | 来源 | 当前可用作用域 |
|---|---|---|
| O | G0原始真值水平中心 | oracle局部上限，不可注入搜索或当已知先验 |
| F | 在真值径向和速度方向加25个误差节点 | 读表诊断；theta/psi没有扫描，不是观测支持 |
| U | 全部可用水平观测定义的多分支集合 | 完整支持未建立 |
| U0 | 原单HLA方位观测定义的多分支集合 | 同样未建立；不能把双节点点估计替代 |

## 已有证据逐项核查

1. P2_RC2_kinematic_boundary/P2_RC2_GPT_SYNC.md：冻结共线基准中r/v对方位无影响；多个条件接受距离铺满45–60 km。它证明几何边界和困难候选来源，不证明H01/H06/H12含噪情形的连续支持覆盖。

2. R4_A1_OFFGRID_BEARING_BOUNDARY：粗网格45–60 km/1 km、theta ±5°/0.5°、v1–3 m/s/0.2 m/s、psi ±15°/1°；门槛为min SSE+13.3 sigma_rad²，无噪声加1e-12。9个off-grid真值在0.02°和0.05°各0/270邻域单元保留，无噪声0/9。SURVIVOR_SET_QUALITY记录多连通分支，例如P01首个nominal实现有4个分量，最多不能只保留最大分支。这是离散真值邻域保持诊断，不是连续空间真实覆盖率。

3. r4_a1_fix_continuous.py:125–151及FIX_METHOD文档：观测导出的Cartesian线性解与32个Sobol起点求 bearing min；含噪情况下4096个(r,v,psi)样本经theta剖面，只保留门槛内点，按固定分辨率去重。沿用13.3 sigma_rad²、无噪声1e-20。解决网格表示不等于穷尽连续似然域，最小值本身也没有全局证书。源码明确“finite…not certified exhaustive confidence region”。不能把有限云导出当U0闭合。

4. R4_A1_FIX2_ACOUSTIC_COVERAGE及SEARCH_TRACTABILITY：前者声学累计18/21，后者通用搜索0/21，都是搜索证据。声学覆盖失败不能反推RC2真实支持不存在，声学找不到也不能用来删掉水平分支。原门槛和原停止结论均冻结，不开放FIX3。

5. R4_AUX_GATE2B_JOINT_ERROR_BUDGET：静态11个离散距离、符号角点和假定误差合同，只给静态range误差要求，不生成含运动多分支水平置信集合；时间和关联未数值验证。

6. R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC：H01/H06/H12来自其冻结12-case off-grid panel。双节点ray初始化、IRLS、soft_l1局部优化输出点估计；局部协方差在r4_a1_new_estimator.py明确NAIVE_LOCAL_DIAGNOSTIC，未纳入系统/导航误差。Anchor A最坏range P95约1.3006%，速度28.9386%；它不是同时四状态的支持集、不是H3完整支持证书，也不能只截取较好的range指标。H3-G0镜像辅助节点没有作为第二复声场孔径使用。

7. R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY：CURRENT_1200S_GEOMETRY_INFORMATION_LIMIT_CONFIRMED_BY_CRLB，所有12个case局部Gaussian速度等效量高于10%。只确认既定模型局部速度瓶颈，不是普适不可行，也不提供新的联合支持。

## 下一次科学Gate的必要内容

U/U0必须由观测定义，登记联合likelihood/系统偏差/导航/DOA模型、物理搜索盒、接受门槛及其统计解释。13.3只是历史继承系数；本轮不重新包装成所有误差模型通用的99%置信阈值。需要独立全域多分支发现/候选覆盖证据与新观测条件下真值保持审计，全部失败和不可界方向保留，不靠真实中心兜底。

深度集合应对全部可接受水平与环境/校准分支联合剖面或求并集。有限网格最小值与求解器成功只是候选上界，不是连续全域下界。无法建立支持覆盖时，输出 NOT_ESTABLISHED，而不是“最优分支测深成功”。本轮不构造新RC2观测、不跑优化器、不收紧网格或门槛。
