# 已关闭路线的适用范围复核

本轮不修改任何旧Gate、结果或失败样本；仅区分实现闭合与方法族的科学边界。下述证据为保存报告和代码核查，本轮未重跑历史实验。

| 路线 | 已冻结判断与证据 | 本轮保留的研究空间 |
|---|---|---|
| 同条件纯方位估计器 | [速度信息审计](../R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY/SPEED_INFO_REPORT.md)：已测L2/soft-L1及观测/真值初始化差异很小 | 关闭继续换损失/初始化/滤波器；新增频谱、传播、空间观测仍可研究 |
| exact-TL A1/FIX1/FIX2及SHGO/DIRECT | tractability未建立；通用observation-only预算内未稳定捕获正确窄区；旧冻结保持 | 不做T4/FIX3/第三优化器；不证明正确解不存在或所有算法不可能 |
| 常数β迁移 | [R3-A1](../R3_A1_fidelity/R3_A1_REPORT.md)，解析二模态机制及E-STD常数β迁移不足 | 未完整复现Cockrell原论文条件；未测试所有局部/分模态条纹。禁止把拟合β称物理不变量 |
| 路径身份/MMAC | [R3-B1](../R3_B1_MMAC_observability/R3_B1_FINAL_REPORT.md)，身份与观测维数未闭合；M=1无独立相对时延 | 聚合相关/无身份相干结构仍需独立原文与接收链；六对时延不等于六个独立观测 |
| CZ保守全域Q2/N0搜索 | [正式closeout](../R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY/B1_CZ_ROUTE_CLOSEOUT.md)：45–60km冷启动、whole-cell外包、注册预算未闭合 | 关闭该认证搜索结构；不否定所有传播表示。不能把局部FIM或采样优化伪称全域证书 |
| Yang普通Fourier mode-ID | [C2.3A](../R3_C2_Yang_SA_depth/R3_C2_3A/R3_C2_3A_REPORT.md)：0.6/1.2/2.4km孔径、密集模态不能唯一mode-ID | 已知理想系数存在深度机制；不证明组能量或所有高分辨方法失败；不得再用δ=0掩盖offset |
| MMAR AR/Hankel桥 | [Hankel桥](../R3_C2_AR_MMAR/R3_C2_4B_3A_HANKEL_AMPLITUDE/R3_C2_4B_3A_REPORT.md)：无噪浅海机制；非Fig.5/6严格复现；D(z)未执行 | 原始幅度/相位链与未知运动仍未闭合；不能把AR峰值当幅度 |
| 条件化深度B1/B1A/B1B | [B1 closeout](../R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT/R4_B1_CLOSEOUT.md)：精确水平oracle有深度签名，但plug-in与保存格点联合剖面未建立工程能力 | 不继续缩小(r,v)步长找通过范围；新增垂向/组统计机制可提出设计，不自动开放depth执行 |

原联合指标A1/R4=0%，速度子指标未通过。以上不是物理不可辨识证明。
新阶段完成的是文献覆盖、未知量分析和两项设计；不计实验通过或硬件验收，不虚增预研百分比。
