# E1 基线比较与贡献边界

[模型原式](E1_ACOUSTIC_MEASUREMENT_MODEL_LOCK.md)与[覆盖表](E1_PRIMARY_SOURCE_COVERAGE.csv)是本表的证据。以下为未来实验合同，本轮无比较数值。

| 标签 | 继承方法/输入 | 对照用途 | 原始条件和当前边界 |
|---|---|---|---|
| B0 | 当前双节点 bearing-only，既有四状态动态模型 | 同1200s信息基线 | 保留旧噪声/几何/运动/全部镜像，不能换panel |
| LIT-F | Sun2024 β+多线接收频率、未知 F_l 状态、随机游走、MFB-AUKF | 最近的真实水声 bearing-frequency baseline | Eq.(4)已锁；完整filter印刷冲突未闭合，严格原始LIT-F禁止执行 |
| LIT-RV | Xu2018 原始复声压→波束域相关→FFT→速度幅值与条件符号→β/q EKF | 层2重要的 extracted-observable baseline | 浅海 v_p/波导符号/相位条件；第一CZ不能直接认定可迁移 |
| LIT-RV-PF | Li2019 同 β/q observable、SIR/VSIR | SAME_OBSERVABLE_ALTERNATIVE_ESTIMATOR_BASELINE | 估计器替换不能宣称新信息；有限初始支持域须登记 |
| UUV-CONTROL | Kita motor载频Doppler、侧带RPM、known vehicle映射、bearing/差分bearing | 展示车辆与频率先验条件收益 | 当前 NOT_TRANSFERABLE_CONTROL，不能使用真位置选择栅瓣 |
| NEW | 双频率未知源频率/漂移及节点频率参考 nuisance 联合剖面 | 候选条件扩展 | 无已证明性能；不得用发射频率/径向速度/谱线身份truth输入 |
| CW-CONTROL | Yang2015 + 2018 erratum已知CW/PLL/Doppler物理控制 | 符号、单位和已知源的机制上界 | 非唯一LIT基线；不是未知源方法 |
| RF-REFERENCE | Gencel2018 FFT/CFO/drift机制 | 频率提取与参考模型参考 | RF，不能用于声学TMA性能比较 |

## 同等资源比较

主比较先做有效信息层：所有分支使用同一目标状态、1200s、节点、bearing与频率采样、同一观测误差定义；将原文方法搬到同条件时标 PAPER_ADAPTED_CONDITION，不能冒称复现其2000s/海试RMSE。单节点 Sun原式作为 PAPER_ORIGINAL_NODE_CONTROL 另列；同样双节点数据但原式与源状态融合的版本标 LIT-F-DUAL-ADAPTED，不能以“多一个节点”当算法突破。两节点共享 F_l 与文献随机游走过程条件须显式列出，与 NEW 无校准/漂移扩展分别计算。

提取层只比较从相同原始输入独立获得的 observable 与其Σ。Sun 的未知 F_l 是既有机制；其随机游走是已有频率变化处理。NEW 是否放宽为显式共享/线特有漂移、独立参考漂移，须做一对一条件报告。为适应误差而降低频率过程噪声、冻结一个clock、使用已知源频率，均算新增先验，不能藏进 baseline。

Xu 的相关链提供实际 q，不必经已知 emitted f0 换算，但要支付相干录音、足够延迟/带宽、v_p 与符号条件。失去符号时保留±分支/NOT_APPLICABLE。Kita航速来自车型映射；第一CZ中不存在已验证的对应PWM强线/转速映射和100m信噪条件，其论文MAE不可直接当本项目指标。

## 冻结研发选路规则

1. 对每个已登记场景，V_B0与V_NEW为剖面掉源/参考/环境 nuisance后的可估 speed 局部方差。输出不可估方向为 UNBOUNDED；不以伪逆0当成功。
2. 仅有限可比方差定义 Δ=1−V_NEW/V_B0。24场景中位Δ≥0.20，且最坏场景有效方差不增大，即 max_s V_NEW,s≤max_s V_B0,s（仅数值等式容差≤10⁻⁹相对）；逐场景恶化仍全量报告才过主选路门槛；任何非可估/非可比场景保留并阻止整体PASS，另报原因。
3. 同时逐场景报告相对最近可迁移真实水声 baseline（优先LIT-F同资源适配）的 V 与 nuisance前后rank、弱方向、协方差条件、源/参考先验、节点/时长成本。LIT-RV不可迁移时必须标条件缺口，不能以 ORACLE_vr 代替实际文献基线。
4. 若文献已经提供类似/更强增益，NEW必须在登记的 source/reference漂移或减少车型/校准先验条件下保留≥20%的B0收益且最坏不恶化，并报告相对LIT-F的表现；仅未知f0不构成新条件。原文 Sun 已具未知f0及随机游走，必须承认继承关系。
5. 达到规则但只复现已知机制：LITERATURE_MECHANISM_REIMPLEMENTED_NO_NEW_EXTENSION。新漂移/参考条件仍保留有效信息：UNKNOWN_SOURCE_FREQUENCY_ROBUSTNESS_EXTENSION，仅预研条件扩展标签；不等于最终创新确认。
6. 仅校准prior使增益成立→CALIBRATION_CONDITIONAL；理想层过而提取未过→EXTRACTION_NOT_ESTABLISHED；零控制不对→IMPLEMENTATION_INVALID。不得回头改panel、阈值、漂移范围或增加预算以获得PASS。

20%为研发路线筛选门槛，不是原速度10% Gate，也非实测P95。当前所有增益、创新标签均 NOT_EVALUATED；严格 LIT-F未冻结使执行仍 NOT_READY。
