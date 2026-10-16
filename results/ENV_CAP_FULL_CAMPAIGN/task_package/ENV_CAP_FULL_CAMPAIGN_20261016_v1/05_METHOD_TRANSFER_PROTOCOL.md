# Yang 2015 / Liang 2018跨环境适用性与多层次机制协议

## 一、原方法身份：不同分层不得混称

**Yang 2015**（JASA 138:1678–1686，DOI 10.1121/1.4929748；2018 Erratum DOI 10.1121/1.5081712）：运动声源相对接收形成径向合成孔径，采用已知/可估径向增量的相干CW接收，利用水平方向模态的空间波数与深度结构判断源深。原论文并不是“单个14m物理HLA同时一次解五维”的算法；如果源频/径向位移/起点delta、模态身份已知，那是能力上限。Erratum的offset相位改正必须进入算法核与验证。

**Liang 2018 MMAR**：HLA/相对运动、AR提供模态波数位置，广义Hankel谱提供对应复幅度、模式顺序匹配与深度评分。原AR峰值不能当成模态幅度。此前在浅海无噪条件观测到模式次序优势；20dB和真实峰选/完整`D(z)`未闭合。若原主文的AR峰选阈值、平滑/带宽、模式数无法从原文锁定，应给`PAPER_IMPLEMENTATION_PARTIAL`，不可把自写算法宣称严格原论文复现。

**基线**：普通Fourier/SAB谱、Yang已知模式oracle、MMAR部分论文机制、自有聚合模态/分辨改进（只有合同先锁且确有新信息时才可作为自有方法）。不要在看到结果后添加一套“刚好通过”的新算法。

## 二、物理与接收层严格区分

层A `IDEAL_MODAL_INFORMATION_UPPER_BOUND`：真实模态波数/深度函数/模式身份（truth-oracle）仅用于机制上限，评价深度区分度和多候选；不能报告可提取精度。

层B `FINITE_RADIAL_APERTURE_SPECTRUM`：已知径向坐标和源CW频率的合成孔径模态分离；普通Fourier的`delta_k≈2π/Lr`仅给谱分离分辨尺度，**不是所有高分辨算法不可逾越的CRLB**。必须保存相邻波数间隔分布、有效模态激励权重、谱峰混叠与身份匹配错误。

层C `NON_ORACLE_COHERENT_EXTRACTION`：生成连续接收相位/谱线，估计端看不到真实f0、source drift、delta、真实模态ID或目标轨迹。考虑source/receiver/clock/navigational共同未知量后计算可用的独立信息和提取失败。若频率/距离增量无法辨识，归档退化而非默认已知。

层D `COHERENT_DEPTH_SCORING_AND_SUPPORT`：只有层C数值/接收可用才比较深度点、置信/候选区间及错误分支；同时剖面水平状态`r theta v psi`和左右舷候选，不能在主实验固定truth水平状态后报告深度工程精度。oracle水平状态仍可作为上限旁栏。

层E `CROSS_ENVIRONMENT_TRANSFER`：跨CZ2环境对照方法提取率/模态密度/观测孔径和不利环境角点；若只在单个CZ2环境通过，禁止写“跨海域稳健”。

## 三、径向合成孔径的正确定义

`Lr = |r(t_end)-r(t_start)|`，注意这不是“拖曳阵14m物理孔径”也不等于`platform_speed × elapsed_time`。相对径向运动受目标运动和平台航向影响；必须由观测侧估计/区间约束，不能在非oracle模型直接使用真实`r(t)`。如果`Lr`很小（横向相对运动或等速追随），不能靠标称1200s骗取大合成孔径。

对每个CZ环境，在预注册轨迹和先验下计算：
- `Lr`真实值仅供oracle上限/成本诊断；`Lr_available`来自可用观测及导航假设。
- `Δk`按可分辨模态对/能量组报告分布，不能只取平均间距或只选特别疏的模式。
- `N_modes_total`, `N_modes_energetic`, `N_modes_Rayleigh_unique`, `N_modes_extracted_correctly`, `mode_group_stats`分别定义。
- 模态能量阈值应按可检测噪声/校准和源谱未知量登记；任意降门槛使结果变好不被允许。
- 已有E-STD 2.4km Fourier `n_Rayleigh_unique=0`只对旧冻结环境/模式集合/数值边界有作用域，不变成对所有首CZ环境的通用结论。

## 四、环境与声压前向可信度硬门

历史T2失败的CHIGH完整性、复波数单精度、模态函数与理想刚性底可能严重影响50–60km相干复压力。本项目的新ENV-CAP不可直接沿用未经认证的旧复杂场当真值。CODEX可以在本阶段E2预算内**另建新提供器**（如适于对应模型的waveguide/PE/ray+validated eigenmodes），但必须事先冻结：底边界、模式/射线完整性、网格/步长收敛、独立对照、误差与观察量尺度关系、原始输出SHA、声音符号、所有失真/漏模。绝不默认“CHIGH越高越对”。

若CZ环境只具有ray/TL数值准入而不具相干复模态准入：允许评估ray几何、传播区间和基本Rayleigh必要资源界；不允许输出Yang/MMAR的声压提取成功率、深度误差、合成孔径相干增益。

## 五、零信息与统计约束

1. 模态完全自由的复增益吸收源深控制（同一源/组随时间的独立自由度须合理）。
2. source frequency unknown、node reference offset/drift与速度投影的nuisance-profiler；按同时发射时刻对齐，不能把接收时刻同步当作发射时刻同步。
3. 单通道未知每频复源谱不许宣称独立绝对距离信息。
4. 左右舷镜像分支保留，单HLA/倾斜x-z平面的同场镜像控制必须显示不可判的结论。
5. 时间衍生观测彼此相关，不能将同一段声压生成的bearing、radial、modal features当相互独立信息累加；缺联合协方差则互斥对照，不联合。
6. 非oracle估计器不可读真值范围窄邻域、已知每个模态的真实排序、生成随机种子；主要统计含拒判INF、假峰/错误舷/错误模态判别。

## 六、条件性能类型

标注`ORACLE_DEPTH_SIGNATURE`, `FOURIER_MODE_ORDER_RESOLVED`, `HIGH_RES_MODE_ORDER_RESOLVED`, `DEPTH_INTERVAL_WITH_COVERAGE`, `DEPTH_POINT_WITH_VALID_ESTIMATOR`, `CROSS_ENV_TRANSFER_STABLE`, `NOT_EVALUATED`。不得只报告误差最小的场景或在未知水平状态时拿known-r评分冒充实际估深。

## 七、如果算法受限

通过预注册困难/控制否定某条具体实现不代表否定整篇Yang或MMAR全部可能变体。记录准确因由，例如`MODE_ID_UNRESOLVED`, `TOO_LITTLE_RADIAL_APERTURE`, `SOURCE_FREQUENCY_GAUGE`, `OFFSET_DELTA_AMBIGUITY`, `PROPAGATION_NUMERIC_NOT_ADMITTED`, `SSP_MISMATCH`, `INPUT_SIGNAL_UNAVAILABLE`，区分机理不存在、机制有而不可提取、提取有而估计不准。若两种方法都不具可用提取链，仍完成边界/资源报告，并停止无根据的性能仿真。
