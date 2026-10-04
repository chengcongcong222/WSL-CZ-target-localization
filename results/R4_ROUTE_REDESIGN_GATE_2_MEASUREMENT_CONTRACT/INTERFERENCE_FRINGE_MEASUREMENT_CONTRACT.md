# Contract B — Interference Fringe（拟议最小被动合同）

状态 `PROPOSED_CONTRACT`；需新增 receiver spectral output 和真实 signal 条件说明，G2 未全部闭合。当前三个频点 **201/235/283 Hz 不是 dense fringe observation**。合同不要求已知 waveform、发射时刻、absolute source phase 或发射端配合。

## 数据与完整链

`target-bearing received intensity spectrum I(f,t) → source/gain trend separation → observable fringe/striation → slope/spacing/curvature or structured ambiguity → joint range/range-rate/environment candidates`。

最低交付物为单 HLA 上目标方向的 **adequately sampled dense intensity spectrogram**：已知 time/frequency axes、频带有效占用与缺测、PSD/强度标度、谱估计/窗口信息、bearing/beam 关联、接收相对频响/AGC/beam gain、噪声与质量标志。采样要能保留目标可观测的条纹变化，不仅把三个离散 line 插值成连续图。不给 Hz 带宽、bin spacing、SNR 或窗口数值门限。

要求足够频带内真实辐射能量，可以来自非合作 broadband machinery-like 内容或其他充分谱占用；这是 REQUIRED_SIGNAL_CONDITION，不是 UUV 必然满足的断言。不强制已知 chirp、爆炸源、active ping 或让目标发射宽带。若实际只有几根孤立线，dense 接收 bins 也不能创造源谱能量；合同不满足，不能用 receiver 升级掩盖 source 缺口。

可由审计过的 beam PSD producer 提供，不强制下游读全部 raw channels；若只有 raw HLA channels，未来上游 beam/谱处理需同步、阵形和相对校准，但本轮不实现。额外 waveform 保存有利于复核而不是 known-source 条件。无需 downstream complex phase/cross-frequency phase；beamforming 内部相对通道 phase 与 absolute source phase 分开。

## Source-spectrum separation

在目标与噪声可分或有噪声模型的条件下，可写 `log I = source log-spectrum + propagation log-intensity + receiver/beam gain + residual`。这只是机制模型；加性噪声并不天然在 log 域成为独立常数，低信号/缺测需质量模型，不能无限补零后解释条纹。

可能处理类别（只作定义，不实现）：frequency high-pass、local detrending、normalized spectral derivative、2D time-frequency slope、cepstral/autocorrelation-like periodicity。它们都要求源谱/接收 trend 与传播条纹有可分离结构；smooth source spectrum 只是假设类别，需要真正的 scale separation，不能凭“平滑”两字设为已知。源自身谐波、机械 comb、moving source-spectrum peaks 也能产生周期/斜纹，应当进入源谱 nuisance 或失效判断。

未知共同 source level 可归一化；任意 S(f,t) 不能自动消除。频率去趋势可能删除 broad propagation 信息，时间 differencing 也依赖 stationary/受约束 source spectrum。需要保留 preprocessing 和被删 trend 的记录，未来检验处理是否去掉了 range 信息。不能先看 truth 再选择频带、删峰或改 trend。

## Motion coordinate 与 absolute-range 机制

平台导航已知不等于目标的 `R(t)` 已知；目标与平台均移动，时域 slope 常含 `range-rate/range`、β-like/模态环境结构。将时间轴乘 truth radial speed 得到 frequency–range map 是 privileged input，不准使用。

未来 joint variable 至少包括初始 range、目标运动变量、source depth、receiver-depth uncertainty、环境/β-like structure、fringe family/order 和 source/gain nuisance。β 的有效结构受 CZ 模态与环境限制，不能自由同时吸收所有 range 缩放，也不能硬设浅海 constant β。旧 constant-β CZ 迁移失败保持冻结。

潜在绝对尺度来自给定条件环境下的 fringe spacing/curvature/组合随**物理距离**的预测关系及跨频/时间一致性；它不是孤立 slope。frequency periodicity 可对应 path/modal delay 组合，环境/深度与不同阶次带来多值；range feature 必须保留所有符合的区段，并与 truth-free bearing/motion 共同检查。若仅剩未知 β×range-rate/range 比值，绝对距离机制没有闭合，不能输出唯一 range。

可以形成 lower-dimensional、ambiguity-structured slope/spacing proposal，而非 exact pressure/TL 全域匹配。但这种 searchability 仍是机制假设；本轮不生成条纹、不提特征、不估 β、不算 coverage 或距离候选。

## G2 与未来证伪（不执行）

G2-1 未闭合：dense output 与源谱现状未知；负责人尚未接受该 proposed extension。G2-2 被动、G2-3 条件信号物理可讨论、G2-4 nuisance 明确、G2-5 条件绝对尺度机制、G2-6 可证伪，均只为结构推理，非实际能力/数值 PASS。

以后另行授权才可设计信号趋势可分性、sparse-versus-dense、未知 motion/β compensation、源深/环境失配与阶次 ambiguity 的检验。如果无可观测 striation、所谓条纹其实是源谱、feature 只给比例而没有 coarse range 信息、所有候选仍覆盖 RC1 或处理必须依赖真值，则 acquisition 假设可 FAIL。本轮不定参数/性能门限，不开启 experiment。
