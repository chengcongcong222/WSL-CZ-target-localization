# R2 Gate-1：接收链与信息链审计

结论：当前没有候选同时闭合 G1-1–G1-6；`candidate_admitted_for_quantitative_design=false`。CZ envelope 与 interference fringe 是优先补齐信息链的条件候选，不是已建立的距离锚。R1 cold start 仍未建立；tracking 只保留 post-acquisition 条件角色。本轮无数值实验、特征提取、模型运行或 estimator 实现。

## 接收能力、归档输入与模型内部量

单平台、单 HLA、first CZ、被动 UUV、短观测和小转向是继承场景。硬件可能采集阵元时序，不等于当前证据已提供同步时序、阵形标定、可用频带和目标信号条件。这里分别审计“物理上可测的统计量”和“当前已归档/冻结的 estimator 输入”。

| 层 | 已有证据 | 本轮允许的解释 |
|---|---|---|
| bearing sequence | A1 API 与冻结 observations | 已有观测输入；不能反推出完整 spatial covariance |
| relative level | `r4_a1_offgrid_bearing.py` 的 `Observations`、`direct_features`、`demean_windows`；development API | 201/235/283 Hz，每频、每窗口去均值的幅度轨迹；synthetic matched upper-bound observation |
| 实际谱线 | R3 working envelope 明确 ideal trackable tones，不是通用 UUV 实测谱 | 当前假定可跟踪；真实线谱稳定性、源谱及接收处理增益仍需信号条件说明 |
| 原始阵元信号 | 当前 A1 输入没有该字段 | 未建立可审计采集接口；不能由 bearing/relative level 重建 raw signal |
| dense-band intensity、到达结构、Doppler frequency track、cross-spectra | 当前 API 未提供 | 不能因为声学系统可能具备就当作已可用 |
| complex modal pressure、mode IDs、ray arrival/phase | 传播 forward 或 oracle control 中存在 | 模型内部量，不可直接给 estimator；真值只准评价 |

A1 的声压先取模、取 level，再去均值；归档信息不包含可供后处理恢复的绝对传播相位。符号模型为 `X_m(f,t)=G_m(f,t) S(f,t) H_m(f;x,z,e)+N_m(f,t)`。它只是信息来源的逻辑说明，不是新仿真或已验证接收模型。

## 五种相位必须分开

| 相位 | 接收链中能否分离 | 当前限制 |
|---|---|---|
| absolute source phase | 未知被动源不能默认已知 | 不准把仿真源相位当测量 |
| propagation phase | 接收复谱中的 source phase 与 propagation phase 叠加 | 未知源不能直接取得 `arg H`；已知 waveform/reference 才可能额外约束，当前未授权 |
| array inter-element phase | 同频、同步、校准、同源 cross-spectrum 可消去共同源相位 | 是空间相对相位，不是绝对传播相位；当前 A1 未归档 cross-spectrum |
| cross-frequency phase | 两频的源相位一般不同，作差不自动消去 | 需要跨频源相干/波形结构和通道标定；少量独立 tones 不满足该假设 |
| temporal phase continuity | 稳定 tone 可条件性跟踪接收相位 | source drift、运动、介质变化及 PLL 影响必须处理；不能等同已知源相位或 exact path phase |

Yang2015 使用接收 CW 的合成孔径、已知原始频率估计径向位移，以及 PLL；这与从模型直接拿 modal phase 是不同的信息链。当前没有把这些条件补入 A1。

## 源级与源谱消除的边界

设 received level `L_f(t)=A_f(t)+B_f(t)+P_f(r(t),z,e)+error`，其中 A 为源谱 level，B 为接收/波束增益。每频每窗口去均值仅消去该窗口内的常数项，源谱漂移和增益变化仍会进入特征。constant per-line unknown level 可以消除，不代表任意未知源过程都被消除。

`L_i-L_j` 中仍有 `A_i-A_j`；未知跨频源谱可改变 frequency ordering、ratio 和 slope。若再对该差作窗口内 temporal centering，可去掉稳定的 per-line offset，但会丢掉静态谱形，也不能补回 dense fringe。源谱恒定、共同幅漂或接收标定均是需说明的条件，不是本轮默认为真。

平滑或差分是已有观测的确定性变换，不能创造独立信息；R2 允许换 representation，是为了使已有距离信息更容易利用。其是否保留足够 range structure、是否改善 searchability 都需要以后验证。不得用已知真值距离平滑数据，或只平滑 forward 的 range axis 却把未作同样处理的时序当观测。

## 六条信息链及断点

每条链的最后两步都只是未来角色：保持多个候选与 nuisance，再 conditional/local RC3；本轮不产生 proposal，也不调用 RC3。

### R2-A：CZ intensity/envelope

`raw beam signal → tracked-line power/level → per-line constant-source/gain removal → observation-domain coarse envelope/edge statistic → range/motion candidate support → conditional local RC3`。

断点：实际稳定源与增益条件未冻结；去均值后 coarse envelope 的绝对距离相关性、短窗口能否看到 CZ edge/peak、与运动/深度/环境的混淆尚未建立。频率平均需要处理源谱与可用频点权重；强度域平均与 dB 平滑不能无说明互换。用已知 r(t) 构造 envelope 是 oracle 路径，不准使用。

无需已知 waveform 或 absolute phase；未知常数源级可去掉，时变源级不可自动去掉。relative normalization 可保留形状，也可能去掉 broad pedestal；平滑可能同时去掉判距结构。km/百米 anchor 只是目标尺度问题，未获证明，也未设阈值。若输出仍只有“第一 CZ、45–60 km”，就只是 RC1 重述，不能算新 acquisition 信息。`CONDITIONAL`。

### R2-B：multipath delay difference

`raw passive signal → observable autocorrelation/cepstral delay structure → source-autocorrelation/noise nuisance separation → unlabeled lag set or aggregate lag signature → joint range/nuisance candidates → conditional local RC3`。

断点最早在 raw signal 与 signal structure：当前只有几条 level 轨迹，不能提供 path-resolved delay。未知 waveform **不必然**阻止 autocorrelation；必须有足够宽带/瞬态、近似稳定信道及可分辨副本。若用 matched filtering，则需要已知/可恢复 waveform。tonal 信号的周期性相关和未知源谱使时延模糊，四条线不能默认等价于宽带到达序列。

逻辑上 `R_x(τ)` 含各 path pair 的 source-autocorrelation shifted copies。无须知道发射时间即可存在 delay difference，但未知源自身的相关峰与 path-pair 重叠仍可能混淆。aggregate lag/power-delay feature 不需要持续命名 C0/C1，因此不与旧 MMAC 的稳定 path labeling 完全同题；它仍需要独立的可提取性与数值稳定性证据。旧 single-depth aggregate autocorr 正因稳定性 Gate 未过而关闭，不能用“不标路径”绕开。

Xu2024 的 autocorrelation/secondary correlation + delay-versus-receiver-depth slope 有原文支撑；它仍依赖垂直采样，单 HLA 不能套其 `1/R` slope。all-eigenray oracle positive 只说明物理信息存在。保留 `CONDITIONAL_PENDING_SIGNAL_STRUCTURE`；不重开 B1 或已关闭 single-depth 支线。

### R2-C：frequency–range fringe

`raw signal → dense-band intensity/time-frequency map → source-spectrum/gain treatment → local striation slope/period/distribution → ambiguous range/range-rate support → conditional local RC3`。

当前可读三个离散 level，不等于可观测连续 fringe。三个稀疏且非均匀频点（旧四频也一样）没有建立可估计的 fringe period 或 slope；本轮不宣称所有离散频设计都不可能。无需 absolute source phase；需要目标在足够频带内有能量、可分离传播结构，并说明稳定源谱、采样和运动距离坐标。

时间条纹 slope 常耦合 `range-rate/range` 与 β；未知目标径向速度、β/环境条件时不能单独反演绝对 range。frequency period 也可能给 multipath delay 组合，不自动给唯一距离。周期、多模态、depth/环境混淆必须保留。幅度图可按局部 slope/period 建立 ambiguity-structured feature，是 searchability 的机制期待；目前未证明比 exact-TL 更平滑或覆盖更好。旧 constant-β 直接迁移被拒，generalized/local fringe 保留 `CONDITIONAL`。

### R2-D：ordered/differential frequency

`raw signal → ordered tracked-line levels → temporal centering / defensible spectrum nuisance treatment → signs/contrast changes/local spectral statistic → coarse multi-candidate proposal → conditional local RC3`。

统计量条件可测，不需要 coherent phase。单次跨频 ordering/ratio 不消除未知 per-line source level；只消除共同标量 gain。可讨论 temporal differential contrast，但它需要源谱足够稳定、通道增益可审计，且依然来自旧幅度信息。稀疏 adjacent line 差不是 dense spectral derivative。

order/sign feature 可降低表示精度和评分维数，却也可能丢掉距离信息，保持周期 alias 或形成大片同码区域。不能把同一 exact-TL residual 作差后仍精确全域匹配称为新 coarse proposal。range relevance、nuisance 后的 retained support 与 proposal capture 均未建立。`CONDITIONAL`，作为 A/C 的辅助 representation，不独立宣称 range anchor。

### R2-E：Doppler/radial velocity

`raw tonal signal → tracked received frequency → unknown emitted f0/drift and clock nuisance → radial-velocity/frequency-change constraint → bearing+kinematic candidate support → conditional local RC3`。

received-frequency 可条件测量，但当前 A1 没有独立 Doppler 时序。未知 emitted f0 可吸收恒定 radial velocity 频移；机械频漂不能默认全部来自运动。已知 f0、多个稳定 tone 的共享 fractional shift 或可验证的变化结构可能补信息，但都须明确条件；nominal 201/235/283 label 不是自动已知 emitted-frequency 证书。

Doppler 提供 range-rate/velocity，不能单独叫 absolute-range anchor。与 bearing、小转向几何联合，理论上可帮助打破部分 r-v compensation；旧 P3 恒定径向速度对的增量近零且已降为 historical simplified result，没有现成 cold-start 成功证据。kinematic/nuisance 表示可能较低维，但短观测下独立信息未建立。`CONDITIONAL`，辅约束，不为第一距离锚。

### R2-F：array spatial

`raw synchronized HLA channels → calibrated cross-spectrum/covariance/beam spectrum → source-common phase/power normalization → curvature or resolvable modal/spatial signature → range-sensitive candidate support → conditional local RC3`。

same-frequency spatial phase 可条件测，不需已知 absolute source phase；其归一化可去同源共同幅度/相位，不能去 sensor mismatch、noise 或不同路径结构。far-field plane-wave 相位主要给 bearing/direction cosine；曲率项随 `aperture²/range` 缩小，不能凭“有相位”假定远程距离可测。这里没有新代数数值、装备孔径阈值或实验。

旧 HLA 8×2 m（14 m）及 8×10 m（70 m，仅诊断）MFP 结果没有窄距离锚；不能把诊断孔径当装备保证。Yang/Liang 属 temporal/modal chain，主目标是 depth；需要 coherence、radial displacement、mode ID 与环境函数，当前尚有 Fourier mode-ID 和 oracle offset 条件。未经闭合便从模型拿 mode phase/ID 是 privileged input。换 covariance 名字后重跑 exact coherent MFP 不满足 G1-5。当前额外 range-anchor variant `REJECT_CURRENT_SCOPE`；不扩孔径、加 VLA 或开放 depth。

参见 [文献迁移](R2_LITERATURE_TRANSFERABILITY.md)、[候选排序](R2_CANDIDATE_RANKING.md)、[机器矩阵](R2_OBSERVABLE_ADMISSIBILITY_MATRIX.csv)。上述断链都不能以 truth、已恢复状态、人工窄先验或另一个 optimizer 补齐。
