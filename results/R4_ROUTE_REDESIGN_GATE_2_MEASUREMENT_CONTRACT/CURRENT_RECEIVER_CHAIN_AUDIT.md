# Gate-2 当前接收链审计

基线：`9e50b0d949523c85fb35c1168ff77c37ddf37aeb`，研究负责人已接受 Gate-1。本轮仅冻结接收/信号合同草案；没有实施采集、处理或实验。

当前合同可确认的范围为 **M0：模型生成的 bearing + 三频、分窗口去均值 relative level**。实际接收系统的 M1 输出与校准合同均 `UNKNOWN / NOT_FROZEN`。两份被动合同可以具体定义，不要求发射端配合；它们是 `PROPOSED_CONTRACT`，需要负责人接受扩展后再决定 quantitative design。目前 Outcome B，定量设计未放行。

## 三个层级与证据身份

| 层级 | 内容 | 本项目状态 |
|---|---|---|
| M0 | A1 已实际消费的 bearing_rad、relative_tl；频点 201/235/283 Hz，按频/窗口去均值 | `CURRENTLY_SUPPORTED` 仅指冻结软件/合成观测接口；不是实测系统能力认证 |
| M1 | 同步阵元波形、目标 beam 波形、PSD/频谱、频率轨迹、自相关等可能接收产品 | 只有原则上可提出；本次没有查到设备输出、实测数据或标定验收。只能 `PROPOSED_MEASUREMENT_CONTRACT` |
| M2 | 已知 waveform、绝对 source phase、任意已知 source spectrum、已知 emitted carrier、active ping、发射相干同步 | `NOT_ADMISSIBLE_FOR_CURRENT_PASSIVE_SCOPE`；两个推荐合同都不用这些条件 |

`CURRENT_RECEIVER_CAPABILITY` 不能由模拟器字段、论文阵列或“拟利用目标声信号”的方案文字推导。当前未确认任何 M1 项为该标签。已有接收场景/导航环境假设也不等于设备日志已经交付。

## sensor → sampling → synchronization → calibration → beam → spectrum → archive

| 检查项 | 已找到的证据 | 可冻结的结论 |
|---|---|---|
| sensor / 单 HLA | 原始阶段性讨论稿 V0.2 及 project_context；实验代码 | 单平台被动单 HLA 是任务约束；实际设备型号、阵元数、阵形产品未冻结 |
| element spacing / aperture | [历史 HLA 对照](../R3_A1_fidelity/HLA_MFP_COMPARE.csv)：8×2 m/14 m；8×10 m/70 m 仅诊断；旧 P3 另有模型配置 | 这些是互不替代的模拟配置，不是装备规格；实际间距/孔径 `UNKNOWN / NOT_FROZEN` |
| sampling rate / acoustic ADC | A1 TIMES 为 0–1200 s、10 s 间隔，121 点 | 这是特征时间轴，不是音频采样率；ADC、抗混叠、位深、连续覆盖均 `UNKNOWN / NOT_FROZEN` |
| 阵元级 raw waveform saved | 扫描仓库信号文件及生成/读取 API，没有发现实测阵元波形产品 | `UNKNOWN / NOT_FROZEN`；不能证明系统从不采，也不能写系统已保存 |
| beamformed waveform saved | 文献/方案提出目标方向增强；旧模型计算 beam 输出 | 实际保存接口、beam 权重/指向、增益与处理版本 `UNKNOWN / NOT_FROZEN` |
| synchronization | 有模拟时间轴；无时钟规格/丢样记录/通道对齐验收 | 实际阵元同步、导航对齐、频率时基 `UNKNOWN / NOT_FROZEN` |
| amplitude calibration | 旧 relative level 消去每频常数；5B-FIX 为注入 stress | 没有实测频响、AGC、beam gain、参考声压或校准记录；实际幅度合同未冻结 |
| phase calibration | 旧模型生成 complex pressure/相位 | 实际通道相位/阵形校准 `UNKNOWN / NOT_FROZEN`；模型相位不可给接收端 |
| beamforming / bearing | A1 消费 bearing；观测由 geometry + noise 生成 | bearing 接口 M0 有证据；不是已验收真实 beamformer。镜像/多峰/质量标志应在拟议产品中保留 |
| frequency band | M0 三个频点；旧四频/宽带与论文频带属于其他条件 | 不存在真实 UUV dense-band 占用或设备可用带宽证书；不新增带宽数值 |
| stable line detector | 5B-FIX 有 tracked-frequency 模型控制；A1 频点预选 | 无实测 detector、关联/失锁/重捕获接口证据；ideal tracking 不等于已实现 detector |
| Doppler track | 旧 P3 frequency 模型，当前 A1 API 无频率轨迹字段 | 实际独立频率流 `UNKNOWN / NOT_FROZEN`；不能从 nominal line label 得到 emitted f0 |
| broadband transient detector | 旧论文宽带/爆炸源，模型 delay control | 本项目真实瞬态、detector 与事件波形 `UNKNOWN / NOT_FROZEN` |
| dense spectrum saved | 当前 NPZ 没有 dense frequency axis/PSD 字段 | `UNKNOWN / NOT_FROZEN`；三频 level 不是 dense spectrogram |
| saved observable | [A1 代码](../../r4_a1_offgrid_bearing.py)、[development 读取](../../r4_a1_search_tractability_development.py)、两个原有 NPZ 的 ZIP/NPY headers | bearing、relative_tl、times_s、case/panel IDs；只检查文件头，不运行科学模块或处理数据 |

`CASE_OBSERVATIONS.npz` 的 bearing 为 `(36,121)`、relative_tl 为 `(36,3,121)`；这些是旧归档形状，不是新观测。归档不含原始波形、源谱、传递函数相位或音频 PSD。其他旧 complex/modal/arrival 文件属于传播模型/论文复现控制，不补实际接收合同。

## 原始项目资料的独立检查

外部资料根：`D:\ccc\博士资料\7.项目\2026-716单阵`。完整路径与原始 SHA 见 [VALIDATION](VALIDATION.json) 的 external_project_bindings；不改这些资料，也不把旧阶段边界当成本次 Git 授权。

- `深海单拖曳阵潜艇目标状态约束方案_阶段性讨论稿_V0.2.docx`：第一节把单 HLA 和目标方向信号作为任务理解；第六节明确阵元数/孔径/频段、历史目标噪声/频谱/录音与观测时长待确认。前文技术路线是计划，不能升级成交付数据证明。
- `docs/questions_for_client.md`：阵元级或 beam–frequency–time 数据、阵形/相位校准、稳定线谱/宽带瞬态及频稳性仍待问证。不能把提出的问题当答案。
- `docs/project_context.md`：导航、SSP、海深、阵深/阵形是场景可用信息假设，未给接口、同步或误差验收。
- `docs/observation_models.md`、`measurement_models.md`、`feature_catalog.md`：明确是暂定候选因子/关系；源频、阵形、相位、环境是 nuisance，不能当工程产品清单。

未将讨论稿内论文性能或旧第五参数目标移植为本轮成功证据；没有拓展 depth、SSP 或传感器数量。

## 输出扩展的最小边界

Contract A 最低希望接收端交付 target-bearing 上的多条可关联谱线的 **时间强度/谱级** 及 time/frequency、gain/质量元数据；Contract B 需要有足够谱占用的 **dense intensity spectrogram**。若这些产品由接收端可审计地产生并保存，不强制把所有阵元 waveform 作为 R2 消费输入；若系统只能给 raw HLA waveform，则未来需在同一单阵内形成这些产品，不能宣称已实现。

波形留存可以提高复核能力，但不是以“没有阵元波形”为理由否定一切幅度合同。M0 的 temporal relative-level proxy 可用于以后另行授权的表示假设检验；它不要求新增物理传感器，却不能证明真实接收链成立。提出 pre-normalization level/PSD 是为了使源/噪声/增益与归一化可核查，不是声称数学上只有绝对声压才能做 CZ envelope。

见 [合同矩阵](R2_MEASUREMENT_CONTRACT_MATRIX.csv)、[信号条件](R2_SIGNAL_CONDITION_MATRIX.csv)、[Contract A](CZ_ENVELOPE_MEASUREMENT_CONTRACT.md)、[Contract B](INTERFERENCE_FRINGE_MEASUREMENT_CONTRACT.md)。
