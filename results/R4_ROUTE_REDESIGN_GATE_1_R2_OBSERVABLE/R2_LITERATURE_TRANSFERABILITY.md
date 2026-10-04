# R2 文献机制到场景迁移

本轮读取项目原文与冻结记录；机制存在不代表本场景已有 observable。所有引文性能都属于论文条件，不是本项目指标。本轮没有做复现、特征提取或传播重算。

## 原文身份与读取范围

有效原始中文资料位于 `D:\ccc\博士资料\7.项目\2026-716单阵\论文\论文-声学学报-径向距`。各 PDF 的原始 SHA、页数和完整路径绑定在 [VALIDATION.json](VALIDATION.json)。本轮读取对应章节及本地 PDF 文本，保留纸面页码与 PDF 页码区分。

| ID | 原文/证据身份 | 本轮重点 |
|---|---|---|
| L1 | 徐嘉璘、郭良浩，应用声学 2024，43(2):237–251，DOI 10.11684/j.issn.1000-310X.2024.02.001；15 页，SHA d18e0bd1…29361e0 | pp.239–241 信号/垂直阵/到达结构；pp.244–247 自相关及二次相关；海试条件 |
| L2 | 徐嘉璘、郭良浩、任云，声学学报 2023，48(4):618–631，DOI 10.15949/j.cnki.0371-0025.2023.04.013；14 页 | 深度维 delay curves、Radon、自相关、宽带海试；[期刊原文](https://www.cpsjournals.cn/article/doi/10.15949/j.cnki.0371-0025.2023.04.013) |
| L3 | 刘与涵等，应用声学 2025，44(1):36–54，DOI 10.11684/j.issn.1000-310X.2025.01.003；19 页 | §§2.1–2.3 几何/孔径/TMA；§4.3 会聚区；§6 场景限制，PDF pp.3–5,9–12 |
| L4 | Yang2015，JASA 138:1678–1686，DOI 10.1121/1.4929748；[项目 PDF](../../literature/yang2015_Yang2015_JASA138_1678.pdf)，10 页 | abstract、§II、§III、appendix；已知源频率、CW、radial aperture、PLL、mode ID |
| L5 | Liang等，Match-Mode Autoregressive Method for Moving Source Depth Estimation in Shallow Water Waveguides，2018，7824671，DOI 10.1155/2018/7824671；[项目 PDF](../../literature/liang2018_mmar_7824671.pdf)，15 页 | HLA beam output、Hankel complex spectrum、AR mode-ID、depth score；与 Yang 不同方法身份 |
| L6 | Cockrell & Schmidt 2010，JASA 127:2780–2789，DOI 10.1121/1.3337223；[作者机构 PDF](https://acoustics.mit.edu/faculty/henrik/LAMSS/Pubs/cockrell_schmidt_jasa_127_p2780-2789_2010.pdf) | broadband intensity striation、水平采样/运动尺度、浅海假设；辅助解释 fringe 所需观测，不替代 CZ 验证 |
| L7 | 徐嘉璘等，利用深海海底声反射区频域干涉结构的声源深度估计方法；项目原始 PDF，12 页 | 首页面机制与阵列条件；幅度干涉可用不等于绝对相位已知，也不等于 CZ 距离锚 |
| L8 | Yang2014，Data-based matched-mode source localization for a moving source；[项目 PDF](../../literature/2014.pdf)，13 页 | 原文身份/摘要与 Yang2015 引用依赖；data-based mode replicas 及垂直覆盖限制 |

仓库未跟踪的 `literature/xu_jialin_2024_cz_tdoa.pdf` 是非 PDF 响应（header 非 `%PDF-`，SHA 0116495b…a6423），不是 L1。它完整保留，不改、不加入本轮提交。有效 L1 原文与旧 [PRIMARY_SOURCE_2024_LOCK.json](../RANGE_UB/CZ_TDOA_SINGLE_DEPTH_ORACLE_FIX/PRIMARY_SOURCE_2024_LOCK.json) 的 d18e… SHA 完全相符。因此旧初次检索“未获取 PDF”状态已被后续全文记录取代；本轮按后续材料审计，不能引用旧 pending 推测代替原文。

## 机制 → observable → 场景 → 缺口

| 文献 | mechanism / required observable | paper scenario | our scenario | transferability / missing condition |
|---|---|---|---|---|
| L1 CZ TDOA | delay-difference 随 receiver depth 的曲线斜率给 range；自相关/二次相关可增强某些 delay 结构 | 仿真 LFM 100–300 Hz、垂直 20–1620 m；海试爆炸源、26 元 VLA 111–1869 m | 单深度附近 HLA、UUV、几条 ideal tones、level-only A1 | 当前 HLA 没有 depth slope 采样；宽带/瞬态和可提取 lag 未建立。MECHANISM_INTERESTING_BUT_OBSERVATION_CHAIN_NOT_ESTABLISHED |
| L2 bottom-bounce TDOA | 4 路形成 delay curves；autocorr 与 Radon 从深度维结构提取几何参数 | 海底反射区，宽带/爆炸源、垂直接收；10–33 km 论文测试 | first CZ 45–60 km 条件域、short observation、无垂直接收孔径 | 原文方法当前 REJECT_CURRENT_SCOPE；不能换区名、把 depth slope 改成未知 UUV 的 time slope。aggregate delay 另列条件问题 |
| L3 mobile-HLA review | 按传播区域选择几何、MFP、多途与干涉特征；远程 TMA / CZ 结构需额外条件 | 多场景综述，包括不同 aperture、receiver depth、signal、maneuver | 当前限制固定，不能挑综述里的最强条件静默移植 | §4.3 的 CZ 边缘穿越方法需较长时间及运动条件；局部 fringe 可给态势，β/模态类/环境不同；CZ order ambiguity 尚需证据。支持问题分解，不提供本项目 acquisition 证书 |
| L4 Yang | 单水听器 CW 接收复信号 → Doppler/radial increment → synthetic aperture → mode spectrum → depth | 已知 emitted tone、phase tracking/PLL、足够 radial aperture；单点需 mode depth functions，VLA data-based 情况另分 | 未知 UUV f0/相位，短时小机动，A1 level-only，deep dense modes | 不能把绝对 source phase 视为已知；工程接收 complex signal 也不同于 model H。原实验 δ=0 oracle 对齐、Fourier mode-ID 限制仍在；非新冷启动距离 observable |
| L5 Liang/MMAR | beam complex output → generalized Hankel/AR wavenumber estimates → mode association/complex amplitudes → depth distribution | shallow-water waveguide，HLA + moving source aperture；阵列/运动/相干条件 | 第一 CZ 密集模态、未知运动、未冻结 raw coherent input | high resolution 有机制价值，但没有深海 cold-start range admission。AR peak height 不能当 modal amplitude，不能用 oracle mode IDs/offset。暂停 depth 路线保持不变 |
| L6 WI | broadband intensity 的二维 striation 提取与结构化 range 估计 | shallow water，单 receiver 向固定源移动，带宽/水平采样条件明确 | moving UUV、deep CZ、三个稀疏点、relative level、未知径向距离轨迹 | 无需 absolute phase 是正面特征；dense-band/运动尺度/β 转移条件缺失，不能直接复用性能或 β≈1 |
| L7 frequency interference depth | 幅度谱随 source depth 与垂直到达角变化，构造 depth-angle feature | bottom bounce、宽带 near-surface source、短 VLA，双弹海试 | HLA、未知深度 nuisance、first CZ | 幅度 interference 不必绝对相位；但其目标是 depth 且需要垂直到达角，不能升级为当前 absolute-range anchor |
| L8 data-based matched mode | 实测 mode replicas 减少环境失配 | moving source 与覆盖目标深度的 receivers，未分辨模态仍受限 | 单 HLA 无相应 vertical coverage，复谱/位移条件未建立 | 不能借“data-based”消除实际阵列与 motion requirements；仅解释 Yang 链的来源，不增加新路线 |

## 项目已关闭的结果不能绕开

[MMAC final](../R3_B1_MMAC_observability/R3_B1_FINAL_REPORT.md)：稳定 continuation 且有 modal energy 支持的 branch 仅 1，delay-only 独立维数为 0；旧 `B1_MMAC_PHYSICS_CONFIRMED` 已撤销。它否定所测 labeled MMAC，不能推出所有 aggregate delay 特征无效。

[range-UB closeout](../R3_FINAL_CLOSEOUT/R3_FINAL_RANGE_UB_STATUS.md) 分开 all-eigenray oracle 信息、single-depth extractability、VLA slope。最新 [single-depth decision](../RANGE_UB/CZ_TDOA_EARLY_LAG_AUTOCORR/RANGE_UB_1C_FIX2_DECISION.json) 的旧预注册 shape consistency 只有 4/12；这不是现实录音成功证据。改为 unlabeled aggregate 并不消除该失败，不能复活已关闭支线或做新 beam sweep。

[WI fidelity report](../R3_A1_fidelity/R3_A1_REPORT.md) 保留 `A1-NOT-DIRECTLY-TRANSFERABLE`：E-STD β_eff 非常数，经典 constant-β 直接迁移不成立，generalized/local fringe 留后备。L3 也按模态/声道/位置区分条纹；不应将不同文献的 β 约定或区域趋势拼成通用常数。本轮不从印刷式猜新 β 数值，不把“fringe useless”或“fringe 已可用”写入判定。

[Yang transfer](../R3_C2_Yang_SA_depth/R3_C2_3A/R3_C2_3A_REPORT.md) 中 L≤2.4 km 的普通 Fourier mode-ID 没有唯一模态；这不是所有高分辨算法不可能。现有 [DZ input contract](../R3_C2_AR_MMAR/R3_C2_4B_3A_HANKEL_AMPLITUDE/DZ_INPUT_CONTRACT.md) 仍有 δ=0 oracle offset，不能供工程 cold start。depth 保持 nuisance，不在本轮发展。
