# CZ envelope 表示候选（设计草案，不运行）

基线 `dcad14337491ac5bf90acacf0be34a57eafac3f1`。负责人接受 Contract A **仅用于 quantitative Gate design**；真实 M1 接收能力与 UUV 信号仍未验证。fringe 未准入。本轮只定义候选，没有对旧观测计算一个 feature。

## 输入与固定坐标

`L_f(t)=a_f(t)+g_f(t)+P_f(R(t;x),z_s,z_r,e)+n_f(t)`。允许输入仅为 received level/intensity、频率/时间身份、bearing 与平台状态。未知源谱、gain、深度/环境作为 nuisance。使用候选状态的 `R(t;x)` 预测幅度是未来 forward 操作，不是使用 truth range；观测表示只能沿时间轴处理，不能按真值距离重采样。complex propagation/source phase、path IDs 不进入 feature。

M0 有三个频点 201/235/283 Hz、121 个特征时间样本。继承代码 `N_W1=61`：W1 indices 0–60（t=0–600 s），W2 indices 61–120（t=610–1200 s），不重叠，不在转向处补值。每窗口以已知 first/last sample 定义 `u=(t-t_first)/(t_last-t_first)`，以 `1/N_w` 均匀样本权重定义 inner product。M0 时间步不是 acoustic ADC rate；这些固定坐标不证明窗口内看得到 CZ edge。

定义 `C_w L=L-mean_w(L)`。已 centered 的 M0 再应用 C_w 应幂等；不会复原被删的常数或跨窗 pedestal。下列三项均可从 M0 relative levels 表达，不需 dense PSD、known waveform/phase/carrier 或跨频已知源谱。未来观测与预测必须用同一个 C_w 和同一固定表示。

## E1 — WINDOW_SHAPE_Q2（推荐草案）

在每窗口中，将 `1, s=2u-1, s²` 按上述 inner product 做符号定义的 Gram–Schmidt，得到与常数正交且单位范数的 `b_w1,b_w2`。feature 为 `c_fwk=<C_w L_f,b_wk>_w`，k=1,2。三个 line × 两窗口 × 两系数，12 个系数；不拟合新的时间节点，不从数据选择多项式阶数。

粗化动作是投影到整窗口的最低两个非恒定形状方向（trend、broad curvature），并丢弃其余时间频率结构；不是 pointwise exact-TL residual。符号 Parseval 说明投影不能创造新信息，未保留的高阶形状可能恰是 range 信息。窗口独立 centering 也可能删除 broad offset。

潜在 range mechanism：环境/收发深度条件下，CZ broad edge/maximum 的位置、宽度和跨线变化，会影响候选 `R(t;x)` 穿过该区段时的 trend/curvature 组合。必须同时保留未知 v/psi/theta 和 source-depth nuisance；线性传播片段可能只给 range–velocity 的组合。如果 coefficient 在不同 range 下经 nuisance 后相同，即 FAIL，不以“有二阶项”保证 acquisition。

优先理由：贴近 M0/Contract A；对 static per-line offsets 严格不变；只保留低维幅度形状；可显式看清 N1 linear drift 吸收的方向；不需 bandwidth/edge detector 新假设。推荐不是最优性能声明。

## E2 — LINE_BLOCK_MEANS_B3

每个窗口按 u 的固定 thirds 分三块：`[0,1/3),[1/3,2/3),[2/3,1]`；边界归右，最后端点包含。feature 为各 line 的 `mean_block(C_w L_f)`。保持实际样本数 `n_wb`；三块满足 `sum_b n_wb m_fwb=0`，只有两个独立对比，不将冗余维度当额外信息。评分用 `n_wb/N_w` 权重。

粗化动作是早/中/晚的 broad-event average。保留候选区段的变化与跨 line 一致性，未声称能定位 peak 的瞬时位置。它不取局部尖峰、不估 fringe slope、不做 delay/Doppler。块内 exact-TL 振荡可能未平均干净；mean shape 也可能被 source drift 吸收。粗块长度来自 observation window 的分段，而不是从毫米局部 section 推出。

## E3 — NORMALIZED_INTENSITY_BLOCKS_B3

由 level 定义正强度代理 `q_fw(t)=10^(C_w L_f(t)/10)/mean_w[10^(C_w L_f/10)]`。它是 received intensity 的归一化代理，不是重建绝对声压或传播相位。每 line 独立归一化后 `q_w(t)=(1/3)sum_f q_fw(t)`，再取同 E2 三块的 `mean_block(q_w)-1`，两窗口各两个独立对比。

粗化动作是跨 line 非相干强度平均 + broad temporal blocks。先逐 line 归一化，static source offsets 在代数上消去，不要求源 line ratios 已知；不能先把未知源级不同的裸 power 相加。非线性变换可能加强尖峰、丢掉有用的跨 line 差异；gain/source time drift 仍保留。潜在机制是多线共同的 CZ broad enhancement/edge，未证明它在当前三频和短窗存在。N1 需在 level 域加入 drift 后完整重做此定义，不能直接加 coefficient drift。

## 参数范围、敏感性与未闭合尺度

所有配置均 `DRAFT_FOR_RESEARCH_LEAD_SELECTION`，不是本轮宣布冻结。负责人下一轮选择 **一项 primary**：建议 E1 q=2；E2/E3 为替代草案，不允许运行三项后取最漂亮者。若预注册 sensitivity，仅 E1 q=1（curvature ablation），或已选 block family B=4 作为唯一块尺度对照；事先选择其中一个，敏感性不得替换 primary verdict。三个 line 等权、两个窗口等权；不新增频率、不挑窗口、不按 truth 改权重。

现有 observation windows 为时间坐标提供依据，但没有证据给出 CZ coarse envelope 的有效平滑时间尺度。对 block 分段能否去掉振荡且保留 range 的物理充分性，明确 **`SCALE_SELECTION_NOT_YET_JUSTIFIED`**。E1 q=2 是最小 trend+curvature 的结构草案，也未验证其物理充分性；下一轮 freeze 必须接受这个可证伪限制，不能将其写成已知合适 bandwidth。

最多一项 primary + 一项事先选择的 sensitivity；不做 continuous smoothing sweep、多个阶数/权重试到通过或把失败支线改为 fringe。当前所有定义均未执行。

参见 [nuisance](CZ_ENVELOPE_NUISANCE_MODEL.md)、[support](CZ_ENVELOPE_SUPPORT_ESTIMATOR_DESIGN.md)。
