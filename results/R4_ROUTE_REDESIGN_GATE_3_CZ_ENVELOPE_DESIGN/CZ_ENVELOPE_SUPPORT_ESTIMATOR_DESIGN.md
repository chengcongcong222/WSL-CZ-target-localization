# Observation-derived support estimator 设计（无实现）

输出为 retained range-motion support 和一组 range intervals，不输出 top-1，不用旧 recovered basin/truth/oracle cache 初始化。以下仅数学与伪代码；没有 search、forward 或 feature 计算。

## 域、观测与 score

沿用水平域 `r∈[45,60] km, theta∈[-5,5] deg, v∈[1,3] m/s, psi∈[-15,15] deg`；source-depth labels 150:5:250 m 为旧 profile nuisance，实际 modal mapping 必须明确。single HLA、三频、旧时间轴和平台导航；fringe/TDOA/Doppler/spatial 新观测一律不混入。

bearing likelihood 使用既有 RC2 rule（noisy `cmin+13.3 sigma²` 与旧 noiseless numerical control），不改 Gate。旧 cases 用 observation-derived 归档 cutoff；fresh 的 observation-only cmin 取得方式/上下界须 pre-run 冻结。若只能得 cmin 上界，可产生更保守 cutoff并注明，不能把局部最小当精确全局 cmin。不得借只覆盖已成功样本的 RC2 point cloud 充当整个 cold-start 域。

记选择的表示为 F，观测为 y_F。预测用候选状态轨迹与既定 forward model 的 level，按相同 F 处理。定义兼容性 `s_F(x)=min_{z in frozen profile,d in frozen N} ||y_F - F(L_pred(x,z)+d)||_F`，S=`{x: bearing-compatible AND s_F(x)<=tau_F}`。E1/E2 的线性情况可用 F(L_pred)+F(d)；E3 不行。tau_F 不依赖本 case 最低 score，也不依赖 truth。

E1 norm 为所有 freq/window 投影能量的均值平方根：`s²=(1/6)sum_fw sum_k residual_c_fwk²`，保留投影 level RMS（dB）。E2 用 `s²=(1/6)sum_fw sum_b(n_wb/N_w) residual_m_fwb²`。E3 是 normalized intensity feature 的 dimensionless block norm。不同 F 的数值和 exact-TL J 都不直接横比。

## 连续域覆盖与搜索复杂度

建议 future search 为**确定性、完整起始域的 cell subdivision**；目的先形成可审计粗 support，不重开 SHGO/DIRECT/DE repair。四维 cell 保留全部 r-v-psi-theta，profile nuisance 取 union；不在 truth motion 的 range line 上求解。

每 cell 必须记录 valid bearing lower bound、feature score lower/upper bounds、depth/drift profiling completeness、状态、预算/call counts。只有可审计的全-cell lower bound 大于兼容 cutoff，才可 reject。一个中心点 score 大、局部解失败或有限点样本不匹配不能排除整 cell。

feature enclosure 的可实现性是 pre-run 必须解决的项，不因 F 低维就默认可认证。可研究 interval/analytic Lipschitz bounds 或有明确范围的 approximation error enclosure；它们不得由观测 truth 或新成功样本倒推。模型内部用于生成预测幅度的运算不变成观测相位；若 fast modal oscillation 或 level near-zero 使 bound 太松，输出更多 unresolved support，不能换算毫米网格作为覆盖证明。

未解决 cell 一律保留在 conservative outer support。可另存已验证兼容 witnesses，但 point cloud 不是 support。预算耗尽也保留 unknown cells，标记 `SEARCH_COVERAGE_UNCLOSED`，不伪造科学 PASS。即使 outer width 恰好很小，未经预注册 closure/stability 检查也不能宣称搜索已收敛；外包络相连可能遮住内部 alias islands，须另报 witness/fine-scale roughness 与分辨率作用域。

建议固定两级 joint-state resolution / hard budget（具体值 NOT_FROZEN）；细级为各轴初始步长减半的 deterministic refinement，轴 bounds 不变、预算各自 strict wrapper计数。先定遍历次序和 split rule（最长归一化边，轴 tie-order 固定），不能见结果后添加 seeds/第三级。bound、终止精度、call cap、wall-clock policy 全部须负责人 pre-run 选择；本轮不猜一个能跑完的四维预算。

## 仅伪代码

```text
receive whitelisted observables and frozen contract metadata
compute the selected observation-domain coarse feature [future only]
cover the full inherited 4D domain; never seed from truth/recovered states
for each cell in deterministic order under an enforced hard cap:
    enclose bearing and coarse-feature compatibility over cell+nuisance
    reject only on a valid whole-cell incompatibility bound
    retain or split by the frozen resolution rule
    if unresolved or capped: retain cell and record incomplete coverage
export conservative joint cells, interval projection, bounds and counters
seal outputs; only then evaluator obtains truth and computes retention
```

## support 拓扑与 handoff

export outer cells 的 range projection 取 closed intervals union，不按 truth 距离排序、不丢小岛、不以最大 cluster 代替全域。报告 interval total length W_union 与 convex hull width W_hull、range interval count、joint-cell component count（固定 face adjacency，不用事后距离阈值连接）。若得到空集，报告 `EMPTY_SUPPORT`、retention=false，不把 width=0 当成功。

未来 Gate 用 W_hull 做主要收缩指标，避免两端小岛总长度很小却仍占 15 km。同时报告未解决体积/区间、component bounds、endpoint boundary censoring 与每候选的 nuisance。boundary 或 finite z profile 影响必须显式。

retention evaluator 分开：truth/reference 水平状态是否落在 joint outer support；truth range 是否在 projection；truth-compatible cell 是否已 resolve；true effective-depth 标签是否包含。单有 projection retention 不足以证明 joint motion 已保留。真值不在旧 RC2 support 时先归因 RC2，不能调整 envelope tolerance 或删除该 case。

未来不以 envelope support 自动承诺 downstream RC3 能捕获正确盆地；毫米 matched section 不等于局部 optimizer capture 域，coarse gate 通过最多说明 matched representation information。
