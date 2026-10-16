# 第一会聚区环境族物理准入协议

## 目标与起点

第一会聚区必须**由环境物理证据确立**；“50–60km”“5km海深”仅供历史对照，不能锁死为首CZ位置。浅海Yang正控制不得混入主CZ。主要环境族可按SSP、声道轴、海底形态、源收深度或地理来源分类，但一族不同扰动档不能虚增为独立环境。

## 候选源池（E0时固定ID，允许未取得资料的负例）

1. `ENV_ESTD_IDEAL`: 仓库现存E-STD Munk-like/5000m/刚性底——**数值历史问题已知**，只是名义理论对照，尚不是具可置信海底的物理环境；第一CZ须重证。
2. `ENV_SCS_WOA_17N118E`: 仓库`m452_build_and_cz.py`、`m453_pipeline.py`、`m454_closure.py`的WOA23 T/S源网格重建+ETOPO2022水深路线。先查D盘缓存、正源/版本/垂向坐标和输出；源码存在不等于数据和结果已认证。
3. `ENV_SCS_WOA_18N114E`: 仓库`B3D3_M4_pack_extracted/README_B3D3_ENVBUILD_PREFLIGHT.md`中18N114E/April气候场、深层约束拼接与Snell积分压力试验。**这是派生诊断场，不是实测SSP或已准入CZ**；缺实际底层资料则只能作低证据模型。
4. `ENV_INDEPENDENT_REGIONAL`: 可获得原始元数据的另一个季节/地理区/文献现场SSP候选，由E0**开始前固定择取规则和资料来源**，不能挑性能最好的环境；若数据拿不到，应标未准入不凑数。
5. `ENV_NO_CZ_NEGATIVE`: 浅海Yang原文、常速或其它物理非第一CZ控制。仅用于检测器伪报和论文算法复现，不参加主CZ有利环境统计。

先以水深/SSP物理意义和目标深度合法性筛选；不能为了让每个族都出现第一CZ而调声道/底质/源收深度。新来源的坐标和资料必须事前登记，源/收深度可采用预注册离散对（例如150/200/250m，供候选筛选，非最终强制配置），目标深度超出水深者直接不可用。

## 公开地球物理数据

- NOAA World Ocean Atlas 2023用于T/S气候平均和标准差/标准误；以TEOS-10/GSW在正确盐度、温度、纬经度、压强单位下生成声速；**WOA是气候分析格点，不是指定时刻CTD实测**。
- NOAA ETOPO 2022可以给海深地形的区域栅格，但不是局地测深或海底声学参数。
- 保留来源、标准深度、时间合成类别、插值点、地理格点、质量掩码、海底全深采样情况。若WOA浅层/深层合成需拼接，保存拼接依据及连接处导数，称`RECONSTRUCTED_SSP_CONDITIONAL`而非`MEASURED_SSP`；不允许“为出现首CZ”而任意拼接。
- 底质/沉积层声速、密度、衰减与粗糙度若未知，名义物理对照与有依据的范围敏感性分开；不能让刚性底成为“真实底质”。

## 数值物理准入两步走

**步骤A：射线/几何机制**：选定range-dependent或range-independent 2D剖面，BELLHOP或独立Snell积分产生ray trajectories/turning points；记录一次源附近传播后的**最早折返聚焦区**，要区分深声道折射回返、表面反射和底反射主导的能量峰。折射线族及聚焦路径拓扑必须可追踪，不以一个Range-TL高峰替代。

**步骤B：空间场证据**：沿距离和接收深度形成声强/TL剖面，观察步骤A首个聚焦区是否确有能量增强，保存峰宽、峰形、声学损失、相邻区域基线和差异。采用E0预冻结的适度改变ray launch counts、beam spacing、step sizes和range/depth grid以审查数值稳定；明确BELLHOP/其他solver输入单位及BOX `[depth_m, range_km]`，复核仓库曾出现的误把70km写成5000km的问题。对非稳定峰标`CZ_NUMERIC_UNCERTIFIED`。

推荐距离搜索范围以E0冻结的10–120km粗网格为起点（不是预设会聚位置），发现候选后只在预注册的局部细化窗口做细格，不能一边观察一边扩大到出现希望峰为止。保留无CZ样本和越界CZ样本。

## 可执行Gate需要包含的字段

`env_id,source_id,source_kind,source_SHA,lat,lon,season,depth_bottom_m,source_z_m,receiver_z_m,frequency_hz,water_speed_min/max,ray_turning_evidence,ray_family,first_CZ_range_km,CZ_width_km,TL_prominence_dB,competing_surface_bottom_paths,numerical_convergence_status,independent_check,classification,missing_inputs,allow_methods`

CZ2的**数值阈值**必须在E0读取现有工具设置与噪声/数值底线后锁定。不得现在假造“一律峰值>3dB=物理CZ”的万能定义。分类不能完全依赖阈值，必须记录射线折返身份和能量空间图。独立对照如果只是同一个solver修改束角，应标`SAME_FAMILY_REFINEMENT`，不能标独立solver。

## 环境敏感性与真实性

对每个已通过CZ2的名义族，另做少量**有依据**的SSP/水深/底质角点（数据自带标准差、季节差异或文献有来源的参数范围）；报告CZ位置移动、丢失、折返身份和传播强度变化。角点测试不能升级为“所有海域鲁棒”。实际没有水深测量或底质，则标`CLIENT_GEOPHYSICAL_DATA_REQUIRED`。

## Yang环境正控制

原Yang约88m海深、350Hz量级、源深4/50m等只允许`NON_CZ_METHOD_POSITIVE_CONTROL`，不得把其复现得出的深度P95和主CZ环境混为一个样本或平均。若某论文条件和主CZ在传播机制不同，须明示条件迁移跨度。

## E1失败时仍能做什么

没有CZ2，则主CZ Yang/MMAR性能实验`NOT_EVALUATED_PHYSICAL_ENV_GATE`。可完成以下独立有价值成果：区域SSP来源及缺口清单、CZ探测必要条件、环境角点导致的焦区迁移诊断、Yang浅海复现合同、Fourier有限孔径必要条件、传感器资源依赖矩阵。不允许在另一个随意合成的“好看CZ”上偷偷替代。
