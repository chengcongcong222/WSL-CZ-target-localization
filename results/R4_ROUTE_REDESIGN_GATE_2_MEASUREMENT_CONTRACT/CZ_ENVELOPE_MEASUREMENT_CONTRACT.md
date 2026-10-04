# Contract A — CZ Envelope（拟议最小被动合同）

状态 `PROPOSED_CONTRACT`；尚未被负责人接受为新增 receiver contract，G2 未全部闭合。冻结的是要求与断链，不是系统已有能力。目标是粗距离候选；不承诺 km/百米精度、单值或 RC3 capture，不设频带/采样/平滑/性能门限。

## 数据选择与完整链

`single-HLA target-bearing receive product → time-indexed tracked-line power/level → per-line source/gain nuisance treatment → observation-domain coarse envelope/edge/shape → joint range/motion/nuisance candidate intervals`。

最低草案采用 **B + C：多个可关联 narrowband components 的相对接收级随时间变化**。不要求 A 的 dense PSD；dense-band averaging 可以作为可选变体，不能成为本轮自动新增 UUV 信号假设。单时刻 A 的 dense PSD 在任意未知源谱下也不能自行解决判距。

最低交付物为：每条 component 的 received frequency 与 timestamp、pre-normalization beam power/level 或可回溯的 relative level + 归一化元数据、目标 bearing/beam 与谱线关联、通道/beam 相对频响或稳定性描述、AGC/beam 权重变化与质量/失锁/缺测标志、噪声背景/不确定性说明。无需已知 absolute source level；相对校准或可建模稳定 gain 足够作为草案要求。若保存 PSD，必须给频率 bin、窗口与标度定义，不能把幅度、功率、PSD 和 dB 均值混用。具体参数留待另行 quantitative design。

R2 消费者不必取得阵元级 waveform 或 complex source phase；上游 target-bearing beam 形成若使用阵元相位，需要真实同步、阵形及相对通道校准。可以接收经过审计的 beam spectral 产品及其元数据，不能默认该上游能力已经存在。平台导航、阵深、环境信息应与时间轴关联；接收深度不是已知目标深度。

## Signal condition 与 nuisance

只冻结 REQUIRED_SIGNAL_CONDITION，不当作事实：谱线在窗口内可关联、有可用强度且不任意跳变；每条源级/beam gain 能用常数或明确受约束的漂移模型描述；目标方向隔离与噪声影响可审计；coarse 传播结构未被短窗或归一化抹去。所有条件当前缺真实 UUV 数据支持，属于条件性物理假设，不要求目标配合发射。

令 `L_f(t)=a_f(t)+g_f(t)+P_f(R(t),z_s,z_r,e)+noise`。`a_f` 是每频未知源谱，`g_f` 是接收/beam gain。若在窗口内为常数，按**同频时间轴**去均值可消去各自常数，甚至其跨频值完全不同；这不需要已知源谱或 line ratios。

这一事实不等于消去任意源谱：时间变化的 `a_f(t)`、线谱更换、方向辐射变化及 gain drift 仍在；不能从去均值后的数据还原绝对谱级或进行无依据跨频平均。单次跨频差/ratio 只能消共同标量源级/增益，保留未知 `a_f-a_g`。若跨频 averaging，则另需 stable ratios、可建模 source trend 或按线独立归一化和权重说明；**smooth spectrum / stationary ratios / slower temporal variation** 都是可选模型条件，不是本项目事实。

“源变化比传播慢”只有可分离时才有用；若 CZ coarse envelope 同样缓慢，detrending 可能把两者一起删掉。任意每频每时刻的源级自由度能解释任何残差，不能用它得到 range anchor。窗口独立去均值也可能抹去跨窗 pedestal/跳变；必须保留归一化记录，不把模型窗均值解释成未知源已被观测。

源级/源谱/时变 gain 为 nuisance；源深、SSP/海深/地声不确定性与接收阵深误差进入候选模型的联合/边缘化变量。不能固定为真值、强行增加窄先验或在本轮发展 SSP/depth。无绝对 source phase；相位不是下游特征。

## Absolute-range 机制及失败出口

CZ envelope/edge/broad maximum 是由传播环境、收发深度及频率决定的有限距离区段结构，提供潜在**有物理距离尺度的模板**。可讨论的推理是：在已知平台导航但未知目标运动下，观测域的多线相对变化、曲率/edge 的组合，逐个与有限 CZ 条件域中 `R(t;x)` 的 coarse 区段/事件假设相容性比较，保留能共同解释时间形状及环境/深度不确定性的距离区间。只使用 coarse feature，不能把原 exact-TL residual 换名。

绝对尺度来自条件传播区段的位置/宽度/频率迁移结构，不来自 unknown source level，也不把观测时间直接当 range coordinate。不同 motion、源深/阵深、SSP 可以给同形或移动区段；edge 两侧、多个 broad peaks 或 CZ order 可多值，必须保留多候选。只可在看到足够 shape/事件且 nuisance 后仍有独立差别时，希望区分第一 CZ 内粗区段；不保证总能收缩 45–60 km。若只识别“第一 CZ 内”，输出就是 RC1 band，没有新锚。

当前尚无证据证明短窗能看见可辨 edge/maximum，或源/环境范围内的 centered envelope 仍有 coarse range 信息。M0 已足以表达某些 temporal proxy，这是未来可研究的子方案；它没有自动闭合真实 receiver 与 signal 合同。

## G2 与未来证伪（不执行）

G2-1 未闭合：原始产品现状未知，草案已定义但负责人未接受扩展。G2-2 被动结构满足；G2-3 物理上条件合理而 UUV 未验证；G2-4 nuisance 角色明确；G2-5 有条件绝对区段机制、没有数值成功声明；G2-6 可证伪。这些 true 仅指合同结构推理。

未来另行授权可设计比较：处理前后可用信息、无 edge/plateau、source drift/beam gain、源深/环境导致同形、motion-scale compensation。若 nuisance 后 feature 与 range 无独立关系、所有 coarse support 仍等于 RC1、平滑删除事件或无法区分源变与传播，则该合同的 acquisition 假设可 FAIL。不给新阈值/具体试验 panel，不运行模拟，也不承诺第三 optimizer 能救回。
