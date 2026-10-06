# E1 修订最小实验合同：设计锁，未准入执行

本文件替代父版本 TOP2 文档中的 E1 基线/频率参考定义，父文件作为历史设计保留。E2 不修改、不执行。所有以下参数是申请阶段预研设计，不是设备验收或硬件性能承诺。硬件 UNKNOWN 非阻塞理由；本轮阻塞来自严格原文基线模型未闭合。

## 范围、输入和公平性

当前1200s、原运动/观测噪声、H01–H12与5km Anchor A、双镜像共24确定性场景保持；不生成新truth，不延长时长，不更改机动，不新增MC。本轮连这些场景的Jacobian也不计算，仅固定结构toy代数。未来科学运行需新的单阶段授权与可审计pre-run代码冻结；本文不构成当前运行许可。

基线为 B0、LIT-F、LIT-RV、UUV-CONTROL、NEW；具体继承/条件见 [comparison](E1_LITERATURE_BASELINE_COMPARISON.md)。Li同observable PF为估计器控制，Yang已知CW为物理控制，RF只为机制参考。严格 LIT-F 因Sun印刷冲突未完成；不可使用自行修正的UKF并仍叫原文复现。

## 层1：MECHANISM_UPPER_BOUND

M0–M3 联合均值Jacobian、Σ及源/参考/环境 nuisance投影是输入；目标函数不靠真值初始化性能评分。使用理想frequency观测可作机制模型输入，明确σ_f=0.001/0.005/0.010Hz是登记假设，而非已提取精度。所有方案使用相同60个20s非重叠频率帧中心10,30,…,1190s；既有bearing采样与频率帧联合时间表固定，不虚构121个独立frequency帧。

NO_CALIBRATION 主结果：未知各线源常数、共享分数线性漂移；未知每节点频率offset/drift必须区分 Hz加性与fractional参考。三线模型 motion rank不按3倍计。完整观察子空间白化并投影；零方向、不可估speed方差必须记UNBOUNDED。WITH_CALIBRATION作为另列条件，来源、精度、成本必填，不能替换主结果。

主fractional-reference设计模型：ρ_j无量纲、κ_j s⁻¹；源共享d_c s⁻¹。源d_c取0、±10⁻⁷/s；参考ρ幅度5×10⁻⁵、κ幅度5×10⁻⁸/s，仅约在200Hz处分别对应0.01Hz、10⁻⁵Hz/s，其他线的Hz误差按频率缩放。父版本±0.01Hz/±10⁻⁵Hz/s仅作为 ADDITIVE_HZ_CONTROL 的B/K，不能移植成主分数漂移已知量。时间偏差δt取0、±0.1s。上述生成范围不是给估计器真实参数，也不是给零误差强校准prior。

参数域/维数合同：主 M3-SHARED-FRACTIONAL 中每线一个 logF_l、一个共享d_c、两节点各ρ/κ/δt；B_j=K_j=0、d_l=0。ADDITIVE_HZ_CONTROL 独立运行B/K模型且ρ/κ=0；BOTH_REFERENCE_TYPES 只作为更宽 nuisance 结构控制，不把两套完全耦合的未知量假装已可估。未知参数不根据生成truth置零；冗余 nuisance 列由rank-revealing投影处理，或使用已登记纯坐标约定，绝不伪装为校准观测。层1有效传播控制取c=1500m/s、τ_j=r_j/c，由候选状态计算；这是与原文同类的已知传播条件，不是未知SSP鲁棒结果。若同时profile自由c或模态未知量，必须另列结构条件，未冻结其原始提取forward实现前不宣称CZ层2可执行。

冻结5个确定性角点，不枚举连续鲁棒保证：
C0所有新增源/参考/时间偏差0；
C1 d_c=+10⁻⁷，节点1/2均(+ρ,+κ,+δt)；
C2为C1全部反号；
C3 d_c=+10⁻⁷，节点1(+ρ,+κ,+δt)，节点2(−ρ,−κ,−δt)；
C4为C3全部反号。
这些是新协议里的有限符号控制，尚未生成信号或计算FIM。主共享漂移与线特有线性漂移要分别标模型，不将共享漂移当独立三线随机噪声。后者先作结构控制，未登记新增性能矩阵前不自动扩大主实验。

必须保留零信息控制：恒定径向/未知源；每节点每时刻自由频率参考；多线共享fractional loading；单节点各线任意源函数；双节点共同源任意漂移与逐节点任意漂移分别测试吸收条件。实现若不遵守已推导结构→IMPLEMENTATION_INVALID，停止科学判定。

## 层2：EXTRACTED_OBSERVABLE 输入合同

接收信号→频谱/波束域互相关→f或q估计→联合covariance/失败mask→TMA。原始录音尚未生成，整个管线尚未运行。

生成端隔离：仅生成器可读truth状态、true F_l、路径/模态身份、时延、真实q和source phase；评分器可读truth作误差评估。估计端只能读接收数据、节点名义导航/时间戳、采样率/单位、登记宽band、误差模型及显式独立校准。接收谱线初值来自观测，不提供true F_l。不得给真实谱线身份、true速度符号、发射时间对齐、模态编号或“正确峰”索引。每一数据键的来源须在后续pre-run代码检查。

原设计两节点相位连续波束输出，fs1024Hz、1200s，单线200Hz或三线170/200/230Hz，保持为生成端条件；给估计端只声明150–250Hz搜索带，不能传三根真值中心/编号。盲峰检测、跨帧轨迹与跨节点同源关联规则须在运行前独立冻结。缺线/换线/误关联保留EXTRACTION_FAILURE，不凭truth纠正。

源模型主共享fractional线性漂移对应相位 φ_l(u)=2π F_l[u+d_c u²/2]的一阶慢漂移形式；若用指数频率模型M3，精确相位是2π F_l(expm1(d_c u)/d_c)，d_c=0取极限。生成/估计必须使用同一登记阶次，不能把二者悄悄混用以制造bias。幅度控制为常数及1+0.2sin(2πu/300)，不是直接给源谱模板。每线旁+0.02Hz、幅度0.1、初相π/3的固定干扰是确定性stress，不是随机SNR/P95验证。

接收参考频率的主模型是时钟重参数化/对应fractional phase evolution，加性LO由独立控制模拟；频率误差参数不能只在输出f上加噪却声称真实ADC录音已验证。时间偏差、参考误差相关性和延迟对齐须在pre-run代码里实现并逐字段校验。当前没有这种实现证据。

频率提取合同：20s非重叠Hann窗、FFT和三点插值候选；峰关联不能使用true中心；报告原始频谱、所有竞争峰、选峰规则、缺失mask、时间支撑、f及联合Σ。该处理是 PROJECT_EXTRACTOR，不叫Sun频率提取复现。Yang PLL只作为另列known-CW控制，需先锁采样/增益单位，不靠truth调参。

LIT-RV合同：复声压/波束输出→γ(f,Δt)→窄带时间或宽带频率FFT→|q|→波导符号条件/±分支→β,q,Σ→原EKF；必须报告v_p、β_w、相位补偿来源和实际孔径。论文窄带延迟至150s，不能用20s帧悄悄替代并称严格复现；必须另列原条件控制或 NOT_APPLICABLE。当前CZ没有已验证的相位补偿/β_w符号，不能输入truth符号。若无法复现，明确 EXTRACTOR_NOT_ESTABLISHED，不改用理想vr称NEW/LIT-RV通过。

covariance合同：同录音的bearing/f、差分bearing、各线、跨时重叠/平滑、两节点共同源波动需保留相关。σ_q理想注入不替代提取Σ；确定性干扰bias不等于高斯噪声。未来须给联合原始噪声→提取Jacobian或独立标定的数据来源；没有则只输出点估计/描述性bias，EXTRACTION_UNCERTAINTY_NOT_CALIBRATED，禁止信息/P95通过。

传播合同：先单有效分量验证符号/量纲，随后若单独获得授权才可用同一E-STD多模态forward；相位因子/传播相位不得重复计入，群时延和谱变不许用true路径辅助估计。源-收对齐必须从候选轨迹/时间模型推断，oracle对齐单列机制控制。当前只有合同，未验证该合成器；不能称原始CZ信号链已可运行。

ORACLE_vr、ORACLE_f0、ORACLE_LINE_ID、ORACLE_ALIGNMENT仅允许另列MECHANISM_UPPER_BOUND，绝不填充NEW方法缺失输入。UUV-CONTROL当前 NOT_TRANSFERABLE_CONTROL，不生成虚拟known车型收益来冲抵主任务。

## 资源、指标和停止

继承原有限搜索proposal：2000个observation-only候选、每local不超过200评估、双观测初始化镜像。实施前冻结严格计数器/单位/场景清单；这些数字不是覆盖证书。本轮没有实施或计时。两路float32约9.8MB为原设计原始体量，不包括复数/STFT/所有中间文件，也不等于已测工程可行性。

层1主门槛：相对B0有效speed方差中位降低≥20%，最坏场景有效方差不增大（max_s V_NEW,s≤max_s V_B0,s），逐例变化全量报告；同时对最近声学文献基线报告公平比较和真实新增条件。条件扩展与仅机制复现的标签见 comparison。任何不可比/奇异不可估方向阻止整体PASS，不删除场景。20%非最终10%历史Gate。

最终创新、速度P95、硬件可用、同步与关联均未验证。理想信息无增量→CURRENT_FREQUENCY_MODEL_NOT_SUPPORTED；理想有益而原始提取不成立→EXTRACTION_NOT_ESTABLISHED；仅外部prior有益→CALIBRATION_CONDITIONAL。不得加时长、改阈值、调精度到通过或换优化器救结果。

本轮准入：E1_PRIMARY_SOURCE_OR_MODEL_LOCK_INCOMPLETE。下一个动作 STOP，科学运行=0，E2=NOT_OPENED，FDSL/POMAP=ORIGINAL_FDSL_POMAP_COMPARATOR_PENDING，R4=0%。
