# GPT 同步：E1 主文与公式锁

父版本0db30f60e781676f9168d380678693a32f1325a1的研究重置设计已独立接受，文献覆盖完整性未一并接受。此次按授权只读取主文、推导公式和修订协议；新MC、新信号、新传播、项目FIM/性能运行全部为0。固定无量纲代数11项PASS/0FAIL不是新科学实验。

正式判定：**E1_PRIMARY_SOURCE_OR_MODEL_LOCK_INCOMPLETE**，E1 execution admission NOT_READY。四篇指定主文完整PDF均取得、关键方法已读；没有把“全文下载成功”当成“严格算法基线已完全锁定”。

- Sun2024 Eq.(4)原始 β+多频率量测已锁，未知每线中心频率与随机游走已确认。完整MFB-AUKF的状态转移、UKF权重、全状态Jacobian与仿真初值存在印刷冲突；原式及推导分列，不静默改成标准UKF。严格LIT-F为PARTIAL_METHOD_ONLY。
- Xu2018确有从波束复声压到互相关谱、速度幅值、条件符号再进入TMA的真实链。链的输入/输出已LOCKED，但浅海相位/平均水平相速度/波导符号不能直接转到第一CZ。其理想噪声敏感性仿真与真实提取结果分开。
- Li2019是同bearing+radial velocity观测上的SIR/VSIR替代估计器，不是新的信息源。冷启动支持域/原始提取细节不能由收敛结果推断。
- Kita使用已知PWM频率与极数，RPM→speed需要车辆映射；原文还用导航真位置选择栅瓣，当前标NOT_TRANSFERABLE_CONTROL。不能隐藏为未知目标的纯声学航速输入。

[coverage](E1_PRIMARY_SOURCE_COVERAGE.csv)逐列记录正文条件、read_scope与formula_lock_status；[模型锁](E1_ACOUSTIC_MEASUREMENT_MODEL_LOCK.md)统一M0–M3；[结构分析](E1_NUISANCE_IDENTIFIABILITY.md)给A–F条件；[比较规则](E1_LITERATURE_BASELINE_COMPARISON.md)与[修订协议](E1_REVISED_MINIMAL_EXPERIMENT_DESIGN.md)冻结B0/LIT-F/LIT-RV/UUV-CONTROL/NEW、20%研发门槛和双层输入合同；[准入JSON](E1_EXECUTION_ADMISSION.json)给逐项阻塞理由。

科学结论限结构：未知f0与恒定径向控制有尺度退化；独立节点仿射reference能吸收仿射Doppler；多线共享fractional漂移每节点/时刻仅一个运动loading；双节点同发射时刻的共同任意源漂移仍可能保留差分，不能一概声称全吸收。校准prior收益须另记外部信息。

NEW未知f0机制已有文献；候选差异在源/节点漂移与双节点条件处理。当前 NEW_EXTENSION_NOT_ESTABLISHED；没有有效信息增益计算，也没有创新确认。20%路线门槛不改历史10% Gate。

停止原因不是缺PDF，也不是硬件UNKNOWN；是严格原文基线存在未决模型/估计器定义，且层2实现/不确定度尚未验证。本次交付已经把未决项定位为具体公式，不以自定义修正版冒充通过。若以后另行决定采用修正式，必须明确标 DERIVED_CORRECTION / PAPER_ADAPTED，先冻结新的基线合同，不自动开始科学计算。

E2未开放；FDSL/POMAP原 comparator pending；历史搜索支线及深度/A2/SSP/P5不开放。原联合A1/R4=0%。提交一个文档/公式锁版本，push后校验HEAD==remote/main，然后STOP等待独立审计。
