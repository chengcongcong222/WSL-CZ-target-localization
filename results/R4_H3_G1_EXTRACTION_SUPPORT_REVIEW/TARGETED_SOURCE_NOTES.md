# 本轮公式依据与引用范围

本轮是针对性接收统计审查，不重开无边界文献覆盖或声学传播实验。

| 来源 | 本轮使用 | 不用于证明 |
|---|---|---|
| [SciPy csd官方文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.csd.html) | 窗平均CSD及conj(X)Y参数顺序 | 实际目标谱、环境鲁棒性、H3精度 |
| [SciPy频谱分析官方教程](https://docs.scipy.org/doc/scipy/tutorial/signal.html#spectral-analysis) | 窗函数、频谱分辨率及能量归一化的定义 | 独立样本数、运动场平稳性 |
| [Paz-Linares等2023原始研究Eq.11](https://pmc.ncbi.nlm.nih.gov/articles/PMC10050575/) | 在该研究明确Gaussian频谱模型下的复Wishart似然 | 海洋传播、确定性线谱或低SNR声学恢复 |
| 既存H3 PRIMARY_SOURCE_PROVENANCE.csv及原文 | 模态源深因子和原方法未知量假设的既有审计 | 本CZ的组增益误差范围或组数选择 |
| G0锁定设计、结果、噪声和原始有限分数 | 已接受的局部机制及本轮读表汇总 | 新的信号提取、连续支持或环境试验 |

推导与来源区分：Rhat元素协方差由proper Gaussian四阶配对自行推导；完整联合CSD的充分性由原样本Gaussian密度自行推导；差分时钟预算由2pi f delta tau自行换算；Hann 1 Hz频点相关由固定2 s周期窗的w²离散傅里叶系数自行推导。这些不是论文声学性能复现。没有将随机源协方差模型与G0未知确定性源均值模型混同。

访问日期2026-10-08。网络引用采用主研究或官方文档；不新增UUV硬件能力事实，不宣称原文支持实际源占有、实装校准或本海区环境不确定度。格式、计算规则和历史结果另由本地哈希及DOCUMENTATION_VALIDATION.json复核。
