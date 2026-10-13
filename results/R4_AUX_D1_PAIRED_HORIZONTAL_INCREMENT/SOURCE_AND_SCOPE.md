# AUX-D1来源与作用域

ABSOLUTE_BEARING_SIGN_CONDITION_ASSUMED。本轮比较U0单节点方位与U5双节点方位，同主节点输入逐字节相同。不是原HLA全部声学方法对比，不是控制观测数量后的纯基线因果。

源文件按9060ca4固定读取，Git blob与工作区原始字节分开记录。旧原始方位／报告轨迹缺失不重搜、不从真值声称恢复。新数据是新独立仿真，12状态为已见开发几何，512校准状态为一次独立Uniform抽样。原三档、R3/H1/H2/RX/SP1和原PENDING审计字段不改。

SciPy least_squares官方文档说明其为局部稳健最小二乘，固定损失并不消除整段偏差：https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html 。本机冻结SciPy1.15.3，显式选项见完整任务与代码，不临时升级。Sobol完整128点再取126点，不宣称截取后仍有完整平衡性质：https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.random_base2.html 。

估计器以报告MAIN_t0为参照，评价真值以真实MAIN_t0为参照，保持导航参照影响。未知深度不在此水平模型内，不称已知或已估计。同步、关联、阵姿态与被动有符号方位为条件，不是硬件实测。没有传播、径向、源频相位、垂向、FIM或真实接收计算。
