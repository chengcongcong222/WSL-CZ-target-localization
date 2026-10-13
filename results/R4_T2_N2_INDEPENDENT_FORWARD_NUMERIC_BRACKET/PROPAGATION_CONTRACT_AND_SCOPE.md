# N2传播合同与执行边界

模型保持原G0的水层SSP、CVWT／Thorp选项、密度、真空海面、刚性5000m海底和频率。K的源／收深度列表保持原G0全部深度不变，包含四资源所需深度，额外历史输出深度不参与新增资源比较。FFP仅使用源深180/190/200/210/220m和四资源精确收深度并集。ENV白名单逐行验证，不按结果改底质、SSP、窗口或网格。

官方KRAKEN手册解释相速度窗对纳入模态的限制：https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html 。SCOOTER的独立谱积分方法及相应ENV输入见：https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node64.html 。FIELDS从GRN生成SHD的关系见：https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node65.html 。现代FIELDS源码及SourceReceiverPositions读取规则已存D盘source_review，按SHA冻结；其RPO／NR／范围端点slash格式与旧手册RDB示例不同。没有替换或编译仓库二进制，也不声称官方当前源码与2018二进制是同一次构建。

FIELDS固定正波数柱面积分和内部多项式谱插值，时间相位exp(−ikR)；接收距离线性插值只使用实际SHD float32坐标，实虚分量分别处理，禁止外推。所有四资源复用同一哨兵的同一FIELDS运行，生成全部5源深／11收深度的输出，没有为某阵型另选窗口。RMAX65km和120秒单FIELDS超时均事前冻结，FFT内部扩充NR记录而不改参数。

KRAKEN按三配置×三频执行一次，任何求解失败停止剩余K调用；已完整配置仍进入只读诊断，独立FFP路径若可用继续。FFP路径首次非0返回、GRN/SHD接口错误或坐标不合法即停止余下FFP，不重跑、不采用FIELD替代，未执行清单全部保留NOT_EVALUATED。

完整场包括全部返回模态。cp>1800或2000的分带仅为贡献诊断，不从主求和删模态。模态形状辅助匹配若出现零存储向量，记未评价和零计数，不赋零形状误差，也不删该模态。复杂度、弱场、量化和源谱系数消除都不覆盖模态窗外真实性。

单数值进程、BLAS1线程，科学总时钟6h、8GiB进程树RSS、16GiB新D盘、512MiB新Git；K/S/F单次超时180/900/120s。所有新求解在A之后，固定源哨兵不是实际观测。本轮无波形、DFT、源噪声抽样、MC、优化、FIM或N1模板重评分。原G0、N1、H3结果不变，A/B性能NOT_EVALUATED，CHIGH全域UNCERTIFIED，R4=0%，B后STOP。