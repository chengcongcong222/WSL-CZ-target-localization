# T2-N3来源、物理假设与结果资格

## 1. 固定来源和Git身份

GitHub：`chengcongcong222/WSL-CZ-target-localization`；待执行父提交 `1ddc329edec9cbf07ffd6f2f5bbd63805104964e`。在代码、ENV、已存模态和成果之间保存SHA256逐字节映射，Git blob SHA与原文件SHA256分别列示。此前 `T2-N2` 只完成KRAKEN 9/9、SCOOTER 1/9、FIELDS首次120s超时；没有取得SCOOTER完整空间场。N2的K窗口全部失败、N1旧窗口B正向证书、G0未准入都不可修改。

固定参考路径：

1. `results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/environments/t2_f201_n160001.env` 及同频235/283；记录5000m SSP、表面`V`、底部`R`、`CVWT`、THORP的原始记录与各自字节身份。
2. `results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/ARRAY_GEOMETRY_TABLE.csv`；四资源 `A0_0,A_P4M,B1_V4,B2_DUP`。
3. `results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/MODE_WINDOW_AND_MATCHING.csv`, `FULL_CATALOGUE_KRAKEN_STABILITY.csv`, `WINDOW_MODAL_CONTRIBUTION.csv`, `NUMERICAL_DECISION.json`, `INDEPENDENT_COLD_AUDIT.json`。
4. `results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/T2_N2_SENTINELS_20261013_v1.csv`（12行）。
5. D盘先前存档中的K0/K1/K2 `.mod`、完整样本/贡献数组必须核对原 SHA，缺失不造假或从真值补造。
6. 仓库中的 `tools/acoustics_toolbox/atWin10/at/bin/kraken.exe` 的实际SHA，以及原环境求解器保存SHA。

外部技术资料（只用于方程／工具选项）：

- KRAKEN 原手册： https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html （顶部`V`=pressure release，底部`R`=perfectly rigid；底部`A`=acoustic halfspace，通常用于真实底质，Thorp/attenuation及CHIGH定义）
- 边界条件： https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node10.html （自由/刚性/声学半空间边界区分）
- 理想波导解析基准独立来源： https://journals.ametsoc.org/view/journals/phoc/45/12/jpo-d-14-0256.1.xml （刚性底/自由表面本征垂向波数）；或手册`The Isovelocity Problem`章节。

## 2. 不可混淆的四种水声模型

- **M0 历史E-STD：**Munk-like SSP、约5000m海深、海面自由、海底理想刚性、存档`CVWT/Thorp`；该配置下N2窗口差异巨大，并不是已认证的真实海底。
- **M1 理想等声速自由表面＋刚性海底：**5000m、c=1500m/s、ρ=1，**无吸收**。可解析的闭式模态基准，用于证明KRAKEN能否按正确边界、相速度窗与源／收位置回出模态和空间场。绝不称第一会聚区模型。
- **M2 理想等声速双压力释放边界：**参数与M1一致，唯底部 `V`。这是另一数学边界极限，不是实际吸收海底。
- **M3 历史E-STD但底部 `V`：**与M0水柱、Thorp和频点相同，唯一底部边界不同。仅用于环境假设的**数学敏感性**，不冒充真实底质/海试验证。真实底质密度、c、吸收系数未由甲方提供；不可选择能让方案好看的数值。

M1/M2与M0/M3物理合同不同；禁止以M1/M2的解析参考直接评价M0实际传播正确性。只有已在同一物理环境做并通过独立算例、数值收敛与边界合理性核验，才可以考虑恢复原T2-G0接收性能试验。该恢复属于**后续独立研究**，不在N3。

## 3. 新的物理推断与适用边界

对给定正常模态的水平相速度 `c_p=ω/Re(k_r)`，当地水柱内垂向波数平方 `k_z(z)^2=ω^2/c(z)^2−(ω/c_p)^2`。若`c_p>max_{z∈[0,H]} c(z)`，则在这个名义水柱中 `k_z(z)^2>0` 对所有深度成立，即**不存在水柱内部几何转向点**，对应底边界敏感的全深传播分支。该条件本身不说明真实海域的海底反射强度、实测吸收或真实高阶模态一定可见。

现有E-STD原SSP从1500m/s声道轴向海底上升，5000m处约1635.1811m/s；新增的1800–2000及2000–2400模态均满足上式。N3第一环节必须从原`.env`而非硬编码提取真实SSP最大声速、c(z)范围、海面/底层选项，对每一已存模态的cp作相同分类并核算声场贡献。显示“高相速度无内部转向”只能得到**机制层**结论，不能因此宣称近200m接收必然观察到底反射的强信号。

## 4. 必须区分的失败／未评估结论

旧G0最大两网格形状差0.7344%，54/54最坏值超登记0.2%；旧N1仅旧1800窗口5无噪声模板比较：A 7项严格正、2项严格负、1项未决，B 10项严格正；旧N2 1800→2000最大99.9902%、2000→2400最大98.7434%，2400双网格最大0.3712%，没有FFP空间证书。不能挑N1的正结果与N2的宽窗口拼出“已认证B”。

N3不会更新R4（仍0%）；不会生成新接收波形、DFT、MC、FIM、速度/航向、真值保留率，也不会把压力释放的深海海底称真实现场参数。即使M1/M2解析基准全通过，M0在第一CZ真实传播仍未认证。仅原单HLA为原设资源，A/B均为**未获批准的可选结构增强**，不能改甲方合同。
