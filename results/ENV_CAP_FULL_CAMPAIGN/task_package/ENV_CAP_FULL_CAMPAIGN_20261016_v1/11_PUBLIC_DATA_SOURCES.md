# 外部公开数据与方法来源（仅作为待审来源，不代表已下载或验证）

## 数据产品

1. **NOAA World Ocean Atlas 2023（WOA23）**。官方网站：`https://www.ncei.noaa.gov/products/world-ocean-atlas`；数据入口：`https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/`。包括不同气候期的温盐度客观分析场和部分误差统计。需要记录版本、气候时段、纬经度、深度标准格点和缺失深度；其地理格点均值**不能当具体测站的CTD实测**。
2. **NOAA ETOPO Global Relief Model 2022**：`https://www.ncei.noaa.gov/products/etopo-global-relief-model`。用于区域底深与地形；需记录栅格、分辨率、海底源ID和沿航线抽样；不能当海底声速/密度/衰减信息。
3. **Acoustics Toolbox/BELLHOP 物理/数值手册**：`https://www.oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node61.html`。BELLHOP为声线/波束追踪及场计算器，输入的深度、距离单位与选项必须逐例审计。它并不自动提供包含底质或模态精度完整性证明的“真值”。
4. **Yang, 2015**：*Source depth estimation based on synthetic aperture beamforming for a moving source*, JASA 138(3),1678–1686, DOI `10.1121/1.4929748`；**Erratum 2018**, JASA 144(6),3075, DOI `10.1121/1.5081712`（`https://pubmed.ncbi.nlm.nih.gov/30599696/`）。使用项目已归档全文为公式来源，勘误offset `delta`不可忽视。
5. **Liang et al., 2018 MMAR**：*Match-Mode Autoregressive Method for Moving Source Depth Estimation in Shallow Water Waveguides*, Mathematical Problems in Engineering, Article 7824671。先从原文确认具体AR、峰值筛选、Hankel幅度与深度函数，不允许以历史代码代替原论文。

## 项目已有候选环境（不是自动准入）

- `m452_build_and_cz.py`, `m453_pipeline.py`, `m454_closure.py`：过去WOA23源网格17N118E和ETOPO路线，部分代码默认WSL `/home/...` 目录，需要重新发现数据路径/合法镜像/状态；README报告中的`FIRST_CZ`输出未在本包独立核验。
- `B3D3_M4_pack_extracted/README_B3D3_ENVBUILD_PREFLIGHT.md`：18N114E/April基准及深部重建，用自写Snell积分给“first-deep-return-like”信息，不等于经成熟传播模型确认的第一会聚区。过去其路径族8环境关联失败和无完整实测SSP应保留。
- E-STD约5km水深Munk-like名义算例：旧KRAKEN rigid底模态截断及输出精度未收敛。可以做理论对照，但如无新增可信传播准入，不进入新精度Demo。

## 数据授权与禁区

Codex可按公开许可合法获取小区域子集，不可绕过访问限制；若官网网络迁移/临时不可用，尝试有来源且可核实的正式镜像，记录全部URL及响应、哈希，而不是私下伪造数据。必须保证数据源的历史身份与本轮源码冻结一致。新下载内容身份核查成功后才进入科学计算；下载失败记`SOURCE_UNAVAILABLE`。

## 研究用语

`MEASURED_IN_SITU`, `CLIMATOLOGICAL_OBJECTIVE_ANALYSIS`, `RECONSTRUCTED_WITH_PUBLISHED_CONSTRAINTS`, `IDEAL_ANALYTIC`, `GEOMETRY_DIAGNOSTIC` 是不同的环境来源标签，不可混用。“具有第一会聚区”需实际物理准入，与这些来源类型彼此独立。
