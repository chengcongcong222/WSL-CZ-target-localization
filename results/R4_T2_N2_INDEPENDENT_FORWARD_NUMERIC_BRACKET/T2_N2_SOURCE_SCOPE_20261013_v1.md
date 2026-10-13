# T2-N2证据来源、模型边界与工程接口

## 读取顺序（全部源只读、锁定输入字节）

GitHub仓库 `chengcongcong222/WSL-CZ-target-localization`，基线commit `eeb06809b1b37852ef41d5d03b6348dcf95310a8`。开始后立即核对remote/main、本地HEAD、dirty状态与大文件存档完整性。

1. G0：`results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/{GPT_SYNC.md,DESIGN_FREEZE.json,PROVIDER_PROVENANCE.json,PROVIDER_MODE_WINDOW_AND_MESH.csv,FORWARD_DECISION.json,FORWARD_CATALOGUE_NUMERICAL_STABILITY.csv,ARRAY_GEOMETRY_TABLE.csv,CATALOGUE_primary.npy,CATALOGUE_fine.npy,LOCAL_ARCHIVE_INDEX.csv}`。
2. 原完整模态：`D:/ProjectStorage/WSL-CZ/T2_TILT_G0/modal/t2_f{201,235,283}_n160001.mod`，**根据G0记录的原始SHA256逐个核对**；旧80001模态与完整目录的原始`.npy`复场同样保持只读。
3. N1：`results/R4_T2_N1_ARCHIVED_FIELD_SCORE_STABILITY/{GPT_SYNC.md,NEXT_RESEARCH_DECISION.md,TOP_ERROR_CANDIDATES.csv,PAIR_ADVANTAGE_NUMERIC_INTERVALS.csv,ANALYSIS_FREEZE.json}`，以及D盘旧复场摘要。原B10/10是确定性无噪模板的**两网格数值符号证书**，不是接收性能；A为7项严格正、2项严格负、1项正而未获严格证书。
4. H3-P1：`results/R4_H3_P1_UNIFIED_WINDOW_SSP/GPT_SYNC.md`、`FAILURE_SCOPE.md`。旧CHIGH和SSP割线失败的作用域均不被重判。

## 外部方法文献（仅用于选择独立数值方法）

- KRAKEN官方手册：`https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html`。CHIGH是上相速度界，增大它通常增加计算模态；局部窗内稳定与物理完整性不同；RMAX影响数值求解设置。软件可能限制可计算的泄漏模态。
- SCOOTER官方手册：`https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node64.html`。SCOOTER是基于谱积分／快速场方法的独立声场求解器，使用与KRAKEN相对应的ENV结构，但算法及适用条件不同。
- FIELDS官方手册：`https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node65.html`。`fields.exe`从SCOOTER的`*.GRN`得到`*.SHD`；绝不能误用`field.exe`对`*.MOD`的输出冒充独立SCOOTER。

在repo目录`tools/acoustics_toolbox/atWin10/at/bin/`的当前树中发现`kraken.exe`、`scooter.exe`、`fields.exe`。但二进制**实际能否运行、接口与数值结果是否有效，必须由本机验证**；GitHub存在文件≠已验证运行成功。

## 不可逾越的物理和统计边界

- 本次T2-N2全部为**冻结名义环境下的静态前向数值验证**，并非受声源未知、有限快拍、噪声和频间增益影响的接收观测；没有原单HLA或T2-A/B的动态四五参精度数据。
- 原始水深5000m、海面真空、海底理想刚性边界、名义Munk类SSP、Thorp选项均为仿真条件。SCOOTER与KRAKEN即使一致，也不验证实际海洋底质、海流、阵形畸变和环境失配。
- 两种KRAKEN网格、两个CHIGH截断窗可能同时遗漏模态。SCOOTER的谱积分是**不同数值求解路线**，但其谱窗、有限元与FFT/坐标采样也需要自身收敛；二者相符不意味着绝对真实物理场认证。
- 对每个阵型相同的物理几何位置、同样的水平距离算法与深度取样；不得修改源参数来使某阵通过，不得根据候选匹配结果选取窗口，且CHIGH=2400/FFP=3000以上未检验。
- 先前G0的`FORWARD_NUMERIC_NOT_ADMITTED`、N1的两网格证据、所有原H3-P1失败与R4=0%原样保留。任何数值诊断不能自动解锁G0的正式校准、接收波形、定位MC或工程阵型选择。
