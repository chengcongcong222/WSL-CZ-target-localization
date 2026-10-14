# T2-N3-PREC1 / Codex研究执行合同

**唯一任务：审计已保存的KRAKEN模态二进制的单精度波数舍入，在理想等声速波导中定量归因N3/P04四阵型的0.2%空间形状硬门失败。** 本次批准一次只读数据试验和独立冷复核。严格禁止任何新增声学求解、定位性能MC、网格/窗口微调或历史结果覆盖。

## 0. 工作区与提交检查

仓库 `chengcongcong222/WSL-CZ-target-localization`。原研究父提交：

`8ee2bd53fe1b7a50725833e6c7b5c5dea35fbe23`

开始前：确保本机仓库工作树干净，`git fetch origin main` 后远端 `main` 与上述父提交完全相同；如果已改变先停止并反馈，不得在不明基线上运行。保持项目过去所有存档目录只读。本轮输出新目录 `results/R4_T2_N3_PREC1_MOD_PRECISION_ATTRIBUTION/`。本任务包另含JSON冻结规格、两份Python程序、源研究简报、当前预分析P04表及纯解析C2控制表。将程序放在仓库根目录，保留它们的原始字节SHA并Commit A冻结，**确认Commit A推送后方准运行一次真实归因**。

任何冻结输入与预登记SHA不一致，立即STOP，不寻找相近文件替代。必须读取实际Git仓库中的 `.mod` 二进制，不得用本任务包的预分析CSV替代原始数据。N3设计中的N3 Gate依旧 FAIL；0.2%为历史阈值、不能在此重定义为新的准入阈值。

## 1. 冻结源

必读并字节校验：

- `results/R4_T2_N3_BOUNDARY_PHYSICS_ANALYTIC_REFERENCE/modal_archive/n3_PR_R_f201_h1800_n160001/n3_PR_R_f201_h1800_n160001.mod` — SHA256 `cf5745862a3b367283c9ea9deef723ff194ca354c7f389fa5b185b992f028d0f`
- `results/R4_T2_N3_BOUNDARY_PHYSICS_ANALYTIC_REFERENCE/modal_archive/n3_PR_R_f201_h2400_n160001/n3_PR_R_f201_h2400_n160001.mod` — SHA256 `d9e476d0ed84c08c30f7551ad090e3823067d36918ed41375eeb3ec720afae42`
- `results/R4_T2_N3_BOUNDARY_PHYSICS_ANALYTIC_REFERENCE/ANALYTIC_KRAKEN_MODE_AND_FIELD_CROSSCHECK.csv` — SHA256 `2de1541cdfcc93c5f3dc76310fb3c6b8b47bf2e382518d8534d5ede2f56860bd`
- `results/R4_T2_N3_BOUNDARY_PHYSICS_ANALYTIC_REFERENCE/package/T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv` — SHA256 `bc2ccaebfa52893ba7bc7f217791ec980c47cb37d743abcd20cae055cf732ffd`
- `results/R4_T2_N3_BOUNDARY_PHYSICS_ANALYTIC_REFERENCE/package/t2_n3_analytic_reference.py` — SHA256 `456468df64afb04dc25ebb28cd96c62bbf2f973e03c9bc92df4fcf71d6902467`（只读理论溯源、绝非新数据模板）

首四份在` t2_n3_prec1.py`程序内硬校验。第五份在运行前由Codex单独验SHA。

## 2. 统计与科学问题

同一201 Hz等声速（1500 m/s）、水深5000 m、自由海面、理想刚性底。只用两份已执行且成功保存的MOD：CHIGH=1800与2400，预期模态数分别741、1046；不补跑原N3剩余的14组配置。固定探针P01、P02、P03、P04、P05；固定资源A0_0、A_P4M、B1_V4、B2_DUP，两个窗口共40个独立的“探针×资源×窗口数值对照”（**不是**40个物理独立目标实验）。

对每个资源-探针按源深和所有阵元收深计算空间复声场：

- C0：MOD原单精度水平波数＋MOD原单精度复模态函数。
- C1：独立解析双精度水平波数＋MOD原复模态函数。
- C2：MOD原水平波数＋解析双精度模态函数（源与接收同时替换）。
- C3：解析双精度水平波数＋解析双精度模态函数。

使用解析波数 \(\mu_n=(n+1/2)\pi/5000\)、\(k_{r,n}=\sqrt{(2\pi f/1500)^2-\mu_n^2}\)、\(\phi_n(z)=\sqrt{2/5000}\sin(\mu_n z)\)。模态号不可借助位置真值事后重排，须验证其与MOD存储顺序一致。四方案一致采用原N3 MOD远场径向因子，计算共同复比例消除后的空间形状差异 \(e=\|u-v(v^Hu)\|_2\)，保留源／收深完整阵列，不允许剔除弱信号和P04。N3的解析Hankel核仅用于原始记录差异容限对照，不和本轮远场核悄然混用。

**主判据：**原N3/P04四项C0全部>0.002，若换解析波数C1后四项都≤0.002且波数1787/1787等于解析float32舍入，则在此理想参考模型内支持“波数存储舍入构成主要可复现超限机制”。C2同时报告，**不能**把改善量除以总误差声称独立贡献百分比。主判若失败也必须交付原始错误和阴性结果，不得自动改程序继续做另一次科学执行。

完整40行和单因素影响必须报告：原空间误差与存档N3每行差绝对值≤2e-6、各MOD波数等于解析double在float32下的表示、其余冷审计逐行差≤5e-10；任何违反即不满足本轮归因合同而STOP。

## 3. 执行命令与冷审计

在仓库根目录（同Python环境需要numpy、scipy）运行：

```powershell
python -B t2_n3_prec1.py --selftest
python -B t2_n3_prec1.py --repo . --out results/R4_T2_N3_PREC1_MOD_PRECISION_ATTRIBUTION
python -B t2_n3_prec1_cold_audit.py --repo . --out results/R4_T2_N3_PREC1_MOD_PRECISION_ATTRIBUTION
```

主程序与冷审计程序必须作为不同脚本运行；冷审计不导入主程序，而用`struct`直接二进制解析、`brentq`独立解根、逐阵元`cmath`求和、QR投影算形状。需保存两者的原始输出CSV与独立审核JSON。`test_synthetic_only.py`只是事前程序健壮性检查：输入是人工生成的MOD，**禁止把其1991项PASS称为项目科学审核通过**。真实原MOD只允许一次主执行＋一次冷审计；审计失败STOP，不优化、重跑或重设Gate。

提交及交付：

- `PREC1_ALL_40_CASES.csv`
- `PREC1_P04_FROZEN_FOUR_RESOURCE.csv`
- `PREC1_MODE_WAVENUMBER_STORAGE.csv`
- `PREC1_GATE_ATTRIBUTION.csv`
- `PREC1_COLD_AUDIT_CHECKS.csv`、`PREC1_COLD_AUDIT.json`
- `PREC1_DECISION.json`
- `GPT_SYNC_PREC1.md`（来源SHA、条件结论与反例、溯源、数值不确定性、独立审计和下一研究决策）
- 一幅P04四方案比较图（图注明确理想等声速基准、仅数值形状、非定位精度）
- 源代码SHA与全输出清单、运行日志和Git Commit A/B SHA

高SNR接收、MC、所有定位表现、速度/航向、E-STD源深导数、海底选型、CHIGH窗口修改均不得出现在“已评价”栏。

## 4. 停止与解释

**成功：**本轮只可报告 `PREC1_K_STORAGE_ROUNDING_MECHANISM_SUPPORTED_CONDITIONALLY`，不改变上轮N3 `T2_N3_PARTIAL_NUMERIC_NOT_ADMITTED`、N2 `NUMERIC_BRACKETING_NOT_ESTABLISHED`、原G0 `FORWARD_NUMERIC_NOT_ADMITTED`、A/B `NOT_EVALUATED`、R4=0%。它不意味着在E-STD波导能得到解析双精度\(k_r\)。

**失败：** `PREC1_ATTRIBUTION_INCOMPLETE`，保留证据、及时停止，交给负责人重新制定下一科学问题。不允许连续追加网格、模态窗口或换物理边界来追求结果。

预算：新增KRAKEN/SCOOTER/FIELD/FIELDS/任何传播求解=0；新增MC、录音波形、DFT及后验候选反演=0；最大墙钟600秒、内存上限2GiB、结果新字节上限1GiB。Commit B提交并推送，核对`remote/main==Commit B`，然后STOP。不得以“性能提升”措辞直接面向甲方发布结论。
