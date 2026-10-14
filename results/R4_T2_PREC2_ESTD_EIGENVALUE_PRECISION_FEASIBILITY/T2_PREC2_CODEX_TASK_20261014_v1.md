# Codex执行合同：T2-PREC2 非均匀SSP水平波数精度可行性验证（P0/P1/P2）

**设计状态：待本地冻结与执行。此前会话的连接器稀疏PRT/MOD核查属于预研，不代替本地正式实验。**

## 0. 任务与父提交

- 仓库：`chengcongcong222/WSL-CZ-target-localization`。
- 父提交（精确SHA）：`ac43a8702282fc6c80663319cb6757a30b084e09`。
- 新结果目录：`results/R4_T2_PREC2_ESTD_EIGENVALUE_PRECISION_FEASIBILITY/`。
- 本地D盘原始缓存建议：`D:/ProjectStorage/WSL-CZ/T2_PREC2`，不要改旧N2/N3/PREC1目录。
- 本轮是新项目实验，不是对历史G0/N3 Gate的恢复、复跑或修改。R4=0%。

## 1. 数学问题、允许的证据级别

当前E-STD：水深5000m；应按ENV逐字节读取251个SSP节点，不自行生成曲线；水体密度1；V海面、R海底；`CVWT` Thorp；频率201/235/283Hz；相速度CLOW1500、CHIGH2400m/s，原K2/K3网格160001/240001。

本轮仅针对 *无损自伴、恒密度* 实特征值算子：

\[
\phi''(z)+\left[\frac{\omega^2}{c(z)^2}-k_r^2\right]\phi(z)=0,\quad \phi(0)=0,\quad\phi'(5000)=0.
\]

这与有吸收的完整复波数不是同一研究对象。KRAKEN将衰减作为模态波数的虚部（在该版本中扰动处理），不得把数值估计的实部替换进旧MOD并称作已正确合成的有损复声场。输出含`C-linear`节点间插值的定义、各源SHA、频率与网格、模态序号、实部和原虚部的独立作用域。

\[
|\delta k_r|\le\delta\phi/60000 = (10^{-4})/60000\approx1.67\times10^{-9}\ \mathrm{m^{-1}}.
\]

该值是本轮人工冻结的**单模态相位预算**，不是2e-3空间形状门槛的解析充分条件，不能直接把模态实部成功转化为阵型性能或全复场精度。

## 2. 输入冻结与不可变历史源

从原父提交读取六套K2/K3：

- K2，201/235/283Hz，MOD模态980/1146/1380，160001网格。
- K3，201/235/283Hz，MOD模态980/1146/1380，240001网格。

每一套读取`environments/{stem}.env`、`modal_archive/{stem}.mod`和同名`.prt`。**全部18份文件的SHA256已写入机器冻结配置** `T2_PREC2_EXPERIMENT_SPEC_20261014_v1.json`，来源为旧N2 `SCIENCE_OUTPUT_FREEZE.json`和`SOLVER_PLAN.json`。

先运行包内`t2_prec2_preflight.py`，确认这些文件在**父提交Git blob中**存在且与工作区字节一致。不得借用本地重生成ENV、PRT或MOD作为身份替代，也不得使用当前工作树的未跟踪声学文件作为输入。核对Git commit ancestry、源码哈希、SciPy/Numpy/Python版本、指定CSV字段，形成`PREC2_PREFLIGHT.json`。任何源字节/包清单/科学代码/父提交不一致都立即STOP，记录原异常。

本次必须特别避免PREC1曾出现的“回复中MD SHA写错”。以ZIP中清单、START登记的**逐字节实际计算SHA256**为唯一执行版本。若任何哈希与用户侧预检登记不符，必须请求复核，不能自行以“相近版本”替代。

## 3. 一次冻结、三段研究

### PREC2-P0：PRT→MOD表示精度核查（全量稀疏打印记录）

运行包内` t2_prec2_readonly.py `（去除命令两侧空格），由独立MOD读取器按二进制结构提取所有返回波数。PRT按实际打印行序读取模态索引，每个PRT索引逐位比较`float32(k_PRT.real)==MOD_k.real`，并单独记录原MOD虚部、PRT打印吸收列、波数差及其60km相位量级。

要求输出全部六套PRT抽样；不以PRT的稀疏行插值补足全部模式。根据先前只读预研，期望每两套网格三频分别31/31/30行，合计184行；**这只是强预期和STOP检查，不是把预研当作本次实验执行**。

### PREC2-P1：独立有限模态实特征值计算

从K2每个频率的PRT打印行中，按固定的0起算列表位置`[0,7,15,22,-1]`选5个模式，不据结果择优或替换。每频五模式，共15条；K3用于网格跨来源对照但不替换已冻结选模。

数值方法：取2阶中心差分自伴三对角算子，顶部Dirichlet；底部Neumann ghost消元后作质量对称化，最后一个上/下对角元素系数`sqrt(2)/h^2`，其余`1/h^2`，对角`(omega/c(z_i))^2-2/h^2`。计算从大到小第j个特征值`lambda_j = k_r^2`，调用`scipy.linalg.eigh_tridiagonal`的指定索引，**不构造160001×160001密集矩阵**。

冻结四网格`N=[20000,40000,80000,160000]`，节点与原20m SSP节点整除对齐；节点之间按照原`C-linear`对声速插值，不引入另一套SSP。对lambda先做`h²`、再做`h⁴`与`h⁶`Richardson组合；输出每层原始本征值、外推后实波数、相邻高阶外推的差、与K2原MOD/PRT的偏差。禁止按对结果最有利的网格定义真值。

`|k_R4_coarse-k_R4_fine| <= 1.6666666667e-9`只称有限网格内部收敛**指标**，不能等同于严格误差界。必须记录不收敛的所有模式、弱数值样本和根排序异常。

### PREC2-P2：独立冷复核与有限连续算子交叉核查

以独立文件` t2_prec2_cold_audit.py `按另一路MOD/PRT解析器复核六套原始文件及184个历史抽样标量；另外固定选取201Hz五个模式的第1/3/5项，用Prüfer相位变量对C-linear速度函数直接积分：

\[
\theta'(z)=s\cos^2\theta+\frac{(\omega/c(z))^2-k_r^2}{s}\sin^2\theta,
\qquad s=\omega/1500,\quad \theta(0)=0,
\]

Neumann边界目标：`theta(5000)=(j-0.5)*pi`。采用`DOP853`，两档(max_step,rtol)=(1m,1e-11)、(0.5m,1e-12)，数值根有界Brent求解。该方法不读主程序Richardson实现，不进口旧模型求解代码；PRT及K2实部只可作为寻根初值，不能作为真实根。

对三个模式同时要求：两档IVP根一致于`1.6666666667e-9 m^-1`以内、IVP细根与FDM R6的偏差同样在该设计预算以内。若失败，记录`INDEPENDENT_LOSSLESS_ROOT_NOT_CORROBORATED`与该模式全部信息，保持失败，不得调整选模、扩大时长或改阈值。总冷核验窗口900s。

**特别注意：即使15条FDM和3条IVP数值稳定，也仅是“原名义SSP某些无损实模态”被交叉检查，远远不能证明全部有损复模态、完整声场或阵型定位性能。**

## 4. 数值稳定性、进度预算与中止

- P0/P1只读科学计算：不调用KRAKEN/KRAKENC、FIELD/FIELDS/SCOOTER/BELLHOP，不生成新波形、录音、噪声、DFT、MC、FIM或波束形成结果。
- 冻结主程序预算1800s、冷审计900s、峰RSS≤4GiB，单科学进程，BLAS线程1，Git新增≤64MiB。另在独立启动器层设置全流程wall≤2700s。D盘安全余量>1GiB。避免运行时擅自修复受保护代码。
- 任一源/哈希不符、PRT模式不一致、损失合同不符、数值非有限、独立冷审计失败均STOP；不得默改包内脚本或重跑正式科学流程。
- 记录所有曾经停止或报错的原日志；若正常程序错误可出具说明，但科学重跑必须另行批准。

## 5. 输出与判定

预检必须在**复制任务科学文件到仓库之前**运行：在D盘独立解压目录执行`python -B t2_prec2_preflight.py --repo <干净独立Git工作树路径> --package <本任务解压目录> --out D:/ProjectStorage/WSL-CZ/T2_PREC2/preflight`；以免未跟踪任务文件妨碍第一次工作树清洁检查。通过后将任务文件复制到该独立工作树、把预检记录复制到新结果目录并提交Commit A。

至少交付：

`PREC2_PREFLIGHT.json`、`PREC2_SOURCE_IDENTITIES.csv`、`PREC2_PRT_MOD_ALL_SAMPLES.csv`、`PREC2_PRT_MOD_SUMMARY.csv`、`PREC2_LOSSLESS_SUBSET_FDM.csv`、`PREC2_PRUFER_COLD_ROOTS.csv`、`PREC2_COLD_CHECKS.csv`、`PREC2_COLD_AUDIT.json`、`PREC2_DECISION.json`、`PREC2_SCOPED_INTERPRETATION.md`、`PREC2_NEXT_RESEARCH_DECISION.md`、`GPT_SYNC_PREC2.md`、原始完整日志、SHA索引、预算与失败审计。

即使全通过，本轮最多判为`PREC2_SPARSE_OUTPUT_LOSS_AND_SUBSET_LOSSLESS_REAL_K_CORROBORATED`。若P0通过但独立实根失败，只能`PREC2_SPARSE_ROUNDING_SUPPORTED_SUBSET_NOT_CORROBORATED`。若源身份失败，`PREC2_INPUT_NOT_ADMITTED`。不允许判定`COMPLETE_COMPLEX_PROPAGATION_CERTIFIED`。

研究总结必须报告旧N2窗口效应99%级、N3单精度解析Gate未通过、PREC1理想波导的归因及本轮变化，说明窗口缺漏、真实底质未知、密度/吸收/边界的适用范围。T2-A/B实际接收和定位性能`NOT_EVALUATED`，R4=0%。原G0/N1/N2/N3/PREC1表格和代码不得修改。只有新目录与R4_MASTER追加审计行可以更改。

## 6. Commit与结束

在干净的独立工作树基于固定父提交操作，历史SHA引用必须明确。先冻结以上任务包全部源文件与预检、配置（Commit A）；再仅运行一次正式主程序，成功后运行一次独立冷审计，出报告（Commit B）；推送并核对remote/main==Commit B后STOP。不得在Commit B之后自动启动全模态高精度侧车、E-STD完整声场、海底模型替换或T2实际性能仿真。

主程序（按项目根目录解析相对路径；`out`为绝对或相对均可）：

```powershell
python -B t2_prec2_readonly.py --repo . --config results/R4_T2_PREC2_ESTD_EIGENVALUE_PRECISION_FEASIBILITY/T2_PREC2_EXPERIMENT_SPEC_20261014_v1.json --out results/R4_T2_PREC2_ESTD_EIGENVALUE_PRECISION_FEASIBILITY
```

冷审计：

```powershell
python -B t2_prec2_cold_audit.py --repo . --config results/R4_T2_PREC2_ESTD_EIGENVALUE_PRECISION_FEASIBILITY/T2_PREC2_EXPERIMENT_SPEC_20261014_v1.json --out results/R4_T2_PREC2_ESTD_EIGENVALUE_PRECISION_FEASIBILITY
```

这两个命令是执行入口示例，须在Commit A后由一次性守卫按顺序调用；若冻结文件哈希失败严禁启动。
