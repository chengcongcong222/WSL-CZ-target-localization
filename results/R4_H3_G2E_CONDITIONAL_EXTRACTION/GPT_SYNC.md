# H3-G2E conditional frequency-statistic synchronization

Parent: 3f01469c4a718413d7368db6c564a269c0ca5ce7.
Design commit: DESIGN_COMMIT_VERIFICATION.json.
Execution commit: the commit containing this report. Remote verification is outside the repository to avoid self-SHA recursion.

Decision: H3_G2E_COMPRESSED_RESPONSE_EXTRACTION_UNRELIABLE. R4 remains 0%; submitted then STOP.

The full Gaussian sample/joint-CSD chain passed 1046809 cold checks, 0 FAIL.
All 31104 registered likelihood witnesses preserve every channel and selected cross-frequency block.
Maximum relative NLL difference: 1.8189894e-12.
Sample CSD is never inverted, ridge is never added and samples are never forced to rank one.
K=1 is a singular sample statistic with a positive-definite MODEL covariance.

- K=1: accepted 749/1728; failed template/frequency blocks 4397; oracle failures 4388.
- K=8: accepted 1470/1728; failed template/frequency blocks 1537; oracle failures 1541.
- K=32: accepted 1488/1728; failed template/frequency blocks 1358; oracle failures 1356.

Overall diagnostic quality accepted 3707/5184; grid-unstable cells 3.
Quality means every template/frequency block exceeds the energy gate and unit-projector error <=0.1.
This is a frozen extraction diagnostic, not a new project performance PASS.
The per-block conservative 1% Markov null bound can reject informative samples.
Failure does not establish physical information absence.
Every actual/oracle rejection, before/after error and background-estimation error remains in the records.

Resource-matched B1/B2 accepted 1236/1728 vs 1233/1728.
Median paired B1-minus-B2 maximum detected projector error -0.00011983946398212882.
A response-stability comparison does not demonstrate depth localization.
EXTRACTION_SUMMARY.csv includes every scene/resource/K/noise/package.
NORMALIZED_RESPONSE_BIAS_VARIANCE.csv labels detected-only moments and retains all 16-realization denominators.

Noise is estimated from 256 independent background samples, using full empirical covariance.
Oracle noise is separate. Hann 249/250-Hz correlation is +1/6.
Likelihood witnesses use synthetic non-truth fixed-gain candidates and template-dependent unknown source powers;
they are statistical algebra checks, NOT physical model fits or recovered source/calibration values.
Full CSD is sufficient within the SAME zero-mean Gaussian covariance family.
No compressed-response Fisher efficiency is claimed without a matched finite-sample model.

Unknown source scalar is removed by normalization; unknown fixed channel/frequency gains remain in the observed direction.
Physical propagation response therefore is not identified without external calibration constraints.
Source/gain gauges, separable-component depth absorption, free-modal zero and C2 zero remain explicit.

New stochastic simulation: 48 independent innovation pools (3 geometries x16 realizations), each containing 3 static template draws.
192 scaled noise conditions, K prefixes, resource/frequency masks and two grids are paired derivatives.
5184 diagnostic cells are NOT 5184 independent experiments.
Saved source/noise/background arrays, locked fields and raw Y/background SHA records reconstruct every observation.
No new propagation, localization-error Monte Carlo, time-series extraction or 64-s stationarity evidence was generated.

Execution 107.19s; execution plus cold audit 325.90s.
Stored bytes at audit 213022091; cap 1,000,000,000.
Independent algorithms are internal; external research-lead acceptance remains pending.

H3-G0/G1 conclusions unchanged. Actual UUV source occupancy, constrained propagation, environment/pose robustness,
full RC2 support and continuous depth performance NOT_ESTABLISHED. Time-domain extraction NOT_OPENED.
No next module is authorized.

Formula references: [SciPy spectral analysis](https://docs.scipy.org/doc/scipy/tutorial/signal.html#spectral-analysis)
and [Gaussian spectral-likelihood research](https://pmc.ncbi.nlm.nih.gov/articles/PMC10050575/).


## 中文同步摘要

正式分类：H3_G2E_COMPRESSED_RESPONSE_EXTRACTION_UNRELIABLE。完整频域样本—联合CSD统计链通过；压缩响应按冻结标准质量通过3707/5184。

|K|质量通过|登记单元|未通过能量门限的频率/模板块|
|---|---:|---:|---:|
|1|749|1728|4397|
|8|1470|1728|1537|
|32|1488|1728|1358|

上述门限是预登记的保守提取诊断，失败不能解释为物理不可辨识或深度信息不存在。B1/B2只能作资源匹配的提取稳定性比较，不能由质量通过数重判G0深度信息。

所有相对/固定噪声底、弱场、K=1失败均保留。未知源标量经归一化去除，但未知阵元频率增益仍留在响应中；未经校准的传播响应恢复未建立。检测条件下的偏差/方差明确保留全16次实现分母。

冻结设计SHA：2e0515359ed26248d89db9f3821855a295cb9702。执行SHA为包含本报告的提交，最终remote/main核对记录放在D盘运维目录。
H3-G0/G1原结论不变；实际UUV源占用、传播未知量和完整RC2深度集合未建立；R4=0%。提交推送后停止。


### 登记噪声条件下的提取诊断

|噪声模型|设计比例|质量通过|登记单元|
|---|---:|---:|---:|
|RELATIVE|0.01|1208|1296|
|RELATIVE|0.05|920|1296|
|ABSOLUTE_FLOOR|0.01|1132|1296|
|ABSOLUTE_FLOOR|0.05|447|1296|

1%/5%相对比例及固定绝对噪声底均为设计假设，不代表真实源级或真实SNR。噪声底提高时的失败不应概括为所有H3观测不可用。保守能量门限也可能拒绝仍含信息的样本。
