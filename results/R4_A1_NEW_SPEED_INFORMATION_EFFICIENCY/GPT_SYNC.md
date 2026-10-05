# 速度信息效率审计

CURRENT_1200S_GEOMETRY_INFORMATION_LIMIT_CONFIRMED_BY_CRLB

DEVELOPMENT / INFORMATION EFFICIENCY AUDIT。复用 12000 个保存场景；E0 为 D3 逐 run 恒等，E1 只将 loss 改为 linear，E2 为 ORACLE_INITIALIZATION_DIAGNOSTIC_ONLY。新增 MC、truth case、bearing noise draw 均为 0；R4-A1=0%，R4=0%，depth NOT_OPENED。

| Anchor | Variant | Worst speed P95 % | PROJECT speed | Worst range P95 % | Bearing P95 deg | Heading P95 deg | Max failure rate |
|---|---|---:|---:|---:|---:|---:|---:|
| A | E0 | 24.594634 | 0/12 | 0.573956 | 0.019953 | 2.177292 | 0.000000 |
| A | E1 | 24.404747 | 0/12 | 0.552253 | 0.020786 | 2.241919 | 0.000000 |
| A | E2 | 24.404746 | 0/12 | 0.552253 | 0.020786 | 2.241919 | 0.000000 |
| B | E0 | 26.723484 | 0/12 | 0.585828 | 0.028535 | 2.532162 | 0.000000 |
| B | E1 | 25.962868 | 0/12 | 0.586871 | 0.027814 | 2.559821 | 0.000000 |
| B | E2 | 25.962868 | 0/12 | 0.586871 | 0.027814 | 2.559822 | 0.000000 |

## 模型与统计口径

四状态模型使用真实节点位置及已移除生成 bias 的同一批随机 bearing noise，共 242 measurements。E1 保留原 ray-intersection + IRLS initializer、状态缩放与解析 Jacobian，TRF、max_nfev=100、ftol/xtol/gtol=1e-10。没有裁点、truth warm start、truth 邻域边界、multistart、额外先验或事后调参。E2 唯一区别是真实 Cartesian 初值，不能形成 estimator claim。

F=HᵀH/sigma_theta²；原状态缩放 measurement Jacobian 的 s/smax>1e-10 为数值秩门限。已知 sigma，不用 fitted residual variance rescaling。physical F/C、奇异值及条件数按 matrix_index 保存在 CRLB_MATRIX_RECORDS.npz 与 CRLB_GEOMETRY_RUN_RESULTS.csv。满秩时缩放 SVD 恢复的 C 为 F 的逆，非满秩保留 pseudoinverse，但不确定度必须报告无穷，不能给有限下界。

CRLB 覆盖每个保存布放角及双镜像。每 case 报告 min/median/P95/max，route 使用 max。速度梯度 [0,0,vx/v,vy/v]，报告 1σ、1.645σ、1.96σ；主要比较为 1.96σ。所有等效值标为 LOCAL_LINEAR_GAUSSIAN_EQUIVALENT / NOT_EMPIRICAL_P95 / NOT_GLOBAL_GUARANTEE。r/theta/psi 使用 delta method，只作交叉诊断，不重定义 Gate。

效率比 eta_v 为平均条件局部 CRLB variance / empirical signed speed variance；sample std 使用 ddof=1，偏差、median、RMSE 同时报告。仅全部 run 有效且方差有限时定义；它是局部无偏方差参考，不是总误差 bound。失败保留 INF 并进入分位数/Gate；有限样本的偏差统计显式标为描述性统计。

E1/E2 最大逐 case P95 差为 0.000007522 个百分点，冻结门限为 2 个百分点。逐 run 状态、速度与 cost 差另存 paired basin diagnostics；P95 接近不等于证明全局 MLE 覆盖。

A: 最坏局部 1.96σ 速度等效值=26.208223%；超过 10% 的 case=12/12，超过 5%=12/12。
A E1 mean_speed_bias_mps: case min..max=-0.013423078..0.023798968。
A E1 median_speed_bias_mps: case min..max=-0.022871136..0.018726962。
A E1 speed_std_mps: case min..max=0.167199048..0.238862268。
A E1 speed_RMSE_mps: case min..max=0.167031928..0.239000526。
A E1 eta_v: case min..max=0.922617015..1.186498140。
B: 最坏局部 1.96σ 速度等效值=29.485364%；超过 10% 的 case=12/12，超过 5%=12/12。
B E1 mean_speed_bias_mps: case min..max=-0.011020030..0.018670130。
B E1 median_speed_bias_mps: case min..max=-0.020849937..0.016219646。
B E1 speed_std_mps: case min..max=0.175402210..0.246917468。
B E1 speed_RMSE_mps: case min..max=0.175518128..0.247375976。
B E1 eta_v: case min..max=0.871675120..1.143692497。

## 信息累积与转向窗口

T=300/600/900/1200 s 只截取同一 trajectory 做 analytic FIM，不运行 prefix optimizer，也不外推 1800/2400 s。STRAIGHT 包含 t<=600 的 61 epochs；POST_TURN 为 t>600 的 60 epochs；FULL 为 121 epochs。边界不重复，F_full=F_straight+F_post，所有窗口沿用同一 x0 和全局 t。

速度有效信息使用 Schur complement Fvv-Fvp Fpp^{-1} Fpv 消去初始位置。full 相对 straight 的增益同时包含新增时长、观测数和现有改变的几何，不是 isolated 15° turn 因果效应；没有额外直航 counterfactual，不能单独量化转向效应。

A PREFIX_300: 最坏局部 1.96σ 速度等效值=203.459483%。
A STRAIGHT_600: 最坏局部 1.96σ 速度等效值=73.538634%。
A PREFIX_900: 最坏局部 1.96σ 速度等效值=40.274655%。
A FULL_1200: 最坏局部 1.96σ 速度等效值=26.208223%。
A full/straight: 局部速度方差下降 min/median/max=0.870509590/0.871903319/0.873534280；Schur trace 倍数=7.759875600/7.791321137/7.834114074。
B PREFIX_300: 最坏局部 1.96σ 速度等效值=228.888730%。
B STRAIGHT_600: 最坏局部 1.96σ 速度等效值=82.735084%。
B PREFIX_900: 最坏局部 1.96σ 速度等效值=45.310030%。
B FULL_1200: 最坏局部 1.96σ 速度等效值=29.485364%。
B full/straight: 局部速度方差下降 min/median/max=0.870503071/0.871799021/0.873491201；Schur trace 倍数=7.765236600/7.797767785/7.841250023。

## 路线与作用域

Outcome II 只在 E1 非全部通过、全部 primary case E1/E2 P95 差<=2pp、至少一个 primary case 的局部 1.96σ 速度等效值>10% 时成立。准确结论是 local Gaussian CRLB equivalent itself above target；不是所有算法的经验 P95 下界，也不是被动双节点的普遍物理不可辨识。B 不能替代 A。

本轮路线判定：CURRENT_1200S_GEOMETRY_INFORMATION_LIMIT_CONFIRMED_BY_CRLB；仅建议下一阶段 OBSERVATION_DESIGN。OBSERVATION_DESIGN 对应 R4_A1_NEW_SPEED_OBSERVATION_DESIGN_GATE，需要另行授权。本轮禁止改变 duration/AUX 相对运动/baseline/turn/RMS/truth，不做 Kalman/smoother/Doppler/range-rate/TDOA/depth/acoustics；push 后 STOP。

既有 range/bearing STRONG、heading PROJECT 结论保留。本轮为 PENDING_RESEARCH_LEAD_INDEPENDENT_AUDIT，不给 application credit。

```json
{
  "stage": "R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY_AUDIT",
  "parent_SHA": "6180e27edadb69aec1f68959b60bdb2036e719e5",
  "parent_audit": "ACCEPTED_CURRENT_1200S_BOTTLENECK",
  "primary_anchor": "A",
  "secondary_anchor": "B",
  "nature": "DEVELOPMENT_INFORMATION_EFFICIENCY_AUDIT_NOT_CONFIRMATION",
  "route_decision": "CURRENT_1200S_GEOMETRY_INFORMATION_LIMIT_CONFIRMED_BY_CRLB",
  "next_recommended": "OBSERVATION_DESIGN",
  "observation_design_stage": "R4_A1_NEW_SPEED_OBSERVATION_DESIGN_GATE",
  "PROJECT_speed_cases": {
    "E0": 0,
    "E1": 0,
    "E2": 0
  },
  "E1_STRONG_speed_cases": 0,
  "E1_E2_max_case_P95_gap": 7.521796008269632e-08,
  "CRLB_cases_above_10": 12,
  "CRLB_cases_above_5": 12,
  "CRLB_scope": "LOCAL_LINEAR_GAUSSIAN_EQUIVALENT;NOT_EMPIRICAL_P95;NOT_GLOBAL_GUARANTEE",
  "information_finding_scope": "At least one saved primary geometry has local Gaussian 1.96-sigma speed equivalent above target; not an empirical P95 lower bound or universal physical impossibility.",
  "new_Monte_Carlo_realizations": 0,
  "new_truth_cases": 0,
  "new_bearing_noise_draws": 0,
  "R4_A1_percent": 0,
  "R4_percent": 0,
  "depth": "NOT_OPENED",
  "automatic_next_stage": false,
  "stop_after_commit_B": true,
  "audit_status": "PENDING_RESEARCH_LEAD_INDEPENDENT_AUDIT",
  "policy_SHA": "ec2a5b32205b4cab444b435f2b993efae80112f2"
}
```
