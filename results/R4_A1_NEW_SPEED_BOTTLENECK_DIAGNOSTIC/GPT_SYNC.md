# 速度瓶颈诊断：证据解释与停止边界

正式开发诊断判定为 `SPEED_INFORMATION_INSUFFICIENT_AT_CURRENT_1200S_APPLICATION_SCENARIO`。本轮复用全部 12000 个保存场景，五分支共 60000 条记录；独立重建 1030883 项检查全部 PASS、0 FAIL，4 项实现测试通过。未生成新随机样本，未增加应用性能信用。

A 的 D0–D4 速度 PROJECT 均为 0/12。D0 最坏速度 P95 为 28.9386%；D1 移除生成偏置后为 29.1250%；D2 使用真实导航位置后为 24.6826%；D3 同时移除偏置和导航误差后为 24.5946%；D4 六状态候选为 29.5088%。D4 STRONG 速度也是 0/12。B 为次级证据，不能替代 A 的主判定。

D1 未恢复速度，说明仅移除静态偏置不足以闭合当前瓶颈。D2/D3 相比 D0 改善，但均未达到全部 12 个 case 的 10% 门限，因此不能宣称导航精度改进已经足以建立速度指标。D3 在只保留原始随机方位噪声后仍失败，支持当前噪声、几何与四状态估计器组合下的残余速度瓶颈。它不是对所有算法的误差下界，也不是物理不可辨识证明。D1–D3 均为 ORACLE_DIAGNOSTIC_ONLY；各分支 P95 改变量不是可相加的方差分解。

A 各 case 的绝对速度误差 P95：D0 为 0.3586–0.5382 m/s，D3 为 0.3275–0.4851 m/s，D4 为 0.3679–0.5540 m/s。相对误差受固定 case 的真速度影响，但绝对误差也未消失。这些区间仅属于本轮已测 case，不是普遍噪声下限。

D4 在 A、B 各 6000 次求解均无失败，全部 case 失败率为 0。A 的其余指标最坏 P95 为：range 2.2050%、bearing 0.1295°、heading 2.6413°，仍达到 PROJECT；速度未通过，所以不准入 fresh validation。

D4 的偏置状态未得到稳定恢复：A 各 case 共模触界率 97.0%–99.2%，差模半差触界率 99.8%–100%；共模估计绝对误差 P95 为 0.09378–0.09617°，差模半差为 0.04689–0.04806°。这提示有界偏置状态与运动状态的分离很弱，不能把触界的 nuisance 估计当作校准成果，也不能将该结果解释为所有偏置估计器都不可能成功。

三个窗口的缩放 Jacobian 在两组均数值满秩 100%，阈值为 s_min/s_max > 1e-10。A 的奇异值比中位数：STRAIGHT 4.2669e-7、POST_TURN 4.2628e-7、FULL 4.4056e-7；FULL 范围为 3.7058e-7–5.2107e-7。FULL 比例只小幅改善，满秩不是速度精度证书；这些诊断使用冻结尺度和原始测量 Jacobian，不能升级为额外机动实验结论。

按照冻结规则，执行提交后 STOP。不开 fresh validation，不做 FIX2，不改变观测时长、转向、几何、噪声或门限，不开放 depth、声学增量或其他支线。既有 range/bearing STRONG、heading PROJECT 结论保持；R4-A1=0%，R4=0%。本轮等待研究负责人独立审计。

---

# A1-new speed bottleneck causal diagnostic

SPEED_INFORMATION_INSUFFICIENT_AT_CURRENT_1200S_APPLICATION_SCENARIO

DEVELOPMENT / CAUSAL DIAGNOSTIC. All12000savedrealizations reused, no new seeds/truth/draws/performanceclaim. D0 is byte-identical saved output; D1 removes generatedbias only; D2 uses true positions only; D3 does both; D4 is the sole observation-only candidate. D1-D3 are ORACLE_DIAGNOSTIC_ONLY and cannot be algorithms or project claims.

| Anchor | Branch | Worst speed P95 % | Absolute speed P95 m/s (case min..max) | PROJECT speed cases | Worst range P95 % | Bearing P95 deg | Heading P95 deg | Max case failure rate |
|---|---|---:|---|---:|---:|---:|---:|---:|
| A | D0 | 28.938635 | 0.358552..0.538160 | 0/12 | 1.300624 | 0.085321 | 2.563957 | 0.000000 |
| A | D1 | 29.125026 | 0.357143..0.536438 | 0/12 | 0.645724 | 0.058261 | 2.571439 | 0.000000 |
| A | D2 | 24.682576 | 0.328140..0.487592 | 0/12 | 1.239928 | 0.065251 | 2.191904 | 0.000000 |
| A | D3 | 24.594634 | 0.327543..0.485050 | 0/12 | 0.573956 | 0.019953 | 2.177292 | 0.000000 |
| A | D4 | 29.508847 | 0.367876..0.553986 | 0/12 | 2.204995 | 0.129502 | 2.641340 | 0.000000 |
| B | D0 | 31.169332 | 0.402858..0.548973 | 0/12 | 1.352744 | 0.132820 | 2.850499 | 0.000000 |
| B | D1 | 30.729516 | 0.398577..0.548149 | 0/12 | 0.665325 | 0.087732 | 2.806767 | 0.000000 |
| B | D2 | 26.695240 | 0.344123..0.484908 | 0/12 | 1.290274 | 0.095219 | 2.549061 | 0.000000 |
| B | D3 | 26.723484 | 0.338725..0.481281 | 0/12 | 0.585828 | 0.028535 | 2.532162 | 0.000000 |
| B | D4 | 31.831692 | 0.408334..0.555362 | 0/12 | 2.372862 | 0.198975 | 2.923265 | 0.000000 |

## Attribution scope

The generator includes fixed realization-wise common/HALF-differential biases while historical4stateestimator has no bias nuisance: MISSPECIFIED_4_STATE_MODEL_UNDER_STATIC_BEARING_BIAS. The mismatch was frozen as a hypothesis, not a predeclared cause. Interventions compare exactly paired savednoise/navigation/deployment; P95 reductions are nonadditive and not a variance decomposition or unique-cause proof. Failure counts remain infinite errors, no sample rejection/replacement. Absolute speed error and relative error share fixed casev, their P95 relation is exact within a case; compare acrosscases without asserting a universal floor.

D4 uses E1 from original bearings/noisy positions and starts both biases at0. Bounded unregularizedTRF, analytic6stateJacobian, soft_l1f_scale1.5, scales50km/2m/s/designbiasbound, max300evaluations. No truthbias/targetinitialization, historicalrunwarmstart, multistart or added penalty. Bounds are design assumptions, not calibrated hardware knowledge. D4 bias estimates are nuisance diagnostics, not project outputs.

The full-rank threshold is scaled s/smax>1e-10, evaluated in61straightepochs/60post-turnepochs/full121. Full-rank is not an accuracy certificate; extremely small ratios and bound-hitting rates may expose weak bias/state separation. Turning comparison reuses the same15degturn, no excitation sweep. Fullwindow rank improvement alone does not prove PROJECTprecision.

D3 failures support the requested SPEED_INFORMATION_INSUFFICIENT_AT_CURRENT_1200S_APPLICATION_SCENARIO label only under this frozen geometry/randomnoise and current4stateestimator; they do not establish a lower bound for allpossibleestimators or physical nonidentifiability. D1passwould demonstrate calibrationrescue onthispanel, not a realizableoraclealgorithm. D4passonlyadmitsfreshvalidation, noA1credit.

New MC=0; R4-A1=0%; R4=0%; depthNOT_OPENED; noacoustic/depthscore. Existing A1submetrics preserved:range/bearingSTRONG,headingPROJECT, speedbottleneck. Stop after executionpush; noFIX2, nofreshconfirmation automatically.

```json
{
  "parent_A1_new_SHA": "e39182ac82b988e71b79d73a06a97076ca76f6bf",
  "A1_new_independent_audit": "ACCEPTED_SPEED_BOTTLENECK",
  "stage": "R4_A1_NEW_SPEED_BOTTLENECK_DIAGNOSTIC",
  "nature": "DEVELOPMENT_CAUSAL_DIAGNOSTIC_NOT_CONFIRMATION",
  "primary_anchor": "A",
  "secondary_anchor": "B",
  "primary_attribution": [
    "RANDOM_BEARING_PLUS_CURRENT_GEOMETRY_AND_4_STATE_ESTIMATOR_LIMIT_SPEED"
  ],
  "attribution_scope": "Frozen panel and current estimators only; oracle intervention findings are not unique-cause proofs. D3 failure is not universal physical nonidentifiability.",
  "speed_route_decision": "SPEED_INFORMATION_INSUFFICIENT_AT_CURRENT_1200S_APPLICATION_SCENARIO",
  "primary_PROJECT_speed_cases": {
    "D0": 0,
    "D1": 0,
    "D2": 0,
    "D3": 0,
    "D4": 0
  },
  "primary_STRONG_speed_cases_D4": 0,
  "D4_route_admitted": false,
  "D4_strong_speed_diagnostic": "NOT_ESTABLISHED",
  "new_random_seeds": 0,
  "new_truth_cases": 0,
  "new_noise_draws": 0,
  "new_Monte_Carlo_realizations": 0,
  "new_application_performance_claim": 0,
  "R4_A1_percent": 0,
  "R4_percent": 0,
  "range": "ESTABLISHED_AT_STRONG_TIER_IN_ACCEPTED_A1_NEW",
  "bearing": "ESTABLISHED_AT_STRONG_TIER_IN_ACCEPTED_A1_NEW",
  "heading": "ESTABLISHED_AT_PROJECT_TIER_NEAR_STRONG_IN_ACCEPTED_A1_NEW",
  "speed": "BOTTLENECK_DIAGNOSED_NO_SCIENTIFIC_CREDIT",
  "next_recommended_stage": "STOP",
  "automatic_fresh_validation": false,
  "automatic_FIX2": false,
  "depth": "NOT_OPENED",
  "new_acoustic_propagation": 0,
  "new_depth_score": 0,
  "audit_status": "PENDING_RESEARCH_LEAD_INDEPENDENT_AUDIT",
  "stop_after_commit_B": true,
  "policy_SHA": "4d7ddf444e2aca93a25aa25547a5a7f8bd6a093b"
}
```

## Jacobian A

STRAIGHT: full-rank rate=1.000000; ratio min/median/max=3.234689e-07/4.266921e-07/5.414271e-07.
POST_TURN: full-rank rate=1.000000; ratio min/median/max=3.193155e-07/4.262757e-07/5.206925e-07.
FULL: full-rank rate=1.000000; ratio min/median/max=3.705793e-07/4.405649e-07/5.210745e-07.

## Jacobian B

STRAIGHT: full-rank rate=1.000000; ratio min/median/max=7.205017e-07/9.608172e-07/1.184677e-06.
POST_TURN: full-rank rate=1.000000; ratio min/median/max=7.339980e-07/9.610101e-07/1.189333e-06.
FULL: full-rank rate=1.000000; ratio min/median/max=8.325745e-07/9.907028e-07/1.145871e-06.
