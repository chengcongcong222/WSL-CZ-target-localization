# AUX Gate2B joint static error budget

AUX1_JOINT_STATIC_REQUIREMENT_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED

Primary target remains T5-family lambda=0.50; T10 uses its own distinct vector and ladder.

| Family | Anchor | lambda | Worst P95 (%) | T5 | T10 |
|---|---|---:|---:|---|---|
| T5_BUDGET_FAMILY | A | 0 | 2.903052 | True | True |
| T5_BUDGET_FAMILY | A | 0.25 | 3.251731 | True | True |
| T5_BUDGET_FAMILY | A | 0.5 | 4.048102 | True | True |
| T5_BUDGET_FAMILY | A | 0.75 | 5.130851 | False | True |
| T5_BUDGET_FAMILY | A | 1.0 | 6.424049 | False | True |
| T5_BUDGET_FAMILY | B | 0 | 3.122284 | True | True |
| T5_BUDGET_FAMILY | B | 0.25 | 3.339383 | True | True |
| T5_BUDGET_FAMILY | B | 0.5 | 3.856566 | True | True |
| T5_BUDGET_FAMILY | B | 0.75 | 4.579416 | True | True |
| T5_BUDGET_FAMILY | B | 1.0 | 5.493994 | False | True |
| T10_BUDGET_FAMILY | A | 0 | 2.903052 | True | True |
| T10_BUDGET_FAMILY | A | 0.25 | 3.947383 | True | True |
| T10_BUDGET_FAMILY | A | 0.5 | 5.948927 | False | True |
| T10_BUDGET_FAMILY | A | 0.75 | 8.356482 | False | True |
| T10_BUDGET_FAMILY | A | 1.0 | 11.196358 | False | False |
| T10_BUDGET_FAMILY | B | 0 | 3.122284 | True | True |
| T10_BUDGET_FAMILY | B | 0.25 | 4.167701 | True | True |
| T10_BUDGET_FAMILY | B | 0.5 | 6.288599 | False | True |
| T10_BUDGET_FAMILY | B | 0.75 | 8.936592 | False | True |
| T10_BUDGET_FAMILY | B | 1.0 | 12.042534 | False | False |

## Interpretation and limitations

True MAIN=(0,0), AUX=B[-sin(beta),side*cos(beta)], target=(R,0). Bearings use actual geometry, independent Gaussian random terms plus fixed b1=common-half_diff, b2=common+half_diff. Estimated node positions have independent per-axis Gaussian offsets. Deployment and navigation are separate factors, both active. Random error is not scaled by lambda. Range=norm(target_hat-MAIN_est), compared with true R; absolute2D position error is reported separately.

Each cell has30000 trials; one new frozen PCG64 seed and common six-column standard-normal table across both families. Mirrors reflect bearing noise and navigation y. Every deterministic bias/deployment sign corner is retained. Shared draws/mirrors are not independent repetitions. All parallel, nonfinite and behind-sensor failures remain infinite errors in nearest-rank quantiles. No outlier rejection. The max over11 ranges, both mirrors and8 signed corners is the joint Gate metric, not pooled statistics.

Lambda0 uses one sign identity per family; independent reproduction of Gate2A saved controls and a preregistered two-sample DKW criterion (familywise alpha=.001 over44 unique controls) check the new seed. Failure means EXECUTION_INVALID. This sampling diagnostic is not hardware assurance.

Numerical joint requirements apply to frozen static model, tested signed corners, scaling values and11 discrete50:1:60km ranges. No continuous hyperbox, intermediate magnitude/range, correlated-error, dynamic/time-skew or association-failure certificate. A maximum tested lambda is a discrete observed frontier, not an interpolated threshold. The half-point remains primary even when smaller lambda passes.

Actual auxiliary array, DOA, navigation/attitude, calibration and clocks remain UNKNOWN. Time/association NOT_NUMERICALLY_VALIDATED. No automatic Gate2C or further numerical architecture work. Next action is client requirement confirmation or a research-lead architecture decision. R4=0%; R4-A1-NEW NOT_OPENED; depth closed. Historical Gate0/1/2A artifacts preserved; client table supersedes the operational presentation without rewriting frozen Gate2A files.

```json
{
  "stage": "R4_AUX_GATE2B_JOINT_STATIC_ERROR_BUDGET",
  "parent_sha": "5d8b9c9f2d8efdb74a95da4ea546aae11457a24f",
  "Gate2A_independent_audit": "ACCEPTED_SINGLE_FACTOR_REQUIREMENTS",
  "primary_target": "HALF_T5_JOINT_ENGINEERING_POINT",
  "primary_decision": "AUX1_HALF_T5_JOINT_BUDGET_ESTABLISHED",
  "Gate2B_decision": "AUX1_JOINT_STATIC_REQUIREMENT_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED",
  "half_point_pass": {
    "T5_BUDGET_FAMILY": {
      "A": true,
      "B": true
    },
    "T10_BUDGET_FAMILY": {
      "A": true,
      "B": true
    }
  },
  "maximum_tested_lambdas": [
    {
      "family": "T5_BUDGET_FAMILY",
      "anchor": "A",
      "maximum_tested_lambda_T5": 0.5,
      "maximum_tested_lambda_T10": 1.0
    },
    {
      "family": "T5_BUDGET_FAMILY",
      "anchor": "B",
      "maximum_tested_lambda_T5": 0.75,
      "maximum_tested_lambda_T10": 1.0
    },
    {
      "family": "T10_BUDGET_FAMILY",
      "anchor": "A",
      "maximum_tested_lambda_T5": 0.25,
      "maximum_tested_lambda_T10": 0.75
    },
    {
      "family": "T10_BUDGET_FAMILY",
      "anchor": "B",
      "maximum_tested_lambda_T5": 0.25,
      "maximum_tested_lambda_T10": 0.75
    }
  ],
  "actual_hardware_capability": "UNKNOWN",
  "actual_auxiliary_array_bearing_capability": "UNKNOWN",
  "actual_navigation_attitude_capability": "UNKNOWN",
  "time_synchronization": "NOT_NUMERICALLY_VALIDATED",
  "target_association": "NOT_NUMERICALLY_VALIDATED",
  "automatic_further_numerical_work": "STOP",
  "AUX_GATE2C": "NOT_OPENED",
  "next_action": "CLIENT_REQUIREMENT_CONFIRMATION_OR_RESEARCH_LEAD_ARCHITECTURE_DECISION",
  "R4_A1_NEW": "NOT_OPENED",
  "R4_progress_percent": 0,
  "depth": "CLOSED",
  "new_acoustic_propagation": 0,
  "new_depth_scores": 0,
  "client_confirmations_required": [
    "auxiliary bearing-producing array / aperture / directed-bearing capability",
    "per-node random DOA / heading / calibration error budget",
    "residual common and differential bearing biases",
    "per-axis navigation position accuracy",
    "deployable5/7km lateral baseline and deployment orientation",
    "target-bearing association and front-back resolution",
    "common time reference / clock synchronization",
    "array calibration and relative reference alignment"
  ],
  "audit_status": "PENDING_RESEARCH_LEAD_AUDIT",
  "scope": "STATIC_SNAPSHOT; TWO_FROZEN_FAMILIES; TESTED_SIGN_CORNERS_AND_11_DISCRETE_RANGES_ONLY; NO_CONTINUOUS_HYPERBOX_GUARANTEE",
  "design_sha": "8c8c0e52df09f4cf1111201396d6a2dec3c5f8b0"
}
```
