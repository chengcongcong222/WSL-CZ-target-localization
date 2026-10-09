# RC2 observation-defined multibranch outer support pilot

Parent: 08da59d7cc1b1b01b1dc27e8c6551257ab7ff530
Design SHA: 26b435d190c86fc2cbd976913bb67f2ebaf79e13
Execution SHA: commit containing this report;final remote record written after push.

## Architecture-specific decisions
{
  "U0": "RC2_OUTER_SUPPORT_VALID_BUT_BROAD",
  "U": "RC2_OBSERVATION_OR_NOISE_CONTRACT_INCOMPLETE"
}

U0 uses all21 pre-existing noisy P/Q case arrays; not H01/H06/H12. U dual-node is input-incomplete:draws,truths and point estimates are present but raw bearings/reported node positions are absent. No old scene(truth,...) calls or new observation draws filled the gap. Same-scene single/auxiliary comparison is unavailable. Hardware resources remain separate.

## Conservative support and statistical scope
Physical prior is r45..60km,theta-5..5deg,v1..3m/s,heading-15..15deg;exact prescribed1200s platform track,zero navigation/system bias and0.1deg121epoch Gaussian noise are inherited only for U0. Real navigation/systematic/temporal correlation and HLA extractability are not certified.

The compatibility set uses residual squared sum<=195sigma^2 and every residual<=5sigma plus registered model rounding slack. A reflection union retains possible left/right branches as a conservative expansion of the stored signed-angle model. The conditional simultaneous coverage lower bound over all21 fixed cases is 0.9619116273, from a self-contained Gaussian-integral Chernoff bound and marginal union bound, not these21 empirical successes. Cross-case independence is unnecessary;within-case iid Gaussian covariance is required by the archived contract.

The interval algorithm receives observations,time,platform only. Directed rounded rational Taylor/Machin enclosures bound Cartesian trajectories and bearings. Only certified incompatibility is discarded. Independent Decimal rejection replay and closed dyadic-tree reconstruction certify exported union inclusion. Finite corners/control tests are supplementary to the geometric proof.

## Widths and budget
{
  "range_km_width_min": 15.0,
  "range_km_width_max": 15.0,
  "speed_mps_width_min": 2.0,
  "speed_mps_width_max": 2.0,
  "theta_deg_width_min": 1.25,
  "theta_deg_width_max": 3.4375,
  "heading_deg_width_min": 30.0,
  "heading_deg_width_max": 30.0,
  "budget_exhausted_cases": 21,
  "retained_boxes": 32419,
  "truth_compatible": 21,
  "truth_retained": 21
}

All budget-unresolved and minimum-resolution leaves are included in SUPPORT_BOXES.csv. These are outer boxes, not proven feasible states. Disconnected branches cannot be lost by the union inclusion argument, but exact connected-component counts remain unresolved. Sign-sector labels are not a global modal count. Exhausting the subdivision budget does not turn a broad outer union into localization failure or physical non-observability.

Finite evaluation truth was read only after all support construction. Truth compatibility,actual retention and any compatible-but-deleted event are separate columns. Small-panel retention is not a sea-trial coverage certificate.

## Independent validation and H3 boundary
231500 checks,231500 PASS,0 FAIL;Decimal rejection certificates 5284. Review scope:full partition coverage and every discarded-box inclusion certificate,analytic co-linear/turn/auxiliary geometry controls and supplemental corner tests.

H3 finite25-node comparison cannot be certified across different panels. Those truth-centered nodes were not candidate inputs. H3-G0 conditional local information remains accepted;G1 unresolved,G2E compressed extraction unreliable,P0/P1 unchanged,C1e NOT_EVALUATED. No acoustic/depth or joint5D result was generated.

## Stop
New KRAKEN/FIELD/MC/noise draws/recordings=0;R4=0%. No automatic SSP,matching,time-domain,depth-MC or joint metric stage. U needs separately authorized archiving of actual bearings and reported tracks before a matching support trial. U0 broad/unfinished supports require a separately reviewed tightening or support/depth-set design, never truth supplementation or deletion of unsettled branches.

## Supplementary geometry and finite-grid scale

Weak0.15deg and full15deg maneuver controls are full rank but have different conditioning;straight motion is rank-deficient. Reflection and circular-seam mathematical controls pass. Actual H3 finite25-node spans are retained in H3_HORIZONTAL_GRID_SCALE_DIAGNOSTIC.csv;they are scale comparisons across different panels,not same-observation nesting. No new noisy observations were generated.
