# AUX client confirmation package and architecture freeze

AUX1_JOINT_STATIC_REQUIREMENT_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED_ACCEPTED

This stage is document closeout only. Historical Gate0/1/2A/2B evidence remains frozen; no numerical work executed.

Package contents: source-bound requirement provenance; editable client matrix with UNKNOWN values; Chinese C1-C9 questionnaire; predeclared branches A/B/C; current architecture state; Gate2B audit acceptance.

[Audit acceptance](AUX_GATE2B_AUDIT_ACCEPTANCE.md)
[Requirement provenance](AUX_REQUIREMENT_PROVENANCE.csv)
[Client matrix](AUX_CLIENT_CONFIRMATION_MATRIX.csv)
[Decision after client](AUX_ARCHITECTURE_DECISION_AFTER_CLIENT.md)
[Current status](R4_CURRENT_ARCHITECTURE_STATUS.json)

Strongest tested T5 vectors: A(lambda.50):5km,.05deg random,.05deg common,.025deg HALF-diff,25m/axis,10deg beta,P95=4.0481%; B(lambda.75):7km,.075deg random,.075deg common,.0375deg HALF-diff,37.5m/axis,15deg beta,P95=4.5794%. T10 family bothlambda.75:common and HALF-diff0.075deg,beta15deg,nav A75m/B150m,P95 A8.3565%/B8.9366%. Every row is a whole vector; no cross-row selection.

STATIC instantaneous geometry; independent zero-mean Gaussian random bearing and per-node2D navigation; frozen fixed bias/deployment sign corners; 11 discrete ranges50:1:60km; both mirrors. No continuous hyperbox/interior/range guarantee; no hardware, time-sync, association or dynamic admission.

Hardware UNKNOWN; time/association NOT_NUMERICALLY_VALIDATED. Package prepared for the research lead to send; NOT_SENT. Numerical architecture work STOPPED; no automatic Gate2C; no A1-NEW or depth reopening; R4=0%. The three client branches require evidence mapping and explicit lead authorization, including separate time/association/chain review.
