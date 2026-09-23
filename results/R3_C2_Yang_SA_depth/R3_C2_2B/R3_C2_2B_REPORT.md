# R3-C2.2B 报告：Yang 2015 主文公式锁定 + 复现协议锁定

UTC：2026-09-23T11:30:59.914058+00:00

## 0. 历史状态

- `C2_2A_PRIMARY_FORMULA_NOT_RECOVERED` -> **SUPERSEDED_BY_USER_PRIMARY_PDF**
- KRAKEN：`C2_1_PARSER_VALIDATED`
- 公式恢复：**PRIMARY_SOURCE_LOCKED**
- `YANG_ROUTE_UNDECIDED`（**不变**）
- 本轮 **禁止** Yang 性能 pass/fail

## 1. PDF 归档

见 `YANG2015_PRIMARY_PDF_NOTE.md`。SHA256 `b34952ead13fdaabdf9f1e7d87728fc2cb72fd09497240aac20ebdfc47f4b5b4`，10 页，JASA 1678-1686。

## 2-3. Eq.(1)-(10) 与 A1-A4

见 `YANG2015_EQUATION_LOCK.md`、`YANG2015_EQUATION_TABLE_LOCKED.csv`、`YANG2015_APPENDIX_EQUATIONS.md`。

全部乘法因子按原页锁定（含 e^{i pi/4}、sqrt(2 pi/(k_m r))、alpha_m、sinh、正则化逆）。

## 4. 代数闭环

见 `YANG2015_EQ1_EQ5_ALGEBRA_CHECK.md`。

**Eq3->Eq4(第一形式)->Eq5：PASS。**
Eq.(4) 第二形式与第一形式不一致（论文内部），实现禁用未校正第二形式。

## 5-8. b_m / Eq.(6) / 正则化 / 归一化

见 `YANG2015_BM_DEFINITION.md`、`YANG2015_EQ6_OBSERVABLE_NOTE.md`。

- b_m=(2 e^{-alpha_m r0})/(alpha_m k_m) sinh[(alpha_m dR)/2] phi_m(z_s)
- Eq.(6) 用 **复谱峰**；**不用**波数数值进深度 score
- phi_bar_inv = phi_m/(phi_m^2+Delta^2)，Delta **EMPIRICAL_NOT_UNIQUE**
- 归一化 D(z)/sum_grid D(z)

## 9. Erratum

`CARRY_FORWARD_FROM_PRIOR_ERRATUM_RECOVERY`；**delta=0** 基线。本轮无 Erratum PDF，未声称重核全文。

## 10-13. 复现配置 / Ref.7 / 峰识别 / 孔径

见 `YANG2015_PAPER_REPRO_CONFIG.md`、`YANG2015_REF7_DEPENDENCY.md`。

- 论文模拟：350 Hz、5 kn、z_s=4/50、z_r=18/70、span~**4990 m**、r1=5010（仅生成）
- **ORACLE_MODE_ID** vs **PRACTICAL_MODE_ID** 分层
- L > 2 pi/|k_n-k_m|，示例 lambda_12~2.3 km（示例，非 E-STD 阈值）
- 环境：**PAPER_REPRO_NEEDS_REF7**

## 14-15. 实现规范与自检

见 `YANG2015_IMPLEMENTATION_SPEC.md`、`YANG2015_IMPLEMENTATION_CHECKLIST.csv`。

FAIL：`Eq4_second_form_consistent`（已规避）。
UNRESOLVED：`appendix_vs_eq3_phase`（已记录，主链按 Eq.(3)）——**未忽略**。

## 16. 判定

### C2_2B_METHOD_SPEC_LOCKED

Yang 2015 主文 PDF 已由用户提供并完成视觉+文本双通道锁定：Eq.(1)-(10) 与附录 A1-A4 全部写入；Eq3->Eq4(第一形式)->Eq5 代数闭环 PASS；Eq.(6) 复谱峰 observable 与正则化逆 phi_bar_inv=phi_m/(phi_m^2+Delta^2) 已锁定，Delta=EMPIRICAL_NOT_UNIQUE；归一化 D/sum D 已锁定。论文复现条件已分 DIRECTLY_SPECIFIED / REFERENCED_TO_REF7；环境完整复现状态 PAPER_REPRO_NEEDS_REF7。Eq.(4) 第二形式与附录空间相位存在论文内部不一致，已记录且实现规避。未运行任何算法。YANG_ROUTE_UNDECIDED。

即使 `PAPER_REPRO_NEEDS_REF7`，**YANG_ROUTE 仍为 UNDECIDED，不得 pass**。

## 17. 停止声明

不跑 Yang / E-STD / z_hat / FWHM / PSL / delta 扫描 / Doppler 误差 / MC / AR / P5；不关闭任何候选。
