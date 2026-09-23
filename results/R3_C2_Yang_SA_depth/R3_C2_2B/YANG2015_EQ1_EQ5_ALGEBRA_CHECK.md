# YANG2015_EQ1_EQ5_ALGEBRA_CHECK

UTC: 2026-09-23T11:30:59.914058+00:00

## 目标

验证 **Eq.(3) -> Eq.(1) -> Eq.(4) ->（k_r=k_m）-> Eq.(5)** 代数闭环，
捕捉 steering 相位、传播相位、alpha_m、sqrt(k_r k_m) 因子错误。

## 代入（S(r)~sqrt(r)）

Eq.(1)xEq.(3)：

p*S*e^{i k_r r} ~ sum_m sqrt(2 pi / k_m) phi_m(z_s) phi_m(z_r) e^{i(k_r-k_m)r - alpha_m r - i pi/4}

再乘 e^{i pi/4}/sqrt(2 pi k_r) 并积分：

g ~ sum_m [phi_s phi_r / sqrt(k_r k_m)] int_{r1}^{r2} e^{i(k_r-k_m)r - alpha_m r} dr

**= Eq.(4) 第一形式。** 闭环 PASS

相位核对：
- steering **e^{+i k_r r}**（Eq.1）
- 传播 **e^{-i k_m r}**（Eq.3）
- 合成 **e^{i(k_r-k_m)r}** -> 峰在 k_r=k_m
- 前因子 e^{i pi/4}*e^{-i pi/4}=1
- sqrt(2 pi / k_m)/sqrt(2 pi k_r)=1/sqrt(k_r k_m)
- S~sqrt(r) 消掉 1/sqrt(r)

## k_r=k_m -> Eq.(5)

int_{r1}^{r2} e^{-alpha_m r} dr = (e^{-alpha_m r1}-e^{-alpha_m r2})/alpha_m

= e^{-alpha_m r0} * 2 sinh(alpha_m dR/2)/alpha_m ,  r0=(r1+r2)/2, dR=r2-r1

=> g(k_m,z_r)= [2 e^{-alpha_m r0} sinh(alpha_m dR/2)/(alpha_m k_m)] phi_s phi_r

= b_m phi_m(z_r),

**b_m = (2 e^{-alpha_m r0})/(alpha_m k_m) * sinh[(alpha_m dR)/2] * phi_m(z_s)**

**= Eq.(5) 印刷形式。** 闭环 PASS

## Eq.(4) 第二形式（印刷）内部问题

印刷：

a_m = [e^{i(k_r-k_m)r2 - alpha_m r2} - e^{i(k_r-k_m)r1 - alpha_m r1}]/(i k_m)

g = sum a_m phi_s phi_r / (k_r - k_m - i alpha_m)

与第一形式逐点比较：

1. 积分闭式分母应为 i(k_r-k_m)-alpha_m = i(k_r-k_m+i alpha_m)，对应 **(k_r-k_m+i alpha_m)**，
   印刷为 **(k_r-k_m-i alpha_m)**（alpha_m 虚部符号相反）。
2. a_m 分子印刷为 r2 项减 r1 项；在 k_r=k_m 处给出 **-b_m**（与 Eq.5 差符号）。
   若分子改为 r1-r2，峰值符号才与 Eq.5 一致。

**结论：**
- Eq3->Eq4(第一形式)->Eq5：**PASS / CLOSED**
- Eq.(4) 第二形式与第一形式：**PAPER_INTERNAL_INCONSISTENT**（疑似排版笔误）
- **实现必须使用第一形式或 Eq.(5)**，不得按第二形式未校正符号直接编码
- 本项 **不是** 我方转录冲突（视觉转录与印刷一致）

## 检查清单映射

| 项 | 结果 |
| --- | --- |
| Eq1 steering 符号 | PASS（e^{+i k_r r}） |
| Eq1 前因子 e^{i pi/4}/sqrt(2 pi k_r) | PASS |
| Eq3 传播相位/衰减 | PASS（e^{-i k_m r - alpha_m r - i pi/4}） |
| S(r) 在积分分子 | PASS |
| sqrt(k_r k_m) 因子 | PASS |
| Eq5 b_m 与峰值积分 | PASS |
| Eq4 第二形式自洽 | FAIL（论文内部，见上） |
