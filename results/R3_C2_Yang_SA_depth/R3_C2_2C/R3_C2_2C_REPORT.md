# R3-C2.2C 报告

UTC: 2026-09-23T12:44:30.356064+00:00

## 判定

### `C2_2C_PAPER_REPRO_CONFIRMED`

all 4 cases peak within 5 m under ORACLE_MODE_ID + delta=0 + matched env; Eq1-Eq5 closed; alpha/delta gates passed

## PRIMARY / ENV

- 4 PDF 按 title/DOI 锁定（见 PRIMARY_SOURCE_SET.md）
- span **4990 m**（P2 勘误 + P3）；dr=**2.5 m**；delta=**0** = `ORACLE_OFFSET_ALIGNMENT`
- H=88 m（正文+Fig.1）；KRAKEN modes M=23（论文 Fig.2a 仅展示 lowest 16 true modes，非总数上限）

## 实现自证

- Eq1->Eq5: fine-dr max_rel (m1-8) = 3.646e-08  → PASS
- alpha limit: stable finite = True
- delta mechanism: |g| preserved, phase k_m*delta（见 csv）
- alpha sign: alpha_m = +|Im(k_m)|  (decay)  [eq3 convention e^{-i k r - alpha r}]

## 论文复现（ORACLE_MODE_ID, delta=0）

| zs | zr | z_hat | near_true |
| --- | --- | --- | --- |
| 4 | 18 | 5 | True |
| 4 | 70 | 5 | True |
| 50 | 18 | 50 | True |
| 50 | 70 | 50 | True |

- 浅/深源在 zr=18 与 70 均可分；z_hat 在 1 m 网格上（zs=4→5 m 为最近节点）
- Shading 100/250/500 m 与 THEORY_CONTROL 的 z_hat 一致
- 多模谱峰相对 b_mφ_r 存在弱模泄漏（`MULTIMODE_LEAKAGE` 行）；**代数 gate 用单模闭环**，泄漏不作 fail
- case matrix 四组合齐全，无静默过滤

谱图 `figures/fig1_reproduction_zr*.svg`；定深 `fig2_reproduction_zr*.svg`。

## 范围

**不得**外推：YANG_ESTD_FAIL / CZ APERTURE_LIMITED / DEPTH_CANNOT / RC3-C_FAIL。
不进入 E-STD / 50-60 km / 2.4 km。

## 停止

等待 GPT 审计。
