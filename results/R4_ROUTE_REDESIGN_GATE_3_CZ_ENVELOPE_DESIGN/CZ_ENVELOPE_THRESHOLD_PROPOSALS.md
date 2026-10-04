# Quantitative Gate threshold proposals（未冻结）

状态 **`PROPOSED_NOT_FROZEN`**，供负责人选择一个 primary tier。表内数值是研究设计建议，不是实验结果，也不是旧 Gate 修改或716验收指标。三档使用同一个已选 representation、nuisance、tau_F、support/search rule；不能每档另调评分阈值。

| 建议档 | joint retention | 每个主案例 W_hull 上限 | C_r=W_hull/15 | range components 上限 | joint broad components 上限 | 两级grid endpoint/width变化上限 |
|---|---|---|---|---|---|---|
| lenient | 100% 有效主案例保留 joint truth/reference，range也保留 | 7.5 km | 0.50 | 4 | 8 | 0.75 km |
| nominal（建议） | 同上，不允许以丢真值换收缩 | 5.0 km | 1/3 | 3 | 6 | 0.50 km |
| strict | 同上 | 3.0 km | 0.20 | 2 | 4 | 0.30 km |

主案例是含噪 bearing、N0 matched-level mechanism cases；无噪声 bearing singleton不计入主成功分母。未来N1须在声明的固定D上单独给相同指标；N0成功不能宣称N1robust。要求逐case联合retention和width，而非平均宽度遮盖失败。invalid/missing样本保留在完整清单，并使整轮不能通过；100%只描述冻结有限panel，不是总体覆盖率或统计置信保证。

理由：lenient要求至少减半，排除trivial 1% shrink；nominal要求15→5 km作为coarse acquisition的明确缩域；strict要求80%缩域但仍只是coarse尺度。以约50km代表值讨论项目“<10%”时，5km只提供一个量级动机，**W_hull不是top-1 error、不是置信区间半宽，也不是五参数验收误差**。不以已知truth distance作逐case阈值归一化。

component caps是有限coarse handoff负担的管理建议，不来自已验证物理别名数。两级稳定性容差暂取对应width cap的10%，只是防数值网格主导的草案，不是已证有足够分辨率。必须同时检查端点双侧移动、width以及预注册adjacency下的broad component拓扑；不因width相同就忽略support整体平移。严格科学范围是外包络和指定分辨率的coarse components，不证明其内部没有细alias。

## 所有档共同的强制条件

1. full-domain observation-only search的保守coverage/bounds可审计；预算不足、bound无效或unresolved主导时不能通过。不得删未知cell后报width。
2. 输出非空，joint truth/reference retained；支持不是只用truth motion的一维curve。
3. 相比同caseB0，envelope确有增量；若B0已窄于所选cap，单靠绝对width不能给CZ信用，标记`NO_CZ_INCREMENT_DEMONSTRATED`。B0退化singleton只regression。
4. 两级raw（非累积保留低预算幸存）outer supports及components稳定；coarse→fine只累计保留不算收敛。
5. source nuisance、z profile与scope明确；符合 [停止规则](CZ_ENVELOPE_STOP_RULES.md)。未要求旧exact-TL的0.001dB Gate在coarse阶段通过。

## 下游与评分阈值的未闭合项

旧RC3 local capture域没有得到可用的全局捕获宽度证书，因此三档 **都不承诺足以启动RC3恢复**。严格档3km也不是毫米matched section的推论。若负责人将来要要求coarse support落入已验证local capture域，应先有独立capture证据；当前不能倒推一个数值。

`tau_F`是observation compatibility band，不能从本case Jmin+Δ定义，不能用truth调至100%retention。noise/numeric allowance、N1 D、grid尺度/bounds/预算尚未冻结；range档选择本身也不授权run。未来主档建议nominal，lenient/strict可作为预先列明的同一run报告，不可失败后降档宣称primary成功。
