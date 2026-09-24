# R3-C2.3A E-STD 迁移（理论上限关）

UTC: see DECISION json

## 条件
S2_KNOWN_STABLE_TONE_UPPER_BOUND; delta=0 = ORACLE_OFFSET_ALIGNMENT
FIELD forward truth; M=9999; c_min=1500; r0=50 km; L=0.6/1.2/2.4 km

## 措辞修正
PEAK_IMPLEMENTATION_NOT_CAUSE -> PEAK_SELECTION_ID_NOT_PRIMARY_CAUSE

## 物理预关
- energetic modes (>1/3/5%) : 数百级
- n_rayleigh_unique @ L=2.4km : **0**
- 谱峰几乎全部 UNRESOLVED_MODE_GROUP（9-22 modes/cell）

## 机制分层
| 层 | 结果 |
| --- | --- |
| EQ5_IDEAL (energetic) | **20/20** z_hat=真值, margin_10m >> 1 |
| TRUE_K / ACTUAL UNIQUE | EMPTY（无 unique mode） |
| Fourier 2.4 km | mode identity 不可唯一确定 |

## 判定
### C2_3A_FOURIER_MODE_IDENTITY_LIMITED

深度签名存在；普通 Fourier 合成孔径在 E-STD 密集模态下无法唯一对应 mode。
**只关闭普通 Fourier mode-ID 能力，不等于 AR/高分辨失败，不关闭 Yang 总路线。**

RANGE_ESCAPE_CHECKED（47.5/55 km）未形成稳定 strict 恢复。

## 下一步（预先固定）
AR / high-resolution modal wavenumber route（R3-C2.3B 的 RC2-δ 耦合暂缓，先解 mode ID）。

YANG_ROUTE_UNDECIDED。
