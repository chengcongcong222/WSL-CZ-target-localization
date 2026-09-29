# RANGE-UB-1B-FIX

**SINGLE_DEPTH_ORACLE_TDOA_RANGE_INFORMATION_CONFIRMED**

修复: arr parser col2=delay, 651点全网格, path-metadata association, rigid bottom.

结果: 真值 rank 全 1, n_valid=651, 别名 J_tau > 1000 ms。
距离脊线被打破, 深度也正确恢复 (best_z_star = z_true)。

对比 TL: r-z 脊线 J<0.05 dB vs TDOA J>1000 ms。

旧 1B: RANGE_UB_1B_INVALIDATED_BY_ARRIVAL_PARSER_AND_INCOMPLETE_Z_PROFILE
