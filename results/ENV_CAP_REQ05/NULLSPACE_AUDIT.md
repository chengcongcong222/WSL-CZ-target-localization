# 小型代数审计

18组矩阵全部是预登记的无量纲反例/条件示例，绝不是六种真实观测的五维Fisher测量值。精确有理数消元独立核对信息矩阵、秩和零空间，124项通过，0失败。独立实现由同一代理编写，未冒称外部科研审核。

|控制|有效秩|不可辨维数|
|---|---:|---:|
|LEVEL_CONSTANT|3|2|
|LEVEL_LINEAR|2|3|
|LEVEL_FREE|0|5|
|LEVEL_FLAT|0|5|
|RELATIVE_KNOWN_GAIN|2|3|
|RELATIVE_CONSTANT_GAIN|2|3|
|RELATIVE_FREE_GAIN|0|5|
|RELATIVE_STATIC_GAIN|0|5|
|RELATIVE_COUPLED|1|4|
|PHASE_PLANE_WAVE|1|4|
|PHASE_FREE_CHANNEL|0|5|
|DELAY_COMMON_TIME|2|3|
|DELAY_FREE_PATH|0|5|
|FREQUENCY_UNKNOWN_OFFSET|1|4|
|FREQUENCY_UNKNOWN_DRIFT|0|5|
|APERTURE_FREE_SOURCE_PHASE|0|5|
|AUX_HORIZONTAL_UPPER_BOUND|4|1|
|AUX_FREE_EPOCH_BIAS|0|5|

完整数值矩阵在ALGEBRA_RESULTS.json，精确零向量在EXACT_NULLSPACES.json，输入在ALGEBRA_SPEC.json。AUX_HORIZONTAL_UPPER_BOUND仅用四个抽象独立水平导数说明z仍为零，不表示实际辅助轨迹已取得四维可靠秩。

共享参考D=[[-1,1,0],[-1,0,1]]作用于单位原始协方差时，DD^T=[[2,1],[1,2]]；忽略非对角项会改变信息。主程序正确协方差下的对比信息与原始空间公共源剖面相等，损失矩阵由有理数验证半正定。平面波和运动方位解析导数与独立差分交叉核对。

这里没有由局部秩推断全域唯一性，没有真实环境Jacobian、实际定位CRLB或接收提取精度；镜像、路径/模式置换、离散别名仍按观测合同保留。
