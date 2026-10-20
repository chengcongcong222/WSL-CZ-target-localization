# 恢复执行接口与代码差分说明

原设计dc91982e、来源修订eab37fe及所有旧STOP保持只读。EFFECTIVE_INPUT_BINDING.json绑定原合同与V3完整路径/SHA；原V1文件不覆盖、不在新面板排除中使用。

生产估计器直接加载原t3_w1_model.py的Model/initialize/estimate，原文件没有复制修订。适配仅quality副本common_bound_deg、half_diff_bound_deg在P1乘2。没有额外起点、warm start、替代损失或重试。P0/P1独立调用同一原函数。生成器调用原t3_w1_generate.observations，以冻结的新种子、面板与成对创新为参数。新内容限面板身份、归档隔离、预算、扩展边界诊断和既定统计。

执行分层：main_run生成私有面板/创新与公开数组后，启动独立estimate_public进程，不传真值。该进程的Python audit hook只允许读取公开观测、冻结源码与必要运行库，写入终点/进度/异常记录；实际路径留档。这是应用层读路径强制许可及运行记录，不宣称操作系统级安全沙箱。生产代码及源库版本先封存，随后才物化面板。

私有归档：D:/ProjectStorage/WSL-CZ/CZ_D2/RESUME_01_DATA/private。
公开归档：同根public；48个npz，各500条记录，数组times(121)、reported(500,121,2,2)、bearing(500,121,2)，float64，米/秒/弧度。
估计器唯一量测输入为这些公开数组和质量/支持合同；不读取面板、创新、真实节点位置、偏置或评价标签。

终点：terminals.npy，48000个固定槽位。顺序Q1/Q2、B0/B3、D01..D12、rep0..499、P0/P1；executed=0表示未执行，不能当作失败或零误差。h/initial为km、deg、m/s、deg、deg、deg；偏置距离单位deg。每个槽位保存原solver_status、nfev、cost、optimality、初值回退、交点数、拟合Jacobian秩/比、残差、6维触界及偏置支持距离。

面板算法、主节点轨迹、随机数顺序、统计和预算继承原DESIGN_FREEZE。双侧辅助布放不反射目标。24000公开记录对48000方法槽位；配对不是48000独立样本。

前置纯合成控制不运行优化、不取新面板、不生成随机创新；独立审计模块不导入主生成/拟合/统计/信息函数。正式主/冷各一把新exclusive-create锁；旧锁不触碰。预算看门程序在独立父进程监视进程树并终止超限，记录采样RSS及累计存储；0.2s轮询不构成连续峰值保证。

本轮一旦冻结生产代码，任何前置身份/公式/面板不一致或正式科学完整性异常均保留数据并停止，不修复后续跑。封存与必要的已完成部分唯一审计仍按原合同处理。
