"""Read-only concise reports and archival status; no new scientific runs."""
from t2_tilt_g0_core import *

def main():
 decision=read(OUT/'DECISION.json');validation=read(OUT/'VALIDATION.json');forward=read(OUT/'FORWARD_DECISION.json');A=read(OUT/'EXECUTION_STARTED.json')['design_SHA'];lines=['# T2-TILT-G0同步报告','',f'设计提交：`{A}`。执行提交由外层回报，避免文档自引用。',f'任务SHA256：`{TASKSHA}`。','',f"前向准入：`{forward['status']}`；传播完整性：`CONDITIONAL_MODAL_WINDOW_LIMIT`。",f"双网格最大去公共源空间形状差异：{forward['max_cross_mesh_profiled_error']:.8%}；FIELD最大差异：{forward['max_FIELD_profiled_error']:.8%}；各冻结门槛0.2%。",'']
 if decision['A']=='NOT_EVALUATED':
  lines+=['A、B距离候选收紧：均为`NOT_EVALUATED`。校准和测试波形、DFT、完整候选推断均未执行。不能把前向准入失败解释为阵型没有距离信息。','',f"失败的目录／频率控制：{forward['failed_catalogue_frequency_controls']}；其他前向控制失败：{forward['failed_other_controls']}。",'全部原始MOD、环境、完整候选副本、形状误差及FIELD结果保留在D盘。没有追加网格、模态删减或阈值变更。']
 else:
  c=decision['primary_geometry_counts'];lines += [f"主A：`{decision['A']}`，{c['A']}/3几何；主B：`{decision['B']}`，{c['B']}/3几何。",f"数值分类：`{decision['numerical_status']}`。",'','主条件固定为z200 m、C0、参考20 dB、主目录、160001网格；主倾角固定+16.60155°。','', '|比较|几何|处理阵有效输出／16|水平标签保留／16|处理阵距离总包络中位km|对照km|信号|','|---|---|---:|---:|---:|---:|---|']
  for r in rows(OUT/'PRIMARY_SIGNAL_BY_GEOMETRY.csv'):
   if r['catalogue']=='primary' and r['mesh']=='160001':lines.append(f"|{r['comparison']}|{r['case_id']}|{r['treated_accepted']}|{r['treated_horizontal_retained']}|{r['treated_median_envelope_km']}|{r['control_median_envelope_km']}|{r['signal']}|")
  lines+=['','全部正负倾角、两档噪声、两种校准、三源深及四种数值目录视图均分开保留。5184配置来自144基础创新、3水平几何；源深、片段、资源和频率不算独立水平案例。','完整候选位图及镜像保留，点估计误差与全集合宽度分别报告。空集无条件误差为INF，有限标签单点不称连续精度。']
 lines+=['',f"新增KRAKEN/FIELD：{forward['new_KRAKEN']}/{forward['new_FIELD']}。冷复核：{validation['checks']}项，FAIL={validation['FAIL']}，接受边界翻转={validation['acceptance_flips']}。",'','本轮为静态全模态空间复谱匹配，不是博士论文WSF全法复现。论文PDF本地未找到；仅使用任务已核对的第三章迁移说明。未迁移海底反射区公式。','单线阵左右舷镜像保留，未知通道增益C1没有被假定已消除。三个频率格点、准确阵型、名义匹配环境均为有利控制条件。','CHIGH未认证、H3-P0/P1源深导数与SSP局部性限制不改变；未产生速度、航向、动态定位或实际海试精度。','', 'R4=0%；研究负责人审计待定；STOP。未授权下一轮。']
 (OUT/'GPT_SYNC.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
 (OUT/'CLIENT_TECHNICAL_UPDATE.md').write_text('# T2-TILT-G0阶段交流（EXPLORATORY）\n\n'+('\n'.join(lines[4:]))+'\n',encoding='utf-8',newline='\n')
 print('reports written; STOP',flush=True)
if __name__=='__main__':main()