"""Read-only N2 reports, complete failures and small raw artifacts; no gate upgrade."""
from t2_n2_common import *
import shutil,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 d=read(OUT/'NUMERICAL_DECISION.json');audit=read(OUT/'INDEPENDENT_COLD_AUDIT.json');design=read(OUT/'EXECUTION_STARTED.json')['design_SHA'];calls=rows(OUT/'PROVIDER_CALLS.csv');full=rows(OUT/'FULL_CATALOGUE_KRAKEN_STABILITY.csv');lines=['# T2-N2研究负责人同步','',f'设计SHA：`{design}`；执行SHA由外层回报，避免本文自引用。',f'任务SHA256：`{TASKSHA}`。', '',f"总体：`{d['overall']}`。", '', '|环节|限定判定|最大空间形状差异|完成／计划|','|---|---|---:|---:|']
 for name,key in [('K窗口','K_WINDOW'),('K网格','K_MESH'),('跨方法','CROSS_METHOD')]:
  r=d[key];v=r['max_difference'];lines.append(f"|{name}|{r['status']}|{float(v):.8%}"+f"|{r['completed']}/{r['planned']}|" if isinstance(v,(int,float)) else f"|{name}|{r['status']}|NOT_EVALUATED|{r['completed']}/{r['planned']}|")
 for name,r in [('FFP网格／谱窗',d['FFP_SELF']['mesh_and_spectral']),('FFP输出插值',d['FFP_SELF']['interpolation'])]:
  v=r['max_difference'];lines.append(f"|{name}|{r['status']}|{float(v):.8%}|{r['completed']}/{r['planned']}|" if isinstance(v,(int,float)) else f"|{name}|{r['status']}|NOT_EVALUATED|{r['completed']}/{r['planned']}|")
 lines+=['',f"FFP自身：`{d['FFP_SELF']['status']}`。只有自身收敛条件成立时，跨方法结果才有数值准入解释；其他情况下保存的差异仅为诊断，不能作为真实误差上界。", '']
 for comp in ['1800_TO_2000','2000_TO_2400','MESH160001_TO240001']:
  rr=[r for r in full if r['comparison']==comp]
  if rr:
   r=max(rr,key=lambda r:float(r['max']));lines.append(f"{comp}最坏：{float(r['max']):.8%}，{r['frequency_hz']} Hz，{r['resource']}/{r['catalogue']}，r={r['worst_r_km']} km、|theta|={r['worst_abs_theta_deg']}°、z={r['worst_z_m']} m。")
  else:lines.append(comp+'：NOT_EVALUATED。')
 lines+=['',f"调用：KRAKEN={d['new_calls']['kraken']}/9，SCOOTER={d['new_calls']['scooter']}/9，FIELDS={d['new_calls']['fields']}/144；没有调用KRAKEN FIELD冒充独立参照。",f"冷审计：PASS={audit['PASS']}，检查行{audit['check_rows']}、FAIL={audit['FAIL']}；K哨兵向量{audit['cold_K_sentinel_vectors']}，FFP向量{audit['cold_FIELDS_resource_vectors']}。零检查项单列NOT_EVALUATED，不称全部通过。", '', '物理：沿用同一5000 m、真空海面、理想刚性海底、CVWT／Thorp选项及常密度SSP；模型与现实海洋不同。K和SCOOTER是不同数值路径，不互称真值。', '窗口1800→2000→2400及FFP2400→3000只提供有限窗口诊断，CHIGH全域完整性未认证；新网格和哨兵对照不恢复原G0。', 'FIELDS使用冻结现代RPO／ReadRcvrRanges参数；正波数、柱面扩散，输出复声压实部／虚部只在实际保存范围轴上作线性插值。请求范围不被当成实际轴，不外推、不逐阵元拟合增益。', '', '原N1优势证书未重算，T2-A/B实际接收与定位性能继续NOT_EVALUATED。新录音、波形、DFT、MC、FIM、模板评分和优化器均0；R4=0%。负责人独立审计待定，STOP。']
 if (OUT/'RUN_FAILURES.json').exists():lines+=['','运行失败及未执行部分：',json.dumps(read(OUT/'RUN_FAILURES.json'),ensure_ascii=False)]
 (OUT/'GPT_SYNC.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
 (OUT/'NUMERIC_SCOPE_AND_FAILURES.md').write_text('# 本轮数值作用域与失败\n\n'+'\n'.join(lines[4:])+'\n',encoding='utf-8',newline='\n')
 # Small modal originals are public; bulky GRN/SHD/field matrices remain locally indexed on D.
 folder=OUT/'modal_archive';folder.mkdir()
 for p in (LOCAL/'K').rglob('*'):
  if p.is_file() and p.suffix.lower() in ['.mod','.prt','.txt']:shutil.copyfile(p,folder/p.name)
 archive=[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p),'git_uploaded_here':bool(p.suffix.lower()=='.mod' and 'K' in p.parts)} for p in LOCAL.rglob('*') if p.is_file()];table(OUT/'LOCAL_ARCHIVE_INDEX.csv',archive)
 dump(OUT/'STOP.json',{'STOP':True,'solver_or_grid_followups_authorized':False,'receiver_performance_followup_authorized':False,'CHIGH_global':'UNCERTIFIED','R4_percent':0,'independent_audit_PASS':audit['PASS'],'raw_GRN_SHD_large_field_storage':'LOCAL_D_ONLY_NOT_CLAIMED_PUBLIC','MOD_copy_public':True,'utc':utc()})
 if full:
  (OUT/'figures').mkdir();fig,axs=plt.subplots(1,3,figsize=(12,4))
  for ax,comp in zip(axs,['1800_TO_2000','2000_TO_2400','MESH160001_TO240001']):
   rr=[r for r in full if r['comparison']==comp]
   if rr:
    val=[100*max(float(r['max']) for r in rr if r['resource']==res) for res in RES];ax.bar(RES,val);ax.axhline(.2,color='red',ls='--');ax.tick_params(axis='x',rotation=30);ax.set_ylabel('Maximum spatial shape difference (%)')
   ax.set_title(comp)
  fig.suptitle('Finite nominal propagation comparison; no receiver performance');fig.tight_layout();fig.savefig(OUT/'figures/K_WINDOW_MESH_DIAGNOSTIC.png',dpi=160);plt.close(fig)
 print('N2 REPORT AND STOP',flush=True)
if __name__=='__main__':main()