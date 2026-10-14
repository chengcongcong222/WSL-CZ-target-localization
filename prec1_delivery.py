"""PREC1 report and provenance from sealed outputs only; no extra physical computations."""
from pathlib import Path
import json,csv,hashlib,subprocess,shutil,io
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results/R4_T2_N3_PREC1_MOD_PRECISION_ATTRIBUTION';LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_N3_PREC1')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def main():
 freeze=read(OUT/'DESIGN_FREEZE.json');budget=read(OUT/'EXECUTION_BUDGET.json');seal=read(OUT/'SCIENTIFIC_OUTPUT_SEAL.json');pre=read(OUT/'PRECHECK_REVALIDATION.json');decision=read(OUT/'PREC1_DECISION.json') if (OUT/'PREC1_DECISION.json').exists() else {};cold=read(OUT/'PREC1_COLD_AUDIT.json') if (OUT/'PREC1_COLD_AUDIT.json').exists() else {};allrows=rows(OUT/'PREC1_ALL_40_CASES.csv') if (OUT/'PREC1_ALL_40_CASES.csv').exists() else [];p04=rows(OUT/'PREC1_P04_FROZEN_FOUR_RESOURCE.csv') if (OUT/'PREC1_P04_FROZEN_FOUR_RESOURCE.csv').exists() else []
 for r in freeze['files']:assert sha(r['path'])==r['sha256']
 for r in seal['files']:assert sha(r['path'])==r['sha256']
 for r in pre['source_checks']:assert sha(ROOT/r['file'])==r['sha256']
 supported=budget['failure'] is None and cold.get('PASS',False) and decision.get('status')=='PREC1_K_STORAGE_ROUNDING_MECHANISM_SUPPORTED_CONDITIONALLY'
 finalstatus='PREC1_K_STORAGE_ROUNDING_MECHANISM_SUPPORTED_CONDITIONALLY' if supported else 'PREC1_ATTRIBUTION_INCOMPLETE'
 final={'status':finalstatus,'original_main_decision_preserved':True,'main_rows':len(allrows),'cold_cases':cold.get('cold_cases',0),'cold_checks':cold.get('check_count',0),'cold_FAIL':cold.get('FAIL','NOT_EVALUATED'),'scientific_source_and_result_identities_unchanged':True,'original_N3_GATE':'FAIL_UNCHANGED','N2':'NUMERIC_BRACKETING_NOT_ESTABLISHED_UNCHANGED','G0':'FORWARD_NUMERIC_NOT_ADMITTED_UNCHANGED','T2_A_B_receiver_performance':'NOT_EVALUATED','E_STD_variable_SSP_high_precision_k':'NOT_EVALUATED','CHIGH_and_real_seabed':'NOT_CERTIFIED','R4_percent':0,'research_lead_audit':'PENDING','STOP':True,'next_stage_authorized':False};write(OUT/'FINAL_ADJUDICATION.json',final)
 designSHA=read(OUT/'FORMAL_PIPELINE_STARTED.json')['design_SHA']
 provenance={'parent_SHA':freeze['parent_SHA'],'design_SHA':designSHA,'task_MD_SHA256':pre['task_MD_sha256'],'ZIP_SHA256':pre['ZIP_sha256'],'original_main_program_SHA256':sha(ROOT/'t2_n3_prec1.py'),'original_cold_program_SHA256':sha(ROOT/'t2_n3_prec1_cold_audit.py'),'source_records':pre['source_checks'],'first_STOP_SHA256':pre['initial_STOP_sha256'],'first_STOP_preserved':sha(LOCAL/'preflight/PREFLIGHT_IDENTITY.json')==pre['initial_STOP_sha256'],'external_hash_registration_corrected_by_user':True,'reference_CSV_and_PNG_in_TASK_PACKAGE':'PREANALYSIS_ONLY_NEVER_MAIN_INPUT','radial_kernel':'FARFIELD_SAME_FOR_C0_C1_C2_C3; old N3 Hankel history compared with preregistered2e-6 tolerance','mode_order':'1787 elementwise exact float32 analytic-rounding comparisons verify stored order; cold roots independent','zero_C3_error':'REFERENCE_SELF_COMPARISON_BY_DEFINITION_NOT_LOCALIZATION_ZERO_ERROR','cold_scope':'ACTUAL_GIT_ARCHIVED_MOD; synthetic selftest separate','formal_budget':budget};write(OUT/'PROVENANCE.json',provenance)
 rawlogs=OUT/'runtime_logs';rawlogs.mkdir(exist_ok=True)
 for p in (LOCAL/'logs').glob('*'):
  if p.is_file():shutil.copyfile(p,rawlogs/p.name)
 lines=['# T2-N3-PREC1 同步报告','',f'设计SHA：`{designSHA}`。执行SHA由外层回传，避免自引用。',f'父提交：`{freeze["parent_SHA"]}`。任务MD SHA256：`{pre["task_MD_sha256"]}`。','',f'最终限定判定：`{finalstatus}`。','', '负责人已更正外部任务哈希登记；原ZIP、MD、JSON、两份科学源代码和全部阈值保持原字节。首次预检STOP原件与本次纠正记录随提交保留。干净独立工作树执行，原工作区11个未跟踪文件不进入研究输入。','', '本轮只读两份真实Git封存MOD：201 Hz、5000 m、1500 m/s等声速、自由海面和理想刚性海底，CHIGH=1800与2400。C0使用存档波数与模态函数；C1仅换解析双精度波数；C2仅换解析源/收模态函数；C3为解析远场参考。四臂径向核一致，不与原N3精确Hankel基准混用。','', '|P04、CHIGH=2400|C0 原存档|C1 仅换波数|C2 仅换模态函数|C3 自参考|','|---|---:|---:|---:|---:|']
 for r in p04:lines.append(f'|{r["resource"]}|{100*float(r["C0_stored_k_stored_phi_vs_C3"]):.6f}%|{100*float(r["C1_analytic_k_stored_phi_vs_C3"]):.6f}%|{100*float(r["C2_stored_k_analytic_phi_vs_C3"]):.6f}%|0（定义）|')
 if allrows:
  c1max=max(float(r['C1_analytic_k_stored_phi_vs_C3']) for r in allrows);c2fails=sum(r['C2_below_0p002']=='False' for r in allrows);c0fails=sum(r['C0_below_0p002']=='False' for r in allrows)
  lines += ['',f'完整40个配对数值对照：C0超历史0.2%门槛{c0fails}/40；C1低于门槛{sum(r["C1_below_0p002"]=="True" for r in allrows)}/40，最大{100*c1max:.6f}%；C2仍超限{c2fails}/40。40行不是40个独立物理目标试验。',f'1787/1787返回波数与解析双精度波数转float32后的表示逐项完全相同。原N3场指标重建最大差{decision["max_history_delta"]:.9g}，低于冻结2e-6容限。','']
 lines += [f'实际MOD独立冷审计：{cold.get("check_count",0)}项、{cold.get("FAIL","NOT_EVALUATED")} FAIL、{cold.get("cold_cases",0)}/40案例，最大形状复核差{cold.get("max_shape_delta","NOT_EVALUATED")}。独立struct解析、brentq根、逐阵元cmath累加和QR投影；未用合成自测替代正式审核。','', '作用域：这支持理想波导内单精度波数存储舍入是主要可复现超限机制；不把改善量除以总误差称独立贡献百分比，也不证明未知SSP中已有解析双精度波数可直接替换。原N3的0.2%硬门仍FAIL，不重新准入。N2巨大窗口依赖与海底合同问题独立保留。','',f'正式主程序1次、冷审计1次；新增KRAKEN/SCOOTER/FIELD/FIELDS、MC、接收波形、DFT和后验候选反演均0。科学执行及冷复核壁钟{budget["wall_s"]:.3f} s，进程树峰RSS {budget["peak_tree_RSS_bytes"]/1024**2:.2f} MiB；冻结上限600 s、2 GiB、1 GiB新结果。','', '下一研究决策只作建议：单独设计非均匀SSP水平波数的高精度输出或存储误差控制，并继续保留模态窗口与真实底质合同缺口。本轮不修改提供器或增加传播求解，不把数值归因称为阵型性能提升。','', 'T2-A/B实际接收与定位性能NOT_EVALUATED；R4=0%。研究负责人审计PENDING。STOP。']
 (OUT/'GPT_SYNC_PREC1.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
 if p04:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;import numpy as np
  fig,ax=plt.subplots(figsize=(10,5));x=np.arange(len(p04));width=.2
  fields=['C0_stored_k_stored_phi_vs_C3','C1_analytic_k_stored_phi_vs_C3','C2_stored_k_analytic_phi_vs_C3','C3_analytic_k_analytic_phi_vs_C3'];labels=['C0 stored k + stored phi','C1 analytic k + stored phi','C2 stored k + analytic phi','C3 reference self (0 by definition)']
  for i,(field,label) in enumerate(zip(fields,labels)):ax.bar(x+(i-1.5)*width,[100*float(r[field]) for r in p04],width,label=label)
  ax.axhline(.2,ls='--',color='black',label='Historical 0.2% gate (unchanged)');ax.set_xticks(x,[r['resource'] for r in p04]);ax.set(ylabel='Array spatial shape discrepancy (%)',title='P04, 201 Hz, CHIGH=2400: ideal isovelocity numeric attribution');ax.legend(fontsize=8);fig.text(.5,.01,'Archived MOD comparison; not localization accuracy, actual ocean certification, or N3 readmission.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.045,1,1]);fig.savefig(OUT/'P04_precision_attribution.png',dpi=170);plt.close(fig)
 ledger=ROOT/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv';original=(LOCAL/'WORKTREE_LEDGER_ORIGINAL.csv').read_bytes();assert ledger.read_bytes()==original
 buf=io.StringIO();csv.writer(buf,lineterminator='\n').writerow(['T2-N3-PREC1','../R4_T2_N3_PREC1_MOD_PRECISION_ATTRIBUTION/FINAL_ADJUDICATION.json',finalstatus,f'40_ARCHIVED_MOD_FOUR_ARM_ROWS;1787_STORAGE_CHECKS;NO_NEW_SOLVERS_RECEIVER_MC;N3_GATE_UNCHANGED;R4_delta=0',f'ACTUAL_MOD_COLD_{cold.get("check_count",0)}_CHECKS_{cold.get("FAIL","NOT_EVALUATED")}_FAIL;RESEARCH_LEAD_PENDING;STOP',0]);ledger.write_bytes(original+(b'' if original.endswith(b'\n') else b'\n')+buf.getvalue().encode('utf-8'))
 write(OUT/'LEDGER_APPEND_PROVENANCE.json',{'working_original_sha256':hashlib.sha256(original).hexdigest(),'new_rows':1,'working_prefix_preserved':ledger.read_bytes().startswith(original),'R4_delta':0})
 files=[p for p in OUT.rglob('*') if p.is_file()];write(OUT/'OUTPUT_MANIFEST.json',{'scope':'snapshot excludes own hash to avoid self-reference','files':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size} for p in files]})
 print('DELIVERY',finalstatus,'STOP',flush=True)
if __name__=='__main__':main()