from pathlib import Path
import subprocess,json,csv
R=Path.cwd();O=R/'results/ENV_CAP_FULL_CAMPAIGN';g=lambda *args:subprocess.check_output(['git','-c',f'safe.directory={R.as_posix()}',*args]).decode().strip()
shas={a:g('rev-parse',a) for a in ['bc4a1fe','9e5db58','4dc7429']}
sync=f'''
# ENV-CAP最终同步

结论：**NO_GO_CURRENT_CZ_METHOD_PERFORMANCE_DEMO**。E0—E5按预注册负结果分支收口；未建立新的第一CZ Yang/MMAR深度或五维定位精度。

|阶段|交付与科学状态|
|---|---|
|E0|任务32/32哈希、历史28项、原始T/S与ETOPO点16项独立核对；合同冻结|
|E1|15个登记单元；主程序9个临时CZ2标签／1来源族；冷审计274项、12FAIL；正式环境准入暂缓，已准入CZ2=0|
|E2|浅海存档12项评分一致；总控制25/26PASS，int32平方溢出导致镜像控制FAIL；不重跑|
|E3|30项方法—环境单元均NOT_EVALUATED；新增MC／接收声压0|
|E4|解析孔径、频率及offset设计尺度＋历史T3有限资源图；新CZ可行性格0；只读冷核对126项0FAIL|
|E5|完整科研报告、7组PNG/PDF、来源和原始SHA、全部失败与未评估表；最终推送后STOP|

17N118E条件重建的临时峰位48.375–51.375km。低能量滑动和／dB地板导致冷重建差异，SCS年平均200m最大差出现在11.375km全零窗，焦区内差约2.84e-14dB；这些诊断不能改判原FAIL。八个SCS单元逐项一致仍作为临时数据保留，不包装成完整环境族认证。E-STD三例亦未通过本轮数值门，常速三例没有登记折射回返。

E2镜像FAIL为新代码整数溢出，不是物理不对称；原始代码、NaN与FAIL完整保留。即使无此错误，新CZ相干完整模态提供器和非oracle提取仍缺失。Yang浅海／MMAR历史条件成果不变，不能借给新CZ精度。

资源尺度例：普通Fourier间隔0.001rad/m对应径向孔径约6.283km；不是所有高分辨算法的下界。历史T3的15.297km、Q1速度最坏场景P95约0.1517m/s、距离约0.7026%可作条件参照；它依赖有符号方位、导航和Uniform偏置支持，不是新CZ提取性能或最低装备。

实际记录科研墙钟约321.94s，63次BELLHOP（含3次前置控制），KRAKEN/KRAKENC/FIELD/FIELDS/SCOOTER均0；新MC0、接收声压0、科学重试0。观测到的峰值约369.05MiB，D盘原始归档约1.03GiB。RSS采样和共享缓存统计作用域见BUDGET_USAGE，不宣称超出测量范围的完整资源证书。

父提交：a9a80191ee4d0f0cecd9db296198606b63a68080。
设计Commit-0：{shas['bc4a1fe']}。
CP1：{shas['9e5db58']}。
CP2：{shas['4dc7429']}。
最终执行SHA：包含本报告的最终提交，推送后由外层回报完整SHA，并写D:/ProjectStorage/WSL-CZ/ENV_CAP_FINAL_GIT_VERIFICATION.json，避免文件自引用。
分支：research/env-cap-20261016；不合并main。main保持父SHA。原T2/FUNC2和全部历史Gate保持，R4=0%。

整体审阅入口：[REVIEW_ENTRYPOINT.md](REVIEW_ENTRYPOINT.md)。提交后STOP，不自动恢复失败单元或增加计算；本轮不请求合同扩展。
'''
(O/'GPT_SYNC_ENV_CAP_FINAL.md').write_text(sync.strip()+'\n',encoding='utf-8',newline='\n')
(O/'REVIEW_ENTRYPOINT.md').write_text('''# ENV-CAP整体审阅入口

建议先读[最终同步](GPT_SYNC_ENV_CAP_FINAL.md)与[最终科研判定](ENV_CAP_FINAL_DECISION.json)，再读[完整科研报告](ENV_CAP_FINAL_SCIENTIFIC_REPORT.md)。本轮是负结果及解析资源分支交付，不是新CZ精度PASS。

关键核查：

1. [E1最终裁定](E1_FINAL_ADJUDICATION.json)与[冷审计](E1_INDEPENDENT_COLD_AUDIT.json)：274项中12FAIL，原9个临时标签保留但正式准入为0；[差异诊断](E1_COLD_DISCREPANCY_DIAGNOSTIC.csv)不修复FAIL。
2. [E2最终裁定](E2_FINAL_ADJUDICATION.json)：整数溢出控制失败；[源与公式身份](E2_PAPER_METHOD_IDENTITY.md)和[勘误](E2_YANG_ERRATUM_FORMULA_AUDIT.md)不借oracle条件。
3. [性能缺失表](E3_PERFORMANCE_BY_CASE.csv)：30单元、零尝试，不填零误差；[资源表](E4_RESOURCE_INVERSION.csv)不将历史分位档当ENV-CAP可行性。
4. [申请可用文字](APPLICATION_USABLE_CLAIMS.md)及[禁用表述](APPLICATION_FORBIDDEN_CLAIMS.md)，检查有没有将研究问题写成已完成能力。
5. [审计范围汇总](INDEPENDENT_AUDIT_SUMMARY.json)、[来源冷核对](E0_SOURCE_INDEPENDENT_COLD_AUDIT.json)、[E2/E4只读冷核对](E2_E4_INDEPENDENT_COLD_AUDIT.json)、[预算](BUDGET_USAGE.json)和[历史台账未改](HISTORICAL_LEDGER_PRESERVATION.json)。只读PASS不改变科学FAIL。

图表入口：[七组科研图](figures/)，PNG供浏览、PDF供排版；[图数据血缘](FIGURE_DATA_PROVENANCE.csv)与[大文件原始索引](raw_local_archive_index.csv)逐一绑定D盘输入。新科学与冷审计代码在scripts/，任务包原字节在task_package/。

科研分支research/env-cap-20261016，最终SHA由外层回报和D盘ENV_CAP_FINAL_GIT_VERIFICATION.json核对。main与R4保持原样，STOP=true。最终审核应判断失败作用域、条件证据和下一研究必要性，不要求当前失败被包装为性能通过。
''',encoding='utf-8',newline='\n')
(O/'SESSION_HANDOFF.md').write_text(f'''# ENV-CAP最终恢复状态：STOP

分支research/env-cap-20261016；父SHA a9a80191ee4d0f0cecd9db296198606b63a68080。设计{shas['bc4a1fe']}、CP1 {shas['9e5db58']}、CP2 {shas['4dc7429']}。最终执行提交为包含本文件的最终提交，远端核对记录单独存D:/ProjectStorage/WSL-CZ/ENV_CAP_FINAL_GIT_VERIFICATION.json。

E1一次主执行、一次冷审计均已完成，E2一次控制、E4一次资源分析、其联合一次冷审计及E0来源一次冷核对均已完成。不得重复所有SCIENCE_STARTED所标单元。E1冷FAIL12、E2实现FAIL1，原始输出不覆盖；E3未评价。原始目录D:/ProjectStorage/WSL-CZ/ENV_CAP，完整索引raw_local_archive_index.csv。

最终合法动作仅为审阅与读取本报告；所有新实验、代码修复后的科学重跑、原T2/FUNC2恢复均需新授权／冻结。旧历史目录和R4总台账仍只读。R4=0%。没有待批准的合同扩展请求。

完整报告、七组图、资源需求及失败表已经交付。下一会话先读RUN_STATE、STOP.json和REVIEW_ENTRYPOINT，不启动科学程序。
''',encoding='utf-8',newline='\n')
with (O/'CHECKPOINT_STATUS.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,['checkpoint','commit','remote_branch','science_scope','requires_approval','STOP']);w.writeheader();w.writerows([{'checkpoint':'CP1','commit':shas['9e5db58'],'remote_branch':'research/env-cap-20261016','science_scope':'E1auditfailurepreserved','requires_approval':False,'STOP':False},{'checkpoint':'CP2','commit':shas['4dc7429'],'remote_branch':'research/env-cap-20261016','science_scope':'E2implementationfailureandE3mask','requires_approval':False,'STOP':False},{'checkpoint':'FINAL','commit':'COMMIT_CONTAINING_THIS_FILE','remote_branch':'research/env-cap-20261016','science_scope':'Fullnegativecloseoutandresourceconditions','requires_approval':False,'STOP':True}])
print('Final sync and review/handoff entrypoints prepared.')
