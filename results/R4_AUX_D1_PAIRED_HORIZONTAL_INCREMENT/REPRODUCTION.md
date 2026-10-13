# AUX-D1登记顺序

依赖冻结现有Python3.10.11、NumPy1.26.4、SciPy1.15.3；所有输出在D盘。核心单计算进程，BLAS线程固定1（上限4）。不并行诊断，不安装或升级科学依赖。

一次执行：先aux_d1_preflight.py与设计提交A；随后aux_d1_generate.py calibration冻结阈值，aux_d1_generate.py test生成封存测试；aux_d1_execute.py单次推断；封存后aux_d1_evaluate.py、aux_d1_audit.py，再aux_d1_report.py只读报告。各生成／执行阶段有START守卫，存在时禁止重做；评价或报告读取错误仅允许另记纯恢复，不重新计算优化或改变选点。

主预算18000 s，科学累计21600 s；每起点100个残差调用，393216终点上限；RSS8GiB，新增本地2GiB，提交512MiB。科学预算已耗时在SCIENCE_CLOCK.json，各对checkpoint，单次残差守卫时间硬限制。RSS/存储每配对检查；达到限制当前已执行部分全部保存，未执行NOT_EVALUATED。所有终点使用无损float64 NPY，其字段schema随结果提供；CSV为17位有效数字，校准／选点使用原float64。

公开方位与报告轨迹在public_inputs，私有生成扰动在evaluation_only。估计模块无生成器／私有真值导入，U0输入只有主轴。校准器可用校准真值但优化0；测试评价真值必须在INFERENCE_SEAL后读取。冷复核不运行优化器，仅重建输入、Cartesian得分、接受／选择／去重与已有统计。连续全域证书、真实接收能力及R4验收不由本轮产生。B提交推送并核对后STOP。
