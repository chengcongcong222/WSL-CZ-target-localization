# T2-N1重现与单次边界

Windows PowerShell先载入D:\ProjectStorage\WSL-CZ\runtime-environment.ps1；Python -B，PYTHONPATH=D:\ProjectStorage\WSL-CZ\HLA_H1\runtime_packages，BLAS线程固定1。Python3.10.11、NumPy1.26.4、SciPy1.15.3与历史D盘psutil；不安装或升级依赖。

设计准备只检查包哈希、36复場及18旧误差数组的哈希／头信息和目录精确同一性，十个小确定性向量公式检查，不提前运行新评分。A后t2_n1_execute.py只启动唯一一次t2_n1_analyze.py，随后独立t2_n1_audit.py及只读报告。启动守卫已存在时不再执行。所有复场np.load mmap_mode=r、allow_pickle=False，旧档案不写入。大型新中间量留D盘T2_N1，仓库只交摘要、精确输入锁和核验结果。

单数值计算进程，最多1800秒、4GiB RSS、2GiB新增D盘、64MiB新增Git。来源身份、公式或预算异常均STOP，不重跑科学、不修改A冻结核心。B只加结果文件；任何所需下一试验必须另行授权。

描述统计中位数采用通常中位数，P90/P95/P99按最近秩。范数等计数十分位稳定按候选ID破同分；距离、绝对方位及源深使用固定边界分层。全部分母保留，不把频率、候选、目录或五个模板当成独立接收几何。