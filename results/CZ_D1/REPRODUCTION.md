# CZ-D1只读复算说明

在D盘父快照的独立工作树中运行，Python3.10及numpy1.26.4；缓存/临时目录均指定D盘。脚本仅读历史文件，在本results/CZ_D1写结果。

顺序：`consolidate.py` → `supplement.py` → `build_reports.py` → `audit.py`。不导入历史主实验程序，不执行估计器、传播、接收生成或Monte Carlo。`consolidate.py`读取CAP06 RAW_ARCHIVE_INDEX所登记的D盘路径，若不可访问会明确报告缺档，不能将其视为科学性能失败。`supplement.py`从保存估计值和真值重新计算误差、重建最终R3标签并补充来源。`audit.py`独立读取结果与源CSV，不调用上述研究函数。

源快照固定ca86d324afc1689af6c96cd35145e03625086bae；设计规则提交在DESIGN_FREEZE和GPT_SYNC中记录。为复现字节，保留UTF8和原来源行尾；身份索引记实际检出字节SHA，固定Git快照与文件最后变化SHA另列。最终提交不把自身SHA写进自身文件；以包含报告的提交及推送核验回传为准。

补充脚本应在consolidate重新生成主表后运行一次，避免对已补充主表重复追加；这只是可重复的只读文档构建顺序，不是重跑历史科研。旧单次执行规则和STOP全部保留。后续CZ-D2没有执行授权。
