# T2-PREC2 执行包

本包旨在以冻结父提交为基础，开展一次小范围、只读的**原E-STD非均匀SSP实特征值精度可行性实验**。可用结果为：六套N2存档的184条PRT/MOD抽样波数已由本次会话经GitHub连接器预研核查；包内另附可在本地仓库运行的主代码、独立冷审计、机器配置、任务MD及合成数据自测。

使用顺序：解压 → 检验MANIFEST内每个文件实际SHA256 → 将START文本粘贴给本地Codex → 在干净独立工作树核对固定父提交、历史18份输入SHA和软件环境 → Commit A → 一次正式主实验+一次冷审计 → Commit B → STOP。

**禁止**直接将`PRELIMINARY_ESTD_PRT_MOD_SPARSE_AUDIT.csv`作为科学程序输入。历史原始数据留在Git，包内不重复分发MODE/PRT二进制。运行成功后仍不表示非均匀SSP全模态、吸收虚部或CHIGH/海底已认证；不解冻G0/定位MC。R4=0%。

留意：任务MD的哈希由本包MANIFEST和START共同准确登记，严禁用先前其它轮次中的哈希代替。
