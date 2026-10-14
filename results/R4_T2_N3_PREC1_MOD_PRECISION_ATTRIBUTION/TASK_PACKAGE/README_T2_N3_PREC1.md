# T2-N3-PREC1 交付索引

当前结果：直接读取提交 `8ee2bd53fe1b7a50725833e6c7b5c5dea35fbe23` 的两份真实201Hz MOD开展了P01–P05、四种阵型、两窗口的四臂数值替换；结果重点见 `T2_N3_PREC1_RESEARCH_NOTE_20261014.md` 和 `CURRENT_PREC1_P04_FOUR_ARM_COMPARISON.csv`。本会话另外运行了独立纯解析C2的40行核对，P04四项与真实MOD读出的C2在`3e-11`内相符；尚未在容器中对真实MOD执行Python正式冷审计。

推荐流程：将此包解压，先粘贴 `T2_N3_PREC1_START_20261014_v1.txt` 给仓库根目录工作的本地Codex。不要把当前目录里的P04结果CSV输入给科学程序，科学程序自己读取仓库中的原始二进制 `.mod`。旧G0/N1/N2/N3均继续冻结，R4仍为0%。

产物身份：`T2_N3_PREC1_CODEX_TASK_20261014_v1.md` 是完整执行合同；`T2_N3_PREC1_EXPERIMENT_SPEC_20261014_v1.json` 是机器冻结配置；`t2_n3_prec1.py` 是分析程序；`t2_n3_prec1_cold_audit.py` 是独立冷审计；`test_synthetic_only.py` 是**仅测试代码**的人工生成MOD单元测试；`T2_N3_PREC1_RESEARCH_NOTE_20261014.md` 是新结果范围说明。

本包不包含两个原始 MOD 二进制，原因是它们已经保存在上述 GitHub 仓库提交中，且执行需在仓库中按SHA256核对和读取原文件。
