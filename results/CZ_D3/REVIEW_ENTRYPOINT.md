# CZ-D3内部材料审阅入口

状态：INTERNAL_DRAFT_NOT_APPROVED_FOR_EXTERNAL_RELEASE。仅用于研究负责人审阅，不向甲方发送。

1. [内部技术报告](CZ_D3_INTERNAL_TECHNICAL_REPORT.md)（[Word可编辑版](CZ_D3_INTERNAL_TECHNICAL_REPORT.docx)）
2. [15页PPT草稿](CZ_D3_INTERNAL_BRIEFING_DRAFT.pptx)（[PDF预览](CZ_D3_INTERNAL_BRIEFING_DRAFT.pdf)）
3. [精度—观测条件表](PRECISION_OBSERVATION_CONDITIONS.csv)及[完整输入条件](OBSERVATION_CONDITIONS.json)
4. [架构与资源矩阵](ARCHITECTURE_AND_RESOURCE_MATRIX.csv)、[实际部署范围](DEPLOYMENT_GEOMETRY.csv)
5. [逐条措辞证据](CLAIM_EVIDENCE_LEDGER.csv)、[来源SHA索引](SOURCE_INDEX.csv)、[封版审阅清单](RELEASE_REVIEW_CHECKLIST.md)
6. [简短同步报告](GPT_SYNC_CZ_D3.md)

四组图件位于figures，分别提供PNG、PDF、SVG。逐侧、逐场景分位数及配对表保留完整分母。PPT每页来源索引与讲解限制见PPT_SLIDE_EVIDENCE.csv和PPT_SPEAKER_NOTES.md。

CZ-D2原完整性Gate仍FAIL：539,763 PASS、1 FAIL。记录未观察到成功私有读取，不构成完整隔离认证；没有整体科研PASS。所有CZ-D2数字仅作为存在完整性例外的条件性仿真结果。单HLA主线、深度NOT_EVALUATED与R4=0%保持。

本轮仅重排和核对保存表格。TABLE_CHECKS与EDITORIAL_VALIDATION是文档检查，不是新的独立科学冷审计。没有重跑传播、估计、Fisher、Monte Carlo或历史Gate。

编制源码build_internal_delivery.py可在安装numpy/matplotlib/python-pptx/python-docx/Pillow的Python环境运行。脚本读取本仓库固定来源，不读取D盘科学原始数据或调用旧科研程序。Windows中文绘图字体默认为微软雅黑；字体与依赖路径为本地文档环境配置，不是科学输入。

PPT已由本地PowerPoint后台实际渲染。Word可编辑稿已检查OOXML结构；后台Word转PDF未完成，未交付技术报告PDF，也不宣称Word已完成逐页渲染检查。PPT PDF和四份科研图PDF均已交付。

提交后STOP，等待负责人对正式对外措辞和例外披露审核。没有外发授权。