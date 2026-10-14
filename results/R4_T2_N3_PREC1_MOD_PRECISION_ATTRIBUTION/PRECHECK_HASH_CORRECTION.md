# PREC1 首次预检哈希冲突更正记录

研究负责人于2026-10-14在当前会话明确确认：使用原ZIP继续；此前聊天中的任务MD哈希属于外部登记错误，不是文件内容或版本变化。不得修改原任务MD、JSON、源代码及科学阈值。

- 错误的外部登记：6d9cabada38585fd00b06eaacdd1477148f1304f254d7fb66ecbce4f72b3bb26。已撤销，不再作为执行依据。
- 唯一有效任务MD SHA256：647a0afcb6546f8d78bbe1ab7c84eaf253f8cd1b711d9c760aa2a10b181cd791。
- 原ZIP SHA256：64e14064a5960575b72f95b0e90aa6b0d2c15356491d5cf483057462cc1f4cff。
- ZIP内任务MD、启动提示词与包内SHA清单使用相同的有效哈希。12份清单文件的长度与SHA、五份原始来源的工作文件和parent Git blob哈希均通过首次只读核验。

首次 STOP 原件：D:/ProjectStorage/WSL-CZ/T2_N3_PREC1/preflight/PREFLIGHT_IDENTITY.json。该原件和当时读取的任务MD保留，不覆盖、不删除。

本次授权只更正外部哈希登记，不豁免任何其它检查。须重新检查全部来源、代码、父提交及干净工作树；出现其它不一致仍立即STOP。正式MOD归因执行和正式冷审计次数在本记录创建时均为0；没有新增声学求解。