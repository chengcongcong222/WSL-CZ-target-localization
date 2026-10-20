# 冷审计失败作用域与封存说明

正式主分析完成48,000项配置。唯一一次冷审计为539,763项PASS、1项FAIL，随后STOP；未重试，未修改冻结科学代码或改判。

失败检查为 `private_read_isolation`。冻结检查要求被拒绝读取列表为空，同时要求成功读取日志中没有private/panel路径。实际拒绝了两项Python标准库字节码缓存路径：

- `D:\Caches\python\Users\ccc\AppData\Local\Programs\Python\Python310\lib\zipfile.cpython-310.pyc`
- `D:\Caches\python\Users\ccc\AppData\Local\Programs\Python\Python310\lib\encodings\cp437.cpython-310.pyc`

成功读取日志包含对应标准库源文件，没有成功读取私有真值或面板的路径。这与缓存探测受阻后回退到源文件导入一致。拒绝记录仅保存路径，未保存完整open模式及flags，因此不能据此升级为新的完整隔离认证。

本次审计实现把任何被拒绝的读取均判为失败，未区分标准库缓存与私有数据。30项预检验证了私有路径拦截，但没有覆盖运行时NPZ导入的缓存访问。这属于本次审计实现的不足。保留原始FAIL，不在看到结果后修正并重跑。

最终保持 `integrity_pass=false`、`IMPLEMENTATION_OR_NUMERIC_GATE_FAIL`。报告中的数值仅为封存的描述性结果，不能作为已通过完整科学准入的性能结论。全部数值检查通过也不抵消此项完整性失败。

本说明仅解释已保存日志，不是第二次冷审计。历史Gate、main、R4=0%不变，不启动CZ-D3。详细日志身份及限制见同名JSON。