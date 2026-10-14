# PREC3：源码与二进制来源审查

|证据对象|结论|范围|
|---|---|---|
|公开源码内部运算|SOURCE_CODE_INDICATES|显式 KIND=8 的本征值、外推和复波数变量|
|公开源码 MOD 写盘|SOURCE_CODE_INDICATES|波数和模态函数经默认或 KIND=4 转换|
|本地 Windows 构建对应关系|NOT_ESTABLISHED|有二进制身份，没有源码、编译参数或构建日志|
|旧档案表示关系|PREC2_ACCEPTED|184 条稀疏打印值转 float32 与原 MOD 实部一致|
|全部内部双精度输出能力|NOT_ESTABLISHED|本轮不运行、重编译或修改 KRAKEN/KRAKENC|

本轮源码固定到上游提交 `8b4682bd719d5d17f5a1c798039d6be8235e3859`，完整字节保存在 D 盘，见 OFFICIAL_SOURCE_INDEX.json。公开源码说明内部变量和写盘路径；它不建立旧 Windows 程序与该版本的构建等价关系。MOD 的实际 complex64 布局来自封存二进制数据核查。

[KrakenMod](https://raw.githubusercontent.com/oalib-acoustics/Acoustics-Toolbox/8b4682bd719d5d17f5a1c798039d6be8235e3859/Kraken/KrakenMod.f90)、[KRAKEN 主程序](https://raw.githubusercontent.com/oalib-acoustics/Acoustics-Toolbox/8b4682bd719d5d17f5a1c798039d6be8235e3859/Kraken/kraken.f90)。

吸收合同：原 ENV 为 CVWT，输入吸收列为 0，T 增加 Thorp。公开 CRCI 把 dB/km 换成 Neper/m，再在 SSP 节点形成正虚部复声速；C-linear 插值作用于该复声速。原传播核采用 exp(-ikR)，负虚部 k 对应随距离衰减。本轮冻结相同符号、常密度、水深和 V/R 边界；不把无损实根与原 MOD 虚部直接拼接。原 Thorp 数据并非真实海域实测参数。

[吸收转换](https://raw.githubusercontent.com/oalib-acoustics/Acoustics-Toolbox/8b4682bd719d5d17f5a1c798039d6be8235e3859/misc/AttenMod.f90)、[SSP 插值](https://raw.githubusercontent.com/oalib-acoustics/Acoustics-Toolbox/8b4682bd719d5d17f5a1c798039d6be8235e3859/misc/sspMod.f90)、[官方 KRAKEN 手册](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html)。

虚部计算必须区分：KRAKEN 的模态损耗扰动与直接复算子求解。本轮新计算来自独立 Python 小样本，并非旧程序高精度全模态侧车。两条小样本路径共享已定义介质合同，但不共享最终复波数；另有三固定模态的差分方法及冷复核。内部收敛与独立一致性均不能扩展为全模态、真实海底或完整会聚区声场认证。

二进制 SHA256：
- kraken.exe：997d533594c89f8cc63b520ae29cadec743ecfbdb64898f4c3294ab0254f4eb7
- krakenc.exe：ddbd5b43f937ca83dcb5093cced621396814720fd726c6a6d69ad522823d552a

新增声学传播程序调用为 0。历史 G0/N2/N3/PREC1/PREC2 不变，R4=0%。
