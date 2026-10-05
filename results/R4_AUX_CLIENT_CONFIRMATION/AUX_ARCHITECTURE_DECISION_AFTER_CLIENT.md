# 客户回复后的辅助节点架构决策

当前 Gate2B 联合静态要求已接受，实际硬件 UNKNOWN。数值架构工作 STOPPED；R4-A1-NEW NOT_OPENED；R4=0%。本表只预定义回复后的决策，不授权任何新实验。

客户证据应记录实际值、单位、误差统计定义、工况、资料出处及审核人。随机方位、固定系统残差、导航协方差和实际布放分别映射；不把波束宽度当DOA精度，不把CEP当每轴σ，不忽略相关误差。未知或口径不兼容的资料不能视为通过，也不能视为明确不满足。

## Branch A — T5 hardware-compatible

客户证据能够支持至少一个完整T5联合向量，且其误差模型与测量条件可映射时，由负责人审核后可记录：

`AUX1_HARDWARE_CONDITIONALLY_ADMITTED_FOR_R4_A1_NEW`

这仍是条件性准入。正式打开R4-A1-NEW前，必须另审时间同步、目标关联和measurement-chain provenance。负责人明确授权后，下一阶段才可针对增强架构开展off-grid、方位精度与动态双节点acquisition验证；本包不启动该阶段。

## Branch B — T10 only

不能支持T5，但能支持至少一个完整T10联合向量时，可记录：

`AUX1_SUPPORTS_PROJECT_LEVEL_10PCT_RANGE_SUBMETRIC_ONLY`

负责人决定是否接受、调整项目指标或将距离子指标设为≤10%。距离子指标不代表多参数联合误差已满足项目目标，不自动打开A1。

## Branch C — hardware insufficient

客户证据明确表明方位、基线、导航或标定等关键项使两个T10联合向量均不兼容时，可记录：

`AUX1_BEARING_ONLY_ARCHITECTURE_NOT_ADMITTED`

某一个anchor不兼容不排除另一个anchor。未提供资料仅保持待确认，不能据此触发Branch C。负责人可再考虑TDOA、外部距离先验、主动/协同测距或项目目标调整；不自动开始任何替代路线。

## 未回复或资料不完整

保持`CLIENT_CONFIRMATION_REQUIRED`；所有未确认项为UNKNOWN，数值架构工作继续停止。同步、目标关联和通信只记录事实，不依据本包生成数值验收门限。负责人审查映射与适用范围后才更改状态。

## 共同边界

Gate2B仅覆盖静态瞬时几何、冻结误差模型、11个离散距离点及符号角点。T5/T10向量与各自预算族保持完整，不能使用Gate2A单因素最大值拼接系统能力；也不能从角点实验推断连续误差超矩形内部全覆盖。实际硬件事实与科学仿真证据分列。

无论进入哪一分支，下一实验均须负责人明确授权。不得自动Gate2C、自动重开A1、自动更换架构或恢复深度工程。冻结R3/A1/FIX1/FIX2/tractability及已有路线结论。
