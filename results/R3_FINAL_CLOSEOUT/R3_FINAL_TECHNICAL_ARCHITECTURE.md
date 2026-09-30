# R3_FINAL_TECHNICAL_ARCHITECTURE

UTC: 2026-09-30

## 主线架构

$$
\boxed{
\text{RC1 会聚区环境先验}
\rightarrow
\text{RC2 连续方位/平台运动维持运动学候选}
\leftrightarrow
\text{RC3 多频相对传播形状约束}
\rightarrow
\text{小转向形成虚拟基线}
\rightarrow
\text{跨频别名互补切断距离脊线}
\rightarrow
\text{后续递归更新}
}
$$

## 状态定义

$$
x = [r, \theta, z, v, \psi]^T
$$

处理结构：

```text
main tracking state = [r, theta, v, psi]
z = PROFILED_NUISANCE_VARIABLE
scientific depth role = MIXED_OR_UNRESOLVED
```

**禁止**写成：`single-array one-shot five-parameter inversion`

## 模块说明

| 模块 | 功能 | 状态 |
|---|---|---|
| RC1 | 会聚区环境先验（E-STD Munk, H=5000m） | 已验证 |
| RC2 | 连续方位 + 平台运动 → 运动学候选云 | 已验证 |
| RC3 | 多频相对 TL 形状 → profiled margin | 已验证 |
| 转向 | 小幅平台转向形成虚拟基线 | 已验证（15°） |
| 多频 | 跨频别名互补切断 r-z 脊线 | 已验证 |

## 处理流程

1. RC2 从 P2 搜索盒（114576 节点）维持候选
2. RC3 对 RC2 候选做多频相对传播约束
3. 平台转向形成虚拟观测基线
4. 多频联合排除 r-z 别名
5. 输出候选集合（非单点估计）
