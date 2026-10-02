# R3_FINAL_UNVERIFIED_AND_FORBIDDEN_CLAIMS

UTC: 2026-09-30

Canonical sync: 2026-10-02。R3已冻结；参数候选集审计通过不改变以下未验证范围。

## 未验证项

1. 真实 HLA 数据中能否稳定提取/跟踪当前 relative-TL feature
2. 真实 UUV 一定存在 201/235/283/338 四条稳定线
3. 幅漂+SSP+频漂联合 corner
4. 真实海洋 SSP 误差统计
5. 在线 SSP 估计/补偿能否恢复距离锚
6. frequency-tracking estimation error
7. 45–60 km 所有真实距离均可锚定
8. 所有转向角
9. 多随机 bearing-noise realizations
10. 实际海试
11. 五参数同时 <10%
12. 连续准确测深
13. 端到端观测→特征提取→状态估计
14. P5 未开启

## 禁止承诺

- 单阵实现五参数一次高精度联合反演
- 已实现五参数 <10%
- 已证明多参数联合误差<10%。
- 50–60 km 全范围高精度测距已验证
- 任意 3 条线即可
- 实际 UUV 必有上述线谱
- 2% 未知多普勒误差均鲁棒
- SSP 失配一定是海洋中最大误差
- 深度已准确估计
- 真实 HLA 特征提取已完成
- RC3/系统端到端已经完成

最终matched control的top-1在六个case均精确落在真值节点，但全survivor-set仍保留4°/5°航向节点，统一<10%集合界未建立。theta_truth=0°时通常的相对百分比误差未定义，z是profile nuisance，甲方联合误差范数也尚未定义。因此20%只能作为 `DIAGNOSTIC_ONLY_NOT_REQUIREMENT_VERIFICATION`，不能写成正式验收通过或不达标。
