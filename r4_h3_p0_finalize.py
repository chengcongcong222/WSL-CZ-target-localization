"""Editorial/final-delivery checks only; does not change frozen scientific methods."""
from pathlib import Path
import json,csv
import r4_h3_p0 as E
def main():
    o=E.O;d=E.read(o/"H3_P0_DECISION.json");v=E.read(o/"VALIDATION.json")
    review=E.read(o/"INITIAL_GATE_REVIEW_VALIDATION.json")
    E.verify()
    assert v["PASS"] and review["PASS"] and d["new_KRAKEN_calls"]==0
    reg=E.G.rows(o/"NOMINAL_REGRESSION_AND_CHAIN_CONTROLS.csv")
    fail=[r for r in reg if r["PASS"]!="True"]
    assert len(fail)==20 and all(r["check"].startswith("source_depth_steps_") for r in fail)
    assert not d["science_information_evaluated"]
    old=E.G.rows(E.G.O/"NUMERICAL_STABILITY.csv")
    whole=max(float(r["depth_step_relative"]) for r in old if r["calibration"]=="C1b" and r["scene_id"].endswith("_M"))
    iv=max(float(r["depth_information_relative"]) for r in old if r["calibration"]=="C1b" and r["scene_id"].endswith("_M"))
    # Source-column gate new to this design; historical G0 Gates are not silently replaced.
    text=f"""# H3-P0：前置数值Gate阻断同步报告

正式判定：H3_P0_NUMERICAL_OR_FORMULA_INCOMPLETE

本轮已完成冻结配置的唯一一次前置检查，并按失败停止规则停止。
SSP正负扰动传播尚未执行；新增KRAKEN调用为0/52。
不能回答SSP未知量消除多少深度信息，不能归入A/B/C正负机制判定。

Parent SHA：{E.PARENT}。
设计SHA：{d["design_SHA"]}。
执行SHA为包含本报告的提交；remote/main最终核对记录保存于D盘运维目录，避免自引用SHA。

## 实际停止原因与独立重建

前置检查共252项：232通过、20未通过。
原名义复场、Jacobian、半步Jacobian和有效信息数值回归126/126通过。
72项水平解析—有限差分链式检查通过，18项实际名义自由模态增益深度吸收控制通过。
20项失败全部来自本轮新增的深度单方向差分稳定性检查：
在消除未知源与固定阵元频响后，比较1m与0.5m深度差分的深度方向本身，
相对差异必须不超过设计提交预注册的2%。
36个登记几何/资源/网格/噪声组合中，20个超过2%；最大2.6638063293%。
H01和H06包含失败，H12没有失败。相同失败模式在80001与160001均出现。
该项门槛未因结果而改变，也未补跑更小深度步长。

独立二进制解析、Cartesian场求和及QR投影复核27607项，0重建FAIL。
补充深度方向冷重算36/36与保存SVD结果一致，最大差异{review["max_QR_SVD_difference"]:.8g}。
它确认20项科学准入失败真实保留，不使科学Gate自动通过。
名义C1g-FREE及C2零控制54项全部通过；C1g源深与C2状态的不确定度为INF，无伪逆零误差。

## 与原H3-G0的区别

2%的“深度方向单列归一化”门槛是我在本轮设计中新增的较严格诊断，不是原G0冻结Gate。
G0原检查对整体状态Jacobian归一化，并检查剖面后有效深度信息稳定性。
不同归一化不能混为同一个回归：
本轮继承的C1b历史整体Jacobian深度步长差最大约{whole:.8g}；
历史网格/步长有效深度信息相对差最大约{iv:.8g}，按原G0合同保留原STABLE状态。
G0条件局部深度信息结论保持不变。本轮并没有证明该旧结论失效。
详细对照见G0_VS_P0_DEPTH_STEP_METRIC.csv；其中历史数据仅为只读说明，
没有追加实验，也没有反向修改旧数据、阈值或已接受判断。

## SSP provider准备与尚未验证的风险

已冻结52份正负偏移环境文件以及140项输入/代码绑定。
所有文件只修改E-STD水体SSP的251个深度点声速列，其余传播输入保持原值。
delta_c=-1/0/+1m/s只是数学敏感性诊断，不是海洋实测环境误差范围。
原相速度窗1500–1800m/s保持不变，负偏移使SSP最低声速降至1499m/s。
[KRAKEN官方说明](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html)
指出非零CLOW会排除更慢模态，所以该窗截断是预先登记的provider风险。

由于本轮在源深差分Gate已停止，没有负/正SSP新模态，
尚未实际检测模态进入/退出、两网格SSP导数稳定性或正负SSP割线局部一致性。
这些只能标NOT_EVALUATED，不能写成“观察到模态截断”或“SSP导数失败”。
不自动移动phase-speed窗口、不增加频率/网格、不减小SSP偏移量、不追加KRAKEN。

## 各场景名义深度信息及尚未建立的SSP信息

PROFILED_DEPTH_INFORMATION.csv保存72个名义C1b结果，
覆盖3几何×3资源×2网格×2噪声×2登记比例。
所有C1e-P0深度信息和保留比例留空，标为NOT_EVALUATED_GATE_FAILED，空值不是零信息。
以下为160001、1%设计噪声、B1的原名义有效源深信息（m^-2）；
全部水平状态、未知源和跨时固定阵元频响均按原模型剖面，无水平真值先验：

|几何|相对噪声|固定噪声底|SSP剖面后信息|
|---|---:|---:|---|
"""
    info=E.G.rows(o/"PROFILED_DEPTH_INFORMATION.csv")
    for sc in ["H01_M","H06_M","H12_M"]:
        vals=[next(float(r["C1b_depth_information"]) for r in info if r["scene"]==sc and r["resource"]=="B1" and int(r["mesh"])==160001 and float(r["sigma"])==.01 and r["noise_model"]==n) for n in E.NOISE]
        text+=f"|{sc}|{vals[0]:.6f}|{vals[1]:.6f}|NOT_EVALUATED|\n"
    text+="""\nRESOURCE_MATCHED_INFORMATION.csv保留B1/B2原名义对照；SSP下的B1/B2对照未评估。
没有借用Fisher trace、给水平状态添加真值先验或给不可辨方向赋零误差。
本轮继续采用G0的确定性未知复源均值模型，不混用G2E随机Gaussian协方差族Fisher。

## G2E审计澄清

旧判定和5184个单元原始结果均不改变。
K=1未通过979：847含未检测块，132为全检测后方向误差。
K=8未通过258、K=32未通过240：均含未检测块，全检测后方向误差失败为0。
K=32的240个未通过全部来自5%固定噪声底：
H01 144/144、H06 48/144、H12 0/144；其余三种噪声条件均432/432。
固定Markov门限不随K降低，检测平台不是物理垂向信息极限。
B1/B2响应提取通过数近似，也不推翻名义深度信息增量。
本轮只重聚合保存记录，不调整旧门限、不修复检测器、不追加随机样本。

## 适用范围、路线建议和停止

本轮只证明新的前置数值合同未闭合。SSP—源深混淆及信息生存性尚未评估。
下一阶段若获单独授权，应先决定并认证源深切向量的数值准入合同，
再处理已登记的模态窗与SSP局部差分风险，之后才有依据讨论更完整的物理环境模型。
不能从当前2%新门槛失败直接推断源深机制不可用，也不能为了进入SSP计算事后放宽本轮门槛。
此处是路线建议，不构成自动执行授权。

H3-G0/G1/G2E原结论不变。
完整RC2水平支持、真实UUV源谱、完整C1e环境鲁棒性均NOT_ESTABLISHED。
R4=0%。时域接收、测深Monte Carlo、E2、A2、SSP全阶段、P5及检测器修复未开启。
提交推送、核对remote/main后停止。
"""
    (o/"GPT_SYNC.md").write_text(text,encoding="utf-8")
    # Precise no-evaluation tables replace only the newly-created reporting placeholders.
    E.write(o/"FORWARD_NUMERICAL_FIDELITY.csv",[{"offset_mps":delta,
        "new_KRAKEN_calls":0,"two_grid_forward_certification":"NOT_EVALUATED_PREPROPAGATION_GATE_FAILED",
        "reason":"SOURCE_DEPTH_COLUMN_STEP_LIMIT","interpretation":"NO_OBSERVED_SSP_MODE_FAILURE"}
        for delta in [-1,1]])
    E.write(o/"DEPTH_SSP_DERIVATIVE_ALIGNMENT.csv",[{"scene":s["scene_id"],"resource":r,
        "noise_model":n,"SSP_derivative":"NOT_EVALUATED","depth_SSP_alignment":"NOT_EVALUATED",
        "reason":"PREPROPAGATION_SOURCE_DEPTH_GATE_FAILED"}
        for s in E.scenes() for r in E.RES for n in E.NOISE])
    E.dump(o/"SUPPLEMENTAL_REVIEW_PROVENANCE.json",{
        "review_code_sha256":E.sha(E.ROOT/"r4_h3_p0_initial_gate_review.py"),
        "review_code_added_after_freeze":"INDEPENDENT_READ_ONLY_RECONSTRUCTION_ONLY",
        "additional_historical_input":str(E.G.O/"NUMERICAL_STABILITY.csv").replace("\\","/"),
        "historical_input_sha256":E.sha(E.G.O/"NUMERICAL_STABILITY.csv"),
        "historical_input_use":"METRIC_SCOPE_EXPLANATION_ONLY_NOT_NEW_SCIENCE_ADMISSION",
        "frozen_scientific_code_or_threshold_changed":False,"new_solver_calls":0})
    master=E.ROOT/"results/R4_MASTER"
    plan=master/"R4_PLAN.md"
    src=plan.read_text(encoding="utf-8")
    title="\n\n## H3-P0 structured propagation sensitivity\n\n"
    start=src.rfind(title)
    assert start>=0
    plan.write_text(src[:start]+title+
        "H3_P0_NUMERICAL_OR_FORMULA_INCOMPLETE. Original nominal replay126/126 PASS; "
        "new pre-registered2% depth-column step Gate20/36 FAIL. "
        "0 new KRAKEN; SSP sensitivity NOT_EVALUATED. G0/G1/G2E unchanged. R4=0%; STOP.\n",encoding="utf-8")
    ledger=master/"R4_EVIDENCE_LEDGER.csv"
    src=ledger.read_text(encoding="utf-8")
    src=src.replace("FINITE_SINGLE_SSP_PARAMETER_DIAGNOSTIC;NO_FULL_ENVIRONMENT_OR_DEPTH_ESTIMATION",
                    "PRE_PROPAGATION_DEPTH_COLUMN_STEP_GATE_FAILED;SSP_NOT_EVALUATED")
    ledger.write_text(src,encoding="utf-8")
    E.dump(master/"R4_H3_P0_PROGRESS.json",d)
    assert not list((o/"modal").glob("*.mod"))
    assert len(list((o/"modal").glob("*.env")))==52
    assert len(info)==72 and all(not r["C1e_P0_depth_information"] and not r["retention_ratio"] for r in info)
    zeros=E.G.rows(o/"UNKNOWN_GAIN_ZERO_CONTROLS.csv")
    assert len(zeros)==54 and all(r["PASS"]=="True" for r in zeros)
    E.dump(o/"FINAL_DELIVERY_VALIDATION.json",{
        "PASS":True,"required_artifacts_complete":True,"C1b_replay_rows":72,
        "SSP_information_rows_evaluated":0,"new_KRAKEN_calls":0,
        "prepared_env_files":52,"cold_nominal_review_PASS":v["PASS"],
        "depth_column_QR_review_PASS":review["PASS"],"zero_controls_PASS":54,
        "new_science_gate_failures_retained":20,"original_G0_status_unchanged":True,
        "frozen_code_input_sha_recheck":True,"no_extra_simulation":True,"R4_percent":0})
    required=["DESIGN_FREEZE.json","SSP_PARAMETER_AND_PROVIDER_CONTRACT.md","FORWARD_NUMERICAL_FIDELITY.csv",
      "DEPTH_SSP_DERIVATIVE_ALIGNMENT.csv","PROFILED_DEPTH_INFORMATION.csv","RESOURCE_MATCHED_INFORMATION.csv",
      "UNKNOWN_GAIN_ZERO_CONTROLS.csv","H3_P0_DECISION.json","VALIDATION.json","GPT_SYNC.md"]
    assert all((o/x).exists() for x in required)
    manifest=[{"path":str(p.relative_to(E.ROOT)).replace("\\","/"),"sha256":E.sha(p),"bytes":p.stat().st_size}
       for p in sorted(o.rglob("*")) if p.is_file() and p.name!="DELIVERABLE_MANIFEST.json"]
    assert sum(x["bytes"] for x in manifest)<256000000
    E.dump(o/"DELIVERABLE_MANIFEST.json",{"artifacts":manifest,"bytes":sum(x["bytes"] for x in manifest)})
    print(json.dumps({"PASS":True,"artifacts":len(manifest),"bytes":sum(x["bytes"] for x in manifest),
      "classification":d["decision"],"new_KRAKEN":0,"G0_whole_J_step_max":whole,
      "G0_effective_Iz_stability_max":iv}))
if __name__=="__main__":main()
