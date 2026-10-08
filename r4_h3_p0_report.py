"""Reporting only; retains failed numerical Gates and historical detector evidence."""
import json,csv
from pathlib import Path
import r4_h3_p0 as E
def report():
    o=E.O;d=E.read(o/"H3_P0_DECISION.json");v=E.read(o/"VALIDATION.json")
    numeric=E.G.rows(o/"NUMERICAL_AND_FORMULA_GATES.csv") if (o/"NUMERICAL_AND_FORMULA_GATES.csv").exists() else []
    failed=[r for r in numeric if r["PASS"]!="True"]
    families={}
    for r in failed:
        name=r["check"]
        family=next((x for x in ["fixed_phase_window_branch_set","SSP_secant_locality",
        "SSP_secant_two_grid","effective_info_two_grid","profiled_rank_two_grid",
        "fidelity","relative_chart_covariance","DPI","cold","free_modal"] if name.startswith(x)),name)
        families[family]=families.get(family,0)+1
    diag=E.G.rows(o/"DEPTH_SSP_DERIVATIVE_ALIGNMENT.csv") if (o/"DEPTH_SSP_DERIVATIVE_ALIGNMENT.csv").exists() else []
    frac=[]
    if diag:
        frac=[float(r["positive_negative_secant_disagreement"]) for r in diag]
    text=f"""# H3-P0 SSP传播敏感性同步报告

正式判定：{d["decision"]}

Parent SHA：{E.PARENT}。
设计SHA：{d["design_SHA"]}。
执行SHA为包含本报告的提交；最终remote/main核对保存于D盘运维目录，避免自身提交SHA循环。

## 本轮实际完成及认证

新增KRAKEN调用：{d["new_KRAKEN_calls"]}/52；provider完成：{d["provider_complete"]}。
名义G0回归：{d["nominal_regression_PASS"]}。
数值及局部导数准入：{d["numerical_and_locality_gate_PASS"]}。
内部独立二进制解析、Cartesian场重建和QR剖面复核：{v["checks"]}项，{v["FAIL"]} FAIL。
执行及冷复核总耗时：{v["execution_plus_cold_audit_s"]:.2f}s，预算5400s。
交付存储截至审计：{v["delivered_bytes"]} bytes，预算256,000,000 bytes。

数值Gate失败与实现重建失败分开。冷重算一致，并不能使未通过的物理/局部导数准入自动通过。
失败检查族：{json.dumps(families,ensure_ascii=False)}。
正负投影后割线分歧范围：{min(frac) if frac else None}—{max(frac) if frac else None}；冻结上限0.2。
所有失败频率、资源、场景和噪声条件保留，不改变偏移量、相速度边界、网格或科学规则。

## 物理参数与局部性边界

仅改变E-STD的水体SSP各深度点声速，delta_c=-1/0/+1 m/s。
这是数学敏感性诊断，不是实测海洋环境不确定度范围。
边界、衰减设置、密度、深度、频率、源深、接收位置以及1500–1800 m/s相速度窗口均保持原值。
KRAKEN非零CLOW会排除更慢的模态；名义SSP最小1500 m/s，负偏移后最小1499 m/s。
模态窗的进入/退出和正负割线不一致会污染局部环境Jacobian认证。
因此两网格一致性单独通过，也不足以证明±1 m/s割线是局部连续SSP切向方向。
本轮不自动移动相速度窗、缩小偏移量或追加求解以修复准入。

[官方KRAKEN参数说明](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html)
支持上述CLOW语义；它不替代本轮数值认证。

## 深度信息与资源对照

PROFILED_DEPTH_INFORMATION.csv保留全部72个名义C1b结果（3几何×3资源×2网格×2噪声×2比例）。
若准入失败，C1e-P0深度信息与保留比例留空并明确NOT_EVALUATED_GATE_FAILED，
不能把空值解释为零信息、传播吞噬或物理不可能。DEPTH_SSP_DERIVATIVE_ALIGNMENT.csv此时仅为有限割线诊断。
C1b数字来自原确定性未知复源均值模型的回归。没有给水平状态加入真值先验；
源深信息仍须剖面全部水平状态和未知跨时固定阵元频响。
B0/B1/B2资源完全继承；B1/B2名义增量不改写，未认证SSP下的增量不报告。
自由模态增益/C2严格零控制保留，无界不确定度写为INF，不用伪逆分配虚假零误差。
不得将本轮一个SSP参数称为完整C1e环境鲁棒性。
"""
    if d["science_information_evaluated"]:
        text+="\n本轮局部认证通过；逐条件保留比例见PROFILED_DEPTH_INFORMATION.csv。该结论仅限一个登记SSP参数。\n"
    else:
        text+="\n本轮尚不能回答深度信息在SSP未知量下存活多少，也不能做A/B正负科学解释。应先单独审查传播provider和局部导数认证，之后再决定是否投入更完整的环境模型；本轮不执行该下一阶段。\n"
    text+="""\n## G2E审计澄清及历史冻结

保留H3_G2E_COMPRESSED_RESPONSE_EXTRACTION_UNRELIABLE和所有冻结结果。
K=1质量失败979：847含未检测块，132为全检测后的方向误差；
K=8失败258及K=32失败240均包含未检测块，全检测后超过0.1方向误差的数量为0。
K=32的240个失败全部来自5%固定噪声底：H01 144/144、H06 48/144、H12 0/144。
其余三种登记噪声条件K=32均432/432。
固定Markov门限不随K下降，不能将K=8到32的检测平台解释成物理深度信息极限。
B1/B2提取通过数近似，也不是深度信息增益的替代检验。
此澄清只重聚合已保存结果，没有修复检测器、调整门限或补跑G2E。

## 冻结状态与停止

H3-G0条件局部深度信息：UNCHANGED。
H3-G1传播/支持缺口：UNCHANGED。
H3-G2E：UNCHANGED。
完整RC2支持、真实UUV源谱、完整C1e环境鲁棒性：NOT_ESTABLISHED。
R4=0%。新随机仿真、测深Monte Carlo、时域接收、检测器修复均未执行。
不自动进入E2、A2、SSP全阶段或P5。提交推送、核对remote HEAD后STOP。
"""
    (o/"GPT_SYNC.md").write_text(text,encoding="utf-8")
    # Required tables still explicitly exist if provider/precheck stopped early.
    for name in ["FORWARD_NUMERICAL_FIDELITY.csv","DEPTH_SSP_DERIVATIVE_ALIGNMENT.csv","UNKNOWN_GAIN_ZERO_CONTROLS.csv"]:
        if not (o/name).exists():E.write(o/name,[{"status":"NOT_EVALUATED_PREVIOUS_GATE_FAILED"}])
    master=E.ROOT/"results/R4_MASTER"
    with (master/"R4_PLAN.md").open("a",encoding="utf-8") as f:
        f.write(f"\n\n## H3-P0 structured propagation sensitivity\n\n{d['decision']}. "
                "Single uniform SSP parameter only; numerical/locality Gate retained. "
                "G0/G1/G2E frozen, no full C1e or RC2 coverage. R4=0%; STOP.\n")
    with (master/"R4_EVIDENCE_LEDGER.csv").open("a",encoding="utf-8",newline="") as f:
        csv.writer(f).writerow(["R4_H3_P0_STRUCTURED_PROPAGATION_SENSITIVITY",
          "../R4_H3_P0_PROPAGATION_SENSITIVITY/H3_P0_DECISION.json",d["decision"],
          "FINITE_SINGLE_SSP_PARAMETER_DIAGNOSTIC;NO_FULL_ENVIRONMENT_OR_DEPTH_ESTIMATION",
          "PENDING_RESEARCH_LEAD_AUDIT;STOP;NO_CREDIT",0])
    E.dump(master/"R4_H3_P0_PROGRESS.json",d)
    manifest=[{"path":str(p.relative_to(E.ROOT)).replace("\\","/"),"sha256":E.sha(p),"bytes":p.stat().st_size}
              for p in sorted(o.rglob("*")) if p.is_file() and p.name!="DELIVERABLE_MANIFEST.json"]
    E.dump(o/"DELIVERABLE_MANIFEST.json",{"artifacts":manifest,"bytes":sum(x["bytes"] for x in manifest)})
    print(json.dumps({"decision":d["decision"],"numerical_failure_families":families,"validation":v},indent=2))
if __name__=="__main__":report()
