"""Post-report independent aggregation/integrity review; no new stochastic draws."""
import csv,json,hashlib,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
ROOT=Path(r"D:\XiaomiMiMoProjects\WSL系统汇聚区目标定位")
sys.path.insert(0,str(ROOT))
import r4_h3_g2e as E
o=E.OUT
read=lambda name:list(csv.DictReader((o/name).open(encoding="utf-8-sig")))
cells=read("EXTRACTION_BY_CELL.csv")
checks=0; failure=[]
def check(name,ok):
    global checks
    checks+=1
    if not ok: failure.append(name)
keys=["scene","resource","pack","K","noise_model","sigma","realization"]
lookup={tuple(r[k] for k in keys):r for r in cells}
check("5184_unique_cells",len(cells)==len(lookup)==5184)
required=["DESIGN_FREEZE.json","RAW_VS_CSD_LIKELIHOOD_AUDIT.csv",
"SOURCE_AND_NOISE_COVARIANCE_CONTROLS.csv","EXTRACTION_BY_CELL.csv","RESOURCE_MATCHED_EXTRACTION.csv",
"LOW_ENERGY_FAILURES.csv","GAUGE_AND_NUISANCE_AUDIT.md","H3_G2E_DECISION.json","VALIDATION.json","GPT_SYNC.md"]
for name in required: check("required_"+name,(o/name).is_file())
E.verify_inputs(json.loads((o/"DESIGN_FREEZE.json").read_text(encoding="utf-8")))
for r in read("SOURCE_AND_NOISE_COVARIANCE_CONTROLS.csv"): check("pre_"+r["control"],r["PASS"]=="True")
validation=json.loads((o/"VALIDATION.json").read_text(encoding="utf-8"))
check("frozen_cold_audit",validation["PASS"] and validation["FAIL"]==0)
groups=defaultdict(list)
for r in cells: groups[tuple(r[k] for k in keys[:-1])].append(r)
for r in read("EXTRACTION_SUMMARY.csv"):
    rr=groups[tuple(r[k] for k in keys[:-1])]
    check("summary_denominator",len(rr)==16 and int(r["registered_realizations"])==16)
    check("summary_quality",int(r["accepted_quality_realizations"])==sum(x["accepted_quality"]=="True" for x in rr))
    check("summary_failures",int(r["failed_blocks"])==sum(int(x["failed_blocks"]) for x in rr))
for r in read("RESOURCE_MATCHED_EXTRACTION.csv"):
    kk=[r[x] for x in ["scene","pack","K","noise_model","sigma","realization"]]
    c1=lookup[(kk[0],"B1",*kk[1:])]
    c2=lookup[(kk[0],"B2",*kk[1:])]
    check("paired_B1",r["B1_quality"]==c1["accepted_quality"] and r["B1_failed_blocks"]==c1["failed_blocks"])
    check("paired_B2",r["B2_quality"]==c2["accepted_quality"] and r["B2_failed_blocks"]==c2["failed_blocks"])
low=read("LOW_ENERGY_FAILURES.csv")
expect={tuple(r[k] for k in keys) for r in cells if int(r["failed_blocks"]) or float(r["max_noise_to_channel_amplitude"])>=.1}
check("low_energy_none_removed",{tuple(r[k] for k in keys) for r in low}==expect)
# Read all archived feature arrays once; independent projector moments.
features={}
for path in (o/"features").glob("*n160001.npz"):
    with np.load(path) as z: features[path.stem]={k:z[k] for k in ["u","detected"]}
fieldcache={s:E.fields(s,160001) for s in ["H01","H06","H12"]}
def projector(u):
    return np.outer(u,u.conj())/np.vdot(u,u).real
for r in read("NORMALIZED_RESPONSE_BIAS_VARIANCE.csv"):
    scene=r["scene"]; res=r["resource"]; kind=r["noise_model"]; sigma=float(r["sigma"])
    ri=list(E.CH).index(res); n=len(E.CH[res]); ki=E.KS.index(int(r["K"]))
    t=[0,600,1200].index(int(r["template_s"])); f=list(E.F).index(int(r["frequency_hz"]))
    vv=[]
    for rep in range(16):
        z=features[f"{scene}_{kind}_s{sigma:.2f}_rep{rep:02d}_n160001"]
        if z["detected"][ri,ki,t,f]: vv.append(projector(z["u"][ri,ki,t,f,:n]))
    check("moment_count",int(r["detected"])==len(vv) and int(r["failed"])==16-len(vv))
    if vv:
        mean=sum(vv)/len(vv)
        p=fieldcache[scene][t,E.CH[res],f]
        bias=np.linalg.norm(mean-projector(p))
        var=sum(np.linalg.norm(x-mean)**2 for x in vv)/len(vv)
        check("moment_bias",abs(float(r["conditional_projector_bias"])-bias)<1e-10)
        check("moment_variance",abs(float(r["conditional_projector_variance"])-var)<1e-10)
    else:
        check("no_direction_no_moment",not r["conditional_projector_bias"] and not r["conditional_projector_variance"])
d=json.loads((o/"H3_G2E_DECISION.json").read_text(encoding="utf-8"))
st=read("NUMERICAL_EXTRACTION_STABILITY.csv")
unstable=sum(int(r["detection_disagreements"])>0 or
             (bool(r["max_accepted_projector_grid_difference"]) and float(r["max_accepted_projector_grid_difference"])>.05) for r in st)
accepted=sum(r["accepted_quality"]=="True" for r in cells)
expected="H3_G2E_COMPRESSED_RESPONSE_EXTRACTION_UNRELIABLE" if (
    1-accepted/5184>.05 or unstable/5184>.05) else "H3_G2E_CONDITIONAL_FREQUENCY_STATISTIC_ESTABLISHED"
check("frozen_classification",d["classification"]==expected and d["quality_accepted_cells"]==accepted)
check("boundaries",d["R4_percent"]==0 and d["TIME_DOMAIN_RECEIVED_SIGNAL_EXTRACTION"]=="NOT_OPENED"
      and d["FULL_RC2_DEPTH_SET"]=="NOT_ESTABLISHED" and d["next_stage"]=="NOT_AUTHORIZED;STOP")
designver=json.loads((o/"DESIGN_COMMIT_VERIFICATION.json").read_text(encoding="utf-8-sig"))
check("design_remote_verified",designver["verified"] and designver["design_SHA"]==designver["remote_main_SHA"])
stats=[]
for k in ["1","8","32"]:
    rr=[r for r in cells if r["K"]==k]
    stats.append({"K":int(k),"accepted":sum(r["accepted_quality"]=="True" for r in rr),
                  "denominator":len(rr),"failed_blocks":sum(int(r["failed_blocks"]) for r in rr)})
for kind in ["RELATIVE","ABSOLUTE_FLOOR"]:
    for sig in ["0.01","0.05"]:
        rr=[r for r in cells if r["noise_model"]==kind and float(r["sigma"])==float(sig)]
        stats.append({"noise_model":kind,"sigma":sig,
                      "accepted":sum(r["accepted_quality"]=="True" for r in rr),"denominator":len(rr)})
E.write_json(o/"REPORT_AND_DELIVERY_VALIDATION.json",{"PASS":not failure,"checks":checks,
    "FAIL":len(failure),"failures":failure,"review_script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "new_draws":0,"new_experiments":0,"independent_aggregation_stats":stats})
# Add a Chinese synchronization synopsis without changing frozen scientific report or decision.
append="\n\n## 中文同步摘要\n\n"
append+=f"正式分类：{d['classification']}。完整频域样本—联合CSD统计链通过；压缩响应按冻结标准质量通过{accepted}/5184。\n\n"
append+="|K|质量通过|登记单元|未通过能量门限的频率/模板块|\n|---|---:|---:|---:|\n"
for s in stats[:3]: append+=f"|{s['K']}|{s['accepted']}|{s['denominator']}|{s['failed_blocks']}|\n"
append+="\n上述门限是预登记的保守提取诊断，失败不能解释为物理不可辨识或深度信息不存在。"
append+="B1/B2只能作资源匹配的提取稳定性比较，不能由质量通过数重判G0深度信息。\n\n"
append+="所有相对/固定噪声底、弱场、K=1失败均保留。未知源标量经归一化去除，但未知阵元频率增益仍留在响应中；"
append+="未经校准的传播响应恢复未建立。检测条件下的偏差/方差明确保留全16次实现分母。\n\n"
append+=f"冻结设计SHA：{designver['design_SHA']}。执行SHA为包含本报告的提交，最终remote/main核对记录放在D盘运维目录。\n"
append+="H3-G0/G1原结论不变；实际UUV源占用、传播未知量和完整RC2深度集合未建立；R4=0%。提交推送后停止。\n"
with (o/"GPT_SYNC.md").open("a",encoding="utf-8") as f: f.write(append)
manifest=[{"path":str(p.relative_to(ROOT)).replace("\\","/"),"sha256":E.sha(p),"bytes":p.stat().st_size}
          for p in sorted(o.rglob("*")) if p.is_file() and p.name!="DELIVERABLE_MANIFEST.json"]
total=sum(x["bytes"] for x in manifest)
check("final_storage_cap",total<1000000000)
E.write_json(o/"DELIVERABLE_MANIFEST.json",{"artifacts":manifest,"total_bytes":total})
print(json.dumps({"PASS":not failure,"checks":checks,"FAIL":len(failure),"bytes":total,"decision":expected},indent=2))
sys.exit(0 if not failure else 2)
