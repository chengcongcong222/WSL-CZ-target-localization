"""Design-only freezer: no propagation execution, no random draws."""
import json,csv,hashlib
from decimal import Decimal
from pathlib import Path
import r4_h3_p0 as E
O=E.O;O.mkdir(exist_ok=True);(O/"modal").mkdir(exist_ok=True)
def main():
    assert E.G.git("rev-parse","HEAD")==E.PARENT
    task=Path(r"C:\Users\ccc\.codex\attachments\73c5914e-3186-49d1-932c-9711e4fbbd9d\已粘贴的文本.txt")
    (O/"RESEARCH_LEAD_TASK.md").write_text(task.read_text(encoding="utf-8"),encoding="utf-8")
    (O/".gitattributes").write_text("*.mod -text -diff\n*.npz -text -diff\n",encoding="utf-8")
    envchecks=[]
    for d in [-1,1]:
        for g in E.GRIDS:
            for f in E.F:
                nominal=E.P.path_for(f,g).with_suffix(".env")
                lines=nominal.read_text(encoding="utf-8").splitlines();out=[];count=0
                for i,line in enumerate(lines):
                    parts=line.split()
                    try:
                        ssp=len(parts)==6 and Decimal(0)<=Decimal(parts[0])<=Decimal(5000) and Decimal(parts[1])>1000
                    except Exception:ssp=False
                    if ssp:
                        old=Decimal(parts[1]);parts[1]=str(old+Decimal(d));out.append(" ".join(parts));count+=1
                        envchecks.append({"frequency_hz":f,"mesh":g,"offset_mps":d,"line":i+1,
                            "depth_m":parts[0],"nominal_c_mps":str(old),"shifted_c_mps":parts[1],
                            "only_c_changed":True})
                    else:out.append(line)
                assert count==251
                E.modpath(f,g,d).with_suffix(".env").write_text("\n".join(out)+"\n",encoding="utf-8")
    E.write(O/"SSP_PARAMETER_DIFF.csv",envchecks)
    old=E.G.O
    paths=[E.ROOT/x for x in ["r4_h3_p0.py","r4_h3_p0_audit.py","r4_h3_p0_report.py","r4_h3_p0_design.py",
      "r4_h3_g0.py","r4_h3_g0_audit.py","r4_e2_g0.py","r4_e2_paired_recovery.py",
      "r4_e2_forward_diagnostic.py","r4_e2_nine_pair_pilot_audit.py"]]
    paths+=[old/"DESIGN_FREEZE.json",old/"VERTICAL_INFORMATION_BY_SCENE.csv",old/"H3_G0_DECISION.json",
      E.ROOT/"results/R4_H3_G1_EXTRACTION_SUPPORT_REVIEW/H3_G1_READINESS_DECISION.json",
      E.ROOT/"results/R4_H3_G2E_CONDITIONAL_EXTRACTION/H3_G2E_DECISION.json",
      E.ROOT/"results/R4_H3_G2E_CONDITIONAL_EXTRACTION/EXTRACTION_BY_CELL.csv",
      O/"RESEARCH_LEAD_TASK.md",E.M.AT/"kraken.exe"]
    paths += [E.P.path_for(f,g) for g in E.GRIDS for f in E.F]
    paths += [E.P.path_for(f,g).with_suffix(".env") for g in E.GRIDS for f in E.F]
    paths += [old/f"INPUT_{s['scene_id']}_{r}_n{g}.npz" for s in E.scenes() for r in E.RES for g in E.GRIDS]
    paths += [E.modpath(f,g,d).with_suffix(".env") for d in [-1,1] for g in E.GRIDS for f in E.F]
    binding=[{"path":str(p.resolve().relative_to(E.ROOT)).replace("\\","/"),
              "sha256":E.sha(p),"bytes":p.stat().st_size,
              "hash_mode":"CANONICAL_LF" if p.suffix in [".py",".md",".env",".csv",".json",".txt"] else "RAW"} for p in paths]
    design={"stage":"R4_H3_P0_STRUCTURED_PROPAGATION_SENSITIVITY","parent_SHA":E.PARENT,
      "source_model":"G0_DETERMINISTIC_UNKNOWN_COMPLEX_MEAN;NOT_G2E_GAUSSIAN_COVARIANCE_FISHER",
      "scenes":[s["scene_id"] for s in E.scenes()],"resources":E.RES,"snapshots_s":[0,600,1200],
      "source_z_m":200,"depth_steps_m":[1,.5],"state":["r","theta","z","v","psi"],
      "state_scales":E.G.SCALE.tolist(),"horizontal_known_truth_prior":False,
      "frequencies_hz":E.F,"grids":E.GRIDS,"main_grid":160001,
      "physical_parameter":"SSP_UNIFORM_SOUND_SPEED_OFFSET","offsets_mps":[-1,0,1],
      "offset_interpretation":"MATHEMATICAL_DIAGNOSTIC_NOT_EMPIRICAL_ERROR_RANGE",
      "environment_change":"ONLY_WATER_SSP_C_COLUMN;251_ROWS_0_TO_5000M",
      "phase_speed_window_mps":[1500,1800],"other_environment_inputs":"EXACTLY_INHERITED",
      "noise_models":E.NOISE,"sigma":E.SIGMAS,"absolute_reference":E.G.REF,
      "noise_contract":"G0_INDEPENDENT_PROPER_COMPLEX_RAW_FIELD_NOISE;DENSE_REFERENCE_CONTRAST_COVARIANCE",
      "calibration":"C1b_COMPLEX_GAIN_EACH_ELEMENT_FREQ_FIXED_ACROSS_THREE_TIMES",
      "source":"UNKNOWN_COMPLEX_SOURCE_EACH_FREQ_TIME",
      "SSP_unknown":"ONE_REAL_CONTINUOUS_PARAMETER;NO_PRIOR",
      "SSP_derivative":"CENTERED_MINUS_PLUS_PRESSURE_SECANT_AT_DELTA1;MUST_PASS_LOCALITY_BEFORE_INFORMATION",
      "C1eP0":"PROJECT_UNKNOWN_SOURCE_FIXED_GAIN_AND_SSP_DIRECTION_THEN_ALL_HORIZONTAL",
      "zero_controls":["ACTUAL_C1g_FREE_MODAL_SOURCE_COEFFICIENT_ABSORPTION","C2_FREE_ALL_RESPONSES"],
      "finiteF":"NOT_EXECUTED;MAIN_LOCAL_SENSITIVITY_ONLY",
      "full_RC2_support":"NOT_ESTABLISHED","C1e_full_environment":"NOT_ESTABLISHED",
      "provider_budget":{"new_KRAKEN_calls":52,"per_call_s":90,"total_execution_plus_audit_s":5400,
                         "delivered_bytes":256000000,"solver_parallelism":1,"rerun":False},
      "provider_failure":"STOP_NEW_CALLS_PRESERVE_OUTPUTS",
      "numerical_preregistration":{
        "nominal_pressure_J_regression":1e-9,"nominal_information_regression_relative":1e-7,
        "analytic_FD_chain_relative":.002,"source_depth_step_relative":.02,
        "mode_count_grid_match":True,"min_phi_grid_correlation":.99,
        "max_k_grid_phase_difference_60km_rad":.02,
        "source_free_response_error_over_noise1":.5,"SSP_projected_secant_grid_relative":.02,
        "effective_depth_info_grid_relative":.1,"profiled_rank_grid_match":True,
        "mode_count_no_branch_entry_exit_over_offsets":True,
        "projected_positive_negative_secant_difference_over_central":.2,
        "source_relative_chart_F_relative":1e-7,"DPI_relative":1e-7,"free_modal_zero":1e-10,
        "SVD_rank_relative":1e-10,"null_depth_loading_tolerance":1e-8,
        "unidentifiable_uncertainty":"INF;NO_ZERO_VARIANCE_PSEUDOINVERSE"},
      "locality_requirement":"FIXED_PROVIDER_PHASE_WINDOW_MUST_PRESERVE_BRANCH_SET;POSITIVE_NEGATIVE_SECANTS_MUST_AGREE",
      "locality_failure_scope":"TESTED_STENCIL_NOT_CERTIFIED_LOCAL_JACOBIAN;NOT_ENVIRONMENT_PHYSICAL_NO_GO",
      "admission_failure":"H3_P0_NUMERICAL_OR_FORMULA_INCOMPLETE;NO_SSP_INFORMATION_SIGN_INTERPRETATION",
      "classification_after_admission":{
         "A":"ALL_MAIN_GRID_REGISTERED_RESOURCE_SCENE_NOISE_SIGMA_RETENTION>=0.5",
         "B":"ALL_MAIN_GRID_REGISTERED_RESOURCE_SCENE_NOISE_SIGMA_RETENTION<=0.1",
         "C":"OTHER_VALID_HETEROGENEOUS_RESULTS",
         "scientific_scope":"ONE_PARAMETER_CONDITIONAL_DIAGNOSTIC_NOT_PROJECT_PERFORMANCE_PASS"},
      "input_and_code_bindings_count":len(binding),"bindings":binding,
      "G0_G1_G2E":"FROZEN_UNCHANGED","new_random_draws":0,"new_depth_MC":0,
      "detector_repair":"NOT_OPENED","time_domain":"NOT_OPENED","actual_UUV_spectrum":"NOT_ESTABLISHED",
      "R4_percent":0,"next":"NOT_AUTHORIZED;STOP"}
    E.dump(O/"DESIGN_FREEZE.json",design)
    (O/"SSP_PARAMETER_AND_PROVIDER_CONTRACT.md").write_text("""# SSP物理参数与传播provider合同

唯一参数为SSP_UNIFORM_SOUND_SPEED_OFFSET。每个E-STD水体SSP深度点的声速加同一个delta_c；
诊断点-1/0/+1 m/s。偏移不是有实测依据的环境误差范围，不代表真实海洋条件。
环境输入直接从每个名义缓存对应.env逐行复制，只改251个六列水体SSP行的第二列。
原密度、衰减列、真空/刚性边界、深度5000m、频率、CVWT设置、源深及13个接收深度都不改变。
输出深度必须精确完整，禁止插值；复波数保留复部，按G0 outgoing exp(-ikr)重建复场。

原phase-speed window 1500–1800 m/s不移动。名义最低SSP声速1500，负偏移后1499。
非零CLOW会排除更慢的模态：
[官方KRAKEN说明](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html)。
这是一项预先识别的provider风险，不能隐瞒或自行改窗。正负偏移或两网格模态计数不一致须保留。
本轮要求差分模板内无模态进入/退出，并要求未知源/固定增益剖面后的正负割线一致。
±1的有限差分不自动等价连续局部导数；只有全部认证通过才解释C1e-P0 Fisher信息。
未通过时停止科学正负解释；不缩小步长、不增网格、不追加52次之外的调用。
若以后单独修复provider或导数认证，需要新的授权阶段，不能追溯改写本轮。

信息模型沿用G0确定性未知复源均值与噪声，绝不混用G2E零均值Gaussian协方差族Fisher。
C1b每阵元/频率未知复增益跨时固定，源系数每频率/快照未知。
C1e-P0额外加入一个实SSP切向量，无先验；深度信息剖面全部水平状态。
相对响应使用共同参考引入的完整协方差，并与同源原始信息做一致性和数据处理检查。
自由模态增益、C2零控制及不可辨方向INF保留。该模型仍不是完整C1e或完整RC2支持证书。

计算预算：最多52次新KRAKEN，每次90s，总执行/复核5400s，交付256,000,000 bytes。
先提交推送设计并核对remote/main，名义回归通过后唯一一次执行。
只求解正负偏移；两套名义缓存只读复用。单个provider失败即停止后续求解。
数值准入失败保留完整已得证据，C1e-P0信息列为未评估，不以空值或伪逆输出零误差。
提交结果、推送并核对remote/main后停止。
""",encoding="utf-8")
    # Audit clarification uses only saved original detector outcomes, without modifying G2E.
    rr=E.G.rows(E.ROOT/"results/R4_H3_G2E_CONDITIONAL_EXTRACTION/EXTRACTION_BY_CELL.csv")
    group=[]
    for k in [1,8,32]:
        ss=[r for r in rr if int(r["K"])==k]
        fail=[r for r in ss if r["accepted_quality"]!="True"]
        undet=sum(int(r["failed_blocks"])>0 for r in fail)
        group.append({"K":k,"quality_failed":len(fail),"contains_undetected":undet,
                      "all_detected_direction_error":len(fail)-undet,"new_draws":0})
    assert [(r["quality_failed"],r["contains_undetected"],r["all_detected_direction_error"]) for r in group]==[(979,847,132),(258,258,0),(240,240,0)]
    E.write(O/"G2E_READ_ONLY_CLARIFICATION.csv",group)
    print(json.dumps({"bindings":len(binding),"prepared_environment_files":52,"SSP_rows_changed":len(envchecks),
                      "new_solver_calls":0,"design_sha256":E.sha(O/"DESIGN_FREEZE.json")}))
if __name__=="__main__":main()
