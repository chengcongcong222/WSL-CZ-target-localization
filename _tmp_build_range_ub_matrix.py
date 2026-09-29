#!/usr/bin/env python3
"""RANGE-UB-1A: build transferability matrix and claim audit CSVs."""
from pathlib import Path
import pandas as pd

OUT = Path(r"C:\Users\ccc\XiaomiMiMoProjects\WSL系统汇聚区目标定位\results\RANGE_UB\CZ_TDOA_PRIMARY_SOURCE_AUDIT")

rows = [
    dict(
        item="signal_type",
        paper_condition="爆炸声源 (exploding source) = impulsive broadband cooperative",
        E_STD="被动 UUV 窄带线谱 201/235/283/338 Hz",
        difference="主动脉冲 vs 被动窄带连续",
        consequence="REQUIRES_BROADBAND_OR_IMPULSIVE_SIGNAL",
        severity="BLOCKING",
    ),
    dict(
        item="receiver_aperture",
        paper_condition="关键词:大接收深度; 时延差曲线随接收深度变化",
        E_STD="拖曳 HLA 近水平阵, z约200 m, 无垂直孔径",
        difference="需垂直深度采样 vs 水平阵",
        consequence="REQUIRES_VERTICAL_APERTURE",
        severity="BLOCKING",
    ),
    dict(
        item="delay_extraction",
        paper_condition="多途时延差=不同时延峰之间的差值",
        E_STD="窄带线谱无脉冲到达时刻",
        difference="时延峰需脉冲/宽带分辨 vs 窄带连续",
        consequence="TDOA_OBSERVABLE_NOT_AVAILABLE",
        severity="BLOCKING",
    ),
    dict(
        item="source_depth",
        paper_condition="近海面声源 (near sea surface)",
        E_STD="UUV 深度约 200 m",
        difference="海面附近 vs 200 m",
        consequence="DEPTH_MISMATCH_MODERATE",
        severity="WARNING",
    ),
    dict(
        item="range_validated",
        paper_condition="42-52 km (南海实测)",
        E_STD="50-60 km 第一会聚区",
        difference="重叠区 50-52 km; 外推区 52-60 km",
        consequence="PAPER_VALIDATED_RANGE_PARTIAL_OVERLAP",
        severity="WARNING",
    ),
    dict(
        item="CZ_order",
        paper_condition="会聚区(未明确区数); 利用时延差-距离关系",
        E_STD="第一会聚区已知 (RC1 先验)",
        difference="论文未声明需知区号 vs 项目有先验",
        consequence="CZ_ORDER_ASSUMPTION_PROJECT_DEFINED",
        severity="INFO",
    ),
    dict(
        item="bistatic_geometry",
        paper_condition="收发距离(爆炸声源, 收发可能异位)",
        E_STD="单站被动 HLA",
        difference="可能双基/合作 vs 单站被动",
        consequence="GEOMETRY_MISMATCH_NEEDS_PDF",
        severity="WARNING",
    ),
    dict(
        item="sea_depth",
        paper_condition="深海(南海) [NEEDS_PDF 具体海深]",
        E_STD="H=5000 m Munk",
        difference="需核验是否同量级",
        consequence="NEEDS_PDF",
        severity="INFO",
    ),
    dict(
        item="estimation_output",
        paper_condition="距离估计(水平距离)",
        E_STD="需 [r,theta,v,psi] 中的 r",
        difference="论文只出 r vs 项目需全状态",
        consequence="PARTIAL_COMPATIBLE",
        severity="INFO",
    ),
    dict(
        item="error_metric",
        paper_condition="0.6%-6.1% @ 42-52 km",
        E_STD="理论 Demo Gate",
        difference="海试实测 vs 理论 Demo",
        consequence="ORACLE_UPPER_BOUND_REFERENCE",
        severity="INFO",
    ),
]
df = pd.DataFrame(rows)
df.to_csv(OUT / "CZ_TDOA_TRANSFERABILITY_MATRIX.csv", index=False)

claims = [
    dict(
        claim="42-52 km 验证距离范围",
        status="CONFIRMED_BY_PRIMARY_ABSTRACT",
        evidence="官方摘要: '当收发距离在42km~52 km时'",
        page_ref="摘要 [ABSTRACT]",
    ),
    dict(
        claim="0.6%-6.1% 估计误差",
        status="CONFIRMED_BY_PRIMARY_ABSTRACT",
        evidence="官方摘要: '估计误差为0.6%~6.1%'",
        page_ref="摘要 [ABSTRACT]",
    ),
    dict(
        claim="基于虚源理论",
        status="CONFIRMED_BY_PRIMARY_ABSTRACT",
        evidence="官方摘要: '基于虚源理论推导了深海会聚区的多途到达时延结构'",
        page_ref="摘要 [ABSTRACT]",
    ),
    dict(
        claim="时延差随接收深度变化",
        status="CONFIRMED_BY_PRIMARY_ABSTRACT",
        evidence="官方摘要: '时延差曲线随接收深度变化的近似表达式'",
        page_ref="摘要 [ABSTRACT]",
    ),
    dict(
        claim="爆炸声源",
        status="CONFIRMED_BY_PRIMARY_ABSTRACT",
        evidence="官方摘要: '会聚区爆炸声源距离估计'",
        page_ref="摘要 [ABSTRACT]",
    ),
    dict(
        claim="大接收深度",
        status="CONFIRMED_BY_PRIMARY_ABSTRACT",
        evidence="官方关键词列表",
        page_ref="摘要 [ABSTRACT]",
    ),
    dict(
        claim="公式细节/多途数/接收结构",
        status="NEEDS_PDF",
        evidence="摘要未提供",
        page_ref="N/A",
    ),
]
cf = pd.DataFrame(claims)
cf.to_csv(OUT / "CLAIM_AUDIT.csv", index=False)

exp = pd.DataFrame(
    [
        dict(
            condition_type="simulation",
            source_type="爆炸声源",
            source_depth="近海面 [NEEDS_PDF]",
            receiver_type="[NEEDS_PDF]",
            receiver_depth="大接收深度 [NEEDS_PDF]",
            range_km="NEEDS_PDF",
            sea="NEEDS_PDF",
            result="原理验证",
        ),
        dict(
            condition_type="sea_trial",
            source_type="爆炸声源",
            source_depth="近海面 [NEEDS_PDF]",
            receiver_type="[NEEDS_PDF]",
            receiver_depth="大接收深度 [NEEDS_PDF]",
            range_km="42-52",
            sea="南海",
            result="误差 0.6%-6.1%",
        ),
    ]
)
exp.to_csv(OUT / "CZ_TDOA_EXPERIMENT_CONDITIONS.csv", index=False)

print("transferability rows:", len(df))
print("claims:", len(cf))
print("experiments:", len(exp))
print("BLOCKING:", list(df[df.severity == "BLOCKING"]["item"]))
