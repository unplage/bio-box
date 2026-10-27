# -*- coding: utf-8 -*-
"""Wet-lab follow-up plan, tailored to the run type."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from .spec import RunResult


def build_wetlab_plan(run: "RunResult", designs: list[dict]) -> list[dict]:
    top = [str(d.get("design_id")) for d in designs[:3] if d.get("design_id")]
    top_text = "、".join(top) if top else "首选候选"
    target = run.target_name

    sections: list[dict] = []

    if run.case == "docking":
        sections.append(
            {
                "title": "先确认结合模式（本流程未提供复合物）",
                "items": [
                    "本批次只固定了靶点与抗体框架的内部几何，抗体-靶点相对朝向由模型采样，"
                    "因此必须先做表位/结合模式验证，不能直接把预测 pose 当作真实结构。",
                    "用已知配体或母本抗体做竞争 ELISA/BLI，确认新设计结合同一表位。",
                    "条件允许时做 HDX-MS 或交联质谱，锁定接触界面。",
                ],
            }
        )
    elif run.case == "apo":
        sections.append(
            {
                "title": "本批次没有靶点（非靶向设计）",
                "items": [
                    "本轮只优化 CDR 序列与局部构象，未针对任何抗原，结果不能直接解读为亲和力提升。",
                    "建议先把靶点胞外域结构或复合物补齐，再跑一次 docking/complex 模式。",
                    "本批候选可用于评估框架可开发性、表达量与 CDR 容忍度。",
                ],
            }
        )

    if len(run.role_to_chain) == 1:
        sections.append(
            {
                "title": "表达与纯化（VHH）",
                "items": [
                    "将 CDR 区按设计序列合成，克隆到 pET-21b/pET-22b（C 端 His6）或 pMESy4 周质表达载体。",
                    "E. coli BL21(DE3) 或 WK6，16–25 °C 诱导过夜；周质抽提后 Ni-NTA 纯化，再经 SEC 去除聚体。",
                    "同时表达母本 VHH 作为对照；每批留样做 SDS-PAGE 与 SEC-HPLC。",
                ],
            }
        )
    else:
        sections.append(
            {
                "title": "表达与纯化（Fab）",
                "items": [
                    "VH/VL 分别克隆到含信号肽的哺乳动物表达载体，HEK293 共转染或稳定池表达。",
                    "培养上清经 Protein L/Protein G 或 His 亲和纯化，SEC 精制；确认 H/L 等摩尔配对。",
                    "母本 Fab 与靶点胞外域同批表达，作为动力学与表位对照。",
                ],
            }
        )

    sections.append(
        {
            "title": "结合与亲和力验证",
            "items": [
                f"优先测试 {top_text}；先做单浓度 ELISA/流式初筛，再做 BLI 或 SPR 测定 kon/koff/KD。",
                f"靶点（{target}）按 His 或 Fc 标签固定，采用单周期动力学；每个候选至少两个独立重复。",
                "设置母本抗体、无关 VHH/Fab、仅缓冲液三组对照；扣除非特异结合。",
                "若靶点为细胞表面蛋白，补做流式细胞术，确认在天然构象下仍结合。",
            ],
        }
    )

    sections.append(
        {
            "title": "表位与功能验证",
            "items": [
                "与母本抗体做竞争结合；若为阻断型抗体，补做功能 assay（配体-受体阻断、信号抑制等）。",
                "对关键 CDR 突变做丙氨酸扫描或回复突变，确认增益来源。",
                "若亲和力提升但功能下降，优先检查是否改变了结合角度而非单纯提高亲和力。",
            ],
        }
    )

    sections.append(
        {
            "title": "稳定性与可开发性",
            "items": [
                "DSF 测 Tm（每候选至少两次），比较母本与设计的热稳定性。",
                "SEC-HPLC 测单体比例；40 °C/4 周、冻融 5 次后复测聚集。",
                "AC-SINS、PSR、HIC 评估非特异结合与疏水性；VHH 额外关注 CDR3 疏水斑块。",
                "若表达量低或聚集明显，优先回退框架/CDR 边界附近的突变，而不是继续叠加亲和力突变。",
            ],
        }
    )

    sections.append(
        {
            "title": "组合与下一轮设计",
            "items": [
                "把排名靠前且突变位点互不重叠的 CDR 突变做组合，构建 2–4 个组合体。",
                "避免同时在 CDR 中引入多个疏水/芳香残基；组合后必须重测稳定性与聚集。",
                "保留至少一个只改 CDR3 的极简版本，便于判断增益是否来自单一热点。",
            ],
        }
    )
    return sections
