# -*- coding: utf-8 -*-
"""Turn a NIM response into CSV/JSON/Markdown/HTML deliverables."""

from __future__ import annotations

import csv
import difflib
import html
from pathlib import Path

from .seqprops import (
    aromaticity,
    gravy,
    isoelectric_point,
    liability_motifs,
    molecular_weight,
    net_charge,
)
from .structio import mean_plddt_from_chains, parse_structure_text
from .util import fmt_float, log, now_utc, safe_filename, warn, write_json
from .wetlab import build_wetlab_plan

ROLE_LABELS = {"A": "VHH", "H": "VH", "L": "VL"}


def _flatten(value, prefix: str = "") -> dict:
    flat: dict = {}
    if isinstance(value, dict):
        for key, item in value.items():
            name = f"{prefix}.{key}" if prefix else str(key)
            if isinstance(item, dict):
                flat.update(_flatten(item, name))
            else:
                flat[name] = item
    return flat


def _pass_stats(flat: dict) -> tuple[int, int]:
    total = 0
    passed = 0
    for key, value in flat.items():
        lower = key.lower()
        if lower.endswith(".pass") or lower.endswith("_pass") or lower == "pass":
            if isinstance(value, bool):
                total += 1
                passed += int(value)
            elif isinstance(value, (int, float)) and not isinstance(value, bool):
                total += 1
                passed += int(bool(value))
    return passed, total


def diff_sequences(parent: str, child: str) -> list[dict]:
    matcher = difflib.SequenceMatcher(None, parent, child, autojunk=False)
    mutations: list[dict] = []
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            continue
        if tag == "replace":
            width = min(i2 - i1, j2 - j1)
            for k in range(width):
                mutations.append(
                    {
                        "pos": i1 + k + 1,
                        "from": parent[i1 + k],
                        "to": child[j1 + k],
                        "type": "substitution",
                    }
                )
            if (i2 - i1) > (j2 - j1):
                mutations.append(
                    {
                        "pos": i1 + width + 1,
                        "from": parent[i1 + width:i2],
                        "to": "-",
                        "type": "deletion",
                    }
                )
            elif (j2 - j1) > (i2 - i1):
                mutations.append(
                    {
                        "pos": i1 + width + 1,
                        "from": "-",
                        "to": child[j1 + width:j2],
                        "type": "insertion",
                    }
                )
        elif tag == "delete":
            mutations.append(
                {
                    "pos": i1 + 1,
                    "from": parent[i1:i2],
                    "to": "-",
                    "type": "deletion",
                }
            )
        elif tag == "insert":
            mutations.append(
                {
                    "pos": i1 + 1,
                    "from": "-",
                    "to": child[j1:j2],
                    "type": "insertion",
                }
            )
    return mutations


def _identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    return difflib.SequenceMatcher(None, a, b, autojunk=False).ratio()


def _split_by_lengths(sequence: str, roles: list[str], lengths: list[int]) -> dict[str, str]:
    out = {}
    cursor = 0
    for role, length in zip(roles, lengths):
        out[role] = sequence[cursor:cursor + length]
        cursor += length
    return out


def _reference_to_child(reference: str, child: str) -> dict[int, int]:
    """Map every reference position to the matching child position.

    ``SequenceMatcher.get_matching_blocks()`` only reports the *equal* runs, so
    the substituted positions would be dropped and the extracted antibody would
    be shorter than the parent.  Walking the opcodes keeps the substitutions.
    """
    matcher = difflib.SequenceMatcher(None, reference, child, autojunk=False)
    mapping: dict[int, int] = {}
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            for offset in range(i2 - i1):
                mapping[i1 + offset] = j1 + offset
        elif tag == "replace":
            for offset in range(min(i2 - i1, j2 - j1)):
                mapping[i1 + offset] = j1 + offset
    return mapping


def _extract_antibody_segments(design_sequence: str, run) -> dict[str, str]:
    roles = list(run.sequences)
    lengths = [len(run.sequences[role]) for role in roles]
    total = sum(lengths)

    # A single antibody chain and no target: the design is the antibody.
    if len(roles) == 1 and not run.target_chain_ids:
        return {roles[0]: design_sequence}

    target_sequences = [
        run.chains[cid].sequence() for cid in run.target_chain_ids if cid in run.chains
    ]
    reference = "".join(run.sequences[role] for role in roles) + "".join(target_sequences)

    if len(design_sequence) == total and not target_sequences:
        # The design only covers the antibody: pick the chain order that best
        # matches the parent sequences.
        candidates = [_split_by_lengths(design_sequence, roles, lengths)]
        if len(roles) == 2:
            candidates.append(
                _split_by_lengths(design_sequence, roles[::-1], lengths[::-1])
            )
        return max(
            candidates,
            key=lambda cand: sum(
                _identity(cand[role], run.sequences[role]) for role in roles
            ),
        )

    ref_to_child = _reference_to_child(reference, design_sequence)
    out: dict[str, str] = {}
    cursor = 0
    for role, length in zip(roles, lengths):
        chars = [
            design_sequence[ref_to_child[idx]]
            for idx in range(cursor, cursor + length)
            if idx in ref_to_child
        ]
        if len(chars) >= 0.6 * length:
            out[role] = "".join(chars)
        elif len(roles) == 1:
            out[role] = design_sequence
        else:
            out[role] = run.sequences[role]
        cursor += length
    return out


def _cdr_ranges_by_role(run) -> dict[str, list[tuple[int, int, str]]]:
    out: dict[str, list[tuple[int, int, str]]] = {}
    for meta in run.cdr_metadata:
        role = meta.get("role")
        if not role:
            continue
        out.setdefault(role, []).append(
            (int(meta["start"]), int(meta["end"]), str(meta.get("name", "CDR")))
        )
    return out


def enrich_designs(response: dict, run, out_dir: Path) -> list[dict]:
    designs = response.get("designs") or []
    structures_dir = out_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)
    cdr_ranges = _cdr_ranges_by_role(run)
    enriched: list[dict] = []

    for index, design in enumerate(designs):
        design_id = str(design.get("design_id") or f"design_{index}")
        sequence = str(design.get("sequence") or "")
        segments = _extract_antibody_segments(sequence, run) if sequence else {}

        mutations_by_role: dict[str, list[dict]] = {}
        cdr_mutations: list[dict] = []
        framework_mutations: list[dict] = []
        for role, parent in run.sequences.items():
            child = segments.get(role) or parent
            mutations = diff_sequences(parent, child)
            for mutation in mutations:
                in_cdr = any(
                    start <= mutation["pos"] <= end
                    for start, end, _ in cdr_ranges.get(role, [])
                )
                mutation["region"] = "CDR" if in_cdr else "framework"
                mutation["role"] = role
                if in_cdr:
                    cdr_mutations.append(mutation)
                else:
                    framework_mutations.append(mutation)
            mutations_by_role[role] = mutations

        # Structure / pLDDT
        struct = design.get("structure") or {}
        content = struct.get("content") if isinstance(struct, dict) else None
        structure_path = None
        plddt = None
        if content and isinstance(content, str) and len(content) > 200:
            fmt = str(struct.get("format") or "mmcif").lower()
            suffix = ".pdb" if fmt == "pdb" else ".cif"
            filename = str(struct.get("filename") or f"{design_id}{suffix}")
            if not filename.lower().endswith((".cif", ".mmcif", ".pdb")):
                filename += suffix
            path = structures_dir / safe_filename(filename)
            path.write_text(content, encoding="utf-8")
            structure_path = str(path)
            try:
                chains = parse_structure_text(content)
                plddt = mean_plddt_from_chains(chains)
            except Exception as exc:  # structure is optional for the report
                warn(f"Could not parse structure for {design_id}: {exc}")

        antibody_seq = "".join(
            str(segments.get(role, run.sequences[role])) for role in run.sequences
        )
        liabilities = liability_motifs(antibody_seq)
        flat = _flatten(design.get("metrics") or {})
        passed, total = _pass_stats(flat)

        enriched.append(
            {
                "design_id": design_id,
                "replicate_index": design.get("replicate_index", index),
                "seed": design.get("seed"),
                "sequence": sequence,
                "antibody_sequence": antibody_seq,
                "antibody_sequences": segments,
                "mutations": mutations_by_role,
                "cdr_mutations": cdr_mutations,
                "framework_mutations": framework_mutations,
                "n_cdr_mutations": len(cdr_mutations),
                "n_framework_mutations": len(framework_mutations),
                "plddt": plddt,
                "structure_path": structure_path,
                "metrics": flat,
                "pass_passed": passed,
                "pass_total": total,
                "liabilities": liabilities,
                "n_liabilities": len(liabilities),
                "biophys": {
                    "length": len(antibody_seq),
                    "mw_kda": round(molecular_weight(antibody_seq) / 1000.0, 2),
                    "pI": round(isoelectric_point(antibody_seq), 2),
                    "net_charge_pH7.4": round(net_charge(antibody_seq), 2),
                    "gravy": round(gravy(antibody_seq), 3),
                    "cysteines": antibody_seq.count("C"),
                    "aromaticity": round(aromaticity(antibody_seq), 3),
                },
                "diagnostics": design.get("diagnostics") or {},
            }
        )

    # Heuristic ranking: pass criteria first, then pLDDT, then fewer liabilities.
    def sort_key(item: dict):
        frac = item["pass_passed"] / item["pass_total"] if item["pass_total"] else 0.0
        return (
            frac,
            item["plddt"] if item["plddt"] is not None else 0.0,
            -item["n_liabilities"],
        )

    enriched.sort(key=sort_key, reverse=True)
    for rank, item in enumerate(enriched, start=1):
        item["rank"] = rank
    return enriched


def _md_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("|" + "|".join(["---"] * len(headers)) + "|")
    for row in rows:
        lines.append("| " + " | ".join(str(cell) for cell in row) + " |")
    return "\n".join(lines)


def _mutation_text(mutations: list[dict]) -> str:
    if not mutations:
        return "无"
    roles = {mutation.get("role") for mutation in mutations}
    show_role = len(roles) > 1
    parts = []
    for mutation in mutations:
        prefix = f"{mutation['role']}:" if show_role else ""
        if mutation["type"] == "substitution":
            parts.append(f"{prefix}{mutation['from']}{mutation['pos']}{mutation['to']}")
        else:
            parts.append(
                f"{prefix}{mutation['pos']}{mutation['type'][:3]}"
                f"({mutation['from']}→{mutation['to']})"
            )
    return ", ".join(parts)


def _case_label(case: str) -> str:
    return {
        "complex": "A. 提供了抗体–靶点复合物（天然朝向）",
        "docking": "B. 有靶点结构但无复合物（相对朝向由模型采样）",
        "apo": "C. 仅抗体、无靶点（非靶向设计）",
    }.get(case, case)


def build_report(run, response: dict, designs: list[dict]) -> tuple[str, str]:
    model = response.get("model") or {}
    runtime = response.get("runtime") or {}
    warnings = response.get("warnings") or []
    aggregated = response.get("aggregated_metrics") or []
    total_residues = sum(len(chain.residues) for chain in run.chains.values())
    top = designs[:3]
    top_names = "、".join(d["design_id"] for d in top) if top else "（无）"

    lines: list[str] = []
    lines.append("# RFD4-Proteina 抗体亲和力成熟报告")
    lines.append("")
    lines.append("## 1. 运行摘要")
    lines.append("")
    lines.append(
        _md_table(
            ["项目", "内容"],
            [
                ["靶点", run.target_name],
                ["设计模式", "VHH / 纳米抗体" if run.mode == "vhh" else "Fab"],
                ["输入情形", _case_label(run.case)],
                ["抗体链", ", ".join(f"{cid}（{ROLE_LABELS.get(role, role)}）" for role, cid in run.role_to_chain.items())],
                ["靶点链", ", ".join(run.target_chain_ids) or "无"],
                ["体系残基数", str(total_residues)],
                ["请求设计数", str(run.request["sampling"].get("num_replicates", "-"))],
                ["返回设计数", str(len(designs))],
                ["模型", str(model.get("model_name") or model.get("checkpoint") or "rfd4-proteina")],
                ["NIM 版本", str(runtime.get("nim_version") or "-")],
                ["GPU", str(runtime.get("gpu") or "-")],
                ["请求 ID", str(response.get("request_id") or "-")],
                ["耗时", f"{response.get('duration_ms', 0) / 1000:.1f} s" if response.get("duration_ms") else "-"],
                ["报告生成时间 (UTC)", now_utc()],
            ],
        )
    )
    lines.append("")
    lines.append("**结论先行：**")
    if designs:
        best = designs[0]
        lines.append(
            f"- 启发式首推 **{best['design_id']}**：CDR 突变 {best['n_cdr_mutations']} 个，"
            f"框架突变 {best['n_framework_mutations']} 个，"
            f"pLDDT {fmt_float(best['plddt'], 1)}，"
            f"预测 pI {best['biophys']['pI']}，GRAVY {best['biophys']['gravy']}。"
        )
        lines.append(
            f"- 前三位为 {top_names}；优先做小规模表达 + BLI/SPR 初筛，再决定是否组合突变。"
        )
        if run.case == "docking":
            lines.append(
                "- 本批没有实验复合物，**相对朝向是模型采样的**；在验证表位之前，"
                "亲和力数值只能作为排序参考。"
            )
        if run.case == "apo":
            lines.append("- 本批没有靶点，结果只反映 CDR/框架可开发性，不能当作亲和力提升证据。")
    else:
        lines.append("- 响应中没有 design，请检查 warnings 与请求参数。")

    lines.append("")
    lines.append("## 2. 设计空间与输入")
    lines.append("")
    lines.append("### 2.1 链")
    lines.append("")
    lines.append(
        _md_table(
            ["链", "角色", "残基数", "序列"],
            [
                [
                    chain["id"],
                    "抗体" if chain["id"] in set(run.role_to_chain.values()) else "靶点",
                    str(chain["n_residues"]),
                    f"`{chain['sequence'][:80]}{'…' if len(chain['sequence']) > 80 else ''}`",
                ]
                for chain in run.chain_table
            ],
        )
    )
    lines.append("")
    lines.append("### 2.2 CDR 选择")
    lines.append("")
    lines.append(
        _md_table(
            ["CDR", "链", "起始", "结束", "长度", "序列", "识别方式"],
            [
                [
                    meta.get("name", "-"),
                    meta.get("role", "-"),
                    str(meta.get("start", "-")),
                    str(meta.get("end", "-")),
                    str(meta.get("length", "-")),
                    f"`{meta.get('sequence', '')}`",
                    meta.get("method", "-"),
                ]
                for meta in run.cdr_metadata
            ],
        )
    )
    lines.append("")
    lines.append(
        "CDR 由 ANARCI 按 IMGT/Kabat/Chothia 编号识别时最可靠；若使用启发式边界，"
        "请在湿实验前用母本结构核对。"
    )
    hotspots = {
        cid: sorted(plan.hotspot_residues)
        for cid, plan in run.plans.items()
        if plan.hotspot_residues
    }
    if hotspots:
        lines.append("")
        lines.append("### 2.3 表位热点（C_HOT）")
        lines.append("")
        for cid, residues in hotspots.items():
            lines.append(f"- 链 {cid}: {', '.join(str(r) for r in residues)}")
    lines.append("")
    lines.append("### 2.4 条件设置")
    lines.append("")
    if run.case == "complex":
        lines.append(
            "- 抗体框架：**主链绝对固定**（C_CRD），侧链自由重排；靶点：**全原子绝对固定**；序列固定。"
        )
    elif run.case == "docking":
        lines.append(
            "- 靶点：**全原子绝对固定**；抗体框架：**内部坐标组 1**（保持框架内部几何，"
            "整体朝向由模型采样）。"
        )
    else:
        lines.append("- 仅抗体框架：**内部坐标组 1**，无靶点约束。")
    lines.append(
        "- 选中的 CDR：**不施加坐标/序列约束**，序列与局部构象同时重设计。"
    )
    if not run.request.get("sampling", {}).get("num_steps"):
        lines.append(
            f"- 采样：`{run.request['sampling'].get('recycle_strategy')}` 预设"
            "（60 步、185 次 trunk pass）。"
        )
    else:
        lines.append(
            f"- 采样：`{run.request['sampling'].get('recycle_strategy')}`，"
            f"{run.request['sampling'].get('num_steps')} 步。"
        )
    lines.append("")

    lines.append("## 3. 候选设计总表")
    lines.append("")
    lines.append(
        _md_table(
            [
                "排名", "设计", "CDR 突变", "框架突变", "pLDDT",
                "pI", "GRAVY", "净电荷", "风险位点", "通过指标",
            ],
            [
                [
                    str(d["rank"]),
                    d["design_id"],
                    str(d["n_cdr_mutations"]),
                    str(d["n_framework_mutations"]),
                    fmt_float(d["plddt"], 1),
                    str(d["biophys"]["pI"]),
                    str(d["biophys"]["gravy"]),
                    str(d["biophys"]["net_charge_pH7.4"]),
                    str(d["n_liabilities"]),
                    f"{d['pass_passed']}/{d['pass_total']}" if d["pass_total"] else "—",
                ]
                for d in designs
            ],
        )
    )
    lines.append("")
    lines.append(
        "> 排名是**启发式**的：先看通过标准，再看 pLDDT，最后看序列风险位点；"
        "它不等于真实亲和力排序。"
    )
    lines.append("")
    lines.append(
        "> pLDDT 一列取自返回结构的 CA 原子 B 因子；若服务端不在 B 因子里写 pLDDT，"
        "该列会显示输入结构的 B 因子均值，仅作参考。"
    )
    lines.append("")

    lines.append("## 4. 逐候选详情")
    lines.append("")
    for design in designs:
        lines.append(f"### 4.{design['rank']} {design['design_id']}")
        lines.append("")
        antibody_seq = design.get("antibody_sequence") or "-"
        rows = [
            ["seed", str(design.get("seed") or "-")],
            ["抗体序列", f"`{antibody_seq}`"],
        ]
        if design.get("sequence") and design["sequence"] != antibody_seq:
            rows.append(["完整设计序列（含靶点链）", f"`{design['sequence']}`"])
        rows.extend(
            [
                ["CDR 突变", _mutation_text(design["cdr_mutations"])],
                ["框架突变", _mutation_text(design["framework_mutations"])],
                ["pLDDT（B 因子均值）", fmt_float(design["plddt"], 1)],
                ["分子量 (kDa)", str(design["biophys"]["mw_kda"])],
                ["pI", str(design["biophys"]["pI"])],
                ["净电荷 (pH 7.4)", str(design["biophys"]["net_charge_pH7.4"])],
                ["GRAVY", str(design["biophys"]["gravy"])],
                ["半胱氨酸数", str(design["biophys"]["cysteines"])],
                ["结构文件", design.get("structure_path") or "（响应未含结构）"],
            ]
        )
        lines.append(_md_table(["属性", "值"], rows))
        if design["liabilities"]:
            lines.append("")
            lines.append("**序列风险位点：**")
            for hit in design["liabilities"][:30]:
                lines.append(
                    f"- {hit['type']} @ {hit['position']} `{hit['motif']}`：{hit['note']}"
                )
        if design["metrics"]:
            lines.append("")
            lines.append("**NIM 指标：**")
            rows = []
            for key, value in sorted(design["metrics"].items()):
                rows.append([key, str(value)])
            lines.append(_md_table(["指标", "值"], rows))
        lines.append("")

    lines.append("## 5. 聚合指标")
    lines.append("")
    if aggregated:
        lines.append(
            _md_table(
                ["指标", "汇总方式", "值", "样本数"],
                [
                    [
                        item.get("metric", "-"),
                        item.get("reduction", "-"),
                        str(item.get("value", "-")),
                        str(item.get("n", "-")),
                    ]
                    for item in aggregated
                ],
            )
        )
    else:
        lines.append("响应未返回聚合指标（可能未请求相应 metric tag）。")
    lines.append("")

    lines.append("## 6. 湿实验建议")
    lines.append("")
    for section in build_wetlab_plan(run, designs):
        lines.append(f"### {section['title']}")
        lines.append("")
        for item in section["items"]:
            lines.append(f"- {item}")
        lines.append("")

    lines.append("## 7. 局限与风险")
    lines.append("")
    lines.append("- 计算预测不等于实测亲和力；必须做 SPR/BLI 与功能验证。")
    lines.append("- 若输入情形为 B/C，抗体–靶点相对朝向或靶点本身缺失，结合模式不确定性最大。")
    lines.append("- RFD4 可能在不经意间改变 CDR 长度或框架边界；本报告已统计框架突变，"
                 "建议只推进框架突变少、CDR 突变可解释的候选。")
    lines.append("- 序列风险位点只基于基序规则，未做结构级聚集/免疫原性预测。")
    if warnings:
        lines.append("")
        lines.append("**服务端 warnings：**")
        for item in warnings:
            lines.append(f"- {item}")
    lines.append("")

    lines.append("## 8. 文件清单")
    lines.append("")
    lines.append("- `request.json`：发送给 NIM 的请求（API key 不落盘）")
    lines.append("- `response.json`：原始响应")
    lines.append("- `input_complex.cif`：带条件的输入结构")
    lines.append("- `designs.csv` / `designs.json`：候选汇总")
    lines.append("- `structures/`：返回的结构")
    lines.append("- `report.md` / `report.html`：本报告")
    lines.append("")

    markdown = "\n".join(lines)
    return markdown, _render_html(markdown, run, designs, response)


def _render_html(markdown: str, run, designs: list[dict], response: dict) -> str:
    """Small, dependency-free HTML renderer for the generated Markdown."""
    body = _markdown_to_html(markdown)
    summary_cards = [
        ("输入情形", _case_label(run.case)),
        ("候选数", str(len(designs))),
        ("首推", designs[0]["design_id"] if designs else "—"),
        (
            "首推 pLDDT",
            fmt_float(designs[0]["plddt"], 1) if designs else "—",
        ),
    ]
    cards = "".join(
        f'<div class="card"><div class="card-label">{html.escape(label)}</div>'
        f'<div class="card-value">{html.escape(value)}</div></div>'
        for label, value in summary_cards
    )
    return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>RFD4-Proteina 抗体亲和力成熟报告</title>
<style>
:root {{ color-scheme: light dark; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Noto Sans CJK SC", sans-serif;
  max-width: 1180px; margin: 0 auto; padding: 24px 28px 80px; line-height: 1.65;
  background: #fbfbfa; color: #1f2328; }}
@media (prefers-color-scheme: dark) {{ body {{ background:#16181d; color:#e6e6e6; }} }}
h1 {{ font-size: 1.7rem; border-bottom: 2px solid #2f81f7; padding-bottom: .4rem; }}
h2 {{ margin-top: 2rem; border-bottom: 1px solid #d0d7de; padding-bottom: .25rem; }}
h3 {{ margin-top: 1.35rem; }}
table {{ border-collapse: collapse; width: 100%; margin: .8rem 0 1.1rem; font-size: .92rem; }}
th, td {{ border: 1px solid #d0d7de; padding: 6px 9px; text-align: left; vertical-align: top; }}
th {{ background: #eef3f8; }}
tr:nth-child(even) td {{ background: rgba(127,127,127,.045); }}
code {{ background: rgba(127,127,127,.14); padding: 1px 5px; border-radius: 4px;
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: .88em; }}
blockquote {{ border-left: 4px solid #d4a72c; margin: 1rem 0; padding: .4rem 1rem;
  background: rgba(212,167,44,.09); }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
  gap: 12px; margin: 1.2rem 0; }}
.card {{ border: 1px solid #d0d7de; border-radius: 10px; padding: 12px 14px;
  background: rgba(47,129,247,.06); }}
.card-label {{ font-size: .78rem; opacity: .7; }}
.card-value {{ font-size: 1.02rem; font-weight: 650; margin-top: 3px; }}
ul {{ padding-left: 1.35rem; }}
li {{ margin: .2rem 0; }}
.meta {{ font-size: .85rem; opacity: .75; }}
</style>
</head>
<body>
<h1>RFD4-Proteina 抗体亲和力成熟报告</h1>
<div class="cards">{cards}</div>
<p class="meta">模型：{html.escape(str((response.get('model') or {}).get('model_name') or 'rfd4-proteina'))}
 · 请求 ID：{html.escape(str(response.get('request_id') or '-'))}
 · 靶点：{html.escape(run.target_name)}</p>
{body}
</body>
</html>
"""


def _markdown_to_html(markdown: str) -> str:
    """Tiny Markdown subset renderer (headings, tables, lists, quotes, code)."""
    out: list[str] = []
    lines = markdown.splitlines()
    index = 0
    while index < len(lines):
        line = lines[index]
        if not line.strip():
            index += 1
            continue
        if line.startswith("# "):
            index += 1
            continue  # title is rendered separately
        if line.startswith("## "):
            out.append(f"<h2>{html.escape(line[3:])}</h2>")
            index += 1
            continue
        if line.startswith("### "):
            out.append(f"<h3>{html.escape(line[4:])}</h3>")
            index += 1
            continue
        if line.startswith("|"):
            rows = []
            while index < len(lines) and lines[index].startswith("|"):
                rows.append(lines[index])
                index += 1
            out.append(_html_table(rows))
            continue
        if line.startswith(">"):
            quote = []
            while index < len(lines) and lines[index].startswith(">"):
                quote.append(lines[index].lstrip("> ").strip())
                index += 1
            out.append(f"<blockquote>{_inline(' '.join(quote))}</blockquote>")
            continue
        if line.startswith("- "):
            items = []
            while index < len(lines) and lines[index].startswith("- "):
                items.append(f"<li>{_inline(lines[index][2:])}</li>")
                index += 1
            out.append("<ul>" + "".join(items) + "</ul>")
            continue
        paragraph = [line]
        index += 1
        while index < len(lines) and lines[index].strip() and not lines[index].startswith(("#", "|", "- ", ">")):
            paragraph.append(lines[index])
            index += 1
        out.append(f"<p>{_inline(' '.join(paragraph))}</p>")
    return "\n".join(out)


def _inline(text: str) -> str:
    escaped = html.escape(text)
    parts = escaped.split("`")
    rendered = []
    for idx, part in enumerate(parts):
        rendered.append(f"<code>{part}</code>" if idx % 2 else part)
    text = "".join(rendered)
    # bold
    text = text.replace("**", "\u0000")
    pieces = text.split("\u0000")
    text = "".join(
        f"<strong>{piece}</strong>" if idx % 2 else piece for idx, piece in enumerate(pieces)
    )
    return text


def _html_table(lines: list[str]) -> str:
    rows = []
    for line in lines:
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        rows.append(cells)
    if len(rows) >= 2 and all(set(cell) <= set("-: ") for cell in rows[1]):
        header = rows[0]
        body = rows[2:]
    else:
        header = []
        body = rows
    out = ["<table>"]
    if header:
        out.append("<thead><tr>" + "".join(f"<th>{_inline(c)}</th>" for c in header) + "</tr></thead>")
    out.append("<tbody>")
    for row in body:
        out.append("<tr>" + "".join(f"<td>{_inline(c)}</td>" for c in row) + "</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def write_outputs(run, response: dict, out_dir: Path, no_html: bool = False) -> list[dict]:
    designs = enrich_designs(response, run, out_dir)
    write_json(out_dir / "designs.json", designs)

    csv_path = out_dir / "designs.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "rank", "design_id", "replicate_index", "seed",
                "cdr_mutations", "framework_mutations", "plddt",
                "pI", "net_charge_pH7.4", "gravy", "liabilities",
                "pass_passed", "pass_total",
                "antibody_sequence", "design_sequence",
            ]
        )
        for design in designs:
            writer.writerow(
                [
                    design["rank"], design["design_id"], design["replicate_index"],
                    design.get("seed"), _mutation_text(design["cdr_mutations"]),
                    _mutation_text(design["framework_mutations"]),
                    fmt_float(design["plddt"], 2), design["biophys"]["pI"],
                    design["biophys"]["net_charge_pH7.4"], design["biophys"]["gravy"],
                    design["n_liabilities"], design["pass_passed"],
                    design["pass_total"], design.get("antibody_sequence"),
                    design.get("sequence"),
                ]
            )

    markdown, html_text = build_report(run, response, designs)
    (out_dir / "report.md").write_text(markdown, encoding="utf-8")
    if not no_html:
        (out_dir / "report.html").write_text(html_text, encoding="utf-8")
    log(f"Wrote {len(designs)} designs, report.md" + ("" if no_html else " and report.html"))
    return designs
