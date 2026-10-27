# -*- coding: utf-8 -*-
"""Command-line entry point."""

from __future__ import annotations

import argparse
import json
import os
import traceback
from pathlib import Path

from . import __version__
from .nim import load_replay, run_design
from .report import write_outputs
from .spec import RunConfig, build_run
from .util import (
    PipelineError,
    log,
    set_log_file,
    set_verbose,
    slugify,
    timestamp_slug,
    warn,
    write_json,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rfd4affinity",
        description=(
            "用 NVIDIA RFD4-Proteina NIM 做 VHH/Fab 抗体亲和力成熟："
            "识别 CDR、构建条件结构、提交设计并生成报告。"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例：\n"
            "  # 有复合物结构（最可靠）\n"
            "  python run.py --complex-structure complex.pdb --antibody-chains H,L \\\n"
            "      --target-chains T --epitope T:557-573,T:591-605 --cdrs HCDR3\n\n"
            "  # 只有 VHH 序列 + 靶点结构\n"
            "  python run.py --vhh-seq EVQL... --target-structure 8UCD.pdb \\\n"
            "      --target-chains A --cdrs all --num-designs 8\n\n"
            "  # 只构建请求，不调用 API\n"
            "  python run.py --complex-structure complex.pdb --antibody-chains H,L \\\n"
            "      --dry-run\n"
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--mode", choices=["vhh", "fab"], default="vhh", help="抗体形式（默认 vhh）")

    seq = parser.add_argument_group("序列输入")
    seq.add_argument("--vhh-seq", help="VHH 氨基酸序列或 FASTA 文件（多记录用 文件#记录名）")
    seq.add_argument("--vh-seq", help="Fab 重链（VH）序列或 FASTA 文件")
    seq.add_argument("--vl-seq", help="Fab 轻链（VL）序列或 FASTA 文件")

    structures = parser.add_argument_group("结构输入")
    structures.add_argument("--complex-structure", help="抗体–靶点复合物（PDB/mmCIF）")
    structures.add_argument("--antibody-structure", help="仅抗体结构（PDB/mmCIF）")
    structures.add_argument("--antibody-chains", help="抗体链，逗号分隔；Fab 为 重,轻")
    structures.add_argument("--target-structure", help="靶点结构路径或 4 位 PDB ID")
    structures.add_argument("--target-chains", help="靶点链，逗号分隔（默认全部）")
    structures.add_argument("--target-name", default="target", help="靶点名称（写入报告）")

    cdr = parser.add_argument_group("CDR 选择")
    cdr.add_argument("--cdrs", default="all", help="要设计的 CDR：all、CDR3、HCDR1,LCDR3 等")
    cdr.add_argument(
        "--cdr-scheme",
        choices=["imgt", "kabat", "chothia", "martin"],
        default="imgt",
        help="ANARCI 编号方案（默认 imgt）",
    )
    cdr.add_argument(
        "--cdr-ranges",
        help="手动指定序列坐标，如 H:99-109,L:89-97（覆盖自动识别）",
    )
    cdr.add_argument(
        "--epitope",
        help="靶点热点残基（C_HOT），如 557-573,591-605 或 T:557-573",
    )

    design = parser.add_argument_group("设计参数")
    design.add_argument("--fold", choices=["auto", "esmfold", "none"], default="auto",
                        help="无抗体结构时是否用 ESMFold 折叠 VHH（默认 auto）")
    design.add_argument("--offset", type=float, default=60.0,
                        help="docking 模式下靶点相对抗体的初始质心距离/Å（默认 60）")
    design.add_argument("--fix-sidechains", action="store_true",
                        help="抗体侧链也固定（靶点始终全原子固定；默认抗体只固定主链、侧链重排）")
    design.add_argument("--num-designs", type=int, default=4, help="设计数（默认 4）")
    design.add_argument(
        "--num-steps",
        type=int,
        default=None,
        help="扩散步数（1–200，服务端上限；默认用 recycle 预设）",
    )
    design.add_argument(
        "--recycle",
        choices=["thinking_low", "thinking_medium", "thinking_high"],
        default="thinking_medium",
        help="recycle 策略（默认 thinking_medium）",
    )
    design.add_argument("--seed", type=int, default=None, help="随机种子")
    design.add_argument("--metrics-tags", help="逗号分隔的 metric tag，覆盖默认组合")
    design.add_argument("--example-id", help="写入 example_id（默认 target_mode）")

    run = parser.add_argument_group("运行")
    run.add_argument("--out-dir", help="输出目录（默认 output/<时间>_<模式>）")
    run.add_argument("--api-key", help="NVIDIA API key；默认读 NVIDIA_API_KEY")
    run.add_argument("--base-url", default="https://health.api.nvidia.com/v1")
    run.add_argument("--timeout", type=float, default=1800.0, help="HTTP 超时/秒")
    run.add_argument("--dry-run", action="store_true", help="只构建并保存请求，不调用 API")
    run.add_argument("--replay", help="复用已有的 response.json，不调用 API")
    run.add_argument("--no-html", action="store_true", help="不生成 report.html")
    run.add_argument("--verbose", action="store_true", help="输出调试日志")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    set_verbose(args.verbose)

    out_dir = Path(
        args.out_dir
        or Path("output") / f"{timestamp_slug()}_{slugify(args.target_name)}_{args.mode}"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    set_log_file(out_dir / "run.log")

    api_key = args.api_key or os.environ.get("NVIDIA_API_KEY")
    if api_key:
        log("NVIDIA API key detected")
    else:
        log("No NVIDIA API key found (set NVIDIA_API_KEY or pass --api-key)")
    cfg = RunConfig(
        mode=args.mode,
        vhh_seq=args.vhh_seq,
        vh_seq=args.vh_seq,
        vl_seq=args.vl_seq,
        antibody_structure=args.antibody_structure,
        complex_structure=args.complex_structure,
        antibody_chains=args.antibody_chains,
        target_structure=args.target_structure,
        target_chains=args.target_chains,
        target_name=args.target_name,
        cdrs=args.cdrs,
        cdr_scheme=args.cdr_scheme,
        cdr_ranges=args.cdr_ranges,
        epitope=args.epitope,
        fold=args.fold,
        offset=args.offset,
        fix_sidechains=args.fix_sidechains,
        num_designs=args.num_designs,
        num_steps=args.num_steps,
        recycle=args.recycle,
        seed=args.seed,
        metrics_tags=[t for t in (args.metrics_tags or "").split(",") if t.strip()] or None,
        example_id=args.example_id,
        base_url=args.base_url,
    )

    try:
        run_result = build_run(cfg, api_key, out_dir)
        request_bytes = len(json.dumps(run_result.request).encode("utf-8"))
        log(f"Request prepared: {request_bytes / 1e3:.1f} kB")

        if args.dry_run:
            log("Dry run: request.json and input_complex.cif written; API not called.")
            print(json.dumps(
                {
                    "out_dir": str(out_dir),
                    "case": run_result.case,
                    "antibody_chains": run_result.role_to_chain,
                    "target_chains": run_result.target_chain_ids,
                    "request_kb": round(request_bytes / 1000, 1),
                    "dry_run": True,
                },
                ensure_ascii=False,
                indent=2,
            ))
            return 0

        if args.replay:
            log(f"Replaying response from {args.replay}")
            response = load_replay(Path(args.replay))
            write_json(out_dir / "response.json", response)
        else:
            if not api_key:
                raise PipelineError(
                    "NVIDIA_API_KEY is required to call the hosted NIM. "
                    "Use --dry-run to only build the request, or --replay with a saved response."
                )
            response = run_design(
                run_result.request,
                api_key=api_key,
                base_url=args.base_url,
                timeout=args.timeout,
            )
            write_json(out_dir / "response.json", response)

        designs = write_outputs(run_result, response, out_dir, no_html=args.no_html)
        print(json.dumps(
            {
                "out_dir": str(out_dir),
                "case": run_result.case,
                "designs": len(designs),
                "top": [d["design_id"] for d in designs[:3]],
                "report": str(out_dir / "report.md"),
                "html": None if args.no_html else str(out_dir / "report.html"),
            },
            ensure_ascii=False,
            indent=2,
        ))
        return 0
    except PipelineError as exc:
        log(str(exc), "ERROR")
        warn("Run aborted; see run.log for details.")
        return 2
    except KeyboardInterrupt:
        log("Interrupted", "ERROR")
        return 130
    except Exception as exc:  # noqa: BLE001 - CLI boundary
        log(f"Unexpected failure: {type(exc).__name__}: {exc}", "ERROR")
        log(traceback.format_exc(), "ERROR")
        warn("Run aborted with an unexpected error; the traceback is above and in run.log.")
        return 3


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
