"""NVIDIA Bio Tools - Unified CLI for ProteinMPNN and Boltz-2."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent))

from proteinmpnn import ProteinMPNNClient, predict_sequences
from boltz2 import Boltz2Client, predict_structure


def build_proteinmpnn_parser(subparsers) -> argparse.ArgumentParser:
    """Build ProteinMPNN subcommand parser."""
    parser = subparsers.add_parser(
        "proteinmpnn",
        help="ProteinMPNN: Predict amino acid sequences from protein backbone",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "pdb",
        help="Path to input PDB file",
    )
    parser.add_argument(
        "--api-key",
        help="NVIDIA API key (or set NVIDIA_API_KEY env var)",
    )
    parser.add_argument(
        "--base-url",
        help="Custom base URL for NIM endpoint",
    )
    parser.add_argument(
        "--self-hosted",
        action="store_true",
        help="Use self-hosted NIM at http://localhost:8000",
    )
    parser.add_argument(
        "--ca-only",
        action="store_true",
        help="Use only CA atoms (faster, less accurate)",
    )
    parser.add_argument(
        "--soluble-model",
        action="store_true",
        help="Use soluble protein model variant",
    )
    parser.add_argument(
        "--sampling-temp",
        type=float,
        nargs="+",
        default=[0.1],
        help="Sampling temperature(s)",
    )
    parser.add_argument(
        "--max-lines",
        type=int,
        default=200,
        help="Max ATOM lines to use from PDB",
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output FASTA file path",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print request without sending",
    )
    parser.set_defaults(func=run_proteinmpnn)
    return parser


def build_boltz2_parser(subparsers) -> argparse.ArgumentParser:
    """Build Boltz-2 subcommand parser."""
    parser = subparsers.add_parser(
        "boltz2",
        help="Boltz-2: Predict protein structures",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input group
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--sequence",
        help="Protein sequence (single chain)",
    )
    input_group.add_argument(
        "--chains",
        help="JSON file with chain definitions for multi-chain/ligand",
    )

    parser.add_argument(
        "--api-key",
        help="NVIDIA API key (or set NVIDIA_API_KEY env var)",
    )
    parser.add_argument(
        "--base-url",
        help="Custom base URL for NIM endpoint",
    )
    parser.add_argument(
        "--self-hosted",
        action="store_true",
        help="Use self-hosted NIM at http://localhost:8000",
    )
    parser.add_argument(
        "--chain-id",
        default="A",
        help="Chain ID for single sequence",
    )
    parser.add_argument(
        "--recycling-steps",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--sampling-steps",
        type=int,
        default=50,
    )
    parser.add_argument(
        "--diffusion-samples",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--step-scale",
        type=float,
        default=1.638,
    )
    parser.add_argument(
        "--output-format",
        choices=["mmcif", "pdb"],
        default="mmcif",
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output directory or file",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print request without sending",
    )
    parser.set_defaults(func=run_boltz2)
    return parser


def run_proteinmpnn(args) -> int:
    """Run ProteinMPNN prediction."""
    try:
        if args.dry_run:
            # Build and print request without sending
            path = Path(args.pdb)
            if not path.exists():
                print(f"Error: PDB file not found: {path}", file=sys.stderr)
                return 1
            content = path.read_text(encoding="utf-8")
            lines = content.strip().split("\n")
            atom_lines = [line for line in lines if line.startswith("ATOM")]
            if len(atom_lines) > args.max_lines:
                kept = []
                count = 0
                for line in content.strip().split("\n"):
                    if line.startswith("ATOM"):
                        if count >= args.max_lines:
                            continue
                        count += 1
                    kept.append(line)
                pdb_content = "\n".join(kept)
            else:
                pdb_content = "\n".join(atom_lines)

            payload = {
                "input_pdb": pdb_content,
                "ca_only": args.ca_only,
                "use_soluble_model": args.soluble_model,
                "sampling_temp": args.sampling_temp,
            }
            print(json.dumps(payload, indent=2))
            return 0

        result = predict_sequences(
            pdb_path=args.pdb,
            api_key=args.api_key,
            base_url=args.base_url,
            use_self_hosted=args.self_hosted,
            ca_only=args.ca_only,
            use_soluble_model=args.soluble_model,
            sampling_temp=args.sampling_temp,
            max_atom_lines=args.max_lines,
        )

        if args.output:
            Path(args.output).write_text(result.get("mfasta", ""), encoding="utf-8")
            print(f"Saved to {args.output}")
        else:
            print(result.get("mfasta", ""))

        return 0

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


def run_boltz2(args) -> int:
    """Run Boltz-2 prediction."""
    try:
        if args.chains:
            # Multi-chain from JSON file
            with open(args.chains) as f:
                chains = json.load(f)
            if args.dry_run:
                payload = {
                    "polymers": chains,
                    "recycling_steps": args.recycling_steps,
                    "sampling_steps": args.sampling_steps,
                    "diffusion_samples": args.diffusion_samples,
                    "step_scale": args.step_scale,
                    "output_format": args.output_format,
                }
                print(json.dumps(payload, indent=2))
                return 0
            result = predict_structure(
                chains=chains,
                api_key=args.api_key,
                base_url=args.base_url,
                use_self_hosted=args.self_hosted,
                recycling_steps=args.recycling_steps,
                sampling_steps=args.sampling_steps,
                diffusion_samples=args.diffusion_samples,
                step_scale=args.step_scale,
                output_format=args.output_format,
            )
        else:
            if not args.sequence:
                print("Error: --sequence or --chains required", file=sys.stderr)
                return 1
            if args.dry_run:
                payload = {
                    "polymers": [{"id": args.chain_id, "molecule_type": "protein", "sequence": args.sequence}],
                    "recycling_steps": args.recycling_steps,
                    "sampling_steps": args.sampling_steps,
                    "diffusion_samples": args.diffusion_samples,
                    "step_scale": args.step_scale,
                    "output_format": args.output_format,
                }
                print(json.dumps(payload, indent=2))
                return 0
            result = predict_structure(
                sequence=args.sequence,
                api_key=args.api_key,
                base_url=args.base_url,
                use_self_hosted=args.self_hosted,
                recycling_steps=args.recycling_steps,
                sampling_steps=args.sampling_steps,
                diffusion_samples=args.diffusion_samples,
                step_scale=args.step_scale,
                output_format=args.output_format,
            )

        if args.output:
            output_path = Path(args.output)
            output_dir = output_path.parent if output_path.suffix else output_path
            client = Boltz2Client(
                api_key=args.api_key,
                base_url=args.base_url,
                use_self_hosted=args.self_hosted,
            )
            client.save_all_structures(result, output_dir)
            print(f"Saved {len(result.get('structures', []))} structure(s) to {output_dir}")
        else:
            structures = result.get("structures", [])
            for i, s in enumerate(structures):
                print(f"Structure {i}: {s.get('format', '?')} - {len(s.get('structure', ''))} chars")

        return 0

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="nvidia-bio",
        description="NVIDIA Bio Tools: ProteinMPNN & Boltz-2 pipelines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version="nvidia-bio 0.1.0")

    subparsers = parser.add_subparsers(dest="command", required=True)
    build_proteinmpnn_parser(subparsers)
    build_boltz2_parser(subparsers)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())