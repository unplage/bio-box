# -*- coding: utf-8 -*-
"""Structure sourcing: local files, RCSB download, ESMFold (optional)."""

from __future__ import annotations

import re
from pathlib import Path

from .net import get_text, post_json
from .structio import Chain, mean_plddt_from_chains, parse_structure_text
from .util import PipelineError, log, warn

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.{fmt}"
ESMFOLD_PATH = "/biology/nvidia/esmfold"


def _looks_like_pdb_id(value: str) -> bool:
    return bool(re.fullmatch(r"[0-9][A-Za-z0-9]{3}", value.strip()))


def download_rcsb(pdb_id: str, out_dir: Path, fmt: str = "cif") -> tuple[dict[str, Chain], str]:
    pdb_id = pdb_id.strip().upper()
    out_dir.mkdir(parents=True, exist_ok=True)
    errors = []
    for candidate in (fmt, "cif", "pdb"):
        url = RCSB_URL.format(pdb_id=pdb_id.lower(), fmt=candidate)
        try:
            log(f"Downloading {pdb_id} from RCSB ({candidate})")
            text = get_text(url, timeout=60)
            path = out_dir / f"{pdb_id}.{candidate}"
            path.write_text(text, encoding="utf-8")
            chains = parse_structure_text(text)
            return chains, str(path)
        except PipelineError as exc:
            errors.append(str(exc))
    raise PipelineError(
        f"Could not download {pdb_id} from RCSB: {'; '.join(errors)}"
    )


def load_target(value: str, out_dir: Path) -> tuple[dict[str, Chain], str]:
    """Accept a local structure path or a 4-character PDB ID."""
    path = Path(value)
    if path.exists():
        return parse_structure_text(path.read_text(encoding="utf-8", errors="replace")), str(path)
    if _looks_like_pdb_id(value):
        return download_rcsb(value, out_dir)
    raise PipelineError(
        f"--target-structure {value!r} is neither an existing file nor a 4-character PDB ID"
    )


def _extract_pdb(result: dict) -> str:
    candidates = []
    for key in ("pdbs", "pdb", "structures", "output"):
        if key in result:
            candidates.append(result[key])
    for candidate in candidates:
        if isinstance(candidate, str) and candidate.lstrip().startswith(("ATOM", "HEADER", "data_")):
            return candidate
        if isinstance(candidate, list) and candidate:
            first = candidate[0]
            if isinstance(first, str):
                return first
            if isinstance(first, dict):
                for inner in ("pdb", "pdbs", "content", "structure"):
                    if isinstance(first.get(inner), str):
                        return first[inner]
        if isinstance(candidate, dict):
            for inner in ("pdbs", "pdb", "content", "structure"):
                value = candidate.get(inner)
                if isinstance(value, str):
                    return value
                if isinstance(value, list) and value and isinstance(value[0], str):
                    return value[0]
    raise PipelineError(
        "ESMFold response did not contain a structure; keys were "
        f"{sorted(result)}"
    )


def esmfold_sequence(
    sequence: str,
    api_key: str,
    base_url: str = "https://health.api.nvidia.com/v1",
    timeout: float = 600.0,
) -> tuple[dict[str, Chain], str]:
    url = base_url.rstrip("/") + ESMFOLD_PATH
    log(f"Folding sequence with ESMFold ({len(sequence)} aa)")
    result = post_json(url, {"sequence": sequence}, api_key=api_key, timeout=timeout)
    pdb_text = _extract_pdb(result)
    chains = parse_structure_text(pdb_text)
    if len(chains) != 1:
        warn(f"ESMFold returned {len(chains)} chains; using the longest one")
        key = max(chains, key=lambda k: len(chains[k].residues))
        chains = {key: chains[key]}
    chain = next(iter(chains.values()))
    if len(chain.residues) != len(sequence):
        if abs(len(chain.residues) - len(sequence)) > 3:
            raise PipelineError(
                f"ESMFold returned {len(chain.residues)} residues for a "
                f"{len(sequence)}-aa sequence; refusing to guess the mapping."
            )
    plddt = mean_plddt_from_chains(chains)
    if plddt is not None:
        log(f"ESMFold mean pLDDT {plddt:.1f}")
    return chains, pdb_text


def fold_vhh(
    sequence: str,
    api_key: str | None,
    base_url: str,
    timeout: float,
) -> tuple[dict[str, Chain], str]:
    if not api_key:
        raise PipelineError(
            "Folding a VHH without a structure needs NVIDIA_API_KEY (ESMFold NIM). "
            "Pass --antibody-structure instead, or use --fold none."
        )
    try:
        return esmfold_sequence(sequence, api_key, base_url, timeout)
    except PipelineError as exc:
        raise PipelineError(
            "ESMFold could not fold the VHH sequence. The public ESMFold NIM is "
            "marked deprecated; supply --antibody-structure instead. "
            f"Original error: {exc}"
        ) from exc


def describe_chains(chains: dict[str, Chain]) -> str:
    parts = []
    for chain in chains.values():
        parts.append(f"{chain.id}:{len(chain.residues)}aa")
    return ", ".join(parts)
