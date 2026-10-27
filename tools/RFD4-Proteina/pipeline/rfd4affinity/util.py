# -*- coding: utf-8 -*-
"""Shared helpers, constants and lightweight logging."""

from __future__ import annotations

import json
import os
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


# 20 standard amino acids plus common structure-only substitutions.
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}

AA1_TO_3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}

STANDARD_AA1 = set("ARNDCQEGHILKMFPSTWYV")


def now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def timestamp_slug() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def slugify(value: str, default: str = "run") -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "-", (value or "").strip())
    value = value.strip("-._")
    return value or default


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(obj, handle, ensure_ascii=False, indent=2, default=str)


def read_json(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def redact(value: str) -> str:
    """Mask a secret so it can safely appear in logs."""
    if not value:
        return "<unset>"
    if len(value) <= 8:
        return "***"
    return value[:4] + "..." + value[-3:]


_LOG_VERBOSE = False
_LOG_FILE: Path | None = None


def set_verbose(flag: bool) -> None:
    global _LOG_VERBOSE
    _LOG_VERBOSE = flag


def set_log_file(path: Path) -> None:
    global _LOG_FILE
    _LOG_FILE = path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("", encoding="utf-8")


def log(message: str, level: str = "INFO") -> None:
    stamp = datetime.now().strftime("%H:%M:%S")
    line = f"[{stamp}] {level:5s} {message}"
    print(line, file=sys.stderr, flush=True)
    if _LOG_FILE is not None:
        with _LOG_FILE.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")


def debug(message: str) -> None:
    if _LOG_VERBOSE:
        log(message, "DEBUG")


def warn(message: str) -> None:
    log(message, "WARN")


class PipelineError(RuntimeError):
    """Expected, user-actionable pipeline failure."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PipelineError(message)


def read_fasta(path: Path) -> dict[str, str]:
    """Minimal FASTA reader supporting one-line or wrapped sequences."""
    records: dict[str, list[str]] = {}
    name: str | None = None
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip().split()[0] or f"seq{len(records) + 1}"
                records.setdefault(name, [])
            elif name is not None:
                records[name].append(line)
    return {key: "".join(parts).upper() for key, parts in records.items()}


def read_sequence_arg(value: str) -> str:
    """Accept a raw sequence or an existing FASTA/plain file path."""
    candidate = Path(value)
    if candidate.exists() and candidate.is_file():
        text = candidate.read_text(encoding="utf-8")
        if text.lstrip().startswith(">"):
            recs = read_fasta(candidate)
            if not recs:
                raise PipelineError(f"FASTA file is empty: {candidate}")
            return next(iter(recs.values()))
        return re.sub(r"\s+", "", text).upper()
    return re.sub(r"\s+", "", value).upper()


def validate_sequence(seq: str, label: str = "sequence") -> str:
    seq = (seq or "").strip().upper()
    if not seq:
        raise PipelineError(f"{label} is empty")
    bad = sorted(set(seq) - STANDARD_AA1)
    if bad:
        raise PipelineError(
            f"{label} contains non-standard residues {bad!r}; only 20 standard "
            "amino acids are supported. Remove expression tags and stop codons."
        )
    return seq


def fmt_float(value, digits: int = 3) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "-"


def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "design"
