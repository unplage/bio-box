# -*- coding: utf-8 -*-
"""Shared helpers, constants and lightweight logging."""

from __future__ import annotations

import json
import re
import sys
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

STANDARD_AA1 = set("ARNDCQEGHILKMFPSTWYV")

# Anything that looks like an NVIDIA API key must never reach a log file.
SECRET_PATTERN = re.compile(r"nvapi-[A-Za-z0-9_\-]{16,}")


def now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def timestamp_slug() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def slugify(value: str, default: str = "run") -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "-", (value or "").strip())
    value = value.strip("-._")
    return value or default


def cif_data_name(value: str, default: str = "design") -> str:
    """Return a valid mmCIF ``data_`` block code.

    CIF block codes must start with a letter, so ``8UCD_vhh`` has to become
    ``d8UCD_vhh``; otherwise strict parsers reject the whole file.
    """
    name = slugify(value, default=default)
    if not name or not name[0].isalpha():
        name = f"d{name}" if name else default
    return name


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
    return "***REDACTED***"


_LOG_VERBOSE = False
_LOG_FILE: Path | None = None


def set_verbose(flag: bool) -> None:
    global _LOG_VERBOSE
    _LOG_VERBOSE = flag


def set_log_file(path: Path | None) -> None:
    global _LOG_FILE
    _LOG_FILE = path
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("", encoding="utf-8")


def log(message: str, level: str = "INFO") -> None:
    stamp = datetime.now().strftime("%H:%M:%S")
    line = f"[{stamp}] {level:5s} {SECRET_PATTERN.sub('***REDACTED***', message)}"
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


def _looks_like_path(value: str) -> bool:
    """Decide whether *value* is a file name rather than a raw sequence.

    A long raw sequence must never be handed to ``Path.is_file()``: the
    kernel rejects it with ``OSError: File name too long``.
    """
    stripped = value.strip()
    if not stripped or len(stripped) > 240 or "\n" in stripped or "\x00" in stripped:
        return False
    if "/" in stripped or "\\" in stripped:
        return True
    if stripped.lower().endswith((".fa", ".fasta", ".faa", ".txt", ".seq")):
        return True
    try:
        return Path(stripped).is_file()
    except OSError:
        return False


def _parse_fasta(text: str) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(parts).upper()
            name = line[1:].strip().split()[0] or f"seq{len(records) + 1}"
            parts = []
        elif name is not None:
            parts.append(line)
    if name is not None:
        records[name] = "".join(parts).upper()
    return {key: val for key, val in records.items() if val}


def _select_record(records: dict[str, str], selector: str, label: str) -> str:
    if not records:
        raise PipelineError(f"{label} contains no sequence")
    if selector:
        lowered = selector.lower()
        for key, value in records.items():
            if key.lower() == lowered:
                return value
        for key, value in records.items():
            if key.lower().startswith(lowered):
                return value
        raise PipelineError(
            f"No record named {selector!r} in {label}; available: "
            f"{', '.join(records)}"
        )
    if len(records) > 1:
        raise PipelineError(
            f"{label} contains {len(records)} records ({', '.join(records)}); "
            "select one with '<file>#<record name>'."
        )
    return next(iter(records.values()))


def read_sequence_arg(value: str) -> str:
    """Accept a raw sequence, an inline FASTA record or a FASTA/plain file.

    ``<file>#<record>`` selects one record from a multi-record FASTA.
    """
    raw, _, selector = value.partition("#")
    raw = raw.strip()
    selector = selector.strip()
    if raw.startswith(">"):
        return _select_record(_parse_fasta(raw), selector, "Sequence")
    if _looks_like_path(raw):
        candidate = Path(raw)
        try:
            exists = candidate.is_file()
        except OSError:
            exists = False
        if not exists:
            raise PipelineError(f"Sequence file not found: {candidate}")
        text = candidate.read_text(encoding="utf-8")
        if text.lstrip().startswith(">"):
            return _select_record(
                _parse_fasta(text), selector, f"FASTA file {candidate}"
            )
        if selector:
            raise PipelineError(
                f"{candidate} holds a single raw sequence; '#{selector}' is not usable."
            )
        return re.sub(r"\s+", "", text).upper()
    if selector:
        raise PipelineError(
            f"Record selector '#{selector}' only works with a FASTA file or record."
        )
    return re.sub(r"\s+", "", raw).upper()


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
