# -*- coding: utf-8 -*-
"""CDR identification.

Primary path: ANARCI (https://github.com/oxpig/ANARCI) with IMGT/Kabat/Chothia
numbering.  ANARCI is an optional dependency; when it is unavailable we fall
back to a conservative motif-based heuristic and say so loudly, because
mislabelled CDRs would silently redesign the wrong loop.
"""

from __future__ import annotations

from dataclasses import dataclass

from .util import PipelineError, validate_sequence, warn

# Scheme position ranges (inclusive) for the three antibody CDRs.
CDR_BOUNDS = {
    "imgt": {1: (27, 38), 2: (56, 65), 3: (105, 117)},
    "kabat": {1: (31, 35), 2: (50, 65), 3: (95, 102)},
    "chothia": {1: (26, 32), 2: (52, 56), 3: (95, 102)},
    "martin": {1: (26, 35), 2: (50, 58), 3: (95, 102)},
}

CDR_NAMES = ("CDR1", "CDR2", "CDR3")


@dataclass
class CDRHit:
    role: str          # A / H / L (as used in the report and structure)
    cdr: int           # 1, 2 or 3
    start: int         # 1-based position in the input sequence
    end: int
    sequence: str
    method: str        # 'anarci:imgt' or 'heuristic'

    @property
    def name(self) -> str:
        prefix = "" if self.role == "A" else self.role
        return f"{prefix}CDR{self.cdr}"

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "role": self.role,
            "cdr": self.cdr,
            "start": self.start,
            "end": self.end,
            "length": self.end - self.start + 1,
            "sequence": self.sequence,
            "method": self.method,
        }


def _anarci_rows(sequence: str, scheme: str):
    """Return numbering rows for *sequence* across ANARCI API versions."""
    try:
        from anarci import run_anarci
    except Exception as exc:  # pragma: no cover - environment dependent
        raise PipelineError(
            "ANARCI is not importable; install it or pass --cdr-ranges manually."
        ) from exc

    result = run_anarci([("query", sequence)], scheme=scheme)
    # ANARCI >= 2026 returns (inputs, numbered, details, hit_table); older
    # versions return (numbered, details, hit_table).
    numbered = result[1] if len(result) == 4 else result[0]
    hits = numbered[0] if numbered else []
    if not hits:
        return None
    hit = hits[0]
    # Newer: (rows, start, end); older: (start, end, rows)
    if isinstance(hit, (tuple, list)) and len(hit) == 3:
        if isinstance(hit[0], (list, tuple)) and hit and isinstance(hit[0][0], tuple):
            return list(hit[0])
        if isinstance(hit[2], (list, tuple)):
            return list(hit[2])
    if isinstance(hit, (tuple, list)) and len(hit) == 2 and isinstance(hit[0], list):
        return list(hit[0])
    return None


def _anarci_cdr_hits(sequence: str, scheme: str, role: str) -> list[CDRHit] | None:
    try:
        rows = _anarci_rows(sequence, scheme)
    except Exception as exc:
        warn(f"ANARCI numbering failed: {exc}")
        return None
    if not rows:
        warn("ANARCI found no immunoglobulin hit; falling back to heuristic CDRs.")
        return None

    scheme_positions: list[tuple[int, int, str]] = []  # (scheme_pos, seq_pos, aa)
    seq_pos = 0
    for position, aa in rows:
        number = position[0] if isinstance(position, (tuple, list)) else int(position)
        if aa in ("-", " ", ".") or aa is None:
            continue
        seq_pos += 1
        if seq_pos > len(sequence):
            break
        scheme_positions.append((int(number), seq_pos, aa))
    if seq_pos != len(sequence):
        warn(
            f"ANARCI numbered {seq_pos} residues but the input has {len(sequence)}; "
            "mapping may be incomplete."
        )

    bounds = CDR_BOUNDS[scheme]
    hits: list[CDRHit] = []
    for cdr in (1, 2, 3):
        lo, hi = bounds[cdr]
        # Filter on the *numbering* position but report the *sequence* positions:
        # the two differ whenever the framework has indels (IMGT insertions),
        # and using the numbering positions would redesign framework residues.
        positions = [
            seq_pos
            for scheme_pos, seq_pos, _ in scheme_positions
            if lo <= scheme_pos <= hi
        ]
        if not positions:
            continue
        start = min(positions)
        end = max(positions)
        sub = sequence[start - 1:end]
        hits.append(CDRHit(role, cdr, start, end, sub, f"anarci:{scheme}"))
    return hits or None


def _heuristic_cdr_hits(sequence: str, role: str) -> list[CDRHit]:
    """Very conservative Kabat-like anchors.

    This is a fallback, not a replacement for ANARCI.  It assumes a canonical
    Ig domain with Cys23, Trp41 and Cys92-like anchors and will be wrong for
    unusual frameworks.  Callers must surface ``method=heuristic``.
    """
    s = sequence.upper()
    n = len(s)
    cys = [i for i, aa in enumerate(s) if aa == "C"]
    cys1_candidates = [i for i in cys if i < 50]
    if not cys1_candidates:
        raise PipelineError(
            "Could not locate the first conserved cysteine for heuristic CDR "
            "detection; use ANARCI or --cdr-ranges."
        )
    cys1 = cys1_candidates[0]
    trp = next((i for i in range(cys1 + 8, n) if s[i] == "W"), None)
    cys2 = next((i for i in cys if i >= cys1 + 55), None)

    # CDR1 (Kabat ~31-35) sits a few residues after the first cysteine.
    c1_start = min(cys1 + 8, n)
    c1_end = min(cys1 + 12, n - 1)
    # CDR2 (Kabat ~50-65) follows the conserved tryptophan.
    if trp is not None:
        c2_start = min(trp + 9, n)
        c2_end = min(trp + 16, n - 1)
    else:
        c2_start = min(c1_end + 16, n)
        c2_end = min(c2_start + 7, n - 1)
    # CDR3 (Kabat ~95-102) follows the second cysteine and ends before WGxG.
    if cys2 is not None:
        c3_start = min(cys2 + 3, n)
        c3_end = n - 1
        for j in range(cys2 + 7, max(cys2 + 8, n - 3)):
            if j + 3 >= n:
                break
            if s[j] == "W" and s[j + 1] == "G" and s[j + 3] == "G":
                c3_end = j - 1
                break
        else:
            c3_end = min(cys2 + 12, n - 1)
    else:
        raise PipelineError(
            "Could not locate the second conserved cysteine for heuristic CDR "
            "detection; use ANARCI or --cdr-ranges."
        )
    spans = [(1, c1_start, c1_end), (2, c2_start, c2_end), (3, c3_start, c3_end)]
    hits = []
    for cdr, start, end in spans:
        if end < start:
            continue
        hits.append(
            CDRHit(role, cdr, start + 1, end + 1, s[start:end + 1], "heuristic")
        )
    return hits


def is_immunoglobulin(sequence: str, scheme: str = "imgt") -> bool | None:
    """True/False when ANARCI can decide, None when ANARCI is unavailable."""
    try:
        rows = _anarci_rows(sequence, scheme)
    except Exception:
        return None
    return bool(rows)


def number_cdrs(
    sequence: str,
    role: str,
    scheme: str = "imgt",
    prefer_anarci: bool = True,
) -> list[CDRHit]:
    """Return all three CDR hits for one antibody chain."""
    sequence = validate_sequence(sequence, f"{role} chain sequence")
    if scheme not in CDR_BOUNDS:
        raise PipelineError(
            f"Unsupported numbering scheme {scheme!r}; choose from {sorted(CDR_BOUNDS)}"
        )
    if prefer_anarci:
        hits = _anarci_cdr_hits(sequence, scheme, role)
        if hits:
            return hits
    warn(
        f"{role}: using heuristic CDR boundaries; verify before running. "
        "For production use ANARCI or pass --cdr-ranges."
    )
    return _heuristic_cdr_hits(sequence, role)


def _parse_cdr_names(mode: str, raw: str) -> list[tuple[str, int]]:
    """Turn a comma list such as ``HCDR3,LCDR1`` into (role, cdr) pairs."""
    roles = ["A"] if mode == "vhh" else ["H", "L"]
    text = (raw or "").strip().upper()
    if not text or text == "ALL":
        return [(role, cdr) for role in roles for cdr in (1, 2, 3)]

    selected: list[tuple[str, int]] = []
    for token in [t.strip() for t in text.split(",") if t.strip()]:
        digits = [char for char in token if char.isdigit()]
        if len(digits) != 1 or digits[0] not in "123":
            raise PipelineError(
                f"Unrecognised CDR name {token!r}; use e.g. CDR3 or HCDR3"
            )
        digit = int(digits[0])
        role = None
        for candidate in roles:
            if token.startswith(candidate + "CDR"):
                role = candidate
                break
        if role is None and token == f"CDR{digit}":
            if len(roles) == 1:
                role = roles[0]
            else:
                raise PipelineError(
                    f"{token!r} is ambiguous for a Fab; write HCDR{digit} or LCDR{digit}."
                )
        if role is None and token == f"{roles[0]}CDR{digit}":
            role = roles[0]
        if role is None:
            raise PipelineError(f"Unrecognised CDR name {token!r} for mode {mode}.")
        if (role, digit) not in selected:
            selected.append((role, digit))
    return selected


def _parse_cdr_ranges(raw: str, roles: list[str]) -> dict[str, list[tuple[int, int]]]:
    """Parse ``H:99-109,L:89-97`` (1-based sequence coordinates)."""
    out: dict[str, list[tuple[int, int]]] = {}
    for token in [t.strip() for t in (raw or "").split(",") if t.strip()]:
        if ":" not in token:
            raise PipelineError(
                f"Invalid CDR range {token!r}; expected CHAIN:START-END (e.g. H:99-109)"
            )
        chain, span = token.split(":", 1)
        chain = chain.strip().upper()
        if chain not in roles:
            if len(roles) == 1 and chain in ("CDR1", "CDR2", "CDR3", "A"):
                chain = roles[0]
            else:
                raise PipelineError(
                    f"Chain {chain!r} not valid for this run; expected one of {roles}"
                )
        span = span.replace("..", "-")
        if "-" in span:
            start_s, end_s = span.split("-", 1)
        else:
            start_s = end_s = span
        try:
            start, end = int(start_s), int(end_s)
        except ValueError:
            raise PipelineError(
                f"Invalid CDR range {span!r} for chain {chain}; expected integers "
                "such as 99-109."
            ) from None
        if start < 1 or end < start:
            raise PipelineError(f"Invalid range {span!r} for chain {chain}")
        out.setdefault(chain, []).append((start, end))
    return out


def resolve_cdr_specs(
    sequences: dict[str, str],
    mode: str,
    cdrs: str | None = None,
    cdr_ranges: str | None = None,
    scheme: str = "imgt",
    prefer_anarci: bool = True,
) -> tuple[dict[str, list[tuple[int, int]]], list[dict], str]:
    """Resolve the user's CDR selection to per-role sequence ranges.

    Returns ``(ranges_by_role, hit_metadata, method_summary)`` where ranges are
    1-based inclusive positions in each role's input sequence.
    """
    roles = list(sequences.keys())
    if cdr_ranges:
        ranges = _parse_cdr_ranges(cdr_ranges, roles)
        metadata = []
        for role, spans in ranges.items():
            for start, end in spans:
                metadata.append(
                    {
                        "name": f"{'' if role == 'A' else role}CDR?",
                        "role": role,
                        "start": start,
                        "end": end,
                        "length": end - start + 1,
                        "sequence": sequences[role][start - 1:end],
                        "method": "manual",
                    }
                )
        return ranges, metadata, "manual ranges"

    wanted = _parse_cdr_names(mode, cdrs or "all")
    ranges: dict[str, list[tuple[int, int]]] = {}
    metadata: list[dict] = []
    methods: set[str] = set()
    for role in roles:
        hits = number_cdrs(sequences[role], role, scheme=scheme, prefer_anarci=prefer_anarci)
        by_index = {hit.cdr: hit for hit in hits}
        for sel_role, cdr in wanted:
            if sel_role != role:
                continue
            hit = by_index.get(cdr)
            if hit is None:
                raise PipelineError(f"Could not identify {role} CDR{cdr}.")
            ranges.setdefault(role, []).append((hit.start, hit.end))
            metadata.append(hit.as_dict())
            methods.add(hit.method)
    if not ranges:
        raise PipelineError("No CDRs selected; nothing to design.")
    return ranges, metadata, ", ".join(sorted(methods))


def merge_ranges(ranges: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping/adjacent inclusive ranges."""
    if not ranges:
        return []
    ordered = sorted(ranges)
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(int(a), int(b)) for a, b in merged]


def expand_ranges(ranges: list[tuple[int, int]]) -> list[int]:
    out: list[int] = []
    for start, end in ranges:
        out.extend(range(start, end + 1))
    return out
