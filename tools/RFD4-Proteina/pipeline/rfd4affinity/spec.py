# -*- coding: utf-8 -*-
"""Turn user inputs into an annotated mmCIF plus a NIM request body."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .annotate import ChainPlan, LoopSegment, write_annotated_mmcif
from .cdr import expand_ranges, is_immunoglobulin, merge_ranges, resolve_cdr_specs
from .fold import describe_chains, fold_vhh, load_target
from .structio import (
    Chain,
    chain_table,
    combine_chains,
    map_sequence_positions,
    read_structure,
)
from .util import (
    PipelineError,
    debug,
    log,
    read_sequence_arg,
    require,
    slugify,
    validate_sequence,
    warn,
    write_json,
)

# Metric tags understood by the RFD4-Proteina NIM.
#
# ``ORACLE_METRIC_TAGS`` need an optional runtime that the hosted trial
# deployment does not ship (PyRosetta or the structure-prediction oracle), so
# the server rejects them with HTTP 422.  They are never enabled by default.
CORE_METRIC_TAGS = {
    "design",
    "hotspot_conditioning",
    "rasa_conditioning",
    "centroid_conditioning",
    "interchain_interactions",
    "3d_coordinate_conditioning",
    "internal_coordinate_conditioning",
    "unindexed",
    "partial_diffusion",
    "small_molecule_only",
    "no_ligand",
    "cms",
    "cath",
    "structure_prediction",
}

ORACLE_METRIC_TAGS = {
    "antibody",
    "developability",
    "unconditional_refolding",
    "ppi_refolding",
    "motif_scaffold_refolding",
    "internal_coordinate_motif_scaffold_refolding",
    "sm_binder_refolding",
    "dna_refolding",
    "apo_refolding",
    "ame_refolding",
    "ame_no_ligand_refolding",
    "disco_refolding",
}

METRIC_TAGS = CORE_METRIC_TAGS | ORACLE_METRIC_TAGS

# The hosted NIM only guarantees the ``design`` tag (that is what NVIDIA's own
# examples request).  Extra conditioning tags are opt-in via --metrics-tags.
DEFAULT_METRIC_TAGS = ["design"]

MAX_REQUEST_BYTES = 32 * 1024 * 1024
MAX_TOKENS = 1024
MAX_NUM_STEPS = 200


@dataclass
class RunConfig:
    mode: str = "vhh"
    vhh_seq: str | None = None
    vh_seq: str | None = None
    vl_seq: str | None = None
    antibody_structure: str | None = None
    complex_structure: str | None = None
    antibody_chains: str | None = None
    target_structure: str | None = None
    target_chains: str | None = None
    target_name: str = "target"
    cdrs: str = "all"
    cdr_scheme: str = "imgt"
    cdr_ranges: str | None = None
    epitope: str | None = None
    fold: str = "auto"
    offset: float = 60.0
    fix_sidechains: bool = False
    num_designs: int = 4
    num_steps: int | None = None
    recycle: str = "thinking_medium"
    seed: int | None = None
    metrics_tags: list[str] | None = None
    example_id: str | None = None
    base_url: str = "https://health.api.nvidia.com/v1"
    cdr_encoding: str = "expand"
    cdr_length_range: str | None = None


@dataclass
class RunResult:
    case: str
    mode: str
    request: dict
    cif_text: str
    chains: dict[str, Chain]
    plans: dict[str, ChainPlan]
    role_to_chain: dict[str, str]
    target_chain_ids: list[str]
    cdr_metadata: list[dict]
    cdr_method: str
    target_name: str
    sequences: dict[str, str]
    warnings: list[str] = field(default_factory=list)
    chain_table: list[dict] = field(default_factory=list)
    # CDR runs per role, in *sequence* coordinates of the parent, with the
    # length bounds that were handed to the model.
    cdr_segments: dict[str, list[dict]] = field(default_factory=dict)
    cdr_encoding: str = "expand"


def _split_chains(raw: str | None) -> list[str]:
    if not raw:
        return []
    tokens = [t.strip() for t in raw.replace(",", " ").split() if t.strip()]
    return tokens


def _parse_length_overrides(raw: str | None) -> dict[str, tuple[str, int, int]]:
    """Parse ``--cdr-length-range`` into ``{selector: (kind, lo, hi)}``.

    Accepted selectors: ``CDR3`` (any chain), ``HCDR3`` / ``LCDR3`` (one
    chain) and ``*`` (every CDR).  Accepted values: ``min-max`` (absolute
    lengths) or ``+n`` (keep the native length, allow up to ``n`` residues
    more).
    """
    if not raw:
        return {}
    overrides: dict[str, tuple[str, int, int]] = {}
    for token in raw.replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            raise PipelineError(
                f"Invalid --cdr-length-range entry {token!r}; expected e.g. "
                "'HCDR3:12-21' or '*:+4'."
            )
        selector, value = token.split(":", 1)
        selector = selector.strip().upper()
        value = value.strip()
        if not selector:
            raise PipelineError(f"Invalid --cdr-length-range entry {token!r}: empty selector.")
        if value.startswith("+"):
            try:
                delta = int(value[1:])
            except ValueError:
                raise PipelineError(
                    f"Invalid --cdr-length-range entry {token!r}; expected '+n' or 'min-max'."
                ) from None
            if delta < 0:
                raise PipelineError(
                    f"Invalid --cdr-length-range entry {token!r}; the delta cannot be negative."
                )
            overrides[selector] = ("rel", 0, delta)
            continue
        if "-" not in value:
            raise PipelineError(
                f"Invalid --cdr-length-range entry {token!r}; expected 'min-max' or '+n'."
            )
        low_s, high_s = value.split("-", 1)
        try:
            low, high = int(low_s), int(high_s)
        except ValueError:
            raise PipelineError(
                f"Invalid --cdr-length-range entry {token!r}; the bounds must be integers."
            ) from None
        if low < 1 or high < low:
            raise PipelineError(f"Invalid --cdr-length-range entry {token!r}: {low}-{high}.")
        overrides[selector] = ("abs", low, high)
    return overrides


def _length_bounds(
    name: str,
    native: int,
    overrides: dict[str, tuple[str, int, int]],
) -> tuple[int, int]:
    """Resolve the length bounds for one CDR."""
    candidates = []
    if name:
        candidates.append(name.upper())
        generic = "CDR" + "".join(ch for ch in name.upper() if ch.isdigit())
        if generic and generic not in candidates:
            candidates.append(generic)
    candidates.append("*")
    for selector in candidates:
        if selector not in overrides:
            continue
        kind, low, high = overrides[selector]
        if kind == "rel":
            return native, native + high
        return low, high
    return native, native


def _build_cdr_segments(
    sequences: dict[str, str],
    cdr_ranges: dict[str, list[tuple[int, int]]],
    cdr_metadata: list[dict],
    overrides: dict[str, tuple[str, int, int]],
) -> dict[str, list[dict]]:
    """Turn the selected CDR ranges into contiguous, length-bounded segments."""
    segments: dict[str, list[dict]] = {}
    for role, ranges in cdr_ranges.items():
        sequence = sequences[role]
        role_meta = [m for m in cdr_metadata if m.get("role") == role]
        out: list[dict] = []
        for start, end in merge_ranges(ranges):
            members = [m for m in role_meta if m["start"] <= end and m["end"] >= start]
            names = sorted({str(m.get("name") or "CDR?") for m in members})
            native = end - start + 1
            min_len, max_len = native, native
            has_absolute = False
            for member in members:
                low, high = _length_bounds(str(member.get("name") or ""), native, overrides)
                # Detect whether this member had an absolute override
                name = str(member.get("name") or "")
                is_absolute = False
                for selector in (name.upper(), "CDR" + "".join(ch for ch in name.upper() if ch.isdigit()), "*"):
                    if selector in overrides and overrides[selector][0] == "abs":
                        is_absolute = True
                        break
                if is_absolute:
                    if not has_absolute:
                        min_len, max_len = low, high
                        has_absolute = True
                    else:
                        min_len, max_len = min(min_len, low), max(max_len, high)
                else:
                    if not has_absolute:
                        min_len = min(min_len, low)
                        max_len = max(max_len, high)
            out.append(
                {
                    "name": "+".join(names) if names else "CDR?",
                    "start": start,
                    "end": end,
                    "native_length": native,
                    "min_len": min_len,
                    "max_len": max_len,
                    "sequence": sequence[start - 1:end],
                }
            )
        segments[role] = out
    return segments


def centroid(chains: dict[str, Chain]) -> tuple[float, float, float]:
    xs = ys = zs = 0.0
    n = 0
    for chain in chains.values():
        for residue in chain.residues:
            for x, y, z in residue.atoms.values():
                xs += x
                ys += y
                zs += z
                n += 1
    if not n:
        return (0.0, 0.0, 0.0)
    return (xs / n, ys / n, zs / n)


def _choose_antibody_chains(
    chains: dict[str, Chain], cfg: RunConfig
) -> tuple[list[str], list[str]]:
    ids = list(chains)
    explicit = _split_chains(cfg.antibody_chains)
    if explicit:
        missing = [cid for cid in explicit if cid not in chains]
        require(
            not missing,
            f"--antibody-chains {missing} not found in the structure ({ids}).",
        )
        require(
            len(set(explicit)) == len(explicit),
            "--antibody-chains contains duplicate chain IDs.",
        )
    if cfg.mode == "vhh":
        if explicit:
            require(len(explicit) == 1, "--antibody-chains must name exactly one VHH chain")
            ab = explicit
        else:
            ig_like = [
                cid for cid in ids
                if is_immunoglobulin(chains[cid].sequence()) is True
            ]
            if len(ig_like) == 1:
                ab = ig_like
            elif len(ids) == 1:
                ab = ids
            elif len(ig_like) > 1:
                raise PipelineError(
                    "More than one chain looks like an immunoglobulin; pass "
                    f"--antibody-chains explicitly (candidates: {ig_like})."
                )
            else:
                raise PipelineError(
                    "Could not identify the VHH chain automatically. Pass "
                    f"--antibody-chains (chains: {ids}) or install ANARCI."
                )
            log(f"Auto-selected antibody chain {ab[0]}")
        target = [cid for cid in ids if cid not in ab]
        if cfg.target_chains:
            target = _split_chains(cfg.target_chains)
        return ab, target

    # Fab
    if explicit:
        require(
            len(explicit) == 2,
            "--antibody-chains for a Fab must name two chains (heavy, light)",
        )
        ab = explicit
        target = [cid for cid in ids if cid not in ab]
        if cfg.target_chains:
            target = _split_chains(cfg.target_chains)
        return ab, target

    if len(ids) == 2:
        ordered = sorted(ids, key=lambda cid: len(chains[cid].residues), reverse=True)
        warn(
            "No --antibody-chains given; treating the longer chain as heavy (H) "
            "and the shorter as light (L). Verify this."
        )
        target = [cid for cid in ids if cid not in ordered]
        if cfg.target_chains:
            target = _split_chains(cfg.target_chains)
        return ordered, target
    raise PipelineError(
        "For a Fab, pass --antibody-chains (two chains, heavy then light); "
        f"the structure has chains {ids}."
    )


def _load_antibody(
    cfg: RunConfig,
    api_key: str | None,
    out_dir: Path,
) -> tuple[dict[str, Chain], str]:
    if cfg.complex_structure:
        path = Path(cfg.complex_structure)
        require(path.exists(), f"--complex-structure not found: {path}")
        log(f"Reading antibody–target complex {path}")
        return read_structure(path), "complex"
    if cfg.antibody_structure:
        path = Path(cfg.antibody_structure)
        require(path.exists(), f"--antibody-structure not found: {path}")
        log(f"Reading antibody structure {path}")
        return read_structure(path), "antibody"
    if cfg.mode == "fab":
        raise PipelineError(
            "A Fab cannot be assembled from sequence alone. Supply "
            "--antibody-structure (VH/VL) or --complex-structure."
        )
    sequence = cfg.vhh_seq
    require(
        bool(sequence),
        "VHH mode needs --vhh-seq (or --antibody-structure / --complex-structure).",
    )
    sequence = validate_sequence(read_sequence_arg(str(sequence)), "VHH")
    if cfg.fold == "none":
        raise PipelineError(
            "No antibody structure supplied and --fold none was requested; "
            "provide --antibody-structure."
        )
    chains, _ = fold_vhh(sequence, api_key, cfg.base_url, timeout=600.0)
    return chains, "folded"


def _target_epitope_residues(
    raw: str,
    target_chains: dict[str, Chain],
) -> dict[str, set[int]]:
    """Parse ``557-573,591-605`` or ``T:557-573`` into resseq sets."""
    spans: dict[str, list[tuple[int, int]]] = {}
    for token in [t.strip() for t in raw.replace(";", ",").split(",") if t.strip()]:
        chain_id = None
        if ":" in token:
            chain_id, token = token.split(":", 1)
            chain_id = chain_id.strip().upper()
        token = token.replace("..", "-")
        try:
            if "-" in token:
                start_s, end_s = token.split("-", 1)
                start, end = int(start_s), int(end_s)
            else:
                start = end = int(token)
        except ValueError:
            raise PipelineError(f"Invalid epitope range {token!r}") from None
        if start < 1 or end < start:
            raise PipelineError(f"Invalid epitope range {token!r}")
        spans.setdefault(chain_id or "*", []).append((start, end))

    out: dict[str, set[int]] = {}
    for chain_id, ranges in spans.items():
        if chain_id and chain_id != "*" and chain_id not in target_chains:
            raise PipelineError(
                f"Epitope chain {chain_id!r} is not a target chain "
                f"({list(target_chains)})."
            )
        targets = (
            [target_chains[chain_id]]
            if chain_id and chain_id != "*"
            else list(target_chains.values())
        )
        for chain in targets:
            resseqs = set(chain.resseqs())
            for start, end in ranges:
                if start in resseqs and end in resseqs:
                    resolved = set(range(start, end + 1))
                elif 1 <= start <= len(chain.residues) and 1 <= end <= len(chain.residues):
                    positions = list(range(start, end + 1))
                    resolved = set(map_sequence_positions(chain.sequence(), chain, positions))
                else:
                    raise PipelineError(
                        f"Epitope range {start}-{end} is outside chain {chain.id}: "
                        f"author residues {min(resseqs)}-{max(resseqs)}, "
                        f"{len(chain.residues)} residues. Check the chain ID and numbering."
                    )
                if not resolved:
                    raise PipelineError(
                        f"Epitope range {start}-{end} mapped to no residues in chain {chain.id}."
                    )
                out.setdefault(chain.id, set()).update(resolved)
    return out


def _validate_tags(tags: list[str]) -> list[str]:
    cleaned = []
    for tag in tags:
        tag = tag.strip()
        if not tag:
            continue
        if tag not in METRIC_TAGS:
            raise PipelineError(
                f"Unknown metric tag {tag!r}. Valid tags: "
                f"{', '.join(sorted(METRIC_TAGS))}."
            )
        if tag in ORACLE_METRIC_TAGS:
            warn(
                f"Metric tag {tag!r} needs an optional runtime (PyRosetta or the "
                "structure-prediction oracle). The hosted NIM will reject it with "
                "HTTP 422 unless the deployment ships that runtime."
            )
        if tag not in cleaned:
            cleaned.append(tag)
    return cleaned


def build_run(
    cfg: RunConfig,
    api_key: str | None,
    out_dir: Path,
) -> RunResult:
    out_dir.mkdir(parents=True, exist_ok=True)
    require(cfg.mode in ("vhh", "fab"), "--mode must be vhh or fab")
    require(cfg.num_designs >= 1, "--num-designs must be >= 1")
    if cfg.num_designs > 32:
        warn(
            f"--num-designs {cfg.num_designs} is large; the NIM batches replicates "
            "into one forward pass and GPU memory grows with the batch size."
        )
    require(
        cfg.recycle in ("thinking_low", "thinking_medium", "thinking_high"),
        "--recycle must be thinking_low, thinking_medium or thinking_high",
    )
    if cfg.num_steps is not None:
        require(
            1 <= cfg.num_steps <= MAX_NUM_STEPS,
            f"--num-steps must be between 1 and {MAX_NUM_STEPS} (server limit).",
        )
    require(cfg.offset >= 0, "--offset must be >= 0")
    if cfg.complex_structure and cfg.target_structure:
        raise PipelineError(
            "--complex-structure already contains the target; do not also pass "
            "--target-structure."
        )

    antibody_chains, source = _load_antibody(cfg, api_key, out_dir)
    antibody_ids, remaining = _choose_antibody_chains(antibody_chains, cfg)

    if source == "folded":
        require(
            len(antibody_ids) == 1,
            "Folded VHH should be a single chain; use --antibody-structure for Fabs.",
        )

    # Target chains come either from the complex or from a separate structure.
    if cfg.target_structure:
        target_chains, target_source = load_target(cfg.target_structure, out_dir)
        if cfg.target_chains:
            target_ids = _split_chains(cfg.target_chains)
        else:
            target_ids = list(target_chains)
        for cid in target_ids:
            require(cid in target_chains, f"Target chain {cid!r} not found in {target_source}")
        target = target_chains
        case = "docking"
    elif cfg.complex_structure:
        target_ids = remaining
        if cfg.target_chains:
            target_ids = _split_chains(cfg.target_chains)
        require(
            bool(target_ids),
            "No target chains were found. Pass --target-chains, or use "
            "--target-structure with --complex-structure omitted.",
        )
        for cid in target_ids:
            require(
                cid in antibody_chains,
                f"Target chain {cid!r} not found in the complex "
                f"({list(antibody_chains)}).",
            )
        overlap = set(target_ids) & set(antibody_ids)
        require(not overlap, f"Chains {sorted(overlap)} cannot be both antibody and target.")
        target = antibody_chains
        case = "complex"
    else:
        if cfg.target_chains:
            warn("--target-chains is ignored when no target structure is supplied.")
        target_ids = []
        target = {}
        case = "apo"

    # Role -> source chain id.
    if cfg.mode == "vhh":
        require(len(antibody_ids) == 1, "VHH mode expects exactly one antibody chain")
        role_to_source = {"A": antibody_ids[0]}
        new_antibody_ids = {"A": "A"}
    else:
        require(
            len(antibody_ids) == 2,
            "Fab mode expects exactly two antibody chains (heavy, light)",
        )
        heavy, light = antibody_ids
        h_seq = antibody_chains[heavy].sequence()
        l_seq = antibody_chains[light].sequence()
        if len(h_seq) < len(l_seq):
            warn("Heavy chain is shorter than light chain; double-check --antibody-chains order.")
        role_to_source = {"H": heavy, "L": light}
        new_antibody_ids = {"H": "H", "L": "L"}

    sequences: dict[str, str] = {}
    user_sequences = {
        "A": cfg.vhh_seq,
        "H": cfg.vh_seq,
        "L": cfg.vl_seq,
    }
    for role, source_id in role_to_source.items():
        chain = antibody_chains[source_id]
        raw = user_sequences.get(role)
        if raw:
            seq = validate_sequence(read_sequence_arg(str(raw)), f"{role} chain")
            if seq != chain.sequence():
                warn(
                    f"Provided {role} sequence differs from structure chain "
                    f"{source_id}; CDR positions will be aligned."
                )
        else:
            seq = validate_sequence(chain.sequence(), f"{role} chain (from structure)")
        sequences[role] = seq

    cdr_ranges, cdr_metadata, cdr_method = resolve_cdr_specs(
        sequences,
        cfg.mode,
        cdrs=cfg.cdrs,
        cdr_ranges=cfg.cdr_ranges,
        scheme=cfg.cdr_scheme,
    )

    # Map CDR sequence ranges onto structure residue numbers.
    free_by_source: dict[str, set[int]] = {}
    for role, ranges in cdr_ranges.items():
        chain = antibody_chains[role_to_source[role]]
        positions = expand_ranges(ranges)
        mapped = map_sequence_positions(sequences[role], chain, positions)
        debug(f"{role}: CDR sequence positions {positions} -> structure {mapped}")
        free_by_source.setdefault(chain.id, set()).update(mapped)

    # CDR runs with the length bounds handed to the model.
    length_overrides = _parse_length_overrides(cfg.cdr_length_range)
    cdr_segments = _build_cdr_segments(
        sequences, cdr_ranges, cdr_metadata, length_overrides
    )
    for role, segments in cdr_segments.items():
        for segment in segments:
            debug(
                f"{role} {segment['name']}: sequence {segment['start']}-{segment['end']} "
                f"({segment['native_length']} aa), model length "
                f"{segment['min_len']}-{segment['max_len']}"
            )

    # Hotspots on the target.
    hotspots_by_source: dict[str, set[int]] = {}
    if cfg.epitope:
        if not target_ids:
            warn("--epitope given but there is no target; ignoring hotspots.")
        else:
            target_subset = {cid: target[cid] for cid in target_ids}
            hotspots_by_source = _target_epitope_residues(cfg.epitope, target_subset)

    # Rename chains into the design namespace.
    antibody_map = {role_to_source[role]: new for role, new in new_antibody_ids.items()}
    target_map: dict[str, str] = {}
    for index, old in enumerate(target_ids):
        target_map[old] = "TUVWXYZABCDEFGHIJKLMNOPQRS"[index % 26]

    target_offset = (0.0, 0.0, 0.0)
    if case == "docking":
        ab_centroid = centroid({cid: antibody_chains[cid] for cid in antibody_ids})
        tgt_centroid = centroid({cid: target[cid] for cid in target_ids})
        target_offset = (
            ab_centroid[0] + cfg.offset - tgt_centroid[0],
            ab_centroid[1] - tgt_centroid[1],
            ab_centroid[2] - tgt_centroid[2],
        )
        log(
            f"Placing target {cfg.offset:.0f} Å from the antibody centroid "
            "(relative orientation is sampled by the model)"
        )

    combined = combine_chains(
        antibody_chains,
        antibody_map,
        target,
        target_map,
        target_offset=target_offset,
    )

    plans: dict[str, ChainPlan] = {}
    expand_encoding = cfg.cdr_encoding == "expand"
    loop_segments_by_source: dict[str, list[LoopSegment]] = {}
    for role, segments in cdr_segments.items():
        source_id = role_to_source[role]
        chain = antibody_chains[source_id]
        built: list[LoopSegment] = []
        for segment in segments:
            mapped = map_sequence_positions(
                sequences[role], chain, [segment["start"], segment["end"]]
            )
            built.append(
                LoopSegment(
                    start=mapped[0],
                    end=mapped[-1],
                    min_len=segment["min_len"],
                    max_len=segment["max_len"],
                    name=segment["name"],
                )
            )
        loop_segments_by_source[source_id] = built

    for old, new in antibody_map.items():
        plans[new] = ChainPlan(
            chain_id=new,
            role="antibody",
            free_residues=set(free_by_source.get(old, set())),
            hotspot_residues=set(),
            token_group=0 if case == "complex" else 1,
            loop_segments=loop_segments_by_source.get(old, []) if expand_encoding else [],
        )
    for old, new in target_map.items():
        plans[new] = ChainPlan(
            chain_id=new,
            role="target",
            free_residues=set(),
            hotspot_residues=set(hotspots_by_source.get(old, set())),
            token_group=0,
        )

    # Sanity checks.
    total_residues = sum(len(chain.residues) for chain in combined.values())
    input_residues = total_residues
    if expand_encoding:
        # The CDR residues are replaced by one expandable segment each.
        for segments in cdr_segments.values():
            for segment in segments:
                input_residues -= segment["native_length"]
        input_residues += sum(len(segments) for segments in cdr_segments.values())
        max_residues = input_residues + sum(
            segment["max_len"] - 1
            for segments in cdr_segments.values()
            for segment in segments
        )
        if max_residues > MAX_TOKENS:
            raise PipelineError(
                f"The design can reach {max_residues} residues but the hosted NIM "
                f"accepts at most {MAX_TOKENS} tokens. Trim the target or shorten "
                "the allowed CDR lengths."
            )
    elif total_residues > MAX_TOKENS:
        raise PipelineError(
            f"The design has {total_residues} residues but the hosted NIM accepts "
            f"at most {MAX_TOKENS} tokens. Trim the target (use a domain, not the "
            "whole ectodomain) or use a self-hosted NIM."
        )
    designed = sum(len(plan.free_residues) for plan in plans.values())
    if designed == 0:
        raise PipelineError(
            "No CDR residues were selected, so nothing would be designed. "
            "Use --cdrs all or name CDRs explicitly."
        )

    example_id = slugify(cfg.example_id or f"{cfg.target_name}_{cfg.mode}", default="design")
    cif_path = out_dir / "input_complex.cif"
    write_annotated_mmcif(
        combined,
        plans,
        cif_path,
        case=case,
        fix_sidechains=cfg.fix_sidechains,
        data_name=example_id,
    )
    cif_text = cif_path.read_text(encoding="utf-8")
    size = len(cif_text.encode("utf-8"))
    if size > MAX_REQUEST_BYTES:
        raise PipelineError(
            f"Annotated mmCIF is {size / 1e6:.1f} MB, over the 32 MB request limit."
        )

    tags = cfg.metrics_tags
    if not tags:
        tags = list(DEFAULT_METRIC_TAGS)
        if hotspots_by_source:
            tags.append("hotspot_conditioning")
        if case in ("complex", "docking"):
            tags.append("interchain_interactions")
        if case != "complex":
            tags.append("internal_coordinate_conditioning")
    tags = _validate_tags(tags)
    debug(f"metric tags: {tags}")

    request: dict = {
        "task": "antibody_design",
        "input": {
            "format": "annotated_cif",
            "example_id": example_id,
            "file": {
                "filename": cif_path.name,
                "content": cif_text,
            },
        },
        "sampling": {
            "num_replicates": cfg.num_designs,
            "recycle_strategy": cfg.recycle,
        },
        "metrics": {"tags": tags},
    }
    if cfg.num_steps is not None:
        request["sampling"]["num_steps"] = cfg.num_steps
    if cfg.seed is not None:
        request["sampling"]["seed"] = cfg.seed

    write_json(out_dir / "request.json", request)

    result = RunResult(
        case=case,
        mode=cfg.mode,
        request=request,
        cif_text=cif_text,
        chains=combined,
        plans=plans,
        role_to_chain={role: new_antibody_ids[role] for role in new_antibody_ids},
        target_chain_ids=list(target_map.values()),
        cdr_metadata=cdr_metadata,
        cdr_method=cdr_method,
        target_name=cfg.target_name,
        sequences=sequences,
        chain_table=chain_table(combined),
        cdr_segments=cdr_segments,
        cdr_encoding=cfg.cdr_encoding,
    )
    log(
        f"Case {case}: antibody {result.role_to_chain}, target {result.target_chain_ids or 'none'}; "
        f"{designed} CDR residues free; {total_residues} residues total"
    )
    if expand_encoding:
        bounds = "; ".join(
            f"{role} {segment['name']} {segment['min_len']}-{segment['max_len']}"
            for role, segments in cdr_segments.items()
            for segment in segments
        )
        log(f"CDRs are expandable segments (lengths: {bounds})")
    log(f"Input chains: {describe_chains(combined)}")
    return result
