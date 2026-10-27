# -*- coding: utf-8 -*-
"""Minimal, dependency-free structure I/O.

The NIM accepts **UTF-8 mmCIF only** for attached structures (PDB is rejected).
This module parses PDB/mmCIF into a small internal model and writes a clean
``_atom_site`` mmCIF.  ``annotate.py`` builds the annotated variant on top.
"""

from __future__ import annotations

import copy
import difflib
import shlex
from dataclasses import dataclass, field
from pathlib import Path

from .util import AA3_TO_1, PipelineError, warn


@dataclass
class Residue:
    name: str
    resseq: int
    icode: str = ""
    atoms: dict[str, tuple[float, float, float]] = field(default_factory=dict)
    bfactors: dict[str, float] = field(default_factory=dict)

    @property
    def one(self) -> str:
        return AA3_TO_1.get(self.name.upper(), "X")

    @property
    def label(self) -> str:
        return f"{self.name}{self.resseq}{self.icode}".strip()


@dataclass
class Chain:
    id: str
    residues: list[Residue] = field(default_factory=list)

    def sequence(self) -> str:
        return "".join(res.one for res in self.residues)

    def resseqs(self) -> list[int]:
        return [res.resseq for res in self.residues]


def _normalise_resname(name: str) -> str:
    name = (name or "UNK").strip().upper()
    if name == "MSE":
        return "MET"
    return name


def _is_standard(name: str) -> bool:
    return _normalise_resname(name) in AA3_TO_1


def _add_atom(
    chains: dict[str, Chain],
    chain_id: str,
    resseq: int,
    icode: str,
    resname: str,
    atom_name: str,
    xyz: tuple[float, float, float],
    bfactor: float = 0.0,
) -> None:
    resname = _normalise_resname(resname)
    atom_name = atom_name.strip().upper()
    if resname == "MET" and atom_name == "SE":
        atom_name = "SD"
    chain_id = chain_id.strip() or "_"
    chain = chains.setdefault(chain_id, Chain(chain_id))
    residue = None
    if chain.residues:
        last = chain.residues[-1]
        if last.resseq == resseq and last.icode == icode and last.name == resname:
            residue = last
    if residue is None:
        for candidate in reversed(chain.residues):
            if candidate.resseq == resseq and candidate.icode == icode and candidate.name == resname:
                residue = candidate
                break
    if residue is None:
        residue = Residue(resname, resseq, icode)
        chain.residues.append(residue)
    residue.atoms[atom_name] = xyz
    residue.bfactors[atom_name] = bfactor


def parse_pdb(text: str) -> dict[str, Chain]:
    chains: dict[str, Chain] = {}
    seen_model = False
    for line in text.splitlines():
        record = line[:6].strip()
        if record == "MODEL":
            if seen_model:
                break
            seen_model = True
            continue
        if record == "ENDMDL":
            break
        if record not in ("ATOM", "HETATM") or len(line) < 54:
            continue
        altloc = line[16]
        if altloc not in (" ", "A", "1"):
            continue
        resname = line[17:20].strip()
        if record == "HETATM" and resname not in AA3_TO_1:
            continue
        if not _is_standard(resname):
            continue
        chain_id = line[21].strip() or "_"
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        icode = line[26].strip()
        atom = line[12:16].strip()
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        try:
            bfactor = float(line[60:66])
        except (ValueError, IndexError):
            bfactor = 0.0
        _add_atom(chains, chain_id, resseq, icode, resname, atom, xyz, bfactor)
    if not chains:
        raise PipelineError("No standard protein atoms were found in the PDB input.")
    return chains


def _cif_tokens(line: str) -> list[str]:
    try:
        return shlex.split(line, comments=False, posix=True)
    except ValueError:
        return line.split()


def parse_cif(text: str) -> dict[str, Chain]:
    lines = text.splitlines()
    start = None
    for idx, line in enumerate(lines):
        if line.strip().startswith("_atom_site."):
            start = idx
            break
    if start is None:
        raise PipelineError("No _atom_site loop found in the mmCIF input.")

    columns: list[str] = []
    idx = start
    while idx < len(lines) and lines[idx].strip().startswith("_atom_site."):
        columns.append(lines[idx].strip()[len("_atom_site."):])
        idx += 1
    if not columns:
        raise PipelineError("Could not read _atom_site column names from mmCIF.")

    rows: list[list[str]] = []
    while idx < len(lines):
        line = lines[idx].strip()
        if not line or line.startswith("#") or line.startswith("_") or line.lower().startswith(("loop_", "data_")):
            break
        rows.append(_cif_tokens(lines[idx]))
        idx += 1

    def get(row: list[str], *names: str, default: str = "") -> str:
        for name in names:
            if name in columns:
                pos = columns.index(name)
                if pos < len(row):
                    return row[pos]
        return default

    chains: dict[str, Chain] = {}
    for row in rows:
        group = get(row, "group_PDB", default="ATOM")
        if group not in ("ATOM", "HETATM"):
            continue
        model = get(row, "pdbx_PDB_model_num", default="1")
        if model not in ("1", ".", "?"):
            continue
        altloc = get(row, "label_alt_id", default=".")
        if altloc not in (".", "?", "", "A", "1"):
            continue
        resname = get(row, "auth_comp_id", "label_comp_id", default="UNK")
        if not _is_standard(resname):
            continue
        chain_id = get(row, "auth_asym_id", "label_asym_id", default="_")
        seq_id = get(row, "auth_seq_id", "label_seq_id", default="0")
        try:
            resseq = int(float(seq_id))
        except ValueError:
            continue
        icode = get(row, "pdbx_PDB_ins_code", default="")
        if icode in (".", "?"):
            icode = ""
        atom = get(row, "auth_atom_id", "label_atom_id", default="")
        try:
            xyz = (
                float(get(row, "Cartn_x")),
                float(get(row, "Cartn_y")),
                float(get(row, "Cartn_z")),
            )
        except ValueError:
            continue
        try:
            bfactor = float(get(row, "B_iso_or_equiv", default="0"))
        except ValueError:
            bfactor = 0.0
        _add_atom(chains, chain_id, resseq, icode, resname, atom, xyz, bfactor)
    if not chains:
        raise PipelineError("No standard protein atoms were found in the mmCIF input.")
    return chains


def parse_structure_text(text: str) -> dict[str, Chain]:
    head = text.lstrip()[:200].lower()
    if head.startswith("data_") or "_atom_site." in text:
        return parse_cif(text)
    return parse_pdb(text)


def read_structure(path: Path) -> dict[str, Chain]:
    return parse_structure_text(path.read_text(encoding="utf-8", errors="replace"))


def chain_label_ids(chain: Chain) -> list[int]:
    """Per-residue ``label_seq_id`` values.

    Preserve author numbering when it is unique (matches the official
    examples); fall back to 1..N when insertion codes would collide.
    """
    resseqs = [res.resseq for res in chain.residues]
    if len(set(resseqs)) == len(resseqs):
        return resseqs
    return list(range(1, len(chain.residues) + 1))


def _element(atom_name: str) -> str:
    cleaned = atom_name.strip().upper()
    while cleaned and cleaned[0].isdigit():
        cleaned = cleaned[1:]
    if cleaned[:2] in ("CL", "BR", "SE"):
        return cleaned[:2]
    return cleaned[:1] or "C"


def write_mmcif(chains: dict[str, Chain], path: Path, data_name: str = "complex") -> None:
    """Write a plain ``_atom_site`` mmCIF (no condition annotations)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"data_{data_name}",
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.occupancy",
        "_atom_site.pdbx_formal_charge",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.pdbx_PDB_model_num",
        "_atom_site.id",
    ]
    atom_id = 0
    for entity_id, (chain_id, chain) in enumerate(chains.items(), start=1):
        label_ids = chain_label_ids(chain)
        for label_seq_id, residue in zip(label_ids, chain.residues):
            for atom_name, (x, y, z) in residue.atoms.items():
                atom_id += 1
                icode = residue.icode or "."
                bfactor = residue.bfactors.get(atom_name, 0.0)
                lines.append(
                    " ".join(
                        [
                            "ATOM",
                            _element(atom_name),
                            atom_name,
                            ".",
                            residue.name,
                            chain_id,
                            str(entity_id),
                            str(label_seq_id),
                            icode,
                            str(residue.resseq),
                            residue.name,
                            chain_id,
                            atom_name,
                            f"{bfactor:.2f}",
                            "1.00",
                            "?",
                            f"{x:.3f}",
                            f"{y:.3f}",
                            f"{z:.3f}",
                            "1",
                            str(atom_id),
                        ]
                    )
                )
    lines.append("#")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_chain(chain: Chain, new_id: str | None = None) -> Chain:
    clone = copy.deepcopy(chain)
    if new_id:
        clone.id = new_id
    return clone


def translate_chain(chain: Chain, offset: tuple[float, float, float]) -> None:
    dx, dy, dz = offset
    for residue in chain.residues:
        for name, (x, y, z) in residue.atoms.items():
            residue.atoms[name] = (x + dx, y + dy, z + dz)


def combine_chains(
    antibody: dict[str, Chain],
    antibody_map: dict[str, str],
    target: dict[str, Chain] | None = None,
    target_map: dict[str, str] | None = None,
    target_offset: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> dict[str, Chain]:
    combined: dict[str, Chain] = {}
    for old, new in antibody_map.items():
        combined[new] = copy_chain(antibody[old], new)
    if target and target_map:
        for old, new in target_map.items():
            chain = copy_chain(target[old], new)
            if target_offset != (0.0, 0.0, 0.0):
                translate_chain(chain, target_offset)
            combined[new] = chain
    return combined


def chain_table(chains: dict[str, Chain]) -> list[dict]:
    table = []
    for chain in chains.values():
        table.append(
            {
                "id": chain.id,
                "n_residues": len(chain.residues),
                "sequence": chain.sequence(),
                "resseq_first": chain.residues[0].resseq if chain.residues else None,
                "resseq_last": chain.residues[-1].resseq if chain.residues else None,
            }
        )
    return table


def mean_plddt_from_chains(chains: dict[str, Chain]) -> float | None:
    values = [
        res.bfactors.get("CA")
        for chain in chains.values()
        for res in chain.residues
        if res.bfactors.get("CA") is not None
    ]
    values = [v for v in values if v and v > 0]
    if not values:
        return None
    return sum(values) / len(values)


def mean_plddt(pdb_text: str) -> tuple[float | None, list[float]]:
    values: list[float] = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            values.append(float(line[60:66]))
        except (ValueError, IndexError):
            continue
    if not values:
        return None, []
    return sum(values) / len(values), values


def map_sequence_positions(
    query_seq: str,
    chain: Chain,
    positions: list[int],
) -> list[int]:
    """Map 1-based positions in *query_seq* to residue numbers in *chain*."""
    chain_seq = chain.sequence()
    if query_seq == chain_seq:
        return [chain.residues[p - 1].resseq for p in positions if 1 <= p <= len(chain_seq)]

    matcher = difflib.SequenceMatcher(None, query_seq, chain_seq, autojunk=False)
    ratio = matcher.ratio()
    if ratio < 0.80:
        raise PipelineError(
            f"Provided sequence does not match structure chain {chain.id} "
            f"(identity {ratio:.0%}); check --vh-seq/--vl-seq/--vhh-seq and "
            "--antibody-chains."
        )
    if ratio < 0.98:
        warn(
            f"Sequence-to-chain {chain.id} identity is {ratio:.1%}; CDR mapping "
            "uses alignment and should be reviewed."
        )
    mapping: dict[int, int] = {}
    for block in matcher.get_matching_blocks():
        for offset in range(block.size):
            mapping[block.a + offset] = block.b + offset
    out: list[int] = []
    for pos in positions:
        idx = pos - 1
        if idx not in mapping:
            raise PipelineError(
                f"Sequence position {pos} could not be mapped onto chain {chain.id}; "
                "provide a matching structure or use --cdr-ranges on the structure sequence."
            )
        out.append(chain.residues[mapping[idx]].resseq)
    return out


def make_offset(distance: float = 60.0) -> tuple[float, float, float]:
    return (float(distance), 0.0, 0.0)
