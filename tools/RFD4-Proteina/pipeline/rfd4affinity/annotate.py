# -*- coding: utf-8 -*-
"""Build the *annotated* mmCIF that RFD4-Proteina consumes.

We use the canonical ``annotated_cif`` dialect (rather than a condition-spec
JSON) because it is fully self-describing: every condition is a per-atom or
per-residue column next to the atom it applies to.

Column semantics (verified against the official antibody example shipped with
the NIM):

``mask_coordinate_1_atom``
    True  -> absolute coordinate condition (C_CRD) on this atom.
``condition_CoordinateGroup_1_atom``
    C_CRDGRP: -1 unknown, 0 absolute, >=1 an internal-coordinate group.
``mask_InternalCoordinate_1_atom``
    True  -> C_INTCRD target prepared for this atom (used for group >= 1).
``mask_sequence_1_residue``
    True  -> sequence is conditioned/fixed (C_SEQ); False -> designed.
``mask_hotspot_1_atom``
    True  -> hotspot atom (C_HOT).
``condition_index-group_1_token``
    Chain/token-level coordinate-group id used by the tokenizer.
``mask_index_1_residue``
    True  -> residue participates in indexing.

The default recipe: the target is held all-atom (absolute), the antibody
framework backbone is held (absolute for a complex, internal-coordinate group
for docking/apo) while its sidechains repack, and the selected CDRs are free in
both sequence and structure.

Two ways of expressing "this CDR is to be designed":

``residues`` (legacy)
    The CDR residues stay in the file with ``mask_sequence_1_residue=False``.
    The model redesigns the residues in place, so the CDR keeps its length.

``expand`` (what NVIDIA's own antibody examples do)
    The CDR residues are *not in the input at all*.  Each CDR is replaced by a
    single ``HETATM UNK`` pseudo-residue whose ``label_comp_id`` is
    ``"<min>-<max>"`` and whose ``annotation_expsegmin/max`` columns carry the
    allowed length bounds.  The model then re-generates the loop and may pick
    a different length for every candidate - this is what the official VHH
    tutorial describes ("the CDRs are not in the input at all. They are
    expandable segments: stretches of the chain whose length the model chooses
    within the bounds recorded in the input").
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .structio import Chain, _element
from .util import cif_data_name

BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

COLUMNS = [
    "group_PDB",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_PDB_ins_code",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "id",
    "B_iso_or_equiv",
    "occupancy",
    "pdbx_formal_charge",
    "mask_index_1_residue",
    "mask_hotspot_1_atom",
    "mask_InternalCoordinate_1_atom",
    "mask_sequence_1_residue",
    "condition_index-group_1_token",
    "condition_CoordinateGroup_1_atom",
    "mask_coordinate_1_atom",
    "annotation_expsegmax_1_residue",
    "mask_expsegmax_1_residue",
    "mask_expsegmin_1_residue",
    "annotation_expsegmin_1_residue",
    "atom_array_index",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "pdbx_PDB_model_num",
]


@dataclass
class LoopSegment:
    """A CDR that is handed to the model as an expandable segment.

    ``start``/``end`` are the residue numbers of the CDR in the parent
    structure; ``min_len``/``max_len`` are the length bounds written into the
    annotated mmCIF.
    """

    start: int
    end: int
    min_len: int
    max_len: int
    name: str = ""


@dataclass
class ChainPlan:
    chain_id: str
    role: str                       # 'antibody' | 'target'
    free_residues: set[int] = field(default_factory=set)
    hotspot_residues: set[int] = field(default_factory=set)
    token_group: int = 1
    loop_segments: list[LoopSegment] = field(default_factory=list)

    def is_free(self, resseq: int) -> bool:
        return resseq in self.free_residues

    def loop_for(self, resseq: int) -> LoopSegment | None:
        for segment in self.loop_segments:
            if segment.start <= resseq <= segment.end:
                return segment
        return None


def _atom_annotations(
    plan: ChainPlan,
    resseq: int,
    atom_name: str,
    case: str,
    fix_sidechains: bool,
) -> tuple[bool, int, bool, bool, bool]:
    """Return (mask_coordinate, group, internal, mask_sequence, hotspot)."""
    backbone = atom_name in BACKBONE_ATOMS
    free_residue = plan.is_free(resseq)
    # Hotspots are sidechain contacts; matching the official example keeps the
    # backbone out of C_HOT.
    hotspot = (resseq in plan.hotspot_residues) and not backbone

    if free_residue:
        return False, -1, False, False, hotspot

    if plan.role == "target":
        # The antigen is never redesigned: keep every target atom absolute so
        # the epitope chemistry cannot drift.
        return True, 0, False, True, hotspot

    # antibody
    if case == "complex":
        if backbone or fix_sidechains:
            return True, 0, False, True, hotspot
        return False, -1, False, True, hotspot

    # docking / apo: framework keeps its internal geometry but floats freely
    if backbone:
        return False, 1, True, True, hotspot
    if fix_sidechains:
        return False, 1, True, True, hotspot
    return False, -1, False, True, hotspot


def _loop_row(
    chain_id: str,
    entity_id: int,
    label_seq_id: int,
    segment: LoopSegment,
    atom_id: int,
    token_group: int,
) -> list[str]:
    """Build the ``HETATM UNK`` row that replaces a CDR in the input.

    The row mirrors the official antibody examples: ``label_comp_id`` is the
    ``"<min>-<max>"`` length range, the coordinates are ``nan`` and the
    expseg columns carry the bounds.
    """
    comp_id = f"{segment.min_len}-{segment.max_len}"
    return [
        "HETATM",
        "X",
        "UNK",
        ".",
        comp_id,
        chain_id,
        str(entity_id),
        str(label_seq_id),
        ".",
        str(segment.start),
        comp_id,
        chain_id,
        "UNK",
        str(atom_id),
        "nan",
        "0.0",
        "?",
        "True",
        "False",
        "False",
        "False",
        str(token_group),
        "-1",
        "False",
        str(segment.max_len),
        "True",
        "True",
        str(segment.min_len),
        str(atom_id - 1),
        "nan",
        "nan",
        "nan",
        "1",
    ]


def write_annotated_mmcif(
    chains: dict[str, Chain],
    plans: dict[str, ChainPlan],
    path: Path,
    case: str,
    fix_sidechains: bool = False,
    data_name: str = "complex",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data_name = cif_data_name(data_name, default="complex")
    lines = [f"data_{data_name}", "#", "loop_"]
    lines.extend(f"_atom_site.{column}" for column in COLUMNS)

    atom_id = 0
    for entity_id, (chain_id, chain) in enumerate(chains.items(), start=1):
        plan = plans.get(chain_id)
        if plan is None:
            raise ValueError(f"No annotation plan for chain {chain_id}")
        label_seq_id = 0
        for residue in chain.residues:
            segment = plan.loop_for(residue.resseq)
            if segment is not None:
                if residue.resseq != segment.start:
                    # The CDR is not part of the input any more.
                    continue
                label_seq_id += 1
                atom_id += 1
                lines.append(
                    " ".join(
                        _loop_row(
                            chain_id,
                            entity_id,
                            label_seq_id,
                            segment,
                            atom_id,
                            plan.token_group,
                        )
                    )
                )
                continue
            label_seq_id += 1
            for atom_name, (x, y, z) in residue.atoms.items():
                atom_id += 1
                mask_coord, group, internal, mask_seq, hotspot = _atom_annotations(
                    plan, residue.resseq, atom_name, case, fix_sidechains
                )
                bfactor = residue.bfactors.get(atom_name, 0.0)
                row = [
                    "ATOM",
                    _element(atom_name),
                    atom_name,
                    ".",
                    residue.name,
                    chain_id,
                    str(entity_id),
                    str(label_seq_id),
                    residue.icode or ".",
                    str(residue.resseq),
                    residue.name,
                    chain_id,
                    atom_name,
                    str(atom_id),
                    f"{bfactor:.2f}",
                    "1.00",
                    "?",
                    "True",
                    str(hotspot),
                    str(internal),
                    str(mask_seq),
                    str(plan.token_group),
                    str(group),
                    str(mask_coord),
                    "-1",
                    "False",
                    "False",
                    "-1",
                    str(atom_id - 1),
                    f"{x:.3f}",
                    f"{y:.3f}",
                    f"{z:.3f}",
                    "1",
                ]
                lines.append(" ".join(row))
    lines.append("#")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
