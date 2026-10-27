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
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .structio import Chain, _element, chain_label_ids
from .util import slugify

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
class ChainPlan:
    chain_id: str
    role: str                       # 'antibody' | 'target'
    free_residues: set[int] = field(default_factory=set)
    hotspot_residues: set[int] = field(default_factory=set)
    token_group: int = 1

    def is_free(self, resseq: int) -> bool:
        return resseq in self.free_residues


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


def write_annotated_mmcif(
    chains: dict[str, Chain],
    plans: dict[str, ChainPlan],
    path: Path,
    case: str,
    fix_sidechains: bool = False,
    data_name: str = "complex",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data_name = slugify(data_name, default="complex")
    lines = [f"data_{data_name}", "#", "loop_"]
    lines.extend(f"_atom_site.{column}" for column in COLUMNS)

    atom_id = 0
    for entity_id, (chain_id, chain) in enumerate(chains.items(), start=1):
        plan = plans.get(chain_id)
        if plan is None:
            raise ValueError(f"No annotation plan for chain {chain_id}")
        label_ids = chain_label_ids(chain)
        for label_seq_id, residue in zip(label_ids, chain.residues):
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
