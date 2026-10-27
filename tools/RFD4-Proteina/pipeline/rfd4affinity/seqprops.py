# -*- coding: utf-8 -*-
"""Sequence-level developability descriptors (no third-party deps)."""

from __future__ import annotations

import re

KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

MW = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886,
    "C": 103.1388, "E": 129.1155, "Q": 128.1307, "G": 57.0519,
    "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
    "M": 131.1926, "F": 147.1766, "P": 97.1167, "S": 87.0782,
    "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}

PKA_SIDE = {"D": 3.65, "E": 4.25, "C": 8.18, "Y": 10.07, "H": 6.00, "K": 10.53, "R": 12.48}
PKA_NTERM = 9.69
PKA_CTERM = 2.34


def gravy(sequence: str) -> float:
    values = [KD.get(aa, 0.0) for aa in sequence]
    return sum(values) / len(values) if values else 0.0


def molecular_weight(sequence: str) -> float:
    return sum(MW.get(aa, 0.0) for aa in sequence) + 18.01524


def net_charge(sequence: str, ph: float = 7.4) -> float:
    charge = 0.0
    if sequence:
        charge += 1.0 / (1.0 + 10 ** (ph - PKA_NTERM))
        charge -= 1.0 / (1.0 + 10 ** (PKA_CTERM - ph))
    for aa in sequence:
        pka = PKA_SIDE.get(aa)
        if pka is None:
            continue
        if aa in ("K", "R", "H"):
            charge += 1.0 / (1.0 + 10 ** (ph - pka))
        else:
            charge -= 1.0 / (1.0 + 10 ** (pka - ph))
    return charge


def isoelectric_point(sequence: str) -> float:
    low, high = 0.0, 14.0
    for _ in range(80):
        mid = (low + high) / 2.0
        charge = net_charge(sequence, mid)
        if charge > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def aromaticity(sequence: str) -> float:
    if not sequence:
        return 0.0
    return sum(1 for aa in sequence if aa in "FWY") / len(sequence)


def liability_motifs(sequence: str) -> list[dict]:
    seq = sequence.upper()
    hits: list[dict] = []

    def add(kind: str, position: int, motif: str, note: str) -> None:
        hits.append(
            {
                "type": kind,
                "position": position,
                "motif": motif,
                "note": note,
            }
        )

    # N-glycosylation sequons first: an Asn inside one is not counted again as
    # a deamidation site, otherwise every sequon is reported twice.
    sequon_starts = set()
    for match in re.finditer(r"N[^P][ST]", seq):
        sequon_starts.add(match.start())
        add("N-糖基化", match.start() + 1, match.group(0), "N-X-S/T，建议突变以降低异质性")
    for match in re.finditer(r"N[GSQTH]", seq):
        if match.start() in sequon_starts:
            continue
        add("脱酰胺", match.start() + 1, match.group(0), "Asn 脱酰胺热点，关注长期稳定性")
    for match in re.finditer(r"D[GP]", seq):
        add("异构化/酸水解", match.start() + 1, match.group(0), "Asp 异构化或 DP 酸性切割热点")
    for match in re.finditer(r"M", seq):
        add("氧化", match.start() + 1, "M", "Met 易氧化，必要时替换为 Leu/Ile")
    cys_positions = [m.start() + 1 for m in re.finditer(r"C", seq)]
    if len(cys_positions) % 2:
        add(
            "半胱氨酸配对异常",
            cys_positions[-1],
            "C",
            f"共 {len(cys_positions)} 个 Cys（奇数），确认是否存在游离半胱氨酸",
        )
    if seq[:1] in ("E", "Q"):
        add("N端焦谷氨酸", 1, seq[:1], "N 端 E/Q 可能环化为焦谷氨酸")
    return hits


def cdr_composition(sequence: str) -> dict:
    return {
        "length": len(sequence),
        "gravy": round(gravy(sequence), 2),
        "net_charge_pH7.4": round(net_charge(sequence), 2),
        "aromaticity": round(aromaticity(sequence), 3),
        "cysteines": sequence.count("C"),
        "glycines": sequence.count("G"),
        "prolines": sequence.count("P"),
    }
