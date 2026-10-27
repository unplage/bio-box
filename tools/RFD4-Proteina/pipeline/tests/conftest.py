import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

VHH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKG"
    "RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTLVTVSS"
)
VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
TARGET = "CQPQNGSVTCFGPEADQCVACAHYKDPPFCVARCPSIWKFPDEEGACQPCPIN"

AA3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def write_pdb(path: Path, chains, bfactor: float = 80.0) -> Path:
    lines = []
    serial = 1
    for chain_id, sequence, x0 in chains:
        for index, aa in enumerate(sequence, start=1):
            for atom, dx, dy, dz in (
                ("N", 0, 0, 0), ("CA", 1, 0, 0), ("C", 2, 0, 0),
                ("O", 3, 0, 0), ("CB", 1, 1, 0),
            ):
                x = x0 + index * 3.0 + dx
                lines.append(
                    f"ATOM  {serial:5d} {atom:<4s} {AA3[aa]:>3s} {chain_id}{index:4d}    "
                    f"{x:8.3f}{0.0:8.3f}{dz:8.3f}  1.00 {bfactor:5.2f}           {atom[0]:>2s}"
                )
                serial += 1
        lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


@pytest.fixture(scope="session")
def sequences():
    return {"vhh": VHH, "vh": VH, "vl": VL, "target": TARGET}


@pytest.fixture()
def complex_pdb(tmp_path):
    return write_pdb(
        tmp_path / "complex.pdb",
        [("A", VHH, 0.0), ("T", TARGET, 40.0)],
    )


@pytest.fixture()
def split_pdbs(tmp_path):
    ab = write_pdb(tmp_path / "ab.pdb", [("A", VHH, 0.0)])
    target = write_pdb(tmp_path / "target.pdb", [("T", TARGET, 40.0)])
    return ab, target


@pytest.fixture()
def fab_pdb(tmp_path):
    return write_pdb(
        tmp_path / "fab.pdb",
        [("H", VH, 0.0), ("L", VL, 20.0), ("T", TARGET, 80.0)],
    )


@pytest.fixture()
def fake_response(tmp_path):
    def _build(cif_text: str, sequence: str, design_id: str = "design_0"):
        mutated = list(sequence)
        for position, aa in ((105, "A"), (110, "W")):
            mutated[position - 1] = aa
        return {
            "request_id": "req_test",
            "task": "antibody_design",
            "model": {"model_name": "rfd4-proteina"},
            "runtime": {"nim_version": "1.0.0", "gpu": "H100 80GB"},
            "designs": [
                {
                    "design_id": design_id,
                    "replicate_index": 0,
                    "seed": 7,
                    "sequence": "".join(mutated),
                    "structure": {
                        "format": "mmcif",
                        "filename": f"{design_id}.cif",
                        "content": cif_text,
                    },
                    "metrics": {"design": {"pass": True, "radius_of_gyration": 12.0}},
                    "diagnostics": {"num_designed_residues": 35},
                }
            ],
            "aggregated_metrics": [
                {"metric": "design.pass", "reduction": "mean", "value": 1.0, "n": 1}
            ],
            "warnings": [],
            "duration_ms": 1000,
        }

    return _build
