import email.message
import io
import urllib.error

import pytest

from rfd4affinity import net
from rfd4affinity.seqprops import liability_motifs
from rfd4affinity.spec import DEFAULT_METRIC_TAGS, ORACLE_METRIC_TAGS, RunConfig, build_run
from rfd4affinity.structio import map_sequence_positions, read_structure
from rfd4affinity.util import (
    PipelineError,
    cif_data_name,
    log,
    read_sequence_arg,
    set_log_file,
)

AA3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def test_default_metric_tags_are_supported_by_the_hosted_nim(complex_pdb, tmp_path):
    cfg = RunConfig(
        mode="vhh",
        complex_structure=str(complex_pdb),
        antibody_chains="A",
        target_chains="T",
        cdrs="all",
        epitope="5-10",
    )
    run = build_run(cfg, api_key=None, out_dir=tmp_path / "run")
    tags = set(run.request["metrics"]["tags"])
    assert tags <= set(DEFAULT_METRIC_TAGS) | {"hotspot_conditioning", "interchain_interactions", "internal_coordinate_conditioning"}
    # These need PyRosetta / the structure-prediction oracle and are rejected
    # by the hosted deployment with HTTP 422.
    assert not tags & ORACLE_METRIC_TAGS


def test_num_steps_respects_the_server_limit(complex_pdb, tmp_path):
    cfg = RunConfig(
        mode="vhh",
        complex_structure=str(complex_pdb),
        antibody_chains="A",
        target_chains="T",
        num_steps=500,
    )
    with pytest.raises(PipelineError):
        build_run(cfg, api_key=None, out_dir=tmp_path / "run")


def test_long_raw_sequence_is_not_treated_as_a_file_name():
    seq = "A" * 400
    assert read_sequence_arg(seq) == seq


def test_sequence_argument_accepts_inline_fasta_and_record_selector(tmp_path):
    inline = ">vh_desc\nEVQLLESGGG\nLVQPGGS\n"
    assert read_sequence_arg(inline) == "EVQLLESGGGLVQPGGS"

    fasta = tmp_path / "two.fasta"
    fasta.write_text(">VH\nEVQLLESGGG\n>VL\nELVMTQSP\n", encoding="utf-8")
    assert read_sequence_arg(f"{fasta}#VH") == "EVQLLESGGG"
    assert read_sequence_arg(f"{fasta}#vl") == "ELVMTQSP"
    with pytest.raises(PipelineError) as exc:
        read_sequence_arg(str(fasta))
    assert "records" in str(exc.value)
    with pytest.raises(PipelineError):
        read_sequence_arg(f"{fasta}#nope")
    with pytest.raises(PipelineError):
        read_sequence_arg(str(tmp_path / "missing.fasta"))


def test_cif_data_name_always_starts_with_a_letter():
    assert cif_data_name("8UCD_vhh").startswith("d")
    assert cif_data_name("8UCD_vhh")[0].isalpha()
    assert cif_data_name("target")[0].isalpha()


def test_label_seq_id_is_contiguous_and_auth_seq_id_keeps_author_numbering(tmp_path, sequences):
    lines = []
    serial = 1
    for index, aa in enumerate(sequences["vhh"], start=26):
        for atom, dx in (("N", 0.0), ("CA", 1.0), ("C", 2.0), ("O", 3.0), ("CB", 1.5)):
            x = index * 3.0 + dx
            lines.append(
                f"ATOM  {serial:5d} {atom:<4s} {AA3[aa]:>3s} A{index:4d}    "
                f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 80.00           {atom[0]:>2s}"
            )
            serial += 1
    lines.append("END")
    pdb = tmp_path / "numbered.pdb"
    pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")

    chains = read_structure(pdb)
    assert chains["A"].residues[0].resseq == 26

    cfg = RunConfig(
        mode="vhh",
        antibody_structure=str(pdb),
        antibody_chains="A",
        cdrs="all",
        cdr_encoding="residues",
    )
    run = build_run(cfg, api_key=None, out_dir=tmp_path / "run")
    label_ids = [
        int(line.split()[7])
        for line in run.cif_text.splitlines()
        if line.startswith("ATOM")
    ]
    assert sorted(set(label_ids)) == list(range(1, len(sequences["vhh"]) + 1))
    auth_ids = {
        int(line.split()[9])
        for line in run.cif_text.splitlines()
        if line.startswith("ATOM")
    }
    assert auth_ids == set(range(26, 26 + len(sequences["vhh"])))


def test_map_sequence_positions_rejects_out_of_range_positions(complex_pdb):
    chains = read_structure(complex_pdb)
    chain = chains["A"]
    with pytest.raises(PipelineError):
        map_sequence_positions(chain.sequence(), chain, [len(chain.residues) + 1])


def test_mutated_positions_survive_the_antibody_segment_extraction(fab_pdb, tmp_path, sequences):
    """Regression: substituted residues must not be dropped from the segment.

    Aligning only the *equal* runs shortened the extracted chain, so every
    substitution was reported as a deletion instead.
    """
    from rfd4affinity.report import write_outputs

    cfg = RunConfig(
        mode="fab",
        complex_structure=str(fab_pdb),
        antibody_chains="H,L",
        target_chains="T",
        cdrs="HCDR3,LCDR1",
        num_designs=1,
    )
    run = build_run(cfg, api_key=None, out_dir=tmp_path / "run")
    vh = list(sequences["vh"])
    vl = list(sequences["vl"])
    vh[100] = "A"
    vl[30] = "W"
    response = {
        "request_id": "x",
        "designs": [
            {
                "design_id": "d0",
                "sequence": "".join(vh) + "".join(vl),
                "structure": {"format": "mmcif", "content": ""},
                "metrics": {},
            }
        ],
        "aggregated_metrics": [],
        "warnings": [],
    }
    designs = write_outputs(run, response, tmp_path / "run")
    design = designs[0]
    assert len(design["antibody_sequence"]) == len(sequences["vh"]) + len(sequences["vl"])
    substitutions = {
        (m["role"], m["pos"], m["from"], m["to"])
        for m in design["cdr_mutations"]
    }
    assert ("H", 101, sequences["vh"][100], "A") in substitutions
    assert ("L", 31, sequences["vl"][30], "W") in substitutions
    assert all(m["type"] == "substitution" for m in design["cdr_mutations"])


def test_n_glycosylation_sequon_is_not_counted_twice():
    hits = liability_motifs("ANGS")
    kinds = [hit["type"] for hit in hits]
    assert kinds.count("N-糖基化") == 1
    assert "脱酰胺" not in kinds
    # A plain NG that is not a sequon is still a deamidation site.
    hits = liability_motifs("ANG")
    assert [hit["type"] for hit in hits] == ["脱酰胺"]


def test_api_key_is_redacted_from_the_log(tmp_path):
    log_file = tmp_path / "run.log"
    set_log_file(log_file)
    try:
        log("using key nvapi-abcdefghijklmnopqrstuvwxyz0123456789 now")
    finally:
        set_log_file(None)
    text = log_file.read_text(encoding="utf-8")
    assert "nvapi-abcdefghijklmnopqrstuvwxyz0123456789" not in text
    assert "REDACTED" in text


def _http_error(code, payload):
    return urllib.error.HTTPError(
        "https://example.invalid/x",
        code,
        "boom",
        email.message.Message(),
        io.BytesIO(payload.encode("utf-8")),
    )


def test_metric_tag_422_gives_an_actionable_message(monkeypatch):
    def fake_urlopen(request, timeout=None):
        raise _http_error(
            422,
            '{"error":{"message":"metric tag(s) [\'developability\'] require PyRosetta"}}',
        )

    monkeypatch.setattr(net.urllib.request, "urlopen", fake_urlopen)
    with pytest.raises(PipelineError) as excinfo:
        net.post_json("https://example.invalid/x", {"a": 1}, api_key="nvapi-x", retries=2)
    message = str(excinfo.value)
    assert "metric tag" in message
    assert "developability" in message


def test_backend_crash_is_not_retried(monkeypatch):
    calls = []

    def fake_urlopen(request, timeout=None):
        calls.append(1)
        raise _http_error(
            500,
            '{"error":{"message":"AcceleratorError: CUDA error: an illegal memory access was encountered"}}',
        )

    monkeypatch.setattr(net.urllib.request, "urlopen", fake_urlopen)
    with pytest.raises(PipelineError) as excinfo:
        net.post_json("https://example.invalid/x", {"a": 1}, api_key="nvapi-x", retries=4)
    assert len(calls) == 1
    assert "server-side" in str(excinfo.value)


def test_expand_cdr_encoding_produces_pseudo_residues(complex_pdb, tmp_path):
    """Expand mode replaces CDR residues with pseudo-residues carrying length bounds."""
    cfg = RunConfig(
        mode="vhh",
        complex_structure=str(complex_pdb),
        antibody_chains="A",
        target_chains="T",
        target_name="HER2",
        cdrs="all",
        cdr_encoding="expand",
    )
    run = build_run(cfg, api_key=None, out_dir=tmp_path / "run")
    lines = run.cif_text.splitlines()
    # Parse columns
    j = next(k for k, l in enumerate(lines) if l.startswith("_atom_site."))
    cols = []
    while lines[j].startswith("_atom_site."):
        cols.append(lines[j][len("_atom_site."):].strip()); j += 1
    idx = {c: k for k, c in enumerate(cols)}
    # Find pseudo-residue rows (HETATM with non-AA comp_id)
    pseudo_rows = []
    for line in lines[j:]:
        if not line.startswith("HETATM"):
            continue
        parts = line.split()
        if "-" in parts[idx["label_comp_id"]]:
            pseudo_rows.append(parts)
    # Expect 3 pseudo-residues (CDR1, CDR2, CDR3)
    assert len(pseudo_rows) == 3
    for parts in pseudo_rows:
        comp_id = parts[idx["label_comp_id"]]
        assert "-" in comp_id, f"Expected min-max comp_id, got {comp_id}"
        assert parts[idx["mask_expsegmin_1_residue"]] == "True"
        assert parts[idx["mask_expsegmax_1_residue"]] == "True"
        assert parts[idx["mask_sequence_1_residue"]] == "False"
        assert parts[idx["mask_coordinate_1_atom"]] == "False"
        assert parts[idx["condition_CoordinateGroup_1_atom"]] == "-1"
    # CDR residues should not appear as regular ATOM with mask_sequence=False
    atom_rows = [l.split() for l in lines[j:] if l.startswith("ATOM")]
    for parts in atom_rows:
        comp_id = parts[idx["label_comp_id"]]
        if comp_id.isdigit():  # framework residue
            continue
        # Only standard amino acids should have mask_sequence=False if they are CDR
        # but in expand mode, CDR residues are not in the output at all
        assert parts[idx["label_comp_id"]].isalpha()


def test_cdr_length_range_absolute_and_relative(complex_pdb, tmp_path):
    """cdr-length-range supports absolute and relative bounds."""
    cfg = RunConfig(
        mode="vhh",
        complex_structure=str(complex_pdb),
        antibody_chains="A",
        target_chains="T",
        target_name="HER2",
        cdrs="all",
        cdr_encoding="expand",
        cdr_length_range="CDR3:12-21,*:+4",
    )
    run = build_run(cfg, api_key=None, out_dir=tmp_path / "run")
    lines = run.cif_text.splitlines()
    j = next(k for k, l in enumerate(lines) if l.startswith("_atom_site."))
    while lines[j].startswith("_atom_site."):
        j += 1
    pseudo_rows = []
    for line in lines[j:]:
        if not line.startswith("HETATM"):
            continue
        parts = line.split()
        if "-" in parts[4]:
            pseudo_rows.append(parts)
    # Should have 3 pseudo-residues
    assert len(pseudo_rows) == 3
    # CDR3 should have 12-21 (absolute), CDR1/2 should have native+4
    # Native lengths from test fixture: CDR1=8, CDR2=8, CDR3=18
    bounds = {}
    for parts in pseudo_rows:
        comp = parts[4]
        lo, hi = map(int, comp.split("-"))
        bounds[comp] = (lo, hi)
    # Find the 12-21 one (CDR3)
    assert "12-21" in bounds
    # The other two should be native+4: 8+4=12
    for b in bounds.values():
        if b != (12, 21):
            assert b == (8, 12), f"Expected 8-12 or 12-21, got {b}"
