import csv

from rfd4affinity.report import diff_sequences, write_outputs
from rfd4affinity.seqprops import gravy, isoelectric_point, net_charge
from rfd4affinity.spec import RunConfig, build_run


def test_diff_sequences():
    mutations = diff_sequences("ABCDE", "ABXDE")
    assert mutations == [{"pos": 3, "from": "C", "to": "X", "type": "substitution"}]


def test_seqprops_smoke(sequences):
    seq = sequences["vhh"]
    assert -1.0 < gravy(seq) < 1.0
    assert 6.0 < isoelectric_point(seq) < 12.0
    assert net_charge(seq) > 0


def test_report_outputs(complex_pdb, tmp_path, fake_response, sequences):
    cfg = RunConfig(
        mode="vhh",
        complex_structure=str(complex_pdb),
        antibody_chains="A",
        target_chains="T",
        cdrs="CDR3",
        num_designs=1,
        epitope="5-10",
    )
    run = build_run(cfg, api_key=None, out_dir=tmp_path / "run")
    response = fake_response(run.cif_text, sequences["vhh"])
    designs = write_outputs(run, response, tmp_path / "run")

    assert len(designs) == 1
    assert designs[0]["n_cdr_mutations"] >= 1
    assert designs[0]["n_framework_mutations"] == 0
    assert designs[0]["plddt"] is not None

    out = tmp_path / "run"
    assert (out / "report.md").is_file()
    assert (out / "report.html").is_file()
    assert (out / "designs.json").is_file()
    assert (out / "structures" / "design_0.cif").is_file()
    text = (out / "report.md").read_text(encoding="utf-8")
    assert "RFD4-Proteina" in text and "design_0" in text
    assert "湿实验建议" in text

    with (out / "designs.csv").open(encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    assert rows and rows[0]["design_id"] == "design_0"
    assert rows[0]["framework_mutations"] == "无"


def test_fab_design_extraction(fab_pdb, tmp_path, sequences):
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
    vh[100] = "A"   # HCDR3 (97-109)
    vl[30] = "W"    # LCDR1 (27-32)
    response = {
        "request_id": "fab",
        "model": {"model_name": "rfd4-proteina"},
        "runtime": {},
        "designs": [
            {
                "design_id": "fab_0",
                "sequence": "".join(vh) + "".join(vl),
                "structure": {"format": "mmcif", "content": ""},
                "metrics": {},
            }
        ],
        "aggregated_metrics": [],
        "warnings": [],
    }
    designs = write_outputs(run, response, tmp_path / "run")
    assert len(designs) == 1
    roles = {mutation["role"] for mutation in designs[0]["cdr_mutations"]}
    assert roles == {"H", "L"}
    assert designs[0]["n_framework_mutations"] == 0
