from rfd4affinity.spec import RunConfig, build_run


def _cif_rows(cif_text):
    lines = cif_text.splitlines()
    columns = [line[len("_atom_site."):] for line in lines if line.startswith("_atom_site.")]
    rows = [dict(zip(columns, line.split())) for line in lines if line.startswith("ATOM")]
    return columns, rows


def test_complex_case(complex_pdb, tmp_path):
    cfg = RunConfig(
        mode="vhh",
        complex_structure=str(complex_pdb),
        antibody_chains="A",
        target_chains="T",
        target_name="HER2",
        cdrs="all",
        num_designs=2,
        epitope="5-10",
    )
    run = build_run(cfg, api_key=None, out_dir=tmp_path / "run")
    assert run.case == "complex"
    assert run.role_to_chain == {"A": "A"}
    assert run.target_chain_ids == ["T"]
    assert run.request["task"] == "antibody_design"
    assert run.request["input"]["format"] == "annotated_cif"
    assert run.request["sampling"]["num_replicates"] == 2
    assert "hotspot_conditioning" in run.request["metrics"]["tags"]

    columns, rows = _cif_rows(run.cif_text)
    assert "mask_coordinate_1_atom" in columns
    assert "condition_CoordinateGroup_1_atom" in columns
    cdr = [row for row in rows if row["auth_asym_id"] == "A" and row["auth_seq_id"] == "27"]
    assert cdr and all(row["mask_sequence_1_residue"] == "False" for row in cdr)
    framework = [row for row in rows if row["auth_asym_id"] == "A" and row["auth_seq_id"] == "1"]
    assert framework and all(row["mask_sequence_1_residue"] == "True" for row in framework)
    target = [row for row in rows if row["auth_asym_id"] == "T" and row["auth_seq_id"] == "5"]
    assert target
    by_atom = {row["auth_atom_id"]: row for row in target}
    assert by_atom["CB"]["mask_hotspot_1_atom"] == "True"
    assert by_atom["CA"]["mask_hotspot_1_atom"] == "False"
    assert all(row["mask_coordinate_1_atom"] == "True" for row in target)


def test_docking_case(split_pdbs, tmp_path):
    ab, target = split_pdbs
    cfg = RunConfig(
        mode="vhh",
        antibody_structure=str(ab),
        antibody_chains="A",
        target_structure=str(target),
        target_chains="T",
        cdrs="CDR3",
    )
    run = build_run(cfg, api_key=None, out_dir=tmp_path / "run")
    assert run.case == "docking"
    _, rows = _cif_rows(run.cif_text)
    framework = [row for row in rows if row["auth_asym_id"] == "A" and row["auth_seq_id"] == "1"]
    backbone = next(row for row in framework if row["auth_atom_id"] == "CA")
    assert backbone["condition_CoordinateGroup_1_atom"] == "1"
    assert backbone["mask_InternalCoordinate_1_atom"] == "True"
    assert backbone["condition_index-group_1_token"] == "1"
    target_ca = next(
        row for row in rows
        if row["auth_asym_id"] == "T" and row["auth_seq_id"] == "1" and row["auth_atom_id"] == "CA"
    )
    assert target_ca["mask_coordinate_1_atom"] == "True"


def test_fab_case(fab_pdb, tmp_path):
    cfg = RunConfig(
        mode="fab",
        complex_structure=str(fab_pdb),
        antibody_chains="H,L",
        target_chains="T",
        cdrs="HCDR3,LCDR1",
    )
    run = build_run(cfg, api_key=None, out_dir=tmp_path / "run")
    assert run.role_to_chain == {"H": "H", "L": "L"}
    assert run.case == "complex"
    names = {meta["name"] for meta in run.cdr_metadata}
    assert names == {"HCDR3", "LCDR1"}
    assert "L" in run.sequences and "H" in run.sequences


def test_missing_cdr_selection_is_an_error(complex_pdb, tmp_path):
    import pytest
    from rfd4affinity.util import PipelineError

    cfg = RunConfig(
        mode="vhh",
        complex_structure=str(complex_pdb),
        antibody_chains="A",
        target_chains="T",
        cdr_ranges="A:999-1000",
    )
    with pytest.raises(PipelineError):
        build_run(cfg, api_key=None, out_dir=tmp_path / "run")


def test_apo_case(split_pdbs, tmp_path):
    ab, _ = split_pdbs
    cfg = RunConfig(
        mode="vhh",
        antibody_structure=str(ab),
        antibody_chains="A",
        cdrs="CDR3",
    )
    run = build_run(cfg, api_key=None, out_dir=tmp_path / "run")
    assert run.case == "apo"
    assert run.target_chain_ids == []
    _, rows = _cif_rows(run.cif_text)
    assert {row["auth_asym_id"] for row in rows} == {"A"}
    framework = [row for row in rows if row["auth_seq_id"] == "1"]
    backbone = next(row for row in framework if row["auth_atom_id"] == "CA")
    assert backbone["condition_CoordinateGroup_1_atom"] == "1"
    assert backbone["mask_InternalCoordinate_1_atom"] == "True"
    assert backbone["condition_index-group_1_token"] == "1"
    cdr = [row for row in rows if row["auth_seq_id"] == "105"]
    assert cdr and all(row["mask_sequence_1_residue"] == "False" for row in cdr)
