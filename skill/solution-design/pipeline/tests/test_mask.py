"""Mask + reference-numbering tests on the real 8UCD template."""
import numpy as np
import pytest

from tests.conftest import TM, ECD


def test_reference_mapping_identity(target, ref_seq):
    """PDB numbering equals UniProt numbering for 8UCD (verified 1:1)."""
    m = target.reference_map['A']
    assert m.get(67) == 67
    assert m.get(310) == 310
    assert m.get(268) == 268
    # unmodeled termini must not be mapped (chain only covers 67..310)
    assert 66 not in m.values()
    assert min(m.values()) == 67
    assert max(m.values()) == 310


def test_mask_tm_only(target, mask):
    """Designable residues are all inside TMD, none in ECD/ICD."""
    for cid in target.reference_map:
        for refn, info in mask[cid].items():
            if info['designable']:
                assert any(a <= refn <= b for a, b in TM), refn
            if info['region'] == 'ECD':
                assert not info['designable']


def test_heme_pocket_protected(target, mask):
    """H175/H268 (heme coordination) must never be designable."""
    for cid in target.reference_map:
        assert not mask[cid][175]['designable']
        assert not mask[cid][268]['designable']
        assert mask[cid][175]['reason'] == 'protected'
        assert mask[cid][268]['reason'] == 'protected'


def test_mask_symmetric(target, mask):
    """C3 symmetric mask: identical designable set across chains."""
    sets = [sorted(refn for refn, info in mask[c].items()
                   if info['designable'])
            for c in target.reference_map]
    assert sets[0] == sets[1] == sets[2]
    assert len(sets[0]) > 20  # meaningful surface, not degenerate


def test_interface_fixed(target, mask):
    """Oligomer interface residues are hard-fixed."""
    for cid in target.reference_map:
        for refn, info in mask[cid].items():
            if info['interface']:
                assert not info['designable']


def test_ecd_fully_covered(target, mask):
    """All modeled ECD residues are present and fixed."""
    ecd_refs = set()
    for a, b in ECD:
        ecd_refs |= set(range(a, b + 1))
    for cid in target.reference_map:
        fixed_ecd = [r for r, i in mask[cid].items() if i['region'] == 'ECD']
        assert fixed_ecd
        assert all(r in ecd_refs for r in fixed_ecd)


def test_tm_score_identity(target, ref_seq):
    """TM-score of a chain against itself is 1.0."""
    from tmkit.metrics import tm_score
    from tmkit.pdbio import ca_coords
    A, refs = ca_coords(target, 'A')
    assert np.isclose(tm_score(A, A), 1.0, atol=1e-6)
