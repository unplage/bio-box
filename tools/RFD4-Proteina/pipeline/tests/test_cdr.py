import pytest

from rfd4affinity.cdr import (
    CDR_BOUNDS,
    _anarci_rows,
    merge_ranges,
    number_cdrs,
    resolve_cdr_specs,
)
from rfd4affinity.util import PipelineError


def test_imgt_vhh_cdrs(sequences):
    hits = number_cdrs(sequences["vhh"], "A", scheme="imgt")
    assert [hit.cdr for hit in hits] == [1, 2, 3]
    # The reported coordinates are *sequence* positions, not IMGT numbers.
    assert (hits[0].start, hits[0].end) == (26, 33)
    assert hits[0].sequence == "GFTFSSYA"
    assert (hits[1].start, hits[1].end) == (51, 58)
    assert (hits[2].start, hits[2].end) == (97, 114)
    assert hits[2].sequence == "AKDRLSITIRPRYYGLDV"
    assert all(hit.method.startswith("anarci") for hit in hits)


def test_cdr_hits_use_sequence_positions_not_numbering_positions(sequences):
    """Regression: CDR spans must never be taken from the numbering index.

    Taking the numbering positions made CDR3 cover the conserved FR4 motif
    (``WGQGTLVTVSS``) instead of the actual HCDR3.
    """
    sequence = sequences["vhh"]
    rows = _anarci_rows(sequence, "imgt")
    assert rows is not None
    lo, hi = CDR_BOUNDS["imgt"][3]
    seq_pos = 0
    scheme_pairs = []
    for position, aa in rows:
        number = position[0] if isinstance(position, (tuple, list)) else int(position)
        if aa in ("-", " ", ".") or aa is None:
            continue
        seq_pos += 1
        scheme_pairs.append((int(number), seq_pos))

    expected = sorted(seq for scheme, seq in scheme_pairs if lo <= scheme <= hi)
    hit = [h for h in number_cdrs(sequence, "A", scheme="imgt") if h.cdr == 3][0]
    assert list(range(hit.start, hit.end + 1)) == expected
    assert "WGQGTLVTVSS" not in hit.sequence
    assert hit.sequence == sequence[hit.start - 1:hit.end]


def test_amg110_light_chain_cdr3_is_not_the_framework():
    """AMG110 VL: IMGT numbers 105-117, sequence positions 95-103.

    Using the numbering positions put the design region at 105-113, i.e. the
    conserved FR4 ``FGAGTKLEIK``.
    """
    light = (
        "ELVMTQSPSSLTVTAGEKVTMSCKSSQSLLNSGNQKNYLTWYQQKPGQPPKLLIYWASTRESGVPDRFTGSGSGTDFTLTISSVQAEDLAVYYCQNDYSYPLTFGAGTKLEIK"
    )
    hits = number_cdrs(light, "L", scheme="imgt")
    cdr3 = [hit for hit in hits if hit.cdr == 3][0]
    assert (cdr3.start, cdr3.end) == (95, 103)
    assert cdr3.sequence == "QNDYSYPLT"
    assert "FGAGTKLEIK" not in cdr3.sequence


def test_kabat_ranges_are_used(sequences):
    hits = number_cdrs(sequences["vhh"], "A", scheme="kabat")
    assert hits[0].method == "anarci:kabat"
    assert hits[0].sequence == sequences["vhh"][hits[0].start - 1:hits[0].end]


def test_manual_ranges(sequences):
    ranges, metadata, method = resolve_cdr_specs(
        {"A": sequences["vhh"]}, "vhh", cdr_ranges="A:105-110"
    )
    assert ranges == {"A": [(105, 110)]}
    assert method == "manual ranges"
    assert metadata[0]["sequence"] == sequences["vhh"][104:110]


def test_manual_ranges_reject_non_numeric(sequences):
    with pytest.raises(PipelineError):
        resolve_cdr_specs({"A": sequences["vhh"]}, "vhh", cdr_ranges="A:xx-yy")


def test_cdr_name_must_have_one_digit():
    with pytest.raises(PipelineError):
        resolve_cdr_specs({"A": "ACDEFGHIKLMNPQRSTVWY"}, "vhh", cdrs="CDR12")


def test_ambiguous_fab_cdr_requires_chain(sequences):
    with pytest.raises(PipelineError):
        resolve_cdr_specs(
            {"H": sequences["vh"], "L": sequences["vl"]}, "fab", cdrs="CDR3"
        )


def test_merge_ranges():
    assert merge_ranges([(1, 5), (6, 9), (12, 15)]) == [(1, 9), (12, 15)]
