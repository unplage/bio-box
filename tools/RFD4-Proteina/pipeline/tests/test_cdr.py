import pytest

from rfd4affinity.cdr import (
    CDR_BOUNDS,
    merge_ranges,
    number_cdrs,
    resolve_cdr_specs,
)
from rfd4affinity.util import PipelineError


def test_imgt_vhh_cdrs(sequences):
    hits = number_cdrs(sequences["vhh"], "A", scheme="imgt")
    assert [hit.cdr for hit in hits] == [1, 2, 3]
    assert hits[0].start == 27 and hits[0].end == 38
    assert hits[0].sequence == "FTFSSYAMSWVR"
    assert hits[2].sequence.startswith("IRPR")
    assert all(hit.method.startswith("anarci") for hit in hits)


def test_kabat_ranges_are_used(sequences):
    hits = number_cdrs(sequences["vhh"], "A", scheme="kabat")
    lo, hi = CDR_BOUNDS["kabat"][1]
    assert lo <= hits[0].start <= hits[0].end <= hi


def test_manual_ranges(sequences):
    ranges, metadata, method = resolve_cdr_specs(
        {"A": sequences["vhh"]}, "vhh", cdr_ranges="A:105-110"
    )
    assert ranges == {"A": [(105, 110)]}
    assert method == "manual ranges"
    assert metadata[0]["sequence"] == sequences["vhh"][104:110]


def test_ambiguous_fab_cdr_requires_chain(sequences):
    with pytest.raises(PipelineError):
        resolve_cdr_specs(
            {"H": sequences["vh"], "L": sequences["vl"]}, "fab", cdrs="CDR3"
        )


def test_merge_ranges():
    assert merge_ranges([(1, 5), (6, 9), (12, 15)]) == [(1, 9), (12, 15)]
